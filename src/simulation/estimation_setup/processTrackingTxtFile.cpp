/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <algorithm>

#include "tudat/simulation/estimation_setup/processTrackingTxtFile.h"

namespace tudat
{
namespace observation_models
{

std::vector<ObservableType> findAvailableObservableTypes(const std::vector<input_output::TrackingDataType> availableDataTypes)
{
  // Initialise container for available types
  std::vector<ObservableType> availableObservableTypes;

  // Loop over map with observables and their required data types. Add observabletype to vector if those data types are present
  for (const auto& pair : observableRequiredDataTypesMap) {
    std::vector<input_output::TrackingDataType> requiredDataTypeSet = pair.second;
    if (utilities::containsAll(availableDataTypes, requiredDataTypeSet))
      availableObservableTypes.push_back(pair.first);
  }
  return availableObservableTypes;
}

void ProcessedTrackingTxtFileContents::updateObservations()
{
  // Update the observableTypes that one can expect to process
  observationMap_.clear();

  for (const ObservableType observableType : observableTypes_) {
    std::vector<double> observableValues;

    // Convert the raw data to required observables
    // TODO: This function could also take into account the metadata
    switch (observableType) {
      case n_way_range: {
        // Conversion function for n-way range
        auto lightTimeRangeConversion = [](double lightTime, double lightTimeDelay) {
          return (lightTime - lightTimeDelay) * physical_constants::SPEED_OF_LIGHT;
        };

        // Extract columns from raw data and convert to observable values
        std::vector<double> lightTimes = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::n_way_light_time);
        std::vector<double> lightTimeDelays = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::light_time_measurement_delay, 0.0);
        observableValues = utilities::convertVectors(lightTimeRangeConversion, lightTimes, lightTimeDelays);
        break;
      }
      case doppler_measured_frequency: {
        // Conversion function for doppler measured frequency
        auto dopplerFrequencyConversion = [](double dopplerFrequency, double dopplerBaseFrequency) {
          return dopplerFrequency + dopplerBaseFrequency;
        };
        
        // Extract columns from raw data and convert to observable values
        std::vector<double> dopplerFrequencies = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::doppler_measured_frequency);
        std::vector<double> dopplerBaseFrequencies = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::doppler_base_frequency);
        
        // Check if any of the base frequencies are zero
        if (std::any_of(dopplerBaseFrequencies.begin(), dopplerBaseFrequencies.end(), [](double baseFrequency) { return baseFrequency == 0.0; })) {
          std::cerr << "Warning when processing doppler_measured_frequency. Doppler base frequency is zero." << std::endl;
        }

        observableValues = utilities::convertVectors(dopplerFrequencyConversion, dopplerFrequencies, dopplerBaseFrequencies);
        break;
      }
      default: {
        throw std::runtime_error("Error while processing tracking txt file. ObservableType conversion not implemented");
      }
    }

    // Store observables
    observationMap_[observableType] = observableValues;
  }

}

void ProcessedTrackingTxtFileContents::updateObservationTimes()
{
  // Clear any previous values
  observationTimes_.clear();

  // Get data map and time representation
  const auto& numDataRows = rawTrackingTxtFileContents_->getNumRows();
  TimeRepresentation timeRepresentation = getTimeRepresentation();

  // Depending on the time representation, convert further to tdb seconds since j2000
  switch (timeRepresentation) {
    case tdb_seconds_j2000: {
      observationTimes_ = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::tdb_time_j2000);
      break;
    }
    case calendar_day_time: {
      // Convert dates to Julian days since J2000
      std::vector<double> observationJulianDaysSinceJ2000 = utilities::convertVectors(
          basic_astrodynamics::convertCalendarDateToJulianDaySinceJ2000<double>,
          rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::year),
          rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::month),
          rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::day),
          rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::hour),
          rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::minute),
          rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::second)
      );
      // Convert to seconds and add to utc times
      std::vector<double> observationTimesUtc;
      for (double julianDaySinceJ2000 : observationJulianDaysSinceJ2000) {
        observationTimesUtc.push_back(julianDaySinceJ2000 * physical_constants::JULIAN_DAY);
      }
      // Convert to TDB
      observationTimes_ = computeObservationTimesTdbFromJ2000(observationTimesUtc);

      break;
    }
      // Throw error if representation not implemented
    default: {
      throw std::runtime_error("Error while processing tracking txt file: Time representation not recognised or implemented.");
    }
  }

  // Get the delays in the time tag (or set to 0.0 if not specified)
  std::vector<double> timeTagDelays = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::time_tag_delay, 0.0);
  for (size_t idx = 0; idx < observationTimes_.size(); ++idx) {
    observationTimes_[idx] -= timeTagDelays[idx];
  }
}

std::vector<double> ProcessedTrackingTxtFileContents::computeObservationTimesTdbFromJ2000(std::vector<double> observationTimesUtc)
{
  // Get the timescale converter
  earth_orientation::TerrestrialTimeScaleConverter timeScaleConverter = earth_orientation::TerrestrialTimeScaleConverter();

  // Check if there is one LinkEnds per observation time
  if (linkEndsVector_.size() != observationTimesUtc.size()) {
    throw std::runtime_error("Error while processing tracking data: vector of linkEnds and observationTimes not of equal size");
  }

  // Ge a vector of ground station positions
  std::vector<Eigen::Vector3d> groundStationPositions;
  for (const auto& linkEnds : linkEndsVector_) {
    std::string currentGroundStation = linkEnds.at(receiver).getStationName(); // TODO: what if transmitter and receiver different?
    groundStationPositions.push_back(earthFixedGroundStationPositions_.at(currentGroundStation));
  }

  // Convert to TDB using the GS positions
  std::vector<double> observationTimesTdb = timeScaleConverter.getCurrentTimes(basic_astrodynamics::utc_scale,
                                                                               basic_astrodynamics::tdb_scale,
                                                                               observationTimesUtc,
                                                                               groundStationPositions);
  return observationTimesTdb;
}

void ProcessedTrackingTxtFileContents::updateLinkEnds()
{
  // Clear any previous values
  linkEndsVector_.clear();

  // Get information from raw data file
  const auto& metaDataStrMap = rawTrackingTxtFileContents_->getMetaDataStrMap();
  const auto& numDataRows = rawTrackingTxtFileContents_->getNumRows();

  // Deduce linkends representation
  LinkEndsRepresentation linkEndsRepresentation = getLinkEndsRepresentation();

  // Create a vector of LinkEnds based on how they are represented
  // This currently only implements the DSN transmitter and receiver
  switch (linkEndsRepresentation) {

    // TODO: make a cleaner implementation to allow adding different ways of providing the link ends easily
    case dsn_transmitting_receiving_station_nr: {
      const auto& dsnTransmitterIds = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::dsn_transmitting_station_nr);
      const auto& dsnReceiverIds = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::dsn_receiving_station_nr);

      for (size_t i = 0; i < numDataRows; ++i) {
        std::string transmitterName = getStationNameFromStationId(dsnTransmitterIds.at(i));
        std::string receiverName = getStationNameFromStationId(dsnReceiverIds.at(i));
        LinkEnds currentLinkEnds{
            {transmitter, LinkEndId("Earth", transmitterName)},
            {reflector, LinkEndId(spacecraftName_, "")},
            {receiver, LinkEndId("Earth", receiverName)},
        };
        linkEndsVector_.push_back(currentLinkEnds);
      }
      break;
    }

    case vlbi_station: {

      // FIXME: This is a temporary solution
      std::string vlbiStationName = rawTrackingTxtFileContents_->getMetaDataStrMap().at(input_output::TrackingDataType::vlbi_station_name);
//      std::string vlbiStationName = rawTrackingTxtFileContents_->getMetaDataStrMap().at(input_output::TrackingDataType::vlbi_station_name);

      for (size_t i = 0; i < numDataRows; ++i) {
        LinkEnds currentLinkEnds{
            {transmitter, LinkEndId("Earth", "NNO")},
            {reflector, LinkEndId(spacecraftName_, "")},
            {receiver, LinkEndId("Earth", vlbiStationName)},
        };
        linkEndsVector_.push_back(currentLinkEnds);
      }
      break;
    }

      // Throw error if representation not implemented
    default: {
      throw std::runtime_error("Error while processing tracking txt file: LinkEnds representation not recognised or implemented.");
    }
  }

  // Creating a set with all the distinct LinkEnds
  linkEndsSet_ = utilities::vectorToSet(linkEndsVector_);
}

ProcessedTrackingTxtFileContents::TimeRepresentation ProcessedTrackingTxtFileContents::getTimeRepresentation()
{

  // Get all the data types from the raw file contents
  auto const& availableDataTypes = rawTrackingTxtFileContents_->getDataColumnTypes();

  // Return representation based on available data types

  if (utilities::containsAll(availableDataTypes, std::vector<input_output::TrackingDataType>{input_output::TrackingDataType::tdb_time_j2000})) {
    return tdb_seconds_j2000;
  }

  if (utilities::containsAll(availableDataTypes,
                             std::vector<input_output::TrackingDataType>{
                                 input_output::TrackingDataType::year,
                                 input_output::TrackingDataType::month,
                                 input_output::TrackingDataType::day,
                                 input_output::TrackingDataType::hour,
                                 input_output::TrackingDataType::minute,
                                 input_output::TrackingDataType::second
                             })) {
    return calendar_day_time;
  }

  // Throw an error if no match is found
  throw std::runtime_error("Error while processing tracking txt file: Time representation not recognised or implemented.");
}

ProcessedTrackingTxtFileContents::LinkEndsRepresentation ProcessedTrackingTxtFileContents::getLinkEndsRepresentation()
{
  // Get all the available data columns
  auto const& availableDataTypes = rawTrackingTxtFileContents_->getAllAvailableDataTypes();

  // Porvide link Ends representation based on available columns

  if (utilities::containsAll(availableDataTypes,
                             std::vector<input_output::TrackingDataType>{
                                 input_output::TrackingDataType::dsn_transmitting_station_nr,
                                 input_output::TrackingDataType::dsn_receiving_station_nr
                             })) {
    return dsn_transmitting_receiving_station_nr;
  }

  if (utilities::containsAll(availableDataTypes, std::vector<input_output::TrackingDataType>{input_output::TrackingDataType::vlbi_station_name})) {
    return vlbi_station;
  }

  // Throw error if no match is found
  throw std::runtime_error("Error while processing tracking txt file: Link Ends representation not recognised or implemented.");
}

} // namespace observation_models
} // namespace tudat
