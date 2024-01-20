/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_READ_GENERIC_TXT_FILE_H
#define TUDAT_READ_GENERIC_TXT_FILE_H

#include <boost/algorithm/string.hpp>
#include <boost/any.hpp>
#include <string>
#include <cstdarg>
#include <fstream>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "tudat/io/fieldType.h"
#include "tudat/astro/basic_astro.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/interface/spice/spiceInterface.h"

// A TrackingDataType is a unique form of data that can be used by Tudat to derive the observables
// A TrackingFileField is an identifier that specifies a specific type of column data in an input file (with a specific format)
// A TrackingFileFieldConverter is an interface that can convert a raw string input field to its associated double value as expected for the TrackingDataType
// For example:
// - You are reading a file with a column that shows the round trip light time in microseconds
// - The corresponding TrackingFileField will be "round_trip_light_time_microseconds
// - Tudat wants to refer to the data by TrackingDataType::n_way_light_time in seconds
// - The converter

namespace tudat
{
namespace input_output
{



//! Enum describing a unique data type that can later be used to process the information. These are in SI units
enum class TrackingDataType
{
  year,
  month,
  day,
  hour,
  minute,
  second,
  observation_time_scale,
  file_name,
  n_way_light_time,
  light_time_measurement_delay,
  light_time_measurement_accuracy,
  dsn_transmitting_station_nr,
  dsn_receiving_station_nr,
  observation_body, // In case observations corrected for body center.
  observed_body, // In case observations corrected for body center.
  spacecraft_id,
  planet_nr,
  tdb_time_j2000,
  tdb_spacecraft_j2000,
  x_planet_frame,
  y_planet_frame,
  z_planet_frame,
  vx_planet_frame,
  vy_planet_frame,
  vz_planet_frame,
  residual_de405,
  spacecraft_transponder_delay,
  uplink_frequency,
  downlink_frequency,
  signal_to_noise,
  spectral_max,
  doppler_measured_frequency,
  doppler_base_frequency,
  doppler_noise,
  doppler_bandwidth,
  vlbi_station_name,
};

//!Simple converter class that can convert a string data field to a double. One can inherit from this and overload the `toDouble()` method to extend the supported formats
class TrackingFileFieldConverter
{
public:
  explicit TrackingFileFieldConverter(TrackingDataType trackingDataType)
      : doubleDataType_(trackingDataType) {}

  virtual ~TrackingFileFieldConverter() = default;

  virtual double toDouble(std::string& rawField) const
  {
    try {
      return std::stod(rawField);
    } catch (std::invalid_argument&) {
      throw std::runtime_error(
          "The tracking file field cannot be converted correctly. Check your columnTypes. Raw field was \"" + rawField
              + "\".\n");
    }
  }
  const TrackingDataType& getTrackingDataType() { return doubleDataType_; }

private:
  TrackingDataType doubleDataType_;
};

//! A converter specifically for month fields in three-letter representation ("jan", "Jan", ...)
class TrackingFileMonthFieldConverter : public TrackingFileFieldConverter
{
public:
  TrackingFileMonthFieldConverter(TrackingDataType trackingDataType) : TrackingFileFieldConverter(trackingDataType) {}
  double toDouble(std::string& rawField) const
  {
    std::map<std::string, double> monthsMap{{"JAN", 1.},
                                            {"FEB", 2.},
                                            {"MAR", 3.},
                                            {"APR", 4.},
                                            {"MAY", 5.},
                                            {"JUN", 6.},
                                            {"JUL", 7.},
                                            {"AUG", 8.},
                                            {"SEP", 9.},
                                            {"OCT", 10.},
                                            {"NOV", 11.},
                                            {"DEC", 12.}};

    return utilities::upperCaseFromMap(rawField, monthsMap);
  }
};

//! Converter that will convert the raw string to double and then apply a scalar multiplier as specified
class TrackingFileFieldMultiplyingConverter : public TrackingFileFieldConverter
{
public:
  TrackingFileFieldMultiplyingConverter(TrackingDataType trackingDataType, double multiplier)
      : TrackingFileFieldConverter(trackingDataType), multiplier_(multiplier) {}
  double toDouble(std::string& rawField) const
  {
    return multiplier_ * TrackingFileFieldConverter::toDouble(rawField);
  }
private:
  double multiplier_;
};

//! Converter that will convert the raw string to double and then apply a scalar multiplier as specified
class TrackingFileFieldUTCTimeConverter : public TrackingFileFieldConverter
{
public:
  TrackingFileFieldUTCTimeConverter(TrackingDataType trackingDataType = TrackingDataType::tdb_time_j2000)
      : TrackingFileFieldConverter(trackingDataType) {}
  double toDouble(std::string& rawField) const
  {
    return spice_interface::convertDateStringToEphemerisTime(rawField);
  }
};

//! Mapping the `TrackingFileField` to the correct converter, including the `TrackingDataType` it will represent
static const std::map<std::string, std::shared_ptr<TrackingFileFieldConverter>> trackingFileFieldConverterMap = {
    {"spacecraft_id", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::spacecraft_id)},
    {"dsn_transmitting_station_nr", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::dsn_transmitting_station_nr)},
    {"dsn_receiving_station_nr", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::dsn_receiving_station_nr)},
    {"year", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::year)},
    {"month", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::month)},
    {"month_three_letter", std::make_shared<TrackingFileMonthFieldConverter>(TrackingDataType::month)},
    {"day", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::day)},
    {"hour", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::hour)},
    {"minute", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::minute)},
    {"second", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::second)},
    {"round_trip_light_time_seconds", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::n_way_light_time)},
    {"round_trip_light_time_microseconds", std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::n_way_light_time, 1.e-6)},
    {
        "light_time_measurement_delay_microseconds",
        std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::light_time_measurement_delay, 1.e-6)
    },
    {
        "light_time_measurement_accuracy_microseconds",
        std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::light_time_measurement_accuracy, 1.e-6)
    },
    {"planet_nr", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::planet_nr)},
    {"tdb_spacecraft_seconds_j2000", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::tdb_spacecraft_j2000)},
    {"x_planet_frame_km", std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::x_planet_frame, 1.e3)},
    {"y_planet_frame_km", std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::y_planet_frame, 1.e3)},
    {"z_planet_frame_km", std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::z_planet_frame, 1.e3)},
    {"vx_planet_frame_kms", std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::vx_planet_frame, 1.e3)},
    {"vy_planet_frame_kms", std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::vy_planet_frame, 1.e3)},
    {"vz_planet_frame_kms", std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::vz_planet_frame, 1.e3)},
    {"residual_de405_microseconds", std::make_shared<TrackingFileFieldMultiplyingConverter>(TrackingDataType::residual_de405, 1.e-6)},
    {"signal_to_noise_ratio", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::signal_to_noise)},
    {"normalised_spectral_max", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::spectral_max)},
    {"doppler_measured_frequency_hz", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::doppler_measured_frequency)},
    {"doppler_noise_hz", std::make_shared<TrackingFileFieldConverter>(TrackingDataType::doppler_noise)},
    {"utc_datetime_string", std::make_shared<TrackingFileFieldUTCTimeConverter>(TrackingDataType::tdb_time_j2000)}
};

//! Class to extract the raw data from a file with the appropriate conversion to doubles
class TrackingTxtFileContents
{
public:
  TrackingTxtFileContents(std::string fileName,
                          std::vector<std::string> columnTypes,
                          char commentSymbol = '#',
                          std::string valueSeparators = ",: \t")
      : fileName_(std::move(fileName)), columnFieldTypes_(std::move(columnTypes)), commentSymbol_(commentSymbol),
        valueSeparators_(std::move(valueSeparators))
  {
    parseData();
  }

  void parseData();

  void readRawDataMap(std::ifstream& dataFile);

  void addLineToRawDataMap(std::string& rawLine);

  void convertDataMap();

  void addMetaData(TrackingDataType dataType, double value) { metaDataMapDouble_[dataType] = value; }

  void addMetaData(TrackingDataType dataType, const std::string& value) { metaDataMapStr_[dataType] = value; }


// Getters
public:
  size_t getNumColumns() const { return columnFieldTypes_.size(); }
  size_t getNumRows() const { return rawDataMap_.at(columnFieldTypes_[0]).size(); }
  const std::vector<std::string>& getRawColumnTypes() { return columnFieldTypes_; }

  const std::vector<TrackingDataType>& getDataColumnTypes()
  {
    columnDataTypes_.clear();
    for (auto& pair : doubleDataMap_) {
      columnDataTypes_.push_back(pair.first);
    }
    return columnDataTypes_;
  }

  const auto& getRawDataMap() { return rawDataMap_; }
  const auto& getDoubleDataMap() { return doubleDataMap_; }
  const std::vector<double>& getDoubleDataColumn(TrackingDataType dataType);
  const auto& getMetaDataDoubleMap() { return metaDataMapDouble_; }
  const auto& getMetaDataStrMap() { return metaDataMapStr_; }
  const std::vector<TrackingDataType> getMetaDataTypes();
  const std::vector<TrackingDataType> getAllAvailableDataTypes();

private:
  std::string fileName_ = "None";
  std::string separators_ = ":, \t";
  std::vector<std::string> columnFieldTypes_;
  std::vector<TrackingDataType> columnDataTypes_;
  char commentSymbol_;
  std::string valueSeparators_;
  std::map<std::string, std::vector<std::string>> rawDataMap_;
  std::map<TrackingDataType, double> metaDataMapDouble_;
  std::map<TrackingDataType, std::string> metaDataMapStr_;
  std::map<TrackingDataType, std::vector<int>> intDataMap_;
  std::map<TrackingDataType, std::vector<double>> doubleDataMap_;

// Utility variables
private:
  std::vector<std::string> currentSplitRawLine_;
};

/*!
 * Function to create a read out a tracking data file to raw contents
 * @param fileName
 * @param columnTypes column types (string). If known to Tudat, this will define a tracking data type, otherwise, it is not processed and kept as raw data.
 * @param commentSymbol lines that start with this symbol are ignored
 * @param valueSeparators String of characters that separate columns. E.g. ",:" means that every , and : in the file will create a new column
 * @return TrackingFileContents
 */
static inline std::unique_ptr<TrackingTxtFileContents> createTrackingTxtFileContents(const std::string& fileName,
                                                                                     std::vector<std::string>& columnTypes,
                                                                                     char commentSymbol = '#',
                                                                                     const std::string& valueSeparators = ",: \t")
{
  return std::make_unique<TrackingTxtFileContents>(fileName, columnTypes, commentSymbol, valueSeparators);
}

} // namespace input_output
} // namespace tudat

#endif // TUDAT_READ_GENERIC_TXT_FILE_H
