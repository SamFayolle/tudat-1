/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
//
//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"
#include "tudat/simulation/estimation_setup.h"

#include "tudat/io/readOdfFile.h"
#include "tudat/io/readTabulatedMediaCorrections.h"
#include "tudat/io/readTabulatedWeatherData.h"
#include "tudat/simulation/estimation_setup/processOdfFile.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include <boost/date_time/gregorian/gregorian.hpp>

#include "tudat/astro/ground_stations/transmittingFrequencies.h"
//
//namespace tudat
//{
//namespace unit_tests
//{

using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::input_output;
using namespace tudat::simulation_setup;
using namespace tudat::unit_tests;
using namespace tudat;


//BOOST_AUTO_TEST_SUITE( test_dsn_n_way_averaged_doppler_observation_model )
//
//// Comparison of computed observables with Verma (2022), "A Python-based tool for constructing observables from the DSN’s
//// closed-loop archival tracking data files", fig. 4 right
//// Currently the computed observables are slightly different from Verma (2022), so the testing is actually done using
//// values generated by this code.
//BOOST_AUTO_TEST_CASE( testDsnNWayAveragedDopplerModel )
//{
//
//    spice_interface::loadStandardSpiceKernels( );
//    // Verma (2022) uses DE438, but here we use the standard Tudat SPICE kernels as the difference produced by the kernels
//    // is way below the level of the current residuals.
//    // spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/planet-spice/de438.bsp" );
//    spice_interface::loadSpiceKernelInTudat( tudat::paths::getTudatTestDataPath( ) + "dsn_n_way_doppler_observation_model/mgs_map1_ipng_mgs95j.bsp" );
//
//    std::string spacecraftName = "MGS";
//    std::string ephemeridesOrigin = "SSB";
//
//    std::vector< std::string > odfFiles = {
//            tudat::paths::getTudatTestDataPath( ) + "dsn_n_way_doppler_observation_model/9068068a.odf",
//            tudat::paths::getTudatTestDataPath( ) + "dsn_n_way_doppler_observation_model/9068071a.odf" };
//
//    // Define bodies to use.
//    std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Mars" };
//
//    // Create bodies settings needed in simulation
//    BodyListSettings bodySettings;
//    bodySettings = getDefaultBodySettings( bodiesToCreate, ephemeridesOrigin, "J2000" );
//
//
//    bodySettings.at( "Earth" )->shapeModelSettings = fromSpiceOblateSphericalBodyShapeSettings( );
//    bodySettings.at( "Earth" )->rotationModelSettings = gcrsToItrsRotationModelSettings(
//            basic_astrodynamics::iau_2006, "J2000" );
//    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );
//
//    // Create spacecraft
//    bodySettings.addSettings( spacecraftName );
//    bodySettings.at( spacecraftName )->ephemerisSettings =
//            std::make_shared< DirectSpiceEphemerisSettings >( ephemeridesOrigin, "J2000" );
//
//    // Create bodies
//    SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );
//
//    // Read and process ODF file data
//    std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector;
//    for ( std::string odfFile : odfFiles )
//        rawOdfDataVector.push_back( std::make_shared< OdfRawFileContents >( odfFile ) );
//
//    std::shared_ptr< ProcessedOdfFileContents > processedOdfFileContents =
//            std::make_shared< ProcessedOdfFileContents >( rawOdfDataVector, spacecraftName );
//
//    // Create observed observation collection
//    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > observedObservationCollection =
//            observation_models::createOdfObservedObservationCollection< long double, Time >(
//                    processedOdfFileContents, { dsn_n_way_averaged_doppler } );
//
//    // Create computed observation collection
//    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;
//
//    std::map < observation_models::ObservableType, std::vector< observation_models::LinkEnds > > linkEndsPerObservable =
//            observedObservationCollection->getLinkEndsPerObservableType( );
//    for ( auto it = linkEndsPerObservable.begin(); it != linkEndsPerObservable.end(); ++it )
//    {
//        for ( unsigned int i = 0; i < it->second.size(); ++i )
//        {
//            if ( it->first == observation_models::dsn_n_way_averaged_doppler )
//            {
//                observationModelSettingsList.push_back(
//                    std::make_shared< observation_models::DsnNWayAveragedDopplerObservationSettings >( it->second.at( i ) ) );
//            }
//        }
//    }
//
//    observation_models::setOdfInformationInBodies( processedOdfFileContents, bodies );
//
//    std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< long double, Time > > >
//            observationSimulators = observation_models::createObservationSimulators< long double, Time >(
//                    observationModelSettingsList, bodies );
//
//    std::vector< std::shared_ptr< ObservationSimulationSettings< Time > > > observationSimulationSettings =
//            getObservationSimulationSettingsFromObservations< long double, Time >(
//                    observedObservationCollection );
//
//    std::cout<<"Pre-simulation "<<observationSimulationSettings.size( )<<std::endl;
//
//    std::shared_ptr< observation_models::ObservationCollection< long double, Time > >
//            simulatedObservationCollection = simulation_setup::simulateObservations< long double, Time >(
//                    observationSimulationSettings, observationSimulators, bodies );
//    std::cout<<"Post-simulation"<<std::endl;
//
//    LinkEnds dss45MgsLinkEnds;
//    dss45MgsLinkEnds[ transmitter ] = LinkEndId( "Earth", "DSS-45" );
//    dss45MgsLinkEnds[ retransmitter ] = LinkEndId( "MGS" );
//    dss45MgsLinkEnds[ receiver ] = LinkEndId( "Earth", "DSS-45" );
//
//     std::vector< std::shared_ptr< SingleObservationSet< long double, Time > > > singleLinkSimulatedObservations =
//             simulatedObservationCollection->getSingleLinkAndTypeObservationSets( dsn_n_way_averaged_doppler, dss45MgsLinkEnds );
//     std::vector< std::shared_ptr< SingleObservationSet< long double, Time > > > singleLinkObservedObservations =
//             observedObservationCollection->getSingleLinkAndTypeObservationSets( dsn_n_way_averaged_doppler, dss45MgsLinkEnds );
//
//     // Compare value of observations
//     // Min/max values from Verma (2022), "A Python-based tool for constructing observables from the DSN’s closed-loop
//     // archival tracking data files", fig. 4 right
//     std::pair< double, double > vermaFirstBlockMinMaxResidual = std::make_pair( -0.04, 0.10 );
//     std::pair< double, double > vermaSecondBlockMinMaxResidual = std::make_pair( -0.06, 0.07 );
//     // Min/max values obtained with current code
//     std::pair< double, double > firstBlockMinMaxResidual = std::make_pair( -0.025, 0.165 );
//     std::pair< double, double > secondBlockMinMaxResidual = std::make_pair( -0.03, 0.12 );
//
//     // Currently the residuals are compared with firstBlockMinMaxResidual and secondBlockMinMaxResidual, because
//     // the residuals are slightly different from what Verma (2022) presents.
//     // Potential reasons for the differences:
//     // - MGS SPICE kernels: unlikely, trying all the available kernels from NAIF doesn't lead to residuals more similar
//     //         to Verma (2022)
//     // - Earth rotation model: using IAU 2000a/2000b/2006 doesn't make any significant difference in the residuals. Using
//     //         the rotation model from SPICE leads to much larger residuals
//     // - Ground stations position: unlikely, considering that Verma works at NASA he probably used the positions according
//     //         to DSN No. 810-005, 301 (same as used here)
//     // - Some error in the time conversion
//     // - Some error in the conversion between frames
//     // - Difference in the used count interval: Verma uses 60s, the ODF file uses 10s
//
//
//     for ( unsigned int i = 0; i < singleLinkSimulatedObservations.size( ); ++i )
//     {
//         std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings =
//                 singleLinkSimulatedObservations.at( i )->getAncilliarySettings( );
//
//         double referenceFrequency = ancillarySettings->getAncilliaryDoubleData( doppler_reference_frequency );
//
//         // Ignoring observations with reference frequencies of approximately 2.115e9.
//         // The observations with f_ref ~ 2.114e9 are reproduced accurately, while the ones with f_ref ~ 2.115e9 have very
//         // large residuals (the residuals are basically a bias).
//         // Comparing contiguous data points with f_ref ~ 2.114e9 and f_ref ~ 2.115e9, which only differ
//         // in the value of the f_ref, the simulated observations for f_ref ~ 2.115e9 appear to be consistent with the
//         // value of f_ref, while the observed observations do not. Considering the that the data points are identical
//         // besides f_ref, it is suggested that the values of the reference frequency f_ref ~ 2.115e9 are wrongly reported
//         // in the ODF file, hence the associated data points are here ignored.
//         if ( referenceFrequency > 2.115e9 )
//         {
//             continue;
//         }
//
//         std::vector< Eigen::Matrix< long double, Eigen::Dynamic, 1 > > simulatedObservations =
//                 singleLinkSimulatedObservations.at( i )->getObservations( );
//         std::vector< Time > simulatedTimes = singleLinkSimulatedObservations.at( i )->getObservationTimes( );
//         std::vector< Eigen::Matrix< long double, Eigen::Dynamic, 1 > > observedObservations =
//                 singleLinkObservedObservations.at( i )->getObservations( );
//         std::vector< Time > observedTimes = singleLinkObservedObservations.at( i )->getObservationTimes( );
//
//         for ( unsigned int j = 0; j < simulatedTimes.size( ); ++j )
//         {
//             // Ignore a few observations potentially associated with invalid cycle counts. These are observations which
//             // have significantly higher residuals than the previous/subsequent observations, and that in the ODF file are
//             // presented directly before/after invalid observations; the value of the (wrong?) observed observation in the ODF
//             // file might have been computed using some invalid cycle counts associated with the previous/subsequent
//             // invalid observations.
//
//             // t=-25540040.8144906
//             // Simulated observation: 16505.3076701937 Hz
//             // Observed observation: 92427.4388 Hz
//
//             // t=-25706110.8145144
//             // Simulated observation: 182303.286088224 Hz
//             // Observed observation: 182305.413800001 Hz
//
//             // t=-25553670.8144926
//             // Simulated observation: 9431.46921279747 Hz
//             // Observed observation: 9416.685600001 Hz
//
//             if ( simulatedTimes.at( j ).getSeconds< int >( ) == -25540040 ||
//                 simulatedTimes.at( j ).getSeconds< int >( ) == -25706110 ||
//                 simulatedTimes.at( j ).getSeconds< int >( ) == -25553670 )
//             {
//                 continue;
//             }
//
//             // Ignore additional point: single point in the middle of an observation arc with significantly higher
//             // residuals. Contrary to the previous points, it isn't associated with any invalid observation, so no idea
//             // where the larger residuals are coming from.
//             if ( simulatedTimes.at( j ).getSeconds< int >( ) == -25543670 )
//             {
//                 continue;
//             }
//
//             double residual = observedObservations.at( j )( 0 ) - simulatedObservations.at( j )( 0 );
//
//             BOOST_CHECK_EQUAL( simulatedTimes.at( j ), observedTimes.at( j ) );
//
//             if ( simulatedTimes.at( j ) < -2.56e7  )
//             {
//                 BOOST_TEST( residual >= firstBlockMinMaxResidual.first );
//                 BOOST_TEST( residual <= firstBlockMinMaxResidual.second );
//             }
//             else
//             {
//                 BOOST_TEST( residual >= secondBlockMinMaxResidual.first );
//                 BOOST_TEST( residual <= secondBlockMinMaxResidual.second );
//             }
//         }
//     }
//}


//BOOST_AUTO_TEST_CASE( testDsnNWayDataCompression )
int main( )
{
    std::map< int, std::string > spacecraftNamePerSpacecraftId;
    spacecraftNamePerSpacecraftId[ 177 ] = "GRAIL-A";
    spacecraftNamePerSpacecraftId[ 181 ] = "GRAIL-B";

    double initialTimeEnvironment = DateTime( 2012, 03, 02, 0, 0, 0.0 ).epoch< double >();
    double finalTimeEnvironment = DateTime( 2012, 05, 29, 0, 0, 0.0 ).epoch< double >();

    // Load spice kernels
    spice_interface::loadStandardSpiceKernels( );
    spice_interface::loadSpiceKernelInTudat( get_spice_kernels_path( ) + + "/moon_de440_200625.tf");
    spice_interface::loadSpiceKernelInTudat( "/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_kernels/grail_v07.tf" );
    spice_interface::loadSpiceKernelInTudat( "/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_kernels/gra_sclkscet_00013.tsc" );
    spice_interface::loadSpiceKernelInTudat( "/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_kernels/gra_sclkscet_00014.tsc" );
    spice_interface::loadSpiceKernelInTudat( "/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_kernels/grail_120301_120529_sci_v02.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_kernels/gra_rec_120402_120408.bc" );
    spice_interface::loadSpiceKernelInTudat( get_spice_kernels_path( ) + + "/moon_pa_de440_200625.bpc");
    // Create settings for default bodies
    std::vector< std::string > bodiesToCreate = { "Earth", "Moon", "Sun" };
    std::string globalFrameOrigin = "Earth";
    std::string globalFrameOrientation = "J2000";
    BodyListSettings bodySettings = getDefaultBodySettings(
            bodiesToCreate, initialTimeEnvironment, finalTimeEnvironment, globalFrameOrigin, globalFrameOrientation, 120.0 );

    // Add high-accuracy Earth settings
    bodySettings.at( "Earth" )->shapeModelSettings = fromSpiceOblateSphericalBodyShapeSettings( );
    bodySettings.at( "Earth" )->rotationModelSettings = gcrsToItrsRotationModelSettings(
            basic_astrodynamics::iau_2006, globalFrameOrientation,
            std::make_shared< interpolators::InterpolatorGenerationSettings< double > >(
                    interpolators::cubicSplineInterpolation( ), initialTimeEnvironment, finalTimeEnvironment, 3600.0 ),
            std::make_shared< interpolators::InterpolatorGenerationSettings< double > >(
                    interpolators::cubicSplineInterpolation( ), initialTimeEnvironment, finalTimeEnvironment, 3600.0 ),
            std::make_shared< interpolators::InterpolatorGenerationSettings< double > >(
                    interpolators::cubicSplineInterpolation( ), initialTimeEnvironment, finalTimeEnvironment, 60.0 ));
    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );
    bodySettings.at( "Moon" )->ephemerisSettings->resetFrameOrigin( "Earth" );

    // Add spacecraft settings
    std::string spacecraftName = "GRAIL-A";
    std::string spacecraftCentralBody = "Moon";
    bodySettings.addSettings( spacecraftName );
    bodySettings.at( spacecraftName )->ephemerisSettings =
            std::make_shared< InterpolatedSpiceEphemerisSettings >(
                    initialTimeEnvironment, finalTimeEnvironment, 10.0, spacecraftCentralBody, globalFrameOrientation );
    bodySettings.at( spacecraftName )->rotationModelSettings = spiceRotationModelSettings(
            globalFrameOrientation, spacecraftName + "_SPACECRAFT", "" );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );
    bodies.at( "GRAIL-A" )->getVehicleSystems( )->setReferencePointPosition(
            "Antenna", ( Eigen::Vector3d( ) << -0.082, 0.152, -0.810 ).finished( ) );

    // Define ODF data paths
    std::string dataDirectory = "/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_data/";
    std::vector< std::string > odfFiles = { "gralugf2012_097_0235smmmv1.odf" };

    // Laod raw ODF data
    std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector;
    for ( std::string odfFile : odfFiles )
    {
        rawOdfDataVector.push_back( std::make_shared< OdfRawFileContents >( dataDirectory + odfFile ) );
    }

    // Process ODF file data
    std::shared_ptr< ProcessedOdfFileContents< Time > > processedOdfFileContents =
            std::make_shared< ProcessedOdfFileContents< Time > >( rawOdfDataVector, spacecraftName, true );
    observation_models::setOdfInformationInBodies( processedOdfFileContents, bodies );

    // Create data structure that handles Observed Data in Tudat
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > fullObservedObservationCollection =
            splitObservationSets( observation_models::createOdfObservedObservationCollection< long double, Time >(
                    processedOdfFileContents, { dsn_n_way_averaged_doppler } ), observationSetSplitter( time_interval_splitter, 60.0, 10 ) );

    std::shared_ptr< observation_models::SingleObservationSet< long double, Time > > singleDopplerSet =
            fullObservedObservationCollection->getSingleObservationSets( observationParser( dsn_n_way_averaged_doppler ) ).at( 0 );
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > observedObservationCollection =
            std::make_shared< observation_models::ObservationCollection< long double, Time > >(
                    std::vector< std::shared_ptr< SingleObservationSet< long double, Time > > > { singleDopplerSet } );

    std::map< LinkEnds, std::vector< std::shared_ptr< observation_models::SingleObservationSet< long double, Time > > > > observedObservationSets =
            observedObservationCollection->getObservationsSets( ).at( dsn_n_way_averaged_doppler );

    // Define light-time correction settings
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;
    std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back( firstOrderRelativisticLightTimeCorrectionSettings( { "Sun" } ) );
    std::vector< std::string > troposphericCorrectionFileNames =
            {"/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_data/grxlugf2012_092_2012_122.tro" };
    std::vector< std::string > ionosphericCorrectionFileNames =
            {"/home/mfayolle/Tudat/tudat-bundle/tudatpy/examples/estimation/grail_data/gralugf2012_092_2012_122.ion" };
    lightTimeCorrectionSettings.push_back( tabulatedTroposphericCorrectionSettings( troposphericCorrectionFileNames ) );
    lightTimeCorrectionSettings.push_back( tabulatedIonosphericCorrectionSettings( ionosphericCorrectionFileNames, spacecraftNamePerSpacecraftId ) );

    // Define observation model settings
    for ( auto it : observedObservationCollection->getLinkEndsPerObservableType( ) )
    {
        for ( unsigned int i = 0; i < it.second.size( ); ++i )
        {
            if ( it.first == observation_models::dsn_n_way_averaged_doppler )
            {
                observationModelSettingsList.push_back(
                        observation_models::dsnNWayAveragedDopplerObservationSettings(
                                it.second.at( i ), lightTimeCorrectionSettings, nullptr,
                                std::make_shared< LightTimeConvergenceCriteria >( true ) ) );
            }
        }
    }

    // Create observation simulators and compute residuals
    std::vector< std::shared_ptr< ObservationSimulatorBase< long double, Time > > > observationSimulators =
            createObservationSimulators< long double, Time >( observationModelSettingsList, bodies );

    // Create observation simulation settings consistent with observations
    std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< Time > > > observationSimulationSettings =
            getObservationSimulationSettingsFromObservations( observedObservationCollection );

    // Compute residuals
    computeAndSetResiduals< long double, Time >( observedObservationCollection, observationSimulators, bodies );
    Eigen::Matrix< long double, Eigen::Dynamic, 1 > residuals = observedObservationCollection->getConcatenatedResiduals( );

    Eigen::VectorXd startTimes;
    Eigen::VectorXd durations;
    Eigen::VectorXd meanValues;
    Eigen::VectorXd rmsValues;

    getResidualStatistics(
            observedObservationCollection,
            startTimes,
            durations,
            meanValues,
            rmsValues );

    std::cout<<"RMS "<<rmsValues.transpose( )<<std::endl;
    std::cout<<"Mean "<<meanValues.transpose( )<<std::endl<<std::endl;

    input_output::writeMatrixToFile( residuals, "grailTestResidualsUncompressed.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_TimeBiasTest/");

    Eigen::VectorXd observationTimes = utilities::convertStlVectorToEigenVector(
            observedObservationCollection->getConcatenatedTimeVector( ) ).template cast< double >( );
    input_output::writeMatrixToFile( observationTimes, "grailTestTimeUncompresseds.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_TimeBiasTest/");


    Eigen::VectorXd numericalTimeBiasPartials = getNumericalObservationTimePartial< long double, Time >(
            observationSimulationSettings, observationSimulators, bodies, 5.0 );
    Eigen::VectorXd correctedResiduals;
    std::vector< double > timeBiases;
    estimateTimeBiasPerSet( observedObservationCollection, numericalTimeBiasPartials, timeBiases, correctedResiduals );

    input_output::writeMatrixToFile( correctedResiduals, "grailTestCorrectedResidualsUncompressed.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_TimeBiasTest/");

    std::cout<<timeBiases.at( 0 )<<std::endl;


    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > compressedObservedObservationCollection =
            createCompressedDopplerCollection( observedObservationCollection, 60.0 );

    // Create observation simulation settings consistent with observations
    std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< Time > > > compressedObservationSimulationSettings =
            getObservationSimulationSettingsFromObservations( compressedObservedObservationCollection );

    // Compute compressed residuals
    computeAndSetResiduals< long double, Time >( compressedObservedObservationCollection, observationSimulators, bodies );
    Eigen::Matrix< long double, Eigen::Dynamic, 1 > compressedResiduals = observedObservationCollection->getConcatenatedResiduals( );

    getResidualStatistics(
            compressedObservedObservationCollection,
            startTimes,
            durations,
            meanValues,
            rmsValues );

    std::cout<<"RMS "<<rmsValues.transpose( )<<std::endl;
    std::cout<<"Mean "<<meanValues.transpose( )<<std::endl<<std::endl;

    input_output::writeMatrixToFile( startTimes, "grailTestStartTimes.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_TimeBiasTest/");
    input_output::writeMatrixToFile( durations, "grailTestDurations.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_TimeBiasTest/");
    input_output::writeMatrixToFile( meanValues, "grailMeanValues.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_TimeBiasTest/");
    input_output::writeMatrixToFile( rmsValues, "grailRmsValues.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_TimeBiasTest/");

    input_output::writeMatrixToFile( compressedResiduals, "grailTestResiduals.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_TimeBiasTest/");

    observationTimes = utilities::convertStlVectorToEigenVector(
            compressedObservedObservationCollection->getConcatenatedTimeVector( ) ).template cast< double >( );
    input_output::writeMatrixToFile( observationTimes, "grailTestTimes.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_TimeBiasTest/");

    Eigen::VectorXd observationLinkEndsIds = utilities::convertStlVectorToEigenVector(
            compressedObservedObservationCollection->getConcatenatedLinkEndIds( ) ).template cast< double >( );
    input_output::writeMatrixToFile(observationLinkEndsIds , "grailTestLinkEnds.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_TimeBiasTest/");


    Eigen::VectorXd compressedNumericalTimeBiasPartials = getNumericalObservationTimePartial< long double, Time >(
            compressedObservationSimulationSettings, observationSimulators, bodies, 5.0 );
    Eigen::VectorXd correctedCompressedResiduals;
    std::vector< double > compressedtimeBiases;
    estimateTimeBiasPerSet(
            compressedObservedObservationCollection, compressedNumericalTimeBiasPartials,
            compressedtimeBiases, correctedCompressedResiduals );

    input_output::writeMatrixToFile( correctedCompressedResiduals, "grailTestCorrectedResiduals.dat", 16, "/home/mfayolle/Tudat/Data/GRAIL_TimeBiasTest/");

    std::cout<<compressedtimeBiases.at( 0 )<<std::endl;



}
//
//BOOST_AUTO_TEST_SUITE_END( )
//
//}
//
//}