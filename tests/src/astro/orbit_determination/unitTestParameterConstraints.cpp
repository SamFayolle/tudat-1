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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <string>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/simulation.h"
#include "tudat/simulation/estimation.h"

namespace tudat {

namespace unit_tests {

//! Using declarations.
using namespace numerical_integrators;
using namespace spice_interface;
using namespace ephemerides;
using namespace simulation_setup;
using namespace basic_astrodynamics;
using namespace orbital_element_conversions;
using namespace propagators;
using namespace estimatable_parameters;
using namespace unit_conversions;

void addGravityFieldVariationsMoons(
        BodyListSettings& bodySettings,
        const std::complex< double > k2Moon,
        const std::complex< double > k2Io )
{
    // (Complex) Love numbers
    std::map< int, std::vector< std::complex< double > > > loveNumbersIo, loveNumbersOtherMoons;

    std::vector< std::complex< double > > loveNumbersIoList = { k2Io, k2Io, k2Io };
    loveNumbersIo[ 2 ] = loveNumbersIoList;

    std::vector< std::complex< double > > loveNumbersOtherMoonsList = { k2Moon, k2Moon, k2Moon };
    loveNumbersOtherMoons[ 2 ] = loveNumbersOtherMoonsList;

    // Gravity field variations
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationsIo, gravityFieldVariationsOtherMoons;

    std::vector< std::string > deformingBody = { "Jupiter" };
    // For Io
    gravityFieldVariationsIo.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >( deformingBody, loveNumbersIo ) );
    // For the other Galilean moons
    gravityFieldVariationsOtherMoons.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >( deformingBody, loveNumbersOtherMoons ) );

    // Define gravity field variations settings
    bodySettings.at( "Io" )->gravityFieldVariationSettings = gravityFieldVariationsIo;
    bodySettings.at( "Europa" )->gravityFieldVariationSettings = gravityFieldVariationsOtherMoons;
    bodySettings.at( "Ganymede" )->gravityFieldVariationSettings = gravityFieldVariationsOtherMoons;
    bodySettings.at( "Callisto" )->gravityFieldVariationSettings = gravityFieldVariationsOtherMoons;

}

double getInvQFromLoveNumber( const std::complex< double > k2 )
{
    return k2.imag( ) / std::sqrt( k2.real( ) * k2.real( ) + k2.imag( ) * k2.imag( ) );
}

double getTidalPeriodMoon( const std::string moon,
                           const double initialEpoch,
                           const SystemOfBodies& bodies )
{
    Eigen::Vector6d initialState = getInitialStatesOfBodies(
            std::vector< std::string >( { moon } ), std::vector< std::string >( { "Jupiter" } ), bodies, initialEpoch ).segment( 0, 6 );
    Eigen::Vector6d keplerianState = orbital_element_conversions::convertCartesianToKeplerianElements(
            initialState, bodies.at( "Jupiter" )->getGravitationalParameter( ) );
    double semiMajorAxis = keplerianState[ semiMajorAxisIndex ];
    double periodTides = 2.0 * mathematical_constants::PI * std::sqrt( std::pow( semiMajorAxis, 3 ) / bodies.at( "Jupiter" )->getGravitationalParameter( ) );

    return periodTides;
}

double getTimeLagMoonFromLoveNumber( const std::complex< double > k2,
                                     const double tidalPeriod )
{
    double invQ = getInvQFromLoveNumber( k2 );
    return tidalPeriod * std::atan( invQ ) / ( 2.0 * mathematical_constants::PI );
}

BOOST_AUTO_TEST_SUITE( test_parameter_constraints )

BOOST_AUTO_TEST_CASE( testConstraintsSetup )
{
    std::cout.precision( 20 );

    loadStandardSpiceKernels( );

    double initialEpoch = 32.0 * physical_constants::JULIAN_YEAR;
    double finalEpoch = initialEpoch + 20.0 * physical_constants::JULIAN_DAY;

    std::string globalFrameOrientation = "J2000";
    std::string globalFrameOrigin = "Jupiter";

    std::complex< double > k2Moon = std::complex< double >( 0.2, 0.01 );
    std::complex< double > k2Io = std::complex< double >( 0.2, 0.01 );
    double k2MoonNorm = std::sqrt( std::norm( k2Moon ) );
    double k2IoNorm = std::sqrt( std::norm( k2Io ) );

    double invQMoon = getInvQFromLoveNumber( k2Moon );
    double invQIo = getInvQFromLoveNumber( k2Io ) * 1.1;


    // Define bodies settings for simulation
    std::vector< std::string > bodiesToCreate = { "Io", "Europa", "Ganymede", "Callisto", "Jupiter" };

    // Create body objects.
    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, initialEpoch, finalEpoch, globalFrameOrigin, globalFrameOrientation );
    for ( unsigned int i = 0 ; i < bodiesToCreate.size( ) - 1 ; i++ )
    {
        bodySettings.at( bodiesToCreate.at( i ) )->rotationModelSettings = std::make_shared< SynchronousRotationModelSettings >(
                "Jupiter", globalFrameOrientation, "IAU_" + bodiesToCreate.at( i ) );
    }

    // Add temporal variations to gravity fields
    addGravityFieldVariationsMoons( bodySettings, k2Moon, k2Io );


    // Crate system of bodies
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "spacecraft" );

    // Define bodies to propagate
    std::vector< std::string > moonsToPropagate = { "Io", "Europa" };
    unsigned int nbMoons = moonsToPropagate.size( );

    std::vector< std::string > centralBodiesMoons;
    for ( unsigned int i = 0 ; i < nbMoons ; i++ )
    {
        centralBodiesMoons.push_back( "Jupiter" );
    }

    // Add spacecraft around Io
    std::vector< std::string > bodiesToPropagate = moonsToPropagate;
    std::vector< std::string > centralBodies = centralBodiesMoons;
    bodiesToPropagate.push_back( "spacecraft" );
    centralBodies.push_back( "Io" );


    // Create map to store tidal periods
    std::map< std::string, double > tidalPeriods;
    for ( auto moon : moonsToPropagate )
    {
        tidalPeriods[ moon ] = getTidalPeriodMoon( moon, initialEpoch, bodies );
        std::cout << "tidal period for " << moon << " = " << tidalPeriods.at( moon ) << "\n\n";
    }


    // Define accelerations

    // NO NEED FOR THIS (???)
    int maximumMoonDegree = 2;
    int maximumMoonOrder = 2;
    int maximumJupiterDegree = 10;
    int maximumJupiterOrder = 10;

    SelectedAccelerationMap accelerationSettings;


    // Accelerations moons
    bool timeVaryingGravityMoonAcc = false;
    bool timeVaryingGravityScAcc = true;
    for ( const auto moon : moonsToPropagate )
    {
        accelerationSettings[ moon ][ "Jupiter" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >(
                maximumJupiterDegree, maximumJupiterOrder, maximumMoonDegree, maximumMoonOrder, 0, 0, timeVaryingGravityMoonAcc, timeVaryingGravityMoonAcc, timeVaryingGravityMoonAcc ) );

        // Tides in the moon
        if ( moon == "Io" )
        {
//            double timeLagIo = getTimeLagMoonFromLoveNumber( k2Io, tidalPeriods.at( moon ) );
            accelerationSettings[ moon ][ "Jupiter" ].push_back( std::make_shared< DirectTidalDissipationAccelerationSettings >(
                    k2IoNorm, invQIo, tidalPeriods.at( moon ), true, false, false ) );
        }
        else
        {
            double timeLagMoon = getTimeLagMoonFromLoveNumber( k2Moon, tidalPeriods.at( moon ) );
//            accelerationSettings[ moon ][ "Jupiter" ].push_back( std::make_shared< DirectTidalDissipationAccelerationSettings >(
//                    k2MoonNorm, timeLagMoon, true, false, false ) );
            accelerationSettings[ moon ][ "Jupiter" ].push_back( std::make_shared< DirectTidalDissipationAccelerationSettings >(
                    k2MoonNorm, invQMoon, tidalPeriods.at( moon ), true, false, false ) );
        }
    }

    // Acceleration spacecraft
    accelerationSettings[ "spacecraft" ][ "Io" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >(
            maximumMoonDegree, maximumMoonOrder, false, timeVaryingGravityScAcc ) );


    basic_astrodynamics::AccelerationMap accelerations = createAccelerationModelsMap( bodies, accelerationSettings, bodiesToPropagate, centralBodies );


    // Define initial state
    Eigen::Vector6d initialKeplerianStateSc = Eigen::Vector6d::Zero( );
    initialKeplerianStateSc[ semiMajorAxisIndex ] = 500.0e3;
    initialKeplerianStateSc[ inclinationIndex ] = convertDegreesToRadians( 20.0 );
    Eigen::Vector6d initialStateSc = orbital_element_conversions::convertKeplerianToCartesianElements(
            initialKeplerianStateSc, bodies.at( "Io" )->getGravitationalParameter( ) );

    Eigen::VectorXd initialStatesMoons = getInitialStatesOfBodies( moonsToPropagate, centralBodiesMoons, bodies, initialEpoch );
    Eigen::VectorXd initialStates = Eigen::VectorXd::Zero( 6 * bodiesToPropagate.size( ) );
        initialStates.segment( 0, nbMoons * 6 ) = initialStatesMoons;
        initialStates.segment( nbMoons * 6, 6 ) = initialStateSc;

    // Define integrator settings
    double timeStep = 50.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
            initialEpoch, timeStep, CoefficientSets::rungeKuttaFehlberg78, timeStep, timeStep, 1.0e3, 1.0e3 );

    // Define propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
            centralBodies, accelerations, bodiesToPropagate, initialStates, initialEpoch, integratorSettings, std::make_shared< PropagationTimeTerminationSettings >( finalEpoch ) );


    unsigned int nbTestCases = 5;
    for ( unsigned int testCase = 0 ; testCase < nbTestCases ; testCase++  )
    {
        std::cout << "test case " << testCase << "\n\n";
        // TEST CASE 0: NO VARIABLE LOVE NUMBER ESTIMATED -> NO CONSTRAINT
        // TEST CASE 1: VARIABLE LOVE NUMBER ESTIMATED (I+E) BUT NOT FOR DEGREE 2 -> NO CONSTRAINT
        // TEST CASE 2: VARIABLE LOVE NUMBER OF DEGREE 2 ESTIMATED (I+E) BUT NO COMPLEX VALUE -> NO CONSTRAINT
        // TEST CASE 3: FULL DEGREE (Io) + SINGLE DEGREE (Europa) ESTIMATED BUT NO ORDER 2 -> NO CONSTRAINT
        // TEST CASE 4: FULL DEGREE BUT SEVERAL ORDER + INCONSISTENT VALUE (IO) -> ONLY ORDER 2 IS CONSTRAINED + INCONSISTENT VALUE OVERWRITTEN
        // & SINGLE DEGREE + INCONSISTENT VALUE (EUROPA) -> INCONSISTENT VALUE OVERWRITTEN

        // Define list of parameters to estimate
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;

        // Initial states
        for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
        {
            parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                    bodiesToPropagate.at( i ), initialStates.segment( i * 6, 6 ), globalFrameOrigin, globalFrameOrientation ) );
        }


        if ( testCase != 0 )
        {
            // Define whether the complex part of k2 is estimated
            bool useComplexValue = true;
            if ( testCase == 2 )
            {
                useComplexValue = false;
            }

            // Define which degree is estimated for the Love numbers
            int degreeLoveNumber = 2;
            if ( testCase == 1 )
            {
                degreeLoveNumber = 3;
            }

            // Define which orders are estimated for Io's k2
            std::vector< int > ordersK2Io = { 0, 1, 2 };
            if ( testCase == 3 )
            {
                ordersK2Io = { 0, 1 };
            }

            // Single degree tidal Love number for Io
            parameterNames.push_back( std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
                    "Io", degreeLoveNumber, ordersK2Io, "Jupiter", useComplexValue ) );

            // Single degree tidal Love number for Europa
            parameterNames.push_back( std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >(
                    "Europa", degreeLoveNumber, "Jupiter", useComplexValue ) );
        }

        // Inv(Q) for Io
        parameterNames.push_back( std::make_shared< InverseTidalQualityFactorEstimatableParameterSettings >( "Io", "Jupiter" ) ) ;
        // Tidal time lag for Europa
        parameterNames.push_back( std::make_shared< DirectTidalTimeLagEstimatableParameterSettings >( "Europa", "Jupiter" ) );

        // Create parameters to estimate
        bool deactivateConstraints = false;
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate;

        bool isExceptionFound = false;
        try
        {
            parametersToEstimate = createParametersToEstimate< double >(
                    parameterNames, bodies, propagatorSettings, std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > >( ), deactivateConstraints );

            // Print parameters and constraints
            printEstimatableParameterEntries( parametersToEstimate );
            printInterParameterConstraints( parametersToEstimate );

            std::cout << "constraints size: " << parametersToEstimate->getConstraintSize( ) << "\n\n";
        }
        catch( std::runtime_error const& )
        {
            isExceptionFound = true;
        }

        if ( testCase == 0 || testCase == 1 || testCase == 2 )
        {
            BOOST_CHECK( parametersToEstimate->getConstraintSize( ) == 0 );
        }

        // Check if runtime error has occurred
        else if ( testCase == 3 )
        {
            BOOST_CHECK( isExceptionFound );
        }

        else if ( testCase == 4 )
        {
            BOOST_CHECK( parametersToEstimate->getConstraintSize( ) == 2 );

            Eigen::VectorXd parametersValues = parametersToEstimate->getFullParameterValues< double >( );


            std::vector< std::pair< int, int > > indexK2Io = parametersToEstimate->getIndicesForParameterType(
                    std::make_pair( single_degree_variable_tidal_love_number, std::make_pair( "Io", "" ) ) );
            std::vector< std::pair< int, int > > indexInvQIo = parametersToEstimate->getIndicesForParameterType(
                    std::make_pair( inverse_tidal_quality_factor, std::make_pair( "Io", "" ) ) );
            std::vector< std::pair< int, int > > indexK2Europa = parametersToEstimate->getIndicesForParameterType(
                    std::make_pair( full_degree_tidal_love_number, std::make_pair( "Europa", "" ) ) );
            std::vector< std::pair< int, int > > indexTimeLagEuropa = parametersToEstimate->getIndicesForParameterType(
                    std::make_pair( direct_dissipation_tidal_time_lag, std::make_pair( "Europa", "" ) ) );

            Eigen::VectorXd k2Io = parametersValues.segment( indexK2Io.at( 0 ).first, indexK2Io.at( 0 ).second );
            Eigen::VectorXd invQIo = parametersValues.segment( indexInvQIo.at( 0 ).first, indexInvQIo.at( 0 ).second );
            Eigen::VectorXd k2Europa = parametersValues.segment( indexK2Europa.at( 0 ).first, indexK2Europa.at( 0 ).second );
            Eigen::VectorXd timeLagEuropa = parametersValues.segment( indexTimeLagEuropa.at( 0 ).first, indexTimeLagEuropa.at( 0 ).second );

            std::cout << "k2Io: " << k2Io << "\n\n";
            std::cout << "invQIo: " << invQIo << "\n\n";
            std::cout << "k2Europa: " << k2Europa << "\n\n";
            std::cout << "timeLagEuropa: " << timeLagEuropa << "\n\n";
        }



    }



//    unsigned int nbParameters = parametersToEstimate->getParameterSetSize( );
//    std::cout << "nb parameters: " << nbParameters << "\n\n";
//
//    // Print parameters and constraints
//    printEstimatableParameterEntries( parametersToEstimate );
//    printInterParameterConstraints( parametersToEstimate );
//
//    Eigen::MatrixXd constraintStateMultiplier;
//    Eigen::VectorXd constraintRightHandSide;
//    parametersToEstimate->getConstraints( constraintStateMultiplier, constraintRightHandSide );
//
//    std::cout << "constraintStateMultiplier: " << std::endl;
//    std::cout << constraintStateMultiplier << "\n\n";
//
//    std::cout << "parameter values: " << std::endl;
//    std::cout << parametersToEstimate->getFullParameterValues< double >( ).transpose( ) << "\n\n";





//    // Define links in simulation.
//    std::vector< LinkDefinition > linkEnds;
//    linkEnds.resize( 3 );
//    linkEnds[ 0 ][ observed_body ] = std::make_pair< std::string, std::string >( "spacecraft", "" );
////    LinkDefinition linkEndDef;
////    std::string moon = "Io";
////    linkEndDef[ observed_body ] = LinkEndId( moon, "" );
//    linkEnds[ 1 ][ observed_body ] = LinkEndId( "Io", "" );
//    linkEnds[ 2 ][ observed_body ] = LinkEndId( "Europa", "" );
////
////    linkEnds[ 1 ][ observed_body ] = std::make_pair< std::string, std::string >( moon, "" );
////    for ( unsigned int i = 0 ; i < moonsToPropagate.size( ) ; i++ )
////    {
////        std::string moon = moonsToPropagate.at( i );
////        linkEnds[ i+1 ][ observed_body ] = std::make_pair< std::string, std::string >( moon, "" );
////    }
//
//    // Create observation settings list
//    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
//    for ( auto linkEnd : linkEnds )
//    {
//        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( position_observable, linkEnd ) );
//    }
//
//    // Create orbit determination object.
//    OrbitDeterminationManager< double, double > orbitDeterminationManager = OrbitDeterminationManager< double, double >(
//            bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );
//
//    std::map< double, Eigen::VectorXd > stateHistory = std::dynamic_pointer_cast< SingleArcVariationalEquationsSolver< > >(
//            orbitDeterminationManager.getVariationalEquationsSolver( ) )->getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
//    input_output::writeDataMapToTextFile( stateHistory, "stateHistory.dat",
//                                          "/home/mfayolle/Documents/PHD/04 - year 4/tudat-bundle/tudat/applications/Dynamics/SimulationAnalysis/" );
//
//    Eigen::VectorXd initialParameterEstimate =
//            parametersToEstimate->template getFullParameterValues< double >( );
//
//    // Define observation simulation times
//    std::vector< double > obsTimesMoon;
//    for ( double time = initialEpoch + 10.0 * timeStep ; time < finalEpoch - 10.0 * timeStep ; time += 3600.0  )
//    {
//        obsTimesMoon.push_back( time );
//    }
//    std::vector< double > obsTimesSpacecraft;
//    for ( double time = initialEpoch + 10.0 * timeStep ; time < initialEpoch + 7.0 * 86400.0 /*finalEpoch - 10.0 * timeStep*/ ; time += 300.0  )
//    {
//        obsTimesSpacecraft.push_back( time );
//    }
//
//    // Create measurement simulation input
//    std::vector< std::shared_ptr< ObservationSimulationSettings< > > > measurementSimulationInput;
//    measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< > >(
//            position_observable, linkEnds[ 0 ], obsTimesSpacecraft, observed_body ) );
//    for ( unsigned int i = 1 ; i < linkEnds.size( ) ; i++ )
//    {
//        measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< > >(
//                position_observable, linkEnds[ i ], obsTimesMoon, observed_body ) );
//    }
////    measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< > >(
////            velocity_observable, linkEnds[ 0 ], obsTimesMoon, observed_body ) );
////    measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< > >(
////            velocity_observable, linkEnds[ 1 ], obsTimes, observed_body ) );
//
////    double positionNoise = 1.0e-4;
////    double velocityNoise = 1.0e-7;
////    addGaussianNoiseToSingleObservationSimulationSettings( measurementSimulationInput[ 0 ], positionNoise );
////    addGaussianNoiseToSingleObservationSimulationSettings( measurementSimulationInput[ 1 ], positionNoise );
////    addGaussianNoiseToSingleObservationSimulationSettings( measurementSimulationInput[ 0 ], velocityNoise );
////    addGaussianNoiseToSingleObservationSimulationSettings( measurementSimulationInput[ 1 ], velocityNoise );
//
//    // Simulate observations
//    std::shared_ptr< ObservationCollection< > >  observationsAndTimes = simulateObservations< double, double >(
//            measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );
//
//
//    // Perturb parameter values
//    Eigen::VectorXd truthParameters = initialParameterEstimate;
////    initialParameterEstimate[ 0 ] *= 1.0001;
//
////    initialParameterEstimate[ truthParameters.size( ) - 7 ] += 0.005;
////    initialParameterEstimate[ truthParameters.size( ) - 11 ] += 0.005;
////    initialParameterEstimate[ truthParameters.size( ) - 14 ] = getInvQFromLoveNumber(
////            std::complex< double >( k2ReMoon, k2ImIo + 0.005 ) );
//
//    double perturbedk2Im = truthParameters[ nbParameters - 1 ] + 0.0005;
//    double perturbedInvQ = getInvQFromLoveNumber( std::complex< double >( k2ReMoon, perturbedk2Im ) );
//    double perturbedTimeLag = getTimeLagMoonFromLoveNumber( std::complex< double >( k2ReMoon, perturbedk2Im ), tidalPeriods.at( "Io" ) );
//    initialParameterEstimate[ nbParameters - 1 ] = perturbedk2Im;
////    initialParameterEstimate[ nbParameters - 3 ] = perturbedk2Im;
////    initialParameterEstimate[ nbParameters - 5 ] = perturbedInvQ;
//    initialParameterEstimate[ nbParameters - 3 ] = perturbedInvQ;
////    initialParameterEstimate[ nbParameters - 2 ] = perturbedInvQ;
//
//    Eigen::VectorXd initialPerturbation = initialParameterEstimate - truthParameters;
//    double perturbationTidalParameter = initialPerturbation[ nbParameters - 3 ];
//    double perturbationImK2 = initialPerturbation[ nbParameters - 1 ];
//    std::cout << "perturbationImK2: " << perturbationImK2 << "\n\n";
//    std::cout << "perturbationTidalParameter: " << perturbationTidalParameter << "\n\n";
//
//    parametersToEstimate->resetParameterValues( initialParameterEstimate );
//
//    // A priori covariance
//    Eigen::MatrixXd invAprioriCovariance = Eigen::MatrixXd::Zero( nbParameters, nbParameters );
//
//    // Estimate initial states and tidal parameters
//    std::shared_ptr< EstimationInput< double, double > > estimationInput = std::make_shared< EstimationInput< double, double > >(
//            observationsAndTimes, invAprioriCovariance );
//    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 1 ) );
//    estimationInput->applyFinalParameterCorrection_ = true;
//    estimationInput->defineEstimationSettings( true, true, true, true, true, true );
//
////    std::map< observation_models::ObservableType, double > weightPerObservable;
////    weightPerObservable[ position_observable ] = 1.0 / ( positionNoise * positionNoise );
////    weightPerObservable[ velocity_observable ] = 1.0 / ( velocityNoise * velocityNoise );
////
////    estimationInput->setConstantPerObservableWeightsMatrix( weightPerObservable );
//
////    std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );
////
////    // Check if parameters are correctly estimated
////    Eigen::VectorXd estimatedParametervalues = estimationOutput->parameterHistory_.at( estimationOutput->parameterHistory_.size( ) - 1 );
////    Eigen::VectorXd residuals = estimationOutput->residuals_;
////
//////    input_output::writeMatrixToFile( residuals, "residuals" + indexFiles + ".dat", 16, getOutputPath( ) + outputFolder );
//////    input_output::writeMatrixToFile( estimationOutput->getNormalizedWeightedDesignMatrix( ), "normalisedWeightedDesignMatrix" + indexFiles + ".dat", 16, getOutputPath( ) + outputFolder );
//////    input_output::writeMatrixToFile( estimationOutput->getCorrelationMatrix( ), "correlations" + indexFiles + ".dat", 16, getOutputPath( ) + outputFolder );
////
////    std::cout << "-----------------------------------" << "\n\n";
////    std::cout << "truthParameters: " << truthParameters.transpose( ) << "\n\n";
////    std::cout << "initialParameterEstimate: " << initialParameterEstimate.transpose( ) << "\n\n";
////    std::cout << "estimatedParametervalues: " << estimatedParametervalues.transpose( ) << "\n\n";
////
////    Eigen::VectorXd parametersUpdate = estimatedParametervalues - initialParameterEstimate;
////    std::cout << "parametersUpdate: " << parametersUpdate.transpose( ) << "\n\n";
////
////    std::cout << "-----------------------------------" << "\n\n";
////    std::cout << "old k2: " << std::complex< double >( k2ReMoon, k2ImIo ) << "\n\n";
////    std::cout << "old invQ: " << invQIo << "\n\n";
////    std::cout << "old time lag: " << getTimeLagMoonFromLoveNumber( std::complex< double >( k2ReMoon, k2ImIo ), tidalPeriods.at( "Io" ) ) << "\n\n";
////
////    std::complex< double > newK2 = std::complex< double >( estimatedParametervalues[ nbParameters - 2 ], estimatedParametervalues[ nbParameters - 1 ] );
////    std::cout << "new k2: " << newK2 << "\n\n";
////    std::cout << "new tidal parameter: " << estimatedParametervalues[ nbParameters - 3 ] << "\n\n";
////
//////    std::cout << "old k2: " << std::complex< double >( k2ReMoon, k2ImIo ) << "\n\n";
//////    std::cout << "norm k2: " << std::sqrt( k2ReMoon * k2ReMoon + k2ImIo * k2ImIo ) << "\n\n";
//////    std::cout << "ratio InvQ / Im(k2): " << 1.0 / std::sqrt( k2ReMoon * k2ReMoon + k2ImIo * k2ImMoon ) << "\n\n";
////
////    double theoreticalInvQ = getInvQFromLoveNumber( newK2 );
////    double theoreticalTimeLag = getTimeLagMoonFromLoveNumber( newK2, tidalPeriods.at( "Io" ) );
////    std::cout << "theoreticalInvQ: " << theoreticalInvQ << "\n\n";
////    std::cout << "theoreticalTimeLag: " << theoreticalTimeLag << "\n\n";
////    std::cout << "estimated tidal parameter: " << estimatedParametervalues[ nbParameters - 3 ] << "\n\n";
////
//////    std::vector< std::shared_ptr< gravitation::GravityFieldVariations > > allVariations = bodies.at( "Io" )->getGravityFieldVariationSet( )->getVariationObjects( );
//////    std::cout << "size field variations: " << allVariations.size( ) << "\n\n";
//////    std::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > solidTidesVariations =
//////            std::dynamic_pointer_cast< gravitation::BasicSolidBodyTideGravityFieldVariations >( allVariations[ 0 ] );
//////    std::map< int, std::vector< std::complex< double > > > loveNumbers = solidTidesVariations->getLoveNumbers( );
//////    for ( auto itr : loveNumbers )
//////    {
//////        std::cout << itr.first << "\n\n";
//////        for ( unsigned int i = 0 ; i < itr.second.size( ) ; i++ )
//////        {
//////            std::cout << i << " - Re: " << itr.second[ i ].real( ) << " & Im: " << itr.second[ i ].imag( ) << "\n\n";
//////        }
//////    }
////
//////    std::cout << "inverse normalised covariance: " << "\n\n";
//////    std::cout << estimationOutput->getNormalizedInverseCovarianceMatrix( ) << "\n\n";
//////    std::cout << "normalisation terms: " << "\n\n";
//////    std::cout << estimationOutput->getNormalizationTerms( ).transpose( ) << "\n\n";
////
//////    input_output::writeMatrixToFile( estimationOutput->getCorrelationMatrix( ), "correlations.dat", 16,
//////                                     "/home/mfayolle/Documents/PHD/04 - year 4/tudat-bundle/tudat/applications/Dynamics/SimulationAnalysis/" );
//////
//////    input_output::writeMatrixToFile( estimationOutput->getNormalizedWeightedDesignMatrix( ), "weightedPartials.dat", 16,
//////                                     "/home/mfayolle/Documents/PHD/04 - year 4/tudat-bundle/tudat/applications/Dynamics/SimulationAnalysis/" );
//////
//////    std::cout << "formal errors: " << "\n\n";
//////    std::cout << estimationOutput->getFormalErrorVector( ).transpose( ) << "\n\n";
////
//////    std::cout << "size formal errors: " << estimationOutput->getFormalErrorVector( ).size( ) << "\n\n";
//////    std::cout << "global constraints size: " << parametersToEstimate->getConstraintSize( ) << "\n\n";
//////    std::cout << "size parameters to estimate: " << parametersToEstimate->template getFullParameterValues< double >( ).size( ) << "\n\n";
////
////
//    // Test custom constraints
//
//    std::function< Eigen::MatrixXd( Eigen::VectorXd ) > constraintFunctionK2 = [](Eigen::VectorXd parametersVector)
//    {
//        Eigen::MatrixXd constraint = Eigen::MatrixXd::Zero( 1, 2 );
//        constraint(0, 1) = 1.0;
//        return constraint;
//    };
//    Eigen::MatrixXd constraintConstantK2 = Eigen::MatrixXd::Zero( 1, 2 );
//    constraintConstantK2(0, 1) = 1.0;
//
//    std::function< Eigen::MatrixXd( Eigen::VectorXd ) > constraintFunctionInvQ = [](Eigen::VectorXd parametersVector)
//    {
//        Eigen::MatrixXd constraint = Eigen::MatrixXd::Zero( 1, 1 );
//        constraint(0, 0) = - 1.0;
//        return constraint;
//    };
//    Eigen::MatrixXd constraintConstantInvQ = Eigen::MatrixXd::Zero( 1, 1 );
//    constraintConstantInvQ(0, 0) = - 1.0;
//
//    std::vector< std::function< Eigen::MatrixXd( Eigen::VectorXd ) > > linearisedFunctionConstraint = { constraintFunctionK2, constraintFunctionInvQ };
//    std::vector< Eigen::MatrixXd > linearisedConstantConstraint = { constraintConstantK2, constraintConstantInvQ };
//
//    std::vector< estimatable_parameters::EstimatebleParametersEnum > parametersInvolved =
//            { single_degree_variable_tidal_love_number, inverse_tidal_quality_factor };
//    std::vector< std::pair< std::string, std::string > > bodiesInvolvedConstraint; // = { std::make_pair( "Io", "" ), std::make_pair( "Io", "" ) };
//
//    bool applyPerBody = true;
//    bool applyToAllBodies = true;
//    createCustomParameterConstraints( parametersToEstimate, linearisedFunctionConstraint, parametersInvolved, bodiesInvolvedConstraint,
//                                      applyPerBody, applyToAllBodies );
//
//    createConstantCustomParameterConstraints( parametersToEstimate, linearisedConstantConstraint, parametersInvolved,
//                                              bodiesInvolvedConstraint, applyPerBody, applyToAllBodies );
//
//    Eigen::MatrixXd constraintStateMultiplier;
//    Eigen::VectorXd constraintRightHandSide;
//    parametersToEstimate->getConstraints( constraintStateMultiplier, constraintRightHandSide );
//
//    std::cout << "constraintStateMultiplier: " << "\n\n";
//    std::cout << constraintStateMultiplier << "\n\n";

//    std::cout << "constraintRightHandSide: " << "\n\n";
//    std::cout << constraintRightHandSide << "\n\n";
//
//    std::cout << "parameters values: " << "\n\n";
//    std::cout << parametersToEstimate->getFullParameterValues< double >( ).transpose( ) << "\n\n";
//
//    printInterParameterConstraints( parametersToEstimate );
}

BOOST_AUTO_TEST_CASE( testConstrainedEstimation )
{
//    std::cout.precision( 20 );
//
//    loadStandardSpiceKernels( );
//
//    double initialEpoch = 32.0 * physical_constants::JULIAN_YEAR;
//    double finalEpoch = initialEpoch + 20.0 * physical_constants::JULIAN_DAY;
//
//    std::string globalFrameOrientation = "J2000";
//    std::string globalFrameOrigin = "Jupiter";
//
//    std::complex< double > k2Moon = std::complex< double >( 0.2, 0.01 );
//    std::complex< double > k2Io = std::complex< double >( 0.2, 0.01 );
//    double k2MoonNorm = std::sqrt( std::norm( k2Moon ) );
//    double k2IoNorm = std::sqrt( std::norm( k2Io ) );
//    std::cout << "k2MoonNorm: " << k2MoonNorm << "\n\n";
//    std::cout << "k2IoNorm: " << k2IoNorm << "\n\n";
//
//    double invQMoon = getInvQFromLoveNumber( k2Moon );
//    double invQIo = getInvQFromLoveNumber( k2Io );
//
//
//    // Define bodies settings for simulation
//    std::vector< std::string > bodiesToCreate = { "Io", "Europa", "Ganymede", "Callisto", "Jupiter" };
//
//    // Create body objects.
//    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, initialEpoch, finalEpoch, globalFrameOrigin, globalFrameOrientation );
//    for ( unsigned int i = 0 ; i < bodiesToCreate.size( ) - 1 ; i++ )
//    {
//        bodySettings.at( bodiesToCreate.at( i ) )->rotationModelSettings = std::make_shared< SynchronousRotationModelSettings >(
//                "Jupiter", globalFrameOrientation, "IAU_" + bodiesToCreate.at( i ) );
//    }
//
//    // Add temporal variations to gravity fields
//    addGravityFieldVariationsMoons( bodySettings, k2Moon, k2Io );
//
//
//    // Crate system of bodies
//    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
//    bodies.createEmptyBody( "spacecraft" );
//
//    // Define bodies to propagate
//    std::vector< std::string > bodiesToPropagate = { "Io", "Europa" };
//    unsigned int nbMoons = bodiesToPropagate.size( );
//    std::vector< std::string > centralBodies;
//    for ( unsigned int i = 0 ; i < nbMoons ; i++ )
//    {
//        centralBodies.push_back( "Jupiter" );
//    }
//
//    // Add spacecraft around Io
//    bodiesToPropagate.push_back( "spacecraft" );
//    centralBodies.push_back( "Io" );
//
//
//    // Create map to store tidal periods
//    std::map< std::string, double > tidalPeriods;
//    for ( unsigned int k = 0 ; k < nbMoons ; k++ )
//    {
//        tidalPeriods[ bodiesToPropagate.at( k ) ] = getTidalPeriodMoon( bodiesToPropagate.at( k ), initialEpoch, bodies );
//    }
//
//




//    // Define accelerations
//    int maximumMoonDegree = 2;
//    int maximumMoonOrder = 2;
//    int maximumJupiterDegree = 10;
//    int maximumJupiterOrder = 10;
//
//    SelectedAccelerationMap accelerationSettings;
//
//
//    // accelerations moons
//    if ( propagateMoons )
//    {
//        for ( const auto moon : moonsToPropagate )
//        {
//            accelerationSettings[ moon ][ "Jupiter" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >(
//                    maximumJupiterDegree, maximumJupiterOrder, maximumMoonDegree, maximumMoonOrder, 0, 0, shVariationsMoon, shVariationsMoon, shVariationsMoon ) );
//
//            if ( includeTidalForceOnMoons )
//            {
////                double periodTidesMoon = getTidalPeriodMoon( moon, initialEpoch, bodies );
//                std::cout << "periodTides for " << moon << ": " << tidalPeriods.at( moon ) << "\n\n";
//
//                // Tides in the moon
//                if ( moon == "Io" )
//                {
//                    double timeLagIo = getTimeLagMoonFromLoveNumber( std::complex< double >( k2ReMoon, k2ImIo ), tidalPeriods.at( moon ) ); // periodTidesMoon * std::atan( invQIo ) / ( 2.0 * mathematical_constants::PI );
//                    accelerationSettings[ moon ][ "Jupiter" ].push_back( std::make_shared< DirectTidalDissipationAccelerationSettings >(
//                            k2FullIo, invQIo, tidalPeriods.at( moon ), true, false, false ) );
////                accelerationSettings[ moon ][ "Jupiter" ].push_back( std::make_shared< DirectTidalDissipationAccelerationSettings >(
////                        k2FullIo, timeLagIo, true, false, false ) );
//                }
//                else
//                {
//                    double timeLagMoon = getTimeLagMoonFromLoveNumber( std::complex< double >( k2ReMoon, k2ImMoon ), tidalPeriods.at( moon ) ); // periodTidesMoon * std::atan( invQMoon ) / ( 2.0 * mathematical_constants::PI );
//                    accelerationSettings[ moon ][ "Jupiter" ].push_back( std::make_shared< DirectTidalDissipationAccelerationSettings >(
//                            k2FullMoon, invQMoon, tidalPeriods.at( moon ), true, false, false ) );
////                accelerationSettings[ moon ][ "Jupiter" ].push_back( std::make_shared< DirectTidalDissipationAccelerationSettings >(
////                        k2FullMoon, timeLagMoon, true, false, false ) );
//                }
//            }
//        }
//    }
//
//
//    // Acceleration spacecraft
//    accelerationSettings[ "spacecraft" ][ "Io" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >(
//            maximumMoonDegree, maximumMoonOrder, false, shVariationsSc ) );
//
////    for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
////    {
////        if ( setup == 0 || setup == 1 || setup == 4 )
////        {
////            accelerationSettings[ bodiesToPropagate.at( i ) ][ "Jupiter" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >(
////                    maximumJupiterDegree, maximumJupiterOrder, maximumSatelliteDegree, maximumSatelliteOrder, false ) );
////        }
////        else
////        {
////            accelerationSettings[ bodiesToPropagate.at( i ) ][ "Jupiter" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >(
////                    maximumJupiterDegree, maximumJupiterOrder, maximumSatelliteDegree, maximumSatelliteOrder, true ) );
////        }
//
//
////        for( unsigned int j = 0; j < bodiesToPropagate.size( ); j++ )
////        {
////            if( bodiesToPropagate.at( j ) != bodiesToPropagate.at( i ) )
////            {
////                if ( setup == 0 || setup == 1 || setup == 3 || setup == 4 )
////                {
////                    accelerationSettings[ bodiesToPropagate.at( i ) ][ bodiesToPropagate.at( j ) ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >(
////                            maximumSatelliteDegree, maximumSatelliteOrder, false, false ) );
////                }
////                else
////                {
////                    accelerationSettings[ bodiesToPropagate.at( i ) ][ bodiesToPropagate.at( j ) ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >(
////                            maximumSatelliteDegree, maximumSatelliteOrder, false, true ) );
////                }
////
////            }
////        }
////    }
//    basic_astrodynamics::AccelerationMap accelerations = createAccelerationModelsMap(
//            bodies, accelerationSettings, bodiesToPropagate, centralBodies );
//
//
//    // Define propagator settings
//    Eigen::Vector6d initialStateSpacecraftKeplerian = Eigen::Vector6d::Zero( );
//    initialStateSpacecraftKeplerian[ semiMajorAxisIndex ] = 500.0e3;
//    initialStateSpacecraftKeplerian[ inclinationIndex ] = 20.0 * mathematical_constants::PI / 180.0;
//    Eigen::Vector6d initialStateSpacecraft = orbital_element_conversions::convertKeplerianToCartesianElements(
//            initialStateSpacecraftKeplerian, bodies.at( "Io" )->getGravitationalParameter( ) );
//
//    Eigen::VectorXd initialStatesMoons = getInitialStatesOfBodies( moonsToPropagate, centralBodiesMoons, bodies, initialEpoch );
//
//    Eigen::VectorXd initialStates = Eigen::VectorXd::Zero( 6 * bodiesToPropagate.size( ) );
//    if ( propagateMoons )
//    {
//        initialStates.segment( 0, nbMoons * 6 ) = initialStatesMoons;
//        initialStates.segment( nbMoons * 6, 6 ) = initialStateSpacecraft;
//    }
//    else
//    {
//        initialStates.segment( 0, 6 ) = initialStateSpacecraft;
//    }
////    std::cout << "initial states: " << initialStates.transpose( ) << "\n\n";
////    initialStates.segment( 0, 6 ) = initialStateIo;
////    initialStates.segment( 6, 6 ) = initialStateSpacecraft;
//
//    // Define integrator settings
//    double timeStep = 50.0; // 50.0;
//    std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
//            std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
//                    initialEpoch, timeStep, CoefficientSets::rungeKuttaFehlberg78, timeStep, timeStep, 1.0e3, 1.0e3 );
//
////    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
////    dependentVariables.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
////            basic_astrodynamics::spherical_harmonic_gravity, "Io", "Europa" ) );
////    dependentVariables.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
////            basic_astrodynamics::spherical_harmonic_gravity, "Io", "Ganymede" ) );
////    dependentVariables.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
////            basic_astrodynamics::spherical_harmonic_gravity, "Io", "Callisto" ) );
//
//    std::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings =
//            std::make_shared< TranslationalStatePropagatorSettings< > >(
//                    centralBodies, accelerations, bodiesToPropagate, initialStates, initialEpoch, integratorSettings,
//                    std::make_shared< PropagationTimeTerminationSettings >( finalEpoch ), cowell/*, dependentVariables*/ );
//
////    // Propagate dynamics
////    SingleArcDynamicsSimulator< > dynamicsSimulator = SingleArcDynamicsSimulator< >( bodies, propagatorSettings, true );
////    std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
////    std::map< double, Eigen::VectorXd > dependentVariablesHistory = dynamicsSimulator.getDependentVariableHistory( );
////
////    input_output::writeDataMapToTextFile( stateHistory, "stateHistory.dat",
////                                          "/home/mfayolle/Documents/PHD/04 - year 4/tudat-bundle/tudat/applications/Dynamics/SimulationAnalysis/" );
////
////    input_output::writeDataMapToTextFile( dependentVariablesHistory, "dependentVariablesHistory.dat",
////                                          "/home/mfayolle/Documents/PHD/04 - year 4/tudat-bundle/tudat/applications/Dynamics/SimulationAnalysis/");
//
//
//    // Create parameters to estimate
//    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
//    for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
//    {
//        parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
//                bodiesToPropagate.at( i ), initialStates.segment( i * 6, 6 ), globalFrameOrigin, globalFrameOrientation ) );
//    }
//
//    for ( auto moon : moonsToPropagate )
//    {
////        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
////                2, 0, maximumMoonDegree, maximumMoonOrder, moon, spherical_harmonics_cosine_coefficient_block ) );
////        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
////                2, 1,  maximumMoonDegree, maximumMoonOrder, moon, spherical_harmonics_sine_coefficient_block ) );
//
////        parameterNames.push_back( std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
////                moon, 2, std::vector< int >{ 2 }, "Jupiter", true ) );
//        if ( shVariationsSc || shVariationsMoon )
//        {
//            if ( moon == "Io" )
//            {
//                bool useComplexValue = true;
//                parameterNames.push_back( std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
//                        moon, 2, std::vector< int >{ 1, 2 }, "Jupiter", useComplexValue ) );
////            parameterNames.push_back( std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( moon, 2, "Jupiter", useComplexValue ) );
//            }
//            if ( moon == "Europa" )
//            {
//                bool useComplexValue = true;
//                parameterNames.push_back( std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
//                        moon, 2, std::vector< int >{ 2 }, "Jupiter", useComplexValue ) );
//            }
//        }
//
//        if ( includeTidalForceOnMoons )
//        {
//            if ( moon == "Io" )
//            {
//                parameterNames.push_back( std::make_shared< InverseTidalQualityFactorEstimatableParameterSettings >( moon, "Jupiter" ) ) ;
////            parameterNames.push_back( std::make_shared< DirectTidalTimeLagEstimatableParameterSettings >( moon, "Jupiter" ) ) ;
//            }
//            if ( moon == "Europa" )
//            {
//                parameterNames.push_back( std::make_shared< DirectTidalTimeLagEstimatableParameterSettings >( moon, "Jupiter" ) );
//            }
//        }
//    }
//
//    bool deactivateConstraints = false;
//    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate = createParametersToEstimate< double >(
//            parameterNames, bodies, propagatorSettings, std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > >( ),
//            deactivateConstraints );
//    printEstimatableParameterEntries( parametersToEstimate );
//    unsigned int nbParameters = parametersToEstimate->getParameterSetSize( );
////    std::cout << "nb parameters: " << nbParameters << "\n\n";
//
//    // Define links in simulation.
//    std::vector< LinkDefinition > linkEnds;
//    linkEnds.resize( 3 );
//    linkEnds[ 0 ][ observed_body ] = std::make_pair< std::string, std::string >( "spacecraft", "" );
////    LinkDefinition linkEndDef;
////    std::string moon = "Io";
////    linkEndDef[ observed_body ] = LinkEndId( moon, "" );
//    linkEnds[ 1 ][ observed_body ] = LinkEndId( "Io", "" );
//    linkEnds[ 2 ][ observed_body ] = LinkEndId( "Europa", "" );
////
////    linkEnds[ 1 ][ observed_body ] = std::make_pair< std::string, std::string >( moon, "" );
////    for ( unsigned int i = 0 ; i < moonsToPropagate.size( ) ; i++ )
////    {
////        std::string moon = moonsToPropagate.at( i );
////        linkEnds[ i+1 ][ observed_body ] = std::make_pair< std::string, std::string >( moon, "" );
////    }
//
//    // Create observation settings list
//    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
//    for ( auto linkEnd : linkEnds )
//    {
//        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( position_observable, linkEnd ) );
//    }
//
//    // Create orbit determination object.
//    OrbitDeterminationManager< double, double > orbitDeterminationManager = OrbitDeterminationManager< double, double >(
//            bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );
//
//    std::map< double, Eigen::VectorXd > stateHistory = std::dynamic_pointer_cast< SingleArcVariationalEquationsSolver< > >(
//            orbitDeterminationManager.getVariationalEquationsSolver( ) )->getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
//    input_output::writeDataMapToTextFile( stateHistory, "stateHistory.dat",
//                                          "/home/mfayolle/Documents/PHD/04 - year 4/tudat-bundle/tudat/applications/Dynamics/SimulationAnalysis/" );
//
//    Eigen::VectorXd initialParameterEstimate =
//            parametersToEstimate->template getFullParameterValues< double >( );
//
//    // Define observation simulation times
//    std::vector< double > obsTimesMoon;
//    for ( double time = initialEpoch + 10.0 * timeStep ; time < finalEpoch - 10.0 * timeStep ; time += 3600.0  )
//    {
//        obsTimesMoon.push_back( time );
//    }
//    std::vector< double > obsTimesSpacecraft;
//    for ( double time = initialEpoch + 10.0 * timeStep ; time < initialEpoch + 7.0 * 86400.0 /*finalEpoch - 10.0 * timeStep*/ ; time += 300.0  )
//    {
//        obsTimesSpacecraft.push_back( time );
//    }
//
//    // Create measurement simulation input
//    std::vector< std::shared_ptr< ObservationSimulationSettings< > > > measurementSimulationInput;
//    measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< > >(
//            position_observable, linkEnds[ 0 ], obsTimesSpacecraft, observed_body ) );
//    for ( unsigned int i = 1 ; i < linkEnds.size( ) ; i++ )
//    {
//        measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< > >(
//                position_observable, linkEnds[ i ], obsTimesMoon, observed_body ) );
//    }
////    measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< > >(
////            velocity_observable, linkEnds[ 0 ], obsTimesMoon, observed_body ) );
////    measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< > >(
////            velocity_observable, linkEnds[ 1 ], obsTimes, observed_body ) );
//
////    double positionNoise = 1.0e-4;
////    double velocityNoise = 1.0e-7;
////    addGaussianNoiseToSingleObservationSimulationSettings( measurementSimulationInput[ 0 ], positionNoise );
////    addGaussianNoiseToSingleObservationSimulationSettings( measurementSimulationInput[ 1 ], positionNoise );
////    addGaussianNoiseToSingleObservationSimulationSettings( measurementSimulationInput[ 0 ], velocityNoise );
////    addGaussianNoiseToSingleObservationSimulationSettings( measurementSimulationInput[ 1 ], velocityNoise );
//
//    // Simulate observations
//    std::shared_ptr< ObservationCollection< > >  observationsAndTimes = simulateObservations< double, double >(
//            measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );
//
//
//    // Perturb parameter values
//    Eigen::VectorXd truthParameters = initialParameterEstimate;
////    initialParameterEstimate[ 0 ] *= 1.0001;
//
////    initialParameterEstimate[ truthParameters.size( ) - 7 ] += 0.005;
////    initialParameterEstimate[ truthParameters.size( ) - 11 ] += 0.005;
////    initialParameterEstimate[ truthParameters.size( ) - 14 ] = getInvQFromLoveNumber(
////            std::complex< double >( k2ReMoon, k2ImIo + 0.005 ) );
//
//    double perturbedk2Im = truthParameters[ nbParameters - 1 ] + 0.0005;
//    double perturbedInvQ = getInvQFromLoveNumber( std::complex< double >( k2ReMoon, perturbedk2Im ) );
//    double perturbedTimeLag = getTimeLagMoonFromLoveNumber( std::complex< double >( k2ReMoon, perturbedk2Im ), tidalPeriods.at( "Io" ) );
//    initialParameterEstimate[ nbParameters - 1 ] = perturbedk2Im;
////    initialParameterEstimate[ nbParameters - 3 ] = perturbedk2Im;
////    initialParameterEstimate[ nbParameters - 5 ] = perturbedInvQ;
//    initialParameterEstimate[ nbParameters - 3 ] = perturbedInvQ;
////    initialParameterEstimate[ nbParameters - 2 ] = perturbedInvQ;
//
//    Eigen::VectorXd initialPerturbation = initialParameterEstimate - truthParameters;
//    double perturbationTidalParameter = initialPerturbation[ nbParameters - 3 ];
//    double perturbationImK2 = initialPerturbation[ nbParameters - 1 ];
//    std::cout << "perturbationImK2: " << perturbationImK2 << "\n\n";
//    std::cout << "perturbationTidalParameter: " << perturbationTidalParameter << "\n\n";
//
//    parametersToEstimate->resetParameterValues( initialParameterEstimate );
//
//    // A priori covariance
//    Eigen::MatrixXd invAprioriCovariance = Eigen::MatrixXd::Zero( nbParameters, nbParameters );
//
//    // Estimate initial states and tidal parameters
//    std::shared_ptr< EstimationInput< double, double > > estimationInput = std::make_shared< EstimationInput< double, double > >(
//            observationsAndTimes, invAprioriCovariance );
//    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 1 ) );
//    estimationInput->applyFinalParameterCorrection_ = true;
//    estimationInput->defineEstimationSettings( true, true, true, true, true, true );
//
////    std::map< observation_models::ObservableType, double > weightPerObservable;
////    weightPerObservable[ position_observable ] = 1.0 / ( positionNoise * positionNoise );
////    weightPerObservable[ velocity_observable ] = 1.0 / ( velocityNoise * velocityNoise );
////
////    estimationInput->setConstantPerObservableWeightsMatrix( weightPerObservable );
//
////    std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );
////
////    // Check if parameters are correctly estimated
////    Eigen::VectorXd estimatedParametervalues = estimationOutput->parameterHistory_.at( estimationOutput->parameterHistory_.size( ) - 1 );
////    Eigen::VectorXd residuals = estimationOutput->residuals_;
////
//////    input_output::writeMatrixToFile( residuals, "residuals" + indexFiles + ".dat", 16, getOutputPath( ) + outputFolder );
//////    input_output::writeMatrixToFile( estimationOutput->getNormalizedWeightedDesignMatrix( ), "normalisedWeightedDesignMatrix" + indexFiles + ".dat", 16, getOutputPath( ) + outputFolder );
//////    input_output::writeMatrixToFile( estimationOutput->getCorrelationMatrix( ), "correlations" + indexFiles + ".dat", 16, getOutputPath( ) + outputFolder );
////
////    std::cout << "-----------------------------------" << "\n\n";
////    std::cout << "truthParameters: " << truthParameters.transpose( ) << "\n\n";
////    std::cout << "initialParameterEstimate: " << initialParameterEstimate.transpose( ) << "\n\n";
////    std::cout << "estimatedParametervalues: " << estimatedParametervalues.transpose( ) << "\n\n";
////
////    Eigen::VectorXd parametersUpdate = estimatedParametervalues - initialParameterEstimate;
////    std::cout << "parametersUpdate: " << parametersUpdate.transpose( ) << "\n\n";
////
////    std::cout << "-----------------------------------" << "\n\n";
////    std::cout << "old k2: " << std::complex< double >( k2ReMoon, k2ImIo ) << "\n\n";
////    std::cout << "old invQ: " << invQIo << "\n\n";
////    std::cout << "old time lag: " << getTimeLagMoonFromLoveNumber( std::complex< double >( k2ReMoon, k2ImIo ), tidalPeriods.at( "Io" ) ) << "\n\n";
////
////    std::complex< double > newK2 = std::complex< double >( estimatedParametervalues[ nbParameters - 2 ], estimatedParametervalues[ nbParameters - 1 ] );
////    std::cout << "new k2: " << newK2 << "\n\n";
////    std::cout << "new tidal parameter: " << estimatedParametervalues[ nbParameters - 3 ] << "\n\n";
////
//////    std::cout << "old k2: " << std::complex< double >( k2ReMoon, k2ImIo ) << "\n\n";
//////    std::cout << "norm k2: " << std::sqrt( k2ReMoon * k2ReMoon + k2ImIo * k2ImIo ) << "\n\n";
//////    std::cout << "ratio InvQ / Im(k2): " << 1.0 / std::sqrt( k2ReMoon * k2ReMoon + k2ImIo * k2ImMoon ) << "\n\n";
////
////    double theoreticalInvQ = getInvQFromLoveNumber( newK2 );
////    double theoreticalTimeLag = getTimeLagMoonFromLoveNumber( newK2, tidalPeriods.at( "Io" ) );
////    std::cout << "theoreticalInvQ: " << theoreticalInvQ << "\n\n";
////    std::cout << "theoreticalTimeLag: " << theoreticalTimeLag << "\n\n";
////    std::cout << "estimated tidal parameter: " << estimatedParametervalues[ nbParameters - 3 ] << "\n\n";
////
//////    std::vector< std::shared_ptr< gravitation::GravityFieldVariations > > allVariations = bodies.at( "Io" )->getGravityFieldVariationSet( )->getVariationObjects( );
//////    std::cout << "size field variations: " << allVariations.size( ) << "\n\n";
//////    std::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > solidTidesVariations =
//////            std::dynamic_pointer_cast< gravitation::BasicSolidBodyTideGravityFieldVariations >( allVariations[ 0 ] );
//////    std::map< int, std::vector< std::complex< double > > > loveNumbers = solidTidesVariations->getLoveNumbers( );
//////    for ( auto itr : loveNumbers )
//////    {
//////        std::cout << itr.first << "\n\n";
//////        for ( unsigned int i = 0 ; i < itr.second.size( ) ; i++ )
//////        {
//////            std::cout << i << " - Re: " << itr.second[ i ].real( ) << " & Im: " << itr.second[ i ].imag( ) << "\n\n";
//////        }
//////    }
////
//////    std::cout << "inverse normalised covariance: " << "\n\n";
//////    std::cout << estimationOutput->getNormalizedInverseCovarianceMatrix( ) << "\n\n";
//////    std::cout << "normalisation terms: " << "\n\n";
//////    std::cout << estimationOutput->getNormalizationTerms( ).transpose( ) << "\n\n";
////
//////    input_output::writeMatrixToFile( estimationOutput->getCorrelationMatrix( ), "correlations.dat", 16,
//////                                     "/home/mfayolle/Documents/PHD/04 - year 4/tudat-bundle/tudat/applications/Dynamics/SimulationAnalysis/" );
//////
//////    input_output::writeMatrixToFile( estimationOutput->getNormalizedWeightedDesignMatrix( ), "weightedPartials.dat", 16,
//////                                     "/home/mfayolle/Documents/PHD/04 - year 4/tudat-bundle/tudat/applications/Dynamics/SimulationAnalysis/" );
//////
//////    std::cout << "formal errors: " << "\n\n";
//////    std::cout << estimationOutput->getFormalErrorVector( ).transpose( ) << "\n\n";
////
//////    std::cout << "size formal errors: " << estimationOutput->getFormalErrorVector( ).size( ) << "\n\n";
//////    std::cout << "global constraints size: " << parametersToEstimate->getConstraintSize( ) << "\n\n";
//////    std::cout << "size parameters to estimate: " << parametersToEstimate->template getFullParameterValues< double >( ).size( ) << "\n\n";
////
////
//    // Test custom constraints
//
//    std::function< Eigen::MatrixXd( Eigen::VectorXd ) > constraintFunctionK2 = [](Eigen::VectorXd parametersVector)
//    {
//        Eigen::MatrixXd constraint = Eigen::MatrixXd::Zero( 1, 2 );
//        constraint(0, 1) = 1.0;
//        return constraint;
//    };
//    Eigen::MatrixXd constraintConstantK2 = Eigen::MatrixXd::Zero( 1, 2 );
//    constraintConstantK2(0, 1) = 1.0;
//
//    std::function< Eigen::MatrixXd( Eigen::VectorXd ) > constraintFunctionInvQ = [](Eigen::VectorXd parametersVector)
//    {
//        Eigen::MatrixXd constraint = Eigen::MatrixXd::Zero( 1, 1 );
//        constraint(0, 0) = - 1.0;
//        return constraint;
//    };
//    Eigen::MatrixXd constraintConstantInvQ = Eigen::MatrixXd::Zero( 1, 1 );
//    constraintConstantInvQ(0, 0) = - 1.0;
//
//    std::vector< std::function< Eigen::MatrixXd( Eigen::VectorXd ) > > linearisedFunctionConstraint = { constraintFunctionK2, constraintFunctionInvQ };
//    std::vector< Eigen::MatrixXd > linearisedConstantConstraint = { constraintConstantK2, constraintConstantInvQ };
//
//    std::vector< estimatable_parameters::EstimatebleParametersEnum > parametersInvolved =
//            { single_degree_variable_tidal_love_number, inverse_tidal_quality_factor };
//    std::vector< std::pair< std::string, std::string > > bodiesInvolvedConstraint; // = { std::make_pair( "Io", "" ), std::make_pair( "Io", "" ) };
//
//    bool applyPerBody = true;
//    bool applyToAllBodies = true;
//    createCustomParameterConstraints( parametersToEstimate, linearisedFunctionConstraint, parametersInvolved, bodiesInvolvedConstraint,
//                                      applyPerBody, applyToAllBodies );
//
//    createConstantCustomParameterConstraints( parametersToEstimate, linearisedConstantConstraint, parametersInvolved,
//                                              bodiesInvolvedConstraint, applyPerBody, applyToAllBodies );
//
//    Eigen::MatrixXd constraintStateMultiplier;
//    Eigen::VectorXd constraintRightHandSide;
//    parametersToEstimate->getConstraints( constraintStateMultiplier, constraintRightHandSide );
//
//    std::cout << "constraintStateMultiplier: " << "\n\n";
//    std::cout << constraintStateMultiplier << "\n\n";
//
//    std::cout << "constraintRightHandSide: " << "\n\n";
//    std::cout << constraintRightHandSide << "\n\n";
//
//    std::cout << "parameters values: " << "\n\n";
//    std::cout << parametersToEstimate->getFullParameterValues< double >( ).transpose( ) << "\n\n";
//
//    printInterParameterConstraints( parametersToEstimate );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}