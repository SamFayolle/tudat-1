/* git    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <string>
#include <thread>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/ephemerides/keplerEphemeris.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/createSystemModel.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"

namespace tudat
{

namespace unit_tests
{

// Using declarations.
using namespace tudat;
using namespace tudat::estimatable_parameters;
using namespace tudat::orbit_determination;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;

BOOST_AUTO_TEST_SUITE( test_deformation_variational_equations )


template< typename TimeType = double, typename StateScalarType = double >
std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
           std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
executeSimulation(
        const double propagationDuration,
        const Eigen::VectorXd initialStateDifference = Eigen::VectorXd::Zero( 6+5 ),
        const Eigen::VectorXd parameterPerturbation = Eigen::VectorXd::Zero( 10 ),
        const bool propagateVariationalEquations = 1,
        const bool includeTranslationalDynamics = true,
        const bool includeRotationalDynamics = false,
        const bool useSynchronousRotation = false )
{
    int numberOfParametersToEstimate = 10;

    // Specify initial time
    double initialTime = 0.0;
    double finalTime = initialTime + propagationDuration; 

    std::string globalFrameOrigin = "Jupiter";
    std::string globalFrameOrientation = "J2000";

    // Define bodies in simulation
    std::vector< std::string > bodiesToCreate = { "Jupiter", "Io" }; 

    // Get body settings.
    BodyListSettings bodySettings = getDefaultBodySettings( 
        bodiesToCreate, initialTime - 86400.0, finalTime + 86400.0, globalFrameOrigin, globalFrameOrientation );

    // // Set Jupiter's ephemeris
    // bodySettings.at( "Jupiter" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ), "SSB", globalFrameOrientation ); 

    // Set Io's rotation model to synchronous
    if ( useSynchronousRotation )
    {
        bodySettings.get( "Io" )->rotationModelSettings = simulation_setup::synchronousRotationModelSettings( "Jupiter", "J2000", "IAU_Io" );    
    }
    
    // Set Io's ephemeris (TBC)
    double muIo = 5959924010272.5136719;
    double muJupiter = 126686534196012800.0;
    double muEffective = 126692494120023072.0;

    double orbitalPeriodIo = 2.0 * mathematical_constants::PI * std::sqrt( 4.2e8 * 4.2e8 * 4.2e8 / muIo );
    double rotationRateIo = std::sqrt( muEffective / ( 4.2e8 * 4.2e8 * 4.2e8 ) );

        
    Eigen::Vector6d initialKeplerianState = ( Eigen::Vector6d( ) << 4.2e8, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( );
//     bodySettings.at( "Io" )->ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
        //     ( Eigen::Vector6d( ) << 1.0 * 421.8E6, 1.0 * 0.004, 0.0, 0.0, 0.0, 0.0 ).finished( ),
        //     0.0,
        //     getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Io" ),
        //     "Jupiter",
        //     "J2000" );


    // Create bodies 
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    
    // Set Jupiter's rotation as constant
    double rightAscensionPole = ( 358.054324066462 - 90.0 ) * mathematical_constants::PI / 180.0;
    double declinationPole = ( 90.0 - 25.5034135739821 ) * mathematical_constants::PI / 180.0;
    double primeMeridian = ( 284.95 ) * mathematical_constants::PI / 180.0;
    double rotationRateJupiter = ( 870.536 * mathematical_constants::PI / 180.0 ) / 86400.0;
    bodies.at("Jupiter")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
                rightAscensionPole, declinationPole, primeMeridian, rotationRateJupiter, initialTime, globalFrameOrientation, "IAU_Jupiter" ) );


    // Set dummy model for Io's rotation (if not synchronous)
    if ( !useSynchronousRotation )
    {
        bodies.at("Io")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
                0.0, mathematical_constants::PI / 180.0, 0.0, rotationRateIo, initialTime, globalFrameOrientation, "IAU_Io" ) );
    }


    // Set Io's mean moment of inertia
    double scaledMeanMomentOfInertia = 0.37685;
    std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )->setScaledMeanMomentOfInertia( scaledMeanMomentOfInertia );


    // Set the bodies' state and rotation to initial time
    bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );
    bodies.at( "Io" )->setStateFromEphemeris<>( initialTime );
    bodies.at( "Jupiter" )->setStateFromEphemeris<>( initialTime );


    // Create integrator settings
    double timeStep = 100.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > > ( 
        initialTime, timeStep, rungeKutta87DormandPrince, 1.0e-5, 1000.0 ); 


    // Create propagator settings
    std::vector< std::string > bodiesToPropagate = { "Io" };
    std::vector< std::string > centralBodies = { "Jupiter" };


    // Create Maxwell deformation model
    double maxwellRelaxationTime = 179103.0;
    double globalRelaxationTime = 24688.0;
    double fluidLoveNumber = 1.5;
    double radiusIo = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )->getReferenceRadius( ); 

    std::vector< std::string > perturbingBody = {"Jupiter"};
    std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::make_shared< MaxwellDeformationSettings >( 
        maxwellRelaxationTime, globalRelaxationTime, fluidLoveNumber, 2, 2, perturbingBody, Eigen::Vector5d::Zero(), true, true );

    std::map< std::string, std::vector< std::shared_ptr< GravityDeformationSettings > > > gravityDeformationModelMap;   
    gravityDeformationModelMap[ "Io" ] = { maxwellDeformationSettings };

    basic_astrodynamics::GravityDeformationModelMap gravityDeformationModels = createGravityDeformationModelsMap(
        bodies, gravityDeformationModelMap ); 

    // Compute equilibrium coefficients to initialise gravity state   
    Eigen::Vector5d equilibriumCoefficients = Eigen::Vector5d::Zero();

    Eigen::Vector6d initialTranslationalState = orbital_element_conversions::convertKeplerianToCartesianElements( initialKeplerianState, muEffective );
    double distanceIo = initialTranslationalState.segment(0, 3).norm();
    double radiusRatioPowerThree = ( radiusIo / distanceIo ) * ( radiusIo / distanceIo ) * ( radiusIo / distanceIo );
    double gravitationalParametersRatio = muJupiter / muIo;

    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRotationState = getInitialRotationalStateOfBody( "Io", "J2000",  bodies, initialTime );
        initialRotationState[6] = rotationRateIo;
    
    equilibriumCoefficients[ 0 ] = fluidLoveNumber / 2.0 * gravitationalParametersRatio * radiusRatioPowerThree 
                * ( 3.0 * std::sin( 0.0 ) * std::sin( 0.0 ) - 1.0 )
                - fluidLoveNumber * rotationRateIo * rotationRateIo * radiusIo * radiusIo * radiusIo / ( 3.0 * muIo );
    equilibriumCoefficients[ 2 ] = fluidLoveNumber / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
                ( 1.0 - std::sin( 0.0 ) * std::sin( 0.0 ) ) * std::cos( 2.0 * 0.0 );
    equilibriumCoefficients[ 4 ] = fluidLoveNumber / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
                ( 1.0 - std::sin( 0.0 ) * std::sin( 0.0 ) ) * std::sin( 2.0 * 0.0 );
        
    equilibriumCoefficients[ 1 ] = - fluidLoveNumber * gravitationalParametersRatio * radiusRatioPowerThree 
                * ( - std::cos( 0.0 ) * std::sin( 0.0 ) ) * std::cos( 0.0 );
    equilibriumCoefficients[ 3 ] = - fluidLoveNumber * gravitationalParametersRatio * radiusRatioPowerThree 
                * ( - std::cos( 0.0 ) * std::sin( 0.0 ) ) * std::sin( 0.0 );

    Eigen::Vector5d initialGravityState = equilibriumCoefficients;    

    std::shared_ptr< GravityDeformationPropagatorSettings< > > gravityPropagatorSettings = 
        std::make_shared< GravityDeformationPropagatorSettings< > >( bodiesToPropagate, gravityDeformationModels, initialGravityState, integratorSettings,
        std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );    


    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > >  propagatorSettingsList;
    if ( includeTranslationalDynamics )
    {
        // Create translational dynamics model 
        SelectedAccelerationMap accelerationSettingsMap;
        // accelerationSettingsMap[ "Io" ][ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
        accelerationSettingsMap[ "Io" ][ "Jupiter" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );

        AccelerationMap accelerationsMap = createAccelerationModelsMap( bodies, accelerationSettingsMap, bodiesToPropagate, centralBodies );
        std::shared_ptr< TranslationalStatePropagatorSettings<  > > translationalPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >( 
        centralBodies, accelerationsMap, bodiesToPropagate, initialTranslationalState, initialTime, integratorSettings, 
        std::make_shared< PropagationTimeTerminationSettings >( finalTime ), cowell );

        propagatorSettingsList.push_back( translationalPropagatorSettings );
    }       
    if ( includeRotationalDynamics )
    {
        // Create rotational dynamics model
        SelectedTorqueMap torqueSettings;
        torqueSettings[ "Io" ][ "Jupiter" ].push_back( std::make_shared< TorqueSettings >( basic_astrodynamics::second_order_gravitational_torque ) );
        // torqueSettings[ "Io" ][ "Jupiter" ].push_back( std::make_shared<SphericalHarmonicTorqueSettings>(2,2) );
        basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap( bodies, torqueSettings, bodiesToPropagate );

        
        // std::cout << "initialRotationState " << initialRotationState << std::endl;

        std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings = std::make_shared< RotationalStatePropagatorSettings< double > >( 
            torqueModelMap, bodiesToPropagate, initialRotationState, initialTime, integratorSettings, std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

        propagatorSettingsList.push_back( rotationalPropagatorSettings );    
    } 


    // Create full propagator settings
    std::shared_ptr< SingleArcPropagatorProcessingSettings > outputSettings =
            std::make_shared< SingleArcPropagatorProcessingSettings >( );
    outputSettings->setIntegratedResult( false );

    propagatorSettingsList.push_back( gravityPropagatorSettings );
    std::shared_ptr< MultiTypePropagatorSettings< > > propagatorSettings = std::make_shared< MultiTypePropagatorSettings< > >(
            propagatorSettingsList, integratorSettings, initialTime, std::make_shared< PropagationTimeTerminationSettings >( finalTime, true ),
            std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ), outputSettings );

    // Get Io's and Jupiter's state
    Eigen::Vector6d stateIo = bodies.at("Io")->getStateInBaseFrameFromEphemeris( initialTime );
    Eigen::Vector6d stateJupiter = bodies.at("Jupiter")->getStateInBaseFrameFromEphemeris( initialTime );
    Eigen::Vector6d inertialState = stateJupiter - stateIo;  
    // std::cout << "inertialState " << inertialState.transpose() << std::endl; 

    // Reset (perturbed) initial state
    Eigen::VectorXd fullInitialState = Eigen::VectorXd::Zero( 6 + 5 ); 
    fullInitialState.segment( 0, 6 ) = initialTranslationalState;  
    fullInitialState.segment( 6, 5 ) = initialGravityState; 
    fullInitialState += initialStateDifference;
    propagatorSettings->resetInitialStates( fullInitialState ); 
    // std::cout << "fullInitialState " << fullInitialState.transpose() << std::endl; 
    // std::cout << "initialStateDifference " << initialStateDifference.transpose() << std::endl;         

    // Define parameters to estimate
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    {
        parameterNames = getInitialStateParameterSettings< StateScalarType >( propagatorSettings, bodies );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Io", gravitational_parameter ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Jupiter", gravitational_parameter ) );
        // parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Io", constant_rotation_rate ) );

        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
            1, 0, 2, 2, "Io", spherical_harmonics_cosine_coefficient_block ) );
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
            1, 1, 2, 2, "Io", spherical_harmonics_sine_coefficient_block ) );
        // parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", radiation_pressure_coefficient ) );
        // parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", constant_drag_coefficient ) );
        // parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Moon", gravitational_parameter ) );

        // parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
        //         3, 0, 3, 3, "Earth", spherical_harmonics_cosine_coefficient_block ) );
        // parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
        //         3, 1, 3, 3, "Earth", spherical_harmonics_sine_coefficient_block ) );
    }

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodies );
    std::cout << "size parameters set " << parametersToEstimate->getParameterSetSize() << std::endl;

    // Perturb parameters.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > parameterVector =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );
    parameterVector.block( 6+5, 0, numberOfParametersToEstimate, 1 ) += parameterPerturbation;
    parametersToEstimate->resetParameterValues( parameterVector );

    // printEstimatableParameterEntries( parametersToEstimate );

    std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
               std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > results;

    {
        // Create dynamics simulator
        SingleArcVariationalEquationsSolver< StateScalarType, TimeType > dynamicsSimulator =
                SingleArcVariationalEquationsSolver< StateScalarType, TimeType >(
                        bodies,
                        integratorSettings,
                        propagatorSettings,
                        parametersToEstimate,
                        1,
                        std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
                        0,
                        0 );

        // Propagate requested equations.
        if( propagateVariationalEquations )
        {
            dynamicsSimulator.integrateVariationalAndDynamicalEquations( propagatorSettings->getInitialStates( ), 1 );
        }
        else
        {
            dynamicsSimulator.integrateDynamicalEquationsOfMotionOnly( propagatorSettings->getInitialStates( ) );
        }

        // Retrieve test data
        double testEpoch = finalTime;
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > testStates = Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( 6+5 );
        testStates.block( 0, 0, 6, 1 ) = bodies.at( "Io" )->getEphemeris( )->getCartesianState( testEpoch );
        
        std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator.getEquationsOfMotionSolution( );
        testStates = stateHistory.rbegin()->second;

        if( propagateVariationalEquations )
        {
            results.first.push_back(
                    dynamicsSimulator.getStateTransitionMatrixInterface( )->getCombinedStateTransitionAndSensitivityMatrix( testEpoch ) );
            Eigen::MatrixXd testMatrixDirect =
                    dynamicsSimulator.getStateTransitionMatrixInterface( )->getCombinedStateTransitionAndSensitivityMatrix( testEpoch );
            Eigen::MatrixXd testMatrixFull =
                    dynamicsSimulator.getStateTransitionMatrixInterface( )->getFullCombinedStateTransitionAndSensitivityMatrix( testEpoch );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testMatrixDirect, testMatrixFull, std::numeric_limits< double >::epsilon( ) );
        }
        results.second.push_back( testStates );
    }
    return results;
}


// BOOST_AUTO_TEST_CASE( testDeformationAndTranslationalVariationalEquations )
// {
//     double propagationDuration = 10.0 * 3600.0; // 20.0 * physical_constants::JULIAN_DAY;

//     // Load spice kernels.
//     spice_interface::loadStandardSpiceKernels( );

//     std::pair< std::vector< Eigen::MatrixXd >, std::vector< Eigen::VectorXd > > currentOutput;

//     // Define variables for numerical differentiation
//     Eigen::Matrix< double, 11, 1 > perturbedState;
//     Eigen::Matrix< double, 11, 1 > statePerturbation;

//     // Define parameter perturbation
//     int numberOfParametersToEstimate = 10;
//     // Eigen::VectorXd perturbedParameter = Eigen::VectorXd::Zero( numberOfParametersToEstimate );


//     double sphericalHarmonicsPerturbation = 1.0E-6;
//     Eigen::Matrix< double, 10, 1 > perturbedParameter;
//     Eigen::Matrix< double, 10, 1 > parameterPerturbation = ( Eigen::Matrix< double, 10, 1 >( ) <<
//                                                              1.0E10,
//                                                              1.0E10,
//                                                              sphericalHarmonicsPerturbation,
//                                                              sphericalHarmonicsPerturbation,
//                                                              sphericalHarmonicsPerturbation,
//                                                              sphericalHarmonicsPerturbation,
//                                                              sphericalHarmonicsPerturbation,
//                                                              sphericalHarmonicsPerturbation,
//                                                              sphericalHarmonicsPerturbation,
//                                                              sphericalHarmonicsPerturbation ).finished( );

//     // Define state perturbation
//     statePerturbation = ( Eigen::Matrix< double, 11, 1 >( ) 
//         << 100.0, 100.0, 100.0, 0.1, 1.0e-5, 1.0e-5, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6 ).finished( );

//     // Compute state transition and sensitivity matrices
//     currentOutput = executeSimulation< double, double >( propagationDuration, Eigen::Matrix< double, 11, 1 >::Zero( ) );
//     Eigen::MatrixXd stateTransitionAndSensitivityMatrixAtEpoch = currentOutput.first.at( 0 );
//     Eigen::MatrixXd manualPartial = Eigen::MatrixXd::Zero( 11, 11 + numberOfParametersToEstimate );

//     // Numerically compute state transition matrix
//     for( unsigned int j = 0; j < 11; j++ )
//     {
//         std::cout << j << " - perturbation " << statePerturbation( j ) << std::endl;
//         Eigen::VectorXd upPerturbedState, downPerturbedState;
//         perturbedState.setZero( );
//         perturbedState( j ) += statePerturbation( j );
//         upPerturbedState =
//                 executeSimulation< double, double >( propagationDuration, perturbedState, Eigen::Matrix< double, 10, 1 >::Zero( ), 0 ).second.at( 0 );

//         perturbedState.setZero( );
//         perturbedState( j ) -= statePerturbation( j );
//         downPerturbedState =
//                 executeSimulation< double, double >( propagationDuration, perturbedState, Eigen::Matrix< double, 10, 1 >::Zero( ), 0 ).second.at( 0 );

//         manualPartial.block( 0, j, 11, 1 ) =
//                 ( upPerturbedState.segment( 0, 11 ) - downPerturbedState.segment( 0, 11 ) ) / ( 2.0 * statePerturbation( j ) );
//     }

//     // Numerically compute sensitivity matrix
//     for( int j = 0; j < numberOfParametersToEstimate; j++ )
//     {
//         Eigen::VectorXd upPerturbedState, downPerturbedState;
//         perturbedState.setZero( );
//         Eigen::Matrix< double, 10, 1 > upPerturbedParameter, downPerturbedParameter;
//         perturbedParameter.setZero( );
//         perturbedParameter( j ) += parameterPerturbation( j );

//         upPerturbedState = executeSimulation< double, double >( propagationDuration, perturbedState, perturbedParameter ).second.at( 0 );

//         perturbedParameter.setZero( );
//         perturbedParameter( j ) -= parameterPerturbation( j );
//         downPerturbedState = executeSimulation< double, double >( propagationDuration, perturbedState, perturbedParameter ).second.at( 0 );

//         manualPartial.block( 0, j + 11, 11, 1 ) =
//                 ( upPerturbedState.segment( 0, 11 ) - downPerturbedState.segment( 0, 11 ) ) / ( 2.0 * parameterPerturbation( j ) );
//     }

//     std::cout << "manualPartial" << std::endl;
//     std::cout << manualPartial.block(0, 0, 11, 11) << std::endl;

//     std::cout << "stateTransitionAndSensitivityMatrixAtEpoch" << std::endl;
//     std::cout << stateTransitionAndSensitivityMatrixAtEpoch.block(0, 0, 11, 11) << std::endl;

    // Eigen::MatrixXd differences = manualPartial.block(0, 0, 11, 11) - stateTransitionAndSensitivityMatrixAtEpoch.block(0, 0, 11, 11);
    // for ( unsigned int j = 0 ; j < 11 ; j++ )
    // {
    //     for ( unsigned int k = 0 ; k < 11 ; k++ )
    //     {
    //         if ( stateTransitionAndSensitivityMatrixAtEpoch( j,k ) != 0.0 )
    //         {
    //             differences( j,k ) /= stateTransitionAndSensitivityMatrixAtEpoch( j,k );
    //         }
    //     }
    // }
    // std::cout << "relative differences" << std::endl;
    // std::cout << differences << std::endl;

//     // std::cout << "manualPartial" << std::endl;
//     // std::cout << manualPartial.block(0, 11, 11, 10) << std::endl;

//     // std::cout << "stateTransitionAndSensitivityMatrixAtEpoch" << std::endl;
//     // std::cout << stateTransitionAndSensitivityMatrixAtEpoch.block(0, 11, 11, 10) << std::endl;

//     // // Check results
//     // TUDAT_CHECK_MATRIX_CLOSE_FRACTION( stateTransitionAndSensitivityMatrixAtEpoch, manualPartial, 1.0e-6 );
// }

template< typename TimeType = double, typename StateScalarType = double >
std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
           std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >, std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >
executeRotationSimulation( const Eigen::Matrix< StateScalarType, 12, 1 > initialStateDifference,
                           Eigen::Matrix< StateScalarType, 12, 1 >& appliedStateDifference,
                           const Eigen::VectorXd parameterPerturbation = Eigen::VectorXd::Zero( 8 ),
                           const bool propagateVariationalEquations = 1 )
{
    Eigen::VectorXd manualInitialState = Eigen::VectorXd::Zero(12);
    manualInitialState << 0.37863607626213874324, -0.13679006662899653723, 0.79433354992246751891, -0.45492572019178501019,
       2.104497477798486668e-07, 3.6973061209642486676e-07, 4.1201088380643943146e-05, -0.00018593809808886761596,
        -4.1627693959500079984e-05, 0.00026945354354527141071, -0.0002500155717872719059, 0.0001783605477865248467;

    double initialTime = 0.0;
    double finalTime = 20.0 * 86400.0; // 10.0 * 3600.0; //86400.0;
    int numberOfParametersToEstimate = 8;

    std::string globalFrameOrigin = "Jupiter";
    std::string globalFrameOrientation = "J2000";

    std::vector< std::string > bodiesToCreate = { "Jupiter", "Io" }; 

    double muJupiter = getBodyGravitationalParameter( "Jupiter" );
    double muIo = getBodyGravitationalParameter( "Io" );
    double muEffective = muJupiter + muIo;

    // Get body settings.
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodiesToCreate, initialTime - 86400.0, finalTime + 86400.0, globalFrameOrigin, globalFrameOrientation );

    // Set Io's rotation to synchronous        
    bodySettings.get( "Io" )->rotationModelSettings = simulation_setup::synchronousRotationModelSettings( "Jupiter", "J2000", "IAU_Io" );

    // Create system of bodies
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // SystemOfBodies bodies = SystemOfBodies( "Jupiter", "ECLIPJ2000" );
    // bodies.createEmptyBody( "Jupiter", false );
    // bodies.at( "Jupiter" )->setEphemeris(
            // std::make_shared< ephemerides::ConstantEphemeris >( [ = ]( ) { return Eigen::Vector6d::Zero( ); }, "SSB", "J2000" ) );
    // bodies.at( "Jupiter" )->setRotationalEphemeris( simulation_setup::createRotationModel(
    //         simulation_setup::getDefaultRotationModelSettings( "Jupiter", initialTime, finalTime ), "Jupiter" ) );
    // bodies.at( "Jupiter" )->setGravityFieldModel( simulation_setup::createGravityFieldModel(
    //         simulation_setup::getDefaultGravityFieldSettings( "Jupiter", initialTime, finalTime ), "Jupiter", bodies ) );

    // bodies.createEmptyBody( "Io" );

//     Eigen::Matrix3d phobosInertiaTensor = Eigen::Matrix3d::Zero( );
//     phobosInertiaTensor( 0, 0 ) = 0.3615;
//     phobosInertiaTensor( 1, 1 ) = 0.4265;
//     phobosInertiaTensor( 2, 2 ) = 0.5024;

//     phobosInertiaTensor *= ( 11.27E3 * 11.27E3 * 1.0659E16 );

//     double phobosGravitationalParameter = 1.0659E16 * physical_constants::GRAVITATIONAL_CONSTANT;
//     double phobosReferenceRadius = 11.27E3;

//     Eigen::MatrixXd phobosCosineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 ),
//                     phobosSineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 );
//     double phobosScaledMeanMomentOfInertia;
//     gravitation::getDegreeTwoSphericalHarmonicCoefficients( phobosInertiaTensor,
//                                                             phobosGravitationalParameter,
//                                                             phobosReferenceRadius,
//                                                             true,
//                                                             phobosCosineGravityFieldCoefficients,
//                                                             phobosSineGravityFieldCoefficients,
//                                                             phobosScaledMeanMomentOfInertia );

    // Set Io's mean moment of inertia
    std::shared_ptr< SphericalHarmonicsGravityField > ioGravityFieldModel = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) );
    double scaledMeanMomentOfInertia = 0.37685;
    ioGravityFieldModel->setScaledMeanMomentOfInertia( scaledMeanMomentOfInertia );

    Eigen::MatrixXd ioCosineCoefs = ioGravityFieldModel->getCosineCoefficients( );
    Eigen::MatrixXd ioSineCoefs = ioGravityFieldModel->getSineCoefficients( );
    ioCosineCoefs(2,1) = 1.0e-3;
    ioSineCoefs(2,1) = 1.0e-3;
    ioSineCoefs(2,2) = 2.0e-3;
    // ioGravityFieldModel->setCosineCoefficients( ioCosineCoefs );
    // ioGravityFieldModel->setSineCoefficients( ioSineCoefs );
    // std::cout << "ioCosineCoefs" << std::endl;
    // std::cout << ioCosineCoefs << std::endl;
    // std::cout << "ioSineCoefs" << std::endl;
    // std::cout << ioSineCoefs << std::endl;

    // std::shared_ptr< SphericalHarmonicsGravityField > newIoGravityFieldModel = std::make_shared< SphericalHarmonicsGravityField >( 
        // ioGravityFieldModel->getGravitationalParameter( ), ioGravityFieldModel->getReferenceRadius( ), ioCosineCoefs, ioSineCoefs, "IAU_Io", scaledMeanMomentOfInertia );
    // bodies.at( "Io" )->setGravityFieldModel( newIoGravityFieldModel );

    // Retrieve Io's static gravity coefficients
    ioCosineCoefs = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )->getCosineCoefficients( );
    ioSineCoefs = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )->getSineCoefficients( );
    Eigen::Vector5d normalisedStaticGravity = Eigen::Vector5d::Zero( );
    normalisedStaticGravity[0] = ioCosineCoefs(2,0);
    // normalisedStaticGravity[1] = 1.0e-4; // ioCosineCoefs(2,1);
    normalisedStaticGravity[2] = ioCosineCoefs(2,2);
    // normalisedStaticGravity[3] = 1.0e-4; // ioSineCoefs(2,1);
    // normalisedStaticGravity[4] = 1.0e-4; // ioSineCoefs(2,2);
    std::cout << "normalisedStaticGravity " << normalisedStaticGravity.transpose() << std::endl;
    bodies.at( "Io" )->setStaticDegreeTwoCoefficients( normalisedStaticGravity );

    Eigen::Vector6d ioKeplerElements = Eigen::Vector6d::Zero( );
    double ioSemiMajorAxis = 4.2e8;
    ioKeplerElements( 0 ) = ioSemiMajorAxis;
    ioKeplerElements( 1 ) = 0.004;
    ioKeplerElements( 2 ) = 0.002 * mathematical_constants::PI / 180.0;

    // bodies.at( "Io" )
    //         ->setEphemeris( tudat::ephemerides::getTabulatedEphemeris(
    //                 std::make_shared< ephemerides::KeplerEphemeris >( ioKeplerElements, 0.0, muEffective, "Jupiter", "J2000" ),
    //                 initialTime - 3600.0,
    //                 finalTime + 3600,
    //                 60.0 ) );

    // // Set Io's parameters - TO BE MODIFIED
    // double muIo = 5959924010272.5136719;
    // double muJupiter = 126686534196012800.0;
    // double muEffective = 126692494120023072.0;


    double rotationRateIo = std::sqrt( muEffective / ( ioSemiMajorAxis * ioSemiMajorAxis * ioSemiMajorAxis ) );

    // Retrieve initial rotational state from ephemeris
    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRotationState = getInitialRotationalStateOfBody( "Io", "J2000",  bodies, initialTime );
    initialRotationState[6] = rotationRateIo;

    // ///// ALTERNATIVE DEFINITION OF IO'S ROTATIONAL STATE
    // Eigen::Quaterniond noRotationQuaternion = Eigen::Quaterniond( Eigen::AngleAxisd( 1.0E-0, Eigen::Vector3d::UnitZ( ) ) *
    //                                                               Eigen::AngleAxisd( 2.0E-0, Eigen::Vector3d::UnitX( ) ) *
    //                                                               Eigen::AngleAxisd( -0.5E-0, Eigen::Vector3d::UnitZ( ) ) );

    // initialRotationState( 0 ) = noRotationQuaternion.w( );
    // initialRotationState( 1 ) = noRotationQuaternion.x( );
    // initialRotationState( 2 ) = noRotationQuaternion.y( );
    // initialRotationState( 3 ) = noRotationQuaternion.z( );
    // initialRotationState( 4 ) = 1.0e-5;
    // initialRotationState( 5 ) = -1.0e-5;
    // initialRotationState( 6 ) = rotationRateIo;
    // ////

    // // initialRotationState = manualInitialState.segment( 0, 7 );

    Eigen::Matrix< double, 7, 1 > originalRotationState = initialRotationState;
    Eigen::Matrix< double, 7, 1 > stateDifferenceToAdd = initialStateDifference.segment( 0, 7 );
    std::cout << "stateDifferenceToAdd " << stateDifferenceToAdd.transpose() << std::endl;

    initialRotationState += stateDifferenceToAdd;
    initialRotationState( 0 ) = originalRotationState( 0 ) / std::fabs( originalRotationState( 0 ) ) *
            std::sqrt( 1.0 - std::pow( initialRotationState.segment( 1, 3 ).norm( ), 2.0 ) );

    appliedStateDifference.segment( 0, 7 ) = initialRotationState - originalRotationState;
    appliedStateDifference.segment( 7, 5 ) = initialStateDifference.segment( 7, 5 );
    std::cout << "appliedStateDifference" << std::endl;
    std::cout << appliedStateDifference.transpose() << std::endl;
    
    
    std::map< double, Eigen::Matrix< double, 7, 1 > > dummyRotationMap;
    dummyRotationMap[ -1.0E100 ] = initialRotationState;
    dummyRotationMap[ 1.0E100 ] = initialRotationState;

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > dummyInterpolator =
            std::make_shared< interpolators::LinearInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( dummyRotationMap );
    bodies.at( "Io" )->setRotationalEphemeris(
        std::make_shared< TabulatedRotationalEphemeris< double, double > >( dummyInterpolator, "J2000", "IAU_Io" ) );

    // Retrieve body objects for Io and Jupiter
    std::shared_ptr< Body > io = bodies.at( "Io" );
    std::shared_ptr< Body > jupiter = bodies.at( "Jupiter" );

    // Update Jupiter and Io to current state
    io->setStateFromEphemeris( initialTime );
    jupiter->setStateFromEphemeris( initialTime );
    io->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );

    // Check translational intial state
    Eigen::Vector6d initialTranslationalState = io->getState( ) - jupiter->getState( );
    std::cout << "initialTranslationalState " << initialTranslationalState.transpose( ) << std::endl;

    std::cout << "initialRotationState " << initialRotationState.transpose( ) << std::endl;    

//     SelectedAccelerationMap accelerationMap;
// //     std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsIo;
//     //    accelerationMap[ "Phobos" ][ "Mars" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
//     accelerationMap[ "Io" ][ "Jupiter" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 2, 2 ) );

//     std::vector< std::string > translationalBodiesToIntegrate;
//     std::vector< std::string > translationalCentralBodies;

//     translationalBodiesToIntegrate.push_back( "Phobos" );
//     translationalCentralBodies.push_back( "Mars" );

//     AccelerationMap accelerationModelMap =
//             createAccelerationModelsMap( bodies, accelerationMap, translationalBodiesToIntegrate, translationalCentralBodies );

    SelectedTorqueMap torqueMap;
    torqueMap[ "Io" ][ "Jupiter" ].push_back( std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );


    // Create integrator settings
    double timeStep = 100.0 / 2.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > > ( 
        initialTime, timeStep, rungeKutta87DormandPrince, 1.0e-5, 1000.0, 1.0e-12, 1.0e-12 ); 

    // Define propagator settings.
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Io" );

    // Create torque models
    basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap( bodies, torqueMap, bodiesToIntegrate );

    std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
            std::make_shared< RotationalStatePropagatorSettings< double > >(
                    torqueModelMap,
                    bodiesToIntegrate,
                    initialRotationState,
                    initialTime,
                    integratorSettings,
                    std::make_shared< PropagationTimeTerminationSettings >( finalTime, true ) );

    // Create gravity deformation model
    double maxwellRelaxationTime = 179103.0;
    double globalRelaxationTime = 24688.0;
    double fluidLoveNumber = 1.5;
    double radiusIo = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )->getReferenceRadius( ); 

    std::vector< std::string > perturbingBody = {"Jupiter"};
    std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::make_shared< MaxwellDeformationSettings >( 
        maxwellRelaxationTime, globalRelaxationTime, fluidLoveNumber, 2, 2, perturbingBody, normalisedStaticGravity, true, true );

    std::map< std::string, std::vector< std::shared_ptr< GravityDeformationSettings > > > gravityDeformationModelMap;   
    gravityDeformationModelMap[ "Io" ] = { maxwellDeformationSettings };

    basic_astrodynamics::GravityDeformationModelMap gravityDeformationModels = createGravityDeformationModelsMap(
        bodies, gravityDeformationModelMap ); 

    // Compute equilibrium coefficients to initialise gravity state   
    Eigen::Vector5d equilibriumCoefficients = Eigen::Vector5d::Zero();

    // Eigen::Vector6d initialTranslationalState = orbital_element_conversions::convertKeplerianToCartesianElements( ioKeplerElements, muEffective );
    double distanceIo = initialTranslationalState.segment(0, 3).norm();
    double radiusRatioPowerThree = ( radiusIo / distanceIo ) * ( radiusIo / distanceIo ) * ( radiusIo / distanceIo );
    double gravitationalParametersRatio = muJupiter / muIo;

    // std::cout << "translational " << initialTranslationalState.transpose() << std::endl;
    std::cout << "rotational " << initialRotationState.transpose() << std::endl;
    
    equilibriumCoefficients[ 0 ] = fluidLoveNumber / 2.0 * gravitationalParametersRatio * radiusRatioPowerThree 
                * ( 3.0 * std::sin( 0.0 ) * std::sin( 0.0 ) - 1.0 );
                // - fluidLoveNumber * rotationRateIo * rotationRateIo * radiusIo * radiusIo * radiusIo / ( 3.0 * muIo );
    equilibriumCoefficients[ 2 ] = fluidLoveNumber / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
                ( 1.0 - std::sin( 0.0 ) * std::sin( 0.0 ) ) * std::cos( 2.0 * 0.0 );
    equilibriumCoefficients[ 4 ] = fluidLoveNumber / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
                ( 1.0 - std::sin( 0.0 ) * std::sin( 0.0 ) ) * std::sin( 2.0 * 0.0 );
        
    equilibriumCoefficients[ 1 ] = - fluidLoveNumber * gravitationalParametersRatio * radiusRatioPowerThree 
                * ( - std::cos( 0.0 ) * std::sin( 0.0 ) ) * std::cos( 0.0 );
    equilibriumCoefficients[ 3 ] = - fluidLoveNumber * gravitationalParametersRatio * radiusRatioPowerThree 
                * ( - std::cos( 0.0 ) * std::sin( 0.0 ) ) * std::sin( 0.0 );

    std::cout << "equilibriumCoefficients " << equilibriumCoefficients.transpose() << std::endl;

    Eigen::Vector5d initialGravityState = equilibriumCoefficients + initialStateDifference.segment(7,5); 
    // initialGravityState = manualInitialState.segment( 7, 5 ) + initialStateDifference.segment(7,5);    
    std::cout << "initialGravityState " << initialGravityState.transpose() << std::endl;    

    std::shared_ptr< GravityDeformationPropagatorSettings< > > gravityPropagatorSettings = 
        std::make_shared< GravityDeformationPropagatorSettings< > >( bodiesToIntegrate, gravityDeformationModels, initialGravityState, integratorSettings,
        std::make_shared< PropagationTimeTerminationSettings >( finalTime, true ) );  
    
    // Eigen::VectorXd initialTranslationalState;
    // initialTranslationalState =
    //         getInitialStatesOfBodies( translationalBodiesToIntegrate, translationalCentralBodies, bodies, initialTime );

//     initialTranslationalState += initialStateDifference.segment( 0, 6 );
//     std::shared_ptr< TranslationalStatePropagatorSettings<> > translationalPropagatorSettings =
//             std::make_shared< TranslationalStatePropagatorSettings<> >( translationalCentralBodies,
//                                                                         accelerationModelMap,
//                                                                         translationalBodiesToIntegrate,
//                                                                         initialTranslationalState,
//                                                                         finalEphemerisTime,
//                                                                         cowell );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
    propagatorSettingsList.push_back( rotationalPropagatorSettings );
    propagatorSettingsList.push_back( gravityPropagatorSettings );

    // Define dependent variables
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back(
            std::make_shared< SingleDependentVariableSaveSettings >( body_fixed_relative_spherical_position, "Jupiter", "Io" ) );

    std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings = std::make_shared< MultiTypePropagatorSettings< double > >(
            propagatorSettingsList, std::make_shared< PropagationTimeTerminationSettings >( finalTime, true ), dependentVariablesList );   
    

    // Define parameters.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    {
        parameterNames = getInitialStateParameterSettings< StateScalarType >( propagatorSettings, bodies );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Io", mean_moment_of_inertia ) );
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                1, 0, 2, 2, "Io", spherical_harmonics_cosine_coefficient_block ) );
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 1, 2, 2, "Io", spherical_harmonics_sine_coefficient_block ) );
    }

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodies );
    // printEstimatableParameterEntries( parametersToEstimate );

    Eigen::MatrixXd constraintStateMultiplier;
    Eigen::VectorXd constraintRightHandSide;
    parametersToEstimate->getConstraints( constraintStateMultiplier, constraintRightHandSide );
    //    std::cout<<"Unit rotation: "<<std::endl<<unitRotationState.transpose( )<<std::endl;
    //    std::cout<<"Constraints: "<<std::endl<<constraintStateMultiplier.transpose( )<<std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( constraintStateMultiplier.block( 0, 0, 1, 4 ) ),
                                       ( initialRotationState.segment( 0, 4 ) ).transpose( ),
                                       std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( constraintStateMultiplier.block( 0, 4, 1, 9 ) ),
                                       ( Eigen::MatrixXd::Zero( 1, 9 ) ),
                                       std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            ( constraintRightHandSide.block( 0, 0, 1, 1 ) ), ( Eigen::MatrixXd::Zero( 1, 1 ) ), std::numeric_limits< double >::epsilon( ) );

    // Perturb parameters.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > parameterVector =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );
    parameterVector.block( 12, 0, numberOfParametersToEstimate, 1 ) += parameterPerturbation;
    //    std::cout<<"Parameter perturbation "<<
    //               ( parametersToEstimate->template getFullParameterValues< StateScalarType >( ) -
    //                 parameterVector ).transpose( )<<std::endl<<
    //               parameterPerturbation.transpose( )<<std::endl;;
    parametersToEstimate->resetParameterValues( parameterVector );

    std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
               std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >, 
                        std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > > results;

    {
        // Create dynamics simulator
        propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType > dynamicsSimulator =
                propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType >(
                        bodies,
                        integratorSettings,
                        propagatorSettings,
                        parametersToEstimate,
                        1,
                        std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
                        0,
                        0 );

        // Propagate requested equations.
        if( propagateVariationalEquations )
        {
            dynamicsSimulator.integrateVariationalAndDynamicalEquations( propagatorSettings->getInitialStates( ), 1 );
        }
        else
        {
            dynamicsSimulator.integrateDynamicalEquationsOfMotionOnly( propagatorSettings->getInitialStates( ) );
        }

            //    tudat::input_output::writeDataMapToTextFile(
            //                dynamicsSimulator.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( ),
            //                "rotPropTest.dat" );

        // Retrieve test data
        double testEpoch = finalTime; // initialTime + 1.0 * 3600.0;
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > testStates = Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( 12 );
        testStates.segment( 0, 7 ) = dynamicsSimulator.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( ).rbegin()->second.segment( 0, 7 ); //bodies.at( "Io" )->getRotationalEphemeris( )->getRotationStateVector( testEpoch );
        testStates.segment( 7, 5 ) = dynamicsSimulator.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( ).rbegin()->second.segment( 7, 5 );
        std::cout << "final state " << testStates.transpose() << std::endl;\

        if( propagateVariationalEquations )
        {
            // results.first.push_back(
                    // dynamicsSimulator.getStateTransitionMatrixInterface( )->getCombinedStateTransitionAndSensitivityMatrix( testEpoch ) );
            results.first = dynamicsSimulator.getStateTransitionMatrixSolution( );
            Eigen::MatrixXd testMatrixDirect =
                    dynamicsSimulator.getStateTransitionMatrixInterface( )->getCombinedStateTransitionAndSensitivityMatrix( testEpoch );
            Eigen::MatrixXd testMatrixFull =
                    dynamicsSimulator.getStateTransitionMatrixInterface( )->getFullCombinedStateTransitionAndSensitivityMatrix( testEpoch );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testMatrixDirect, testMatrixFull, std::numeric_limits< double >::epsilon( ) );
        }
        // results.second.push_back( testStates );
        results.second.first = dynamicsSimulator.getEquationsOfMotionSolution( );

        // Retrieve dependent variables history
        std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > simulator = dynamicsSimulator.getDynamicsSimulator( );
        std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > dependentVariablesHistory = simulator->getDependentVariableHistory();
        results.second.second = dependentVariablesHistory;

    }
    return results;
}

BOOST_AUTO_TEST_CASE( testDeformationAndRotationalVariationalEquations )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    std::pair< std::map< double, Eigen::MatrixXd >, std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > > currentOutput;

    // Define variables for numerical differentiation
    Eigen::Matrix< double, 12, 1 > perturbedState;
    Eigen::Matrix< double, 12, 1 > statePerturbation;

    // Define parameter perturbation
    int numberOfParametersToEstimate = 8;
    double sphericalHarmonicsPerturbation = 1.0E-4;
    Eigen::Matrix< double, 8, 1 > perturbedParameter;
    Eigen::Matrix< double, 8, 1 > parameterPerturbation;
    parameterPerturbation = Eigen::Matrix< double, 8, 1 >::Constant( sphericalHarmonicsPerturbation );
    parameterPerturbation( 0 ) = 1.0E-4;

    // Compute state transition and sensitivity matrices
    Eigen::Matrix< double, 12, 1 > appliedStateDifference;
    currentOutput = executeRotationSimulation< double, double >( Eigen::Matrix< double, 12, 1 >::Zero( ), appliedStateDifference );
    Eigen::MatrixXd stateTransitionAndSensitivityMatrixAtEpoch = currentOutput.first.rbegin()->second; // currentOutput.first.at( 0 );
    Eigen::VectorXd nominalState = currentOutput.second.first.rbegin()->second; // currentOutput.second.at( 0 );

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > nominalStateInterpolator =
                std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( currentOutput.second.first, 6 );

    std::vector< double > refEpochs;
    for ( auto it : currentOutput.first )
    {
        refEpochs.push_back( it.first );
    }

    std::map< double, Eigen::MatrixXd > analyticalStmHistory = currentOutput.first;
    std::map< double, Eigen::MatrixXd > modifiedAnalyticalStmHistory = analyticalStmHistory;


    //    std::cout<<"Nominal "<<std::endl<<std::endl<<
    //               stateTransitionAndSensitivityMatrixAtEpoch<<std::endl;
    // Define state perturbation
    statePerturbation = ( Eigen::Matrix< double, 12, 1 >( ) << 
                          1.0E-8,
                          1.0E-8,
                          1.0E-8,
                          1.0E-8,
                          1.0E-10,
                          1.0E-10,
                          1.0E-9, 
                          1.0E-6, 
                          1.0E-6,
                          1.0E-6,
                          1.0E-6,
                          1.0E-6 )
                                .finished( );

    Eigen::MatrixXd manualPartial = Eigen::MatrixXd::Zero( 12, 12 + numberOfParametersToEstimate );

    // Numerically compute state transition matrix
    std::map< double, Eigen::MatrixXd > numericalPartialsHistory, interpolatedUpPerturbedStateHistory, interpolatedDownPerturbedStateHistory;
    for ( unsigned int k = 0 ; k < refEpochs.size() ; k++ )
    {
        numericalPartialsHistory[ refEpochs[k] ] = Eigen::MatrixXd::Zero( 12, 12 );
        interpolatedUpPerturbedStateHistory[ refEpochs[k] ] = Eigen::MatrixXd::Zero( 12, 12 );
        interpolatedDownPerturbedStateHistory[ refEpochs[k] ] = Eigen::MatrixXd::Zero( 12, 12 );
    }

    for( unsigned int test = 0; test < 1; test++ )
    {
        double perturbationMultiplier = ( test == 0 ? 1.0 : 1.0E-3 );
        for( unsigned int j = 0; j < 4; j++ )
        {
            Eigen::Matrix< double, 12, 1 > appliedStateDifferenceUp, appliedStateDifferenceDown;

            std::map< double, Eigen::VectorXd > upPerturbedStateHistory, downPerturbedStateHistory;
            Eigen::VectorXd upPerturbedState, downPerturbedState;
            perturbedState.setZero( );
            perturbedState( j ) += perturbationMultiplier * statePerturbation( j );
            upPerturbedStateHistory = executeRotationSimulation< double, double >(
                                       perturbedState, appliedStateDifferenceUp, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 ).second.first;
                                    //    .second.at( 0 );
            upPerturbedState = upPerturbedStateHistory.rbegin()->second;

             std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > upPerturbedsStateInterpolator =
                std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( upPerturbedStateHistory, 6 );


            // std::cout << "upPerturbedState " << upPerturbedState.transpose() << std::endl;

            Eigen::VectorXd stateDifferenceUp = upPerturbedState - nominalState;
            Eigen::VectorXd stateDifferenceStm = stateTransitionAndSensitivityMatrixAtEpoch * appliedStateDifferenceUp;

                       std::cout<<"Test output "<<test<<" "<<j<<"stateDifferenceUp"<<std::endl<<
                                  stateDifferenceUp<<std::endl<<std::endl<<
                                //   "stateTransitionAndSensitivityMatrixAtEpoch"<<std::endl<<
                                                        //  stateTransitionAndSensitivityMatrixAtEpoch<<std::endl<<std::endl<<
                                //   "appliedStateDifferenceUp"<<std::endl<<
                                                        //  appliedStateDifferenceUp<<std::endl<<std::endl<<
                                  "( stateTransitionAndSensitivityMatrixAtEpoch * appliedStateDifferenceUp )"<<std::endl<<
                                                         stateDifferenceStm <<std::endl<<std::endl;

            Eigen::VectorXd relativeErrors = Eigen::VectorXd::Zero(12);
            for ( unsigned int k = 0 ; k < 12 ; k++ )
            {
                relativeErrors[k] = ( stateDifferenceUp[k] - stateDifferenceStm[k] ) / stateDifferenceUp[k];
            }

            // std::cout << "nominalState " << nominalState.transpose() << std::endl;
            std::cout << "relativeErrors" << std::endl;
            std::cout << relativeErrors << std::endl;


            for ( unsigned int k = 0 ; k < refEpochs.size() ; k++ )
            {
                double currentTime = refEpochs[k];
                Eigen::VectorXd currentUpState = upPerturbedsStateInterpolator->interpolate(currentTime).segment( 0, 12 );
                Eigen::VectorXd currentNominalState = nominalStateInterpolator->interpolate(currentTime).segment( 0, 12 );
                Eigen::MatrixXd currentPartial = ( currentUpState - currentNominalState );

                // interpolatedUpPerturbedStateHistory[ currentTime ] = currentUpState;
                // interpolatedDownPerturbedStateHistory[ currentTime ] = currentDownState;
                numericalPartialsHistory.at( currentTime ).block( 0, j, 12, 1 ) = currentPartial;
                modifiedAnalyticalStmHistory.at( currentTime ).block( 0, j, 12, 1 ) = analyticalStmHistory.at( currentTime ) * appliedStateDifferenceUp;

                interpolatedUpPerturbedStateHistory.at( currentTime ).block( 0, j, 12, 1 ) = currentUpState;
                interpolatedDownPerturbedStateHistory.at( currentTime ).block( 0, j, 12, 1 ) = currentNominalState;
            }


        //     if( test == 0 )
        //     {
        //         Eigen::VectorXd testMatrix =
        //                 ( stateTransitionAndSensitivityMatrixAtEpoch.block( 0, 0, 12, 12 ) * appliedStateDifferenceUp );
        //         TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( testMatrix.segment( 0, 6 ) ), ( stateDifferenceUp.segment( 0, 6 ) ), 1.5E-3 );
        //     }
        //     else
        //     {
        //         Eigen::VectorXd testMatrix =
        //                 ( stateTransitionAndSensitivityMatrixAtEpoch.block( 0, 0, 12, 12 ) * appliedStateDifferenceUp );
        //         TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( testMatrix.segment( 6, 7 ) ), ( stateDifferenceUp.segment( 6, 7 ) ), 1.0E-5 );
        //     }
        }
    }

    
    // PARTIALS WRT GRAVITY DEFORMATION
    for( unsigned int j = 7; j < 12; j++ )
    {
        // std::map< double, Eigen::MatrixXd > partialsWrtGravityDeformationHistory, interpolatedUpPerturbedStateHistory, interpolatedDownPerturbedStateHistory;

        Eigen::Matrix< double, 12, 1 > appliedStateDifferenceUp, appliedStateDifferenceDown;

        std::map< double, Eigen::VectorXd > upPerturbedStateHistory, downPerturbedStateHistory;
        Eigen::VectorXd upPerturbedState, downPerturbedState;
        perturbedState.setZero( );
        perturbedState( j ) += statePerturbation( j );
        upPerturbedStateHistory = executeRotationSimulation< double, double >(
                                   perturbedState, appliedStateDifferenceUp, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 ).second.first; //.rbegin()->second;
                                //    .second.at( 0 );
        upPerturbedState = upPerturbedStateHistory.rbegin()->second;

        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > upPerturbedStateInterpolator =
                std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( upPerturbedStateHistory, 6 );

        perturbedState.setZero( );
        perturbedState( j ) -= statePerturbation( j );
        downPerturbedStateHistory = executeRotationSimulation< double, double >(
                                     perturbedState, appliedStateDifferenceDown, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 ).second.first; //.rbegin()->second;
                                    //  .second.at( 0 );
        downPerturbedState = downPerturbedStateHistory.rbegin()->second;

        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > downPerturbedStateInterpolator =
                std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( downPerturbedStateHistory, 6 );

        for ( unsigned int k = 0 ; k < refEpochs.size() ; k++ )
        {
            double currentTime = refEpochs[k];
            Eigen::VectorXd currentUpState = upPerturbedStateInterpolator->interpolate(currentTime).segment( 0, 12 );
            Eigen::VectorXd currentDownState = downPerturbedStateInterpolator->interpolate(currentTime).segment( 0, 12 );
            Eigen::MatrixXd currentPartial = ( currentUpState - currentDownState ) / ( 2.0 * statePerturbation( j ) );

            // interpolatedUpPerturbedStateHistory[ currentTime ] = currentUpState;
            // interpolatedDownPerturbedStateHistory[ currentTime ] = currentDownState;
            numericalPartialsHistory.at( currentTime ).block( 0, j, 12, 1 ) = currentPartial;

            interpolatedUpPerturbedStateHistory.at( currentTime ).block( 0, j, 12, 1 ) = currentUpState;
            interpolatedDownPerturbedStateHistory.at( currentTime ).block( 0, j, 12, 1 ) = currentDownState;
        }

        manualPartial.block( 0, j, 12, 1 ) =
                ( upPerturbedState.segment( 0, 12 ) - downPerturbedState.segment( 0, 12 ) ) / ( 2.0 * statePerturbation( j ) );
    }

//     TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//             ( manualPartial.block( 0, 0, 13, 6 ) ), ( stateTransitionAndSensitivityMatrixAtEpoch.block( 0, 0, 13, 6 ) ), 1.0E-4 );


    // PARTIALS WRT ANGULAR VELOCITY VECTOR
    for( unsigned int j = 4; j < 7 ; j++ )
    {
        Eigen::Matrix< double, 12, 1 > appliedStateDifferenceUp, appliedStateDifferenceDown;

        std::map< double, Eigen::VectorXd > upPerturbedStateHistory, downPerturbedStateHistory;
        Eigen::VectorXd upPerturbedState, downPerturbedState;
        perturbedState.setZero( );
        perturbedState( j ) += statePerturbation( j );
        upPerturbedStateHistory = executeRotationSimulation< double, double >(
                                   perturbedState, appliedStateDifferenceUp, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 ).second.first; //.rbegin()->second;
                                //    .second.at( 0 );
        upPerturbedState = upPerturbedStateHistory.rbegin()->second;

        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > upPerturbedStateInterpolator =
                std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( upPerturbedStateHistory, 6 );

        perturbedState.setZero( );
        perturbedState( j ) -= statePerturbation( j );
        downPerturbedStateHistory = executeRotationSimulation< double, double >(
                                     perturbedState, appliedStateDifferenceDown, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 ).second.first; //.rbegin()->second;
                                    //  .second.at( 0 );
        downPerturbedState = downPerturbedStateHistory.rbegin()->second;

        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > downPerturbedStateInterpolator =
                std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( downPerturbedStateHistory, 6 );

        manualPartial.block( 0, j, 12, 1 ) =
                ( upPerturbedState.segment( 0, 12 ) - downPerturbedState.segment( 0, 12 ) ) / ( 2.0 * statePerturbation( j ) );

        for ( unsigned int k = 0 ; k < refEpochs.size() ; k++ )
        {
            double currentTime = refEpochs[k];
            Eigen::VectorXd currentUpState = upPerturbedStateInterpolator->interpolate(currentTime).segment( 0, 12 );
            Eigen::VectorXd currentDownState = downPerturbedStateInterpolator->interpolate(currentTime).segment( 0, 12 );
            Eigen::MatrixXd currentPartial = ( currentUpState - currentDownState ) / ( 2.0 * statePerturbation( j ) );

            // interpolatedUpPerturbedStateHistory[ currentTime ] = currentUpState;
            // interpolatedDownPerturbedStateHistory[ currentTime ] = currentDownState;
            numericalPartialsHistory.at( currentTime ).block( 0, j, 12, 1 ) = currentPartial;

            interpolatedUpPerturbedStateHistory.at( currentTime ).block( 0, j, 12, 1 ) = currentUpState;
            interpolatedDownPerturbedStateHistory.at( currentTime ).block( 0, j, 12, 1 ) = currentDownState;
        }
    }

    tudat::input_output::writeDataMapToTextFile(
        numericalPartialsHistory,
        "manualStateTransitionMatrixHistory.dat", "/Users/sam.fayolle/Downloads/" );
    tudat::input_output::writeDataMapToTextFile(
        interpolatedUpPerturbedStateHistory,
        "upPerturbedStateHistory.dat", "/Users/sam.fayolle/Downloads/" );
    tudat::input_output::writeDataMapToTextFile(
        interpolatedDownPerturbedStateHistory,
        "downPerturbedStateHistory.dat", "/Users/sam.fayolle/Downloads/" );

    tudat::input_output::writeDataMapToTextFile(
        modifiedAnalyticalStmHistory,
        "stateTransitionMatrixHistory.dat", "/Users/sam.fayolle/Downloads/" );
    tudat::input_output::writeDataMapToTextFile(
        currentOutput.second.first,
        "refStateHistory.dat", "/Users/sam.fayolle/Downloads/" );
    tudat::input_output::writeDataMapToTextFile(
        currentOutput.second.second,
        "refDepVarHistory.dat", "/Users/sam.fayolle/Downloads/" );

    std::cout << " ----------------------------------------------- " << std::endl;
    std::cout << " PARTIALS WRT ANGULAR VELOCITY VECTOR " << std::endl;
    std::cout << "manualPartial " << std::endl;
    std::cout << manualPartial.block(0, 4, 12, 3) << std::endl;
    std::cout << "state transition matrix " << std::endl;
    std::cout << stateTransitionAndSensitivityMatrixAtEpoch.block(0, 4, 12, 3) << std::endl;

    Eigen::MatrixXd differences = manualPartial.block(0, 4, 12, 3) - stateTransitionAndSensitivityMatrixAtEpoch.block(0, 4, 12, 3);
    for ( unsigned int j = 0 ; j < 12 ; j++ )
    {
        for ( unsigned int k = 4 ; k < 4+3 ; k++ )
        {
            if ( stateTransitionAndSensitivityMatrixAtEpoch( j,k ) != 0.0 )
            {
                differences( j,k-4 ) /= stateTransitionAndSensitivityMatrixAtEpoch( j,k );
            }
        }
    }
    std::cout << "relative differences" << std::endl;
    std::cout << differences << std::endl;

    std::cout << " ----------------------------------------------- " << std::endl;
    std::cout << " PARTIALS WRT GRAVITY STATE DEFORMATION " << std::endl;
    std::cout << "manualPartial " << std::endl;
    std::cout << manualPartial.block(0, 7, 12, 5) << std::endl;
    std::cout << "state transition matrix " << std::endl;
    std::cout << stateTransitionAndSensitivityMatrixAtEpoch.block(0, 7, 12, 5) << std::endl;

    Eigen::MatrixXd differences2 = manualPartial.block(0, 7, 12, 5) - stateTransitionAndSensitivityMatrixAtEpoch.block(0, 7, 12, 5);
    for ( unsigned int j = 0 ; j < 12 ; j++ )
    {
        for ( unsigned int k = 7 ; k < 7+5 ; k++ )
        {
            if ( stateTransitionAndSensitivityMatrixAtEpoch( j,k ) != 0.0 )
            {
                differences2( j,k-7 ) /= stateTransitionAndSensitivityMatrixAtEpoch( j,k );
            }
        }
    }
    std::cout << "relative differences" << std::endl;
    std::cout << differences2 << std::endl;


//     // Check element separately
//     BOOST_CHECK_SMALL( std::fabs( manualPartial( 3, 1 + 10 ) - stateTransitionAndSensitivityMatrixAtEpoch( 3, 1 + 10 ) ), 1.0E-2 );
//     manualPartial( 3, 1 + 10 ) = stateTransitionAndSensitivityMatrixAtEpoch( 3, 1 + 10 );

//     TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//             ( manualPartial.block( 0, 10, 6, 3 ) ), ( stateTransitionAndSensitivityMatrixAtEpoch.block( 0, 10, 6, 3 ) ), 1.0E-3 );
//     TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//             ( manualPartial.block( 6, 10, 7, 3 ) ), ( stateTransitionAndSensitivityMatrixAtEpoch.block( 6, 10, 7, 3 ) ), 1.0E-5 );

//     // Numerically compute sensitivity matrix
//     for( int j = 0; j < numberOfParametersToEstimate; j++ )
//     {
//         Eigen::Matrix< double, 13, 1 > appliedStateDifference;

//         Eigen::VectorXd upPerturbedState, downPerturbedState;
//         perturbedState.setZero( );
//         perturbedParameter.setZero( );
//         perturbedParameter( j ) += parameterPerturbation( j );

//         //        std::cout<<"Test "<<j<<" "<<perturbedParameter.transpose( )<<std::endl;

//         upPerturbedState =
//                 executePhobosRotationSimulation< double, double >( perturbedState, appliedStateDifference, perturbedParameter, 0 )
//                         .second.at( 0 );

//         perturbedParameter.setZero( );
//         perturbedParameter( j ) -= parameterPerturbation( j );
//         downPerturbedState =
//                 executePhobosRotationSimulation< double, double >( perturbedState, appliedStateDifference, perturbedParameter, 0 )
//                         .second.at( 0 );

//         manualPartial.block( 0, j + 13, 13, 1 ) =
//                 ( upPerturbedState.segment( 0, 13 ) - downPerturbedState.segment( 0, 13 ) ) / ( 2.0 * parameterPerturbation( j ) );
//     }
//     //    std::cout<<manualPartial<<std::endl<<std::endl
//     //            <<stateTransitionAndSensitivityMatrixAtEpoch<<std::endl<<std::endl<<
//     //              ( manualPartial - stateTransitionAndSensitivityMatrixAtEpoch ).cwiseQuotient(
//     //                  stateTransitionAndSensitivityMatrixAtEpoch )<<std::endl;

//     // Check three values separately: could not find perturbations for which all partials are sufficiently within the linear regime

//     BOOST_CHECK_SMALL( std::fabs( manualPartial( 4, 5 + 13 ) - stateTransitionAndSensitivityMatrixAtEpoch( 4, 5 + 13 ) ), 1.0E-4 );
//     BOOST_CHECK_SMALL( std::fabs( manualPartial( 11, 4 + 13 ) - stateTransitionAndSensitivityMatrixAtEpoch( 11, 4 + 13 ) ), 1.0E-5 );
//     BOOST_CHECK_SMALL( std::fabs( manualPartial( 12, 3 + 13 ) - stateTransitionAndSensitivityMatrixAtEpoch( 12, 3 + 13 ) ), 1.0E-2 );

//     manualPartial( 4, 5 + 13 ) = stateTransitionAndSensitivityMatrixAtEpoch( 4, 5 + 13 );
//     manualPartial( 11, 4 + 13 ) = stateTransitionAndSensitivityMatrixAtEpoch( 11, 4 + 13 );
//     manualPartial( 12, 3 + 13 ) = stateTransitionAndSensitivityMatrixAtEpoch( 12, 3 + 13 );
//     std::cout << manualPartial << std::endl
//               << std::endl
//               << ( manualPartial - stateTransitionAndSensitivityMatrixAtEpoch ).cwiseQuotient( stateTransitionAndSensitivityMatrixAtEpoch )
//               << std::endl;
//     TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//             ( manualPartial.block( 0, 13, 13, 8 ) ), ( stateTransitionAndSensitivityMatrixAtEpoch.block( 0, 13, 13, 8 ) ), 7.5E-3 );
}


// template< typename TimeType = double, typename StateScalarType = double >
// std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
//            std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >, std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >
// executeRotationOnlySimulation( 
//     const Eigen::Matrix< StateScalarType, 7, 1 > initialStateDifference,
//     Eigen::Matrix< StateScalarType, 7, 1 >& appliedStateDifference,
//     const Eigen::VectorXd parameterPerturbation = Eigen::VectorXd::Zero( 8 ),
//     const bool propagateVariationalEquations = 1 )
// {

//     double initialTime = 0.0;
//     double finalTime = 20.0 * 86400.0; // 10.0 * 3600.0; //86400.0;
//     int numberOfParametersToEstimate = 8;

//     std::string globalFrameOrigin = "Jupiter";
//     std::string globalFrameOrientation = "J2000";

//     std::vector< std::string > bodiesToCreate = { "Jupiter", "Io" }; 

//     double muJupiter = getBodyGravitationalParameter( "Jupiter" );
//     double muIo = getBodyGravitationalParameter( "Io" );
//     double muEffective = muJupiter + muIo;

//     // Get body settings.
//     BodyListSettings bodySettings =
//             getDefaultBodySettings( bodiesToCreate, initialTime - 86400.0, finalTime + 86400.0, globalFrameOrigin, globalFrameOrientation );

//     // Set Io's rotation to synchronous        
//     bodySettings.get( "Io" )->rotationModelSettings = simulation_setup::synchronousRotationModelSettings( "Jupiter", "J2000", "IAU_Io" );

//     // Create system of bodies
//     SystemOfBodies bodies = createSystemOfBodies( bodySettings );


//     // Set Io's mean moment of inertia
//     std::shared_ptr< SphericalHarmonicsGravityField > ioGravityFieldModel = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) );
//     double scaledMeanMomentOfInertia = 0.37685;
//     ioGravityFieldModel->setScaledMeanMomentOfInertia( scaledMeanMomentOfInertia );

//     Eigen::MatrixXd ioCosineCoefs = ioGravityFieldModel->getCosineCoefficients( );
//     Eigen::MatrixXd ioSineCoefs = ioGravityFieldModel->getSineCoefficients( );
//     // ioCosineCoefs(2,1) = 1.0e-5;
//     // ioSineCoefs(2,1) = 1.0e-5;
//     // ioSineCoefs(2,2) = 2.0e-5;
//     // ioGravityFieldModel->setCosineCoefficients( ioCosineCoefs );
//     // ioGravityFieldModel->setSineCoefficients( ioSineCoefs );
//     // std::cout << "ioCosineCoefs" << std::endl;
//     // std::cout << ioCosineCoefs << std::endl;
//     // std::cout << "ioSineCoefs" << std::endl;
//     // std::cout << ioSineCoefs << std::endl;

//     // std::shared_ptr< SphericalHarmonicsGravityField > newIoGravityFieldModel = std::make_shared< SphericalHarmonicsGravityField >( 
//         // ioGravityFieldModel->getGravitationalParameter( ), ioGravityFieldModel->getReferenceRadius( ), ioCosineCoefs, ioSineCoefs, "IAU_Io", scaledMeanMomentOfInertia );
//     // bodies.at( "Io" )->setGravityFieldModel( newIoGravityFieldModel );

//     // // Retrieve Io's static gravity coefficients
//     // ioCosineCoefs = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )->getCosineCoefficients( );
//     // ioSineCoefs = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )->getSineCoefficients( );
//     // Eigen::Vector5d normalisedStaticGravity = Eigen::Vector5d::Zero( );
//     // normalisedStaticGravity[0] = ioCosineCoefs(2,0);
//     // // normalisedStaticGravity[1] = 1.0e-5; // ioCosineCoefs(2,1);
//     // normalisedStaticGravity[2] = ioCosineCoefs(2,2);
//     // // normalisedStaticGravity[3] = 1.0e-5; // ioSineCoefs(2,1);
//     // // normalisedStaticGravity[4] = 1.0e-4; // ioSineCoefs(2,2);
//     // std::cout << "normalisedStaticGravity " << normalisedStaticGravity.transpose() << std::endl;
//     // bodies.at( "Io" )->setStaticDegreeTwoCoefficients( normalisedStaticGravity );

//     Eigen::Vector6d ioKeplerElements = Eigen::Vector6d::Zero( );
//     double ioSemiMajorAxis = 4.2e8;
//     ioKeplerElements( 0 ) = ioSemiMajorAxis;
//     ioKeplerElements( 1 ) = 0.004;
//     ioKeplerElements( 2 ) = 0.002 * mathematical_constants::PI / 180.0;


//     double rotationRateIo = std::sqrt( muEffective / ( ioSemiMajorAxis * ioSemiMajorAxis * ioSemiMajorAxis ) );

//     // Retrieve initial rotational state from ephemeris
//     Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRotationState = getInitialRotationalStateOfBody( "Io", "J2000",  bodies, initialTime );
//     initialRotationState[6] = rotationRateIo;

//     // ///// ALTERNATIVE DEFINITION OF IO'S ROTATIONAL STATE
//     // Eigen::Quaterniond noRotationQuaternion = Eigen::Quaterniond( Eigen::AngleAxisd( 1.0E-0, Eigen::Vector3d::UnitZ( ) ) *
//     //                                                               Eigen::AngleAxisd( 2.0E-0, Eigen::Vector3d::UnitX( ) ) *
//     //                                                               Eigen::AngleAxisd( -0.5E-0, Eigen::Vector3d::UnitZ( ) ) );

//     // initialRotationState( 0 ) = noRotationQuaternion.w( );
//     // initialRotationState( 1 ) = noRotationQuaternion.x( );
//     // initialRotationState( 2 ) = noRotationQuaternion.y( );
//     // initialRotationState( 3 ) = noRotationQuaternion.z( );
//     // initialRotationState( 4 ) = 1.0e-5;
//     // initialRotationState( 5 ) = -1.0e-5;
//     // initialRotationState( 6 ) = rotationRateIo;
//     // ////

//     // // initialRotationState = manualInitialState.segment( 0, 7 );

//     Eigen::Matrix< double, 7, 1 > originalRotationState = initialRotationState;
//     Eigen::Matrix< double, 7, 1 > stateDifferenceToAdd = initialStateDifference.segment( 0, 7 );
//     std::cout << "stateDifferenceToAdd " << stateDifferenceToAdd.transpose() << std::endl;

//     initialRotationState += stateDifferenceToAdd;
//     initialRotationState( 0 ) = originalRotationState( 0 ) / std::fabs( originalRotationState( 0 ) ) *
//             std::sqrt( 1.0 - std::pow( initialRotationState.segment( 1, 3 ).norm( ), 2.0 ) );

//     appliedStateDifference.segment( 0, 7 ) = initialRotationState - originalRotationState;
//     // std::cout << "appliedStateDifference" << std::endl;
//     // std::cout << appliedStateDifference.transpose() << std::endl;
    
//     std::map< double, Eigen::Matrix< double, 7, 1 > > dummyRotationMap;
//     dummyRotationMap[ -1.0E100 ] = initialRotationState;
//     dummyRotationMap[ 1.0E100 ] = initialRotationState;

//     std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > dummyInterpolator =
//             std::make_shared< interpolators::LinearInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( dummyRotationMap );
//     bodies.at( "Io" )->setRotationalEphemeris(
//         std::make_shared< TabulatedRotationalEphemeris< double, double > >( dummyInterpolator, "J2000", "IAU_Io" ) );

//     // Retrieve body objects for Io and Jupiter
//     std::shared_ptr< Body > io = bodies.at( "Io" );
//     std::shared_ptr< Body > jupiter = bodies.at( "Jupiter" );

//     // Update Jupiter and Io to current state
//     io->setStateFromEphemeris( initialTime );
//     jupiter->setStateFromEphemeris( initialTime );
//     io->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );

//     // Check translational intial state
//     Eigen::Vector6d initialTranslationalState = io->getState( ) - jupiter->getState( );
//     // std::cout << "initialTranslationalState " << initialTranslationalState.transpose( ) << std::endl;

//     // std::cout << "initialRotationState " << initialRotationState.transpose( ) << std::endl;    

//     SelectedTorqueMap torqueMap;
//     // torqueMap[ "Io" ][ "Jupiter" ].push_back( std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );


//     // Create integrator settings
//     double timeStep = 100.0 / 2.0;
//     std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > > ( 
//         initialTime, timeStep, rungeKutta87DormandPrince, 1.0e-5, 1000.0 ); 

//     // Define propagator settings.
//     std::vector< std::string > bodiesToIntegrate;
//     bodiesToIntegrate.push_back( "Io" );

//     // Create torque models
//     basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap( bodies, torqueMap, bodiesToIntegrate );

//     std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
//             std::make_shared< RotationalStatePropagatorSettings< double > >(
//                     torqueModelMap,
//                     bodiesToIntegrate,
//                     initialRotationState,
//                     initialTime,
//                     integratorSettings,
//                     std::make_shared< PropagationTimeTerminationSettings >( finalTime, true ) );


//     std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
//     propagatorSettingsList.push_back( rotationalPropagatorSettings );

//     // Define dependent variables
//     std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
//     dependentVariablesList.push_back(
//             std::make_shared< SingleDependentVariableSaveSettings >( body_fixed_relative_spherical_position, "Jupiter", "Io" ) );

//     std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings = std::make_shared< MultiTypePropagatorSettings< double > >(
//             propagatorSettingsList, std::make_shared< PropagationTimeTerminationSettings >( finalTime, true ), dependentVariablesList );   
    

//     // Define parameters.
//     std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
//     {
//         parameterNames = getInitialStateParameterSettings< StateScalarType >( propagatorSettings, bodies );
//         parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Io", mean_moment_of_inertia ) );
//         parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
//                 1, 0, 2, 2, "Io", spherical_harmonics_cosine_coefficient_block ) );
//         parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
//                 2, 1, 2, 2, "Io", spherical_harmonics_sine_coefficient_block ) );
//     }

//     // Create parameters
//     std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
//             createParametersToEstimate( parameterNames, bodies );
//     // printEstimatableParameterEntries( parametersToEstimate );

//     Eigen::MatrixXd constraintStateMultiplier;
//     Eigen::VectorXd constraintRightHandSide;
//     parametersToEstimate->getConstraints( constraintStateMultiplier, constraintRightHandSide );
//     //    std::cout<<"Unit rotation: "<<std::endl<<unitRotationState.transpose( )<<std::endl;
//     //    std::cout<<"Constraints: "<<std::endl<<constraintStateMultiplier.transpose( )<<std::endl;

//     TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( constraintStateMultiplier.block( 0, 0, 1, 4 ) ),
//                                        ( initialRotationState.segment( 0, 4 ) ).transpose( ),
//                                        std::numeric_limits< double >::epsilon( ) );
//     TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( constraintStateMultiplier.block( 0, 4, 1, 9 ) ),
//                                        ( Eigen::MatrixXd::Zero( 1, 9 ) ),
//                                        std::numeric_limits< double >::epsilon( ) );
//     TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//             ( constraintRightHandSide.block( 0, 0, 1, 1 ) ), ( Eigen::MatrixXd::Zero( 1, 1 ) ), std::numeric_limits< double >::epsilon( ) );

//     // Perturb parameters.
//     Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > parameterVector =
//             parametersToEstimate->template getFullParameterValues< StateScalarType >( );
//     parameterVector.block( 7, 0, numberOfParametersToEstimate, 1 ) += parameterPerturbation;
//     //    std::cout<<"Parameter perturbation "<<
//     //               ( parametersToEstimate->template getFullParameterValues< StateScalarType >( ) -
//     //                 parameterVector ).transpose( )<<std::endl<<
//     //               parameterPerturbation.transpose( )<<std::endl;;
//     parametersToEstimate->resetParameterValues( parameterVector );

//     std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
//                std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >, 
//                         std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > > results;

//     {
//         // Create dynamics simulator
//         propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType > dynamicsSimulator =
//                 propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType >(
//                         bodies,
//                         integratorSettings,
//                         propagatorSettings,
//                         parametersToEstimate,
//                         1,
//                         std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
//                         0,
//                         0 );

//         // Propagate requested equations.
//         if( propagateVariationalEquations )
//         {
//             dynamicsSimulator.integrateVariationalAndDynamicalEquations( propagatorSettings->getInitialStates( ), 1 );
//         }
//         else
//         {
//             dynamicsSimulator.integrateDynamicalEquationsOfMotionOnly( propagatorSettings->getInitialStates( ) );
//         }

//         // Retrieve test data
//         double testEpoch = finalTime; // initialTime + 1.0 * 3600.0;
//         Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > testStates = Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( 7 );
//         testStates.segment( 0, 7 ) = dynamicsSimulator.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( ).rbegin()->second.segment( 0, 7 ); //bodies.at( "Io" )->getRotationalEphemeris( )->getRotationStateVector( testEpoch );
//         std::cout << "final state " << testStates.transpose() << std::endl;\

//         if( propagateVariationalEquations )
//         {
//             // results.first.push_back(
//                     // dynamicsSimulator.getStateTransitionMatrixInterface( )->getCombinedStateTransitionAndSensitivityMatrix( testEpoch ) );
//             results.first = dynamicsSimulator.getStateTransitionMatrixSolution( );
//             Eigen::MatrixXd testMatrixDirect =
//                     dynamicsSimulator.getStateTransitionMatrixInterface( )->getCombinedStateTransitionAndSensitivityMatrix( testEpoch );
//             Eigen::MatrixXd testMatrixFull =
//                     dynamicsSimulator.getStateTransitionMatrixInterface( )->getFullCombinedStateTransitionAndSensitivityMatrix( testEpoch );
//             TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testMatrixDirect, testMatrixFull, std::numeric_limits< double >::epsilon( ) );
//         }
//         // results.second.push_back( testStates );
//         results.second.first = dynamicsSimulator.getEquationsOfMotionSolution( );

//         // Retrieve dependent variables history
//         std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > simulator = dynamicsSimulator.getDynamicsSimulator( );
//         std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > dependentVariablesHistory = simulator->getDependentVariableHistory();
//         results.second.second = dependentVariablesHistory;

//     }
//     return results;
// }

// BOOST_AUTO_TEST_CASE( testRotationalVariationalEquations )
// {
//     // Load spice kernels.
//     spice_interface::loadStandardSpiceKernels( );

//     std::pair< std::map< double, Eigen::MatrixXd >, std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > > currentOutput;

//     // Define variables for numerical differentiation
//     Eigen::Matrix< double, 7, 1 > perturbedState;
//     Eigen::Matrix< double, 7, 1 > statePerturbation;

//     // Define parameter perturbation
//     int numberOfParametersToEstimate = 8;
//     double sphericalHarmonicsPerturbation = 1.0E-4;
//     Eigen::Matrix< double, 8, 1 > perturbedParameter;
//     Eigen::Matrix< double, 8, 1 > parameterPerturbation;
//     parameterPerturbation = Eigen::Matrix< double, 8, 1 >::Constant( sphericalHarmonicsPerturbation );
//     parameterPerturbation( 0 ) = 1.0E-4;

//     // Compute state transition and sensitivity matrices
//     Eigen::Matrix< double, 7, 1 > appliedStateDifference;
//     currentOutput = executeRotationOnlySimulation< double, double >( Eigen::Matrix< double, 7, 1 >::Zero( ), appliedStateDifference );
//     Eigen::MatrixXd stateTransitionAndSensitivityMatrixAtEpoch = currentOutput.first.rbegin()->second; // currentOutput.first.at( 0 );
//     Eigen::VectorXd nominalState = currentOutput.second.first.rbegin()->second; // currentOutput.second.at( 0 );

//     std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > nominalStateInterpolator =
//                 std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( currentOutput.second.first, 6 );

//     std::vector< double > refEpochs;
//     for ( auto it : currentOutput.first )
//     {
//         refEpochs.push_back( it.first );
//     }

//     std::map< double, Eigen::MatrixXd > analyticalStmHistory = currentOutput.first;
//     std::map< double, Eigen::MatrixXd > modifiedAnalyticalStmHistory = analyticalStmHistory;


//     // Define state perturbation
//     statePerturbation = ( Eigen::Matrix< double, 7, 1 >( ) << 
//                           1.0E-5,
//                           1.0E-5,
//                           1.0E-5,
//                           1.0E-5,
//                           1.0E-9,
//                           1.0E-9,
//                           1.0E-9 )
//                                 .finished( );

//     Eigen::MatrixXd manualPartial = Eigen::MatrixXd::Zero( 7, 7 + numberOfParametersToEstimate );

//     // Numerically compute state transition matrix
//     std::map< double, Eigen::MatrixXd > numericalPartialsHistory, interpolatedUpPerturbedStateHistory, interpolatedDownPerturbedStateHistory;
//     for ( unsigned int k = 0 ; k < refEpochs.size() ; k++ )
//     {
//         numericalPartialsHistory[ refEpochs[k] ] = Eigen::MatrixXd::Zero( 7, 7 );
//         interpolatedUpPerturbedStateHistory[ refEpochs[k] ] = Eigen::MatrixXd::Zero( 7, 7 );
//         interpolatedDownPerturbedStateHistory[ refEpochs[k] ] = Eigen::MatrixXd::Zero( 7, 7 );
//     }

//     for( unsigned int test = 0; test < 1; test++ )
//     {
//         double perturbationMultiplier = ( test == 0 ? 1.0 : 1.0E-3 );
//         for( unsigned int j = 0; j < 4; j++ )
//         {
//             Eigen::Matrix< double, 7, 1 > appliedStateDifferenceUp, appliedStateDifferenceDown;

//             std::map< double, Eigen::VectorXd > upPerturbedStateHistory, downPerturbedStateHistory;
//             Eigen::VectorXd upPerturbedState, downPerturbedState;
//             perturbedState.setZero( );
//             perturbedState( j ) += perturbationMultiplier * statePerturbation( j );
//             upPerturbedStateHistory = executeRotationOnlySimulation< double, double >(
//                                        perturbedState, appliedStateDifferenceUp, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 ).second.first;
//                                     //    .second.at( 0 );
//             upPerturbedState = upPerturbedStateHistory.rbegin()->second;

//              std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > upPerturbedsStateInterpolator =
//                 std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( upPerturbedStateHistory, 6 );


//             // std::cout << "upPerturbedState " << upPerturbedState.transpose() << std::endl;

//             Eigen::VectorXd stateDifferenceUp = upPerturbedState - nominalState;
//             Eigen::VectorXd stateDifferenceStm = stateTransitionAndSensitivityMatrixAtEpoch * appliedStateDifferenceUp;

//                        std::cout<<"Test output "<<test<<" "<<j<<"stateDifferenceUp"<<std::endl<<
//                                   stateDifferenceUp<<std::endl<<std::endl<<
//                                 //   "stateTransitionAndSensitivityMatrixAtEpoch"<<std::endl<<
//                                                         //  stateTransitionAndSensitivityMatrixAtEpoch<<std::endl<<std::endl<<
//                                 //   "appliedStateDifferenceUp"<<std::endl<<
//                                                         //  appliedStateDifferenceUp<<std::endl<<std::endl<<
//                                   "( stateTransitionAndSensitivityMatrixAtEpoch * appliedStateDifferenceUp )"<<std::endl<<
//                                                          stateDifferenceStm <<std::endl<<std::endl;

//             Eigen::VectorXd relativeErrors = Eigen::VectorXd::Zero(7);
//             for ( unsigned int k = 0 ; k < 7 ; k++ )
//             {
//                 relativeErrors[k] = ( stateDifferenceUp[k] - stateDifferenceStm[k] ) / stateDifferenceUp[k];
//             }

//             // std::cout << "nominalState " << nominalState.transpose() << std::endl;
//             std::cout << "relativeErrors" << std::endl;
//             std::cout << relativeErrors << std::endl;


//             for ( unsigned int k = 0 ; k < refEpochs.size() ; k++ )
//             {
//                 double currentTime = refEpochs[k];
//                 Eigen::VectorXd currentUpState = upPerturbedsStateInterpolator->interpolate(currentTime).segment( 0, 7 );
//                 Eigen::VectorXd currentNominalState = nominalStateInterpolator->interpolate(currentTime).segment( 0, 7 );
//                 Eigen::MatrixXd currentPartial = ( currentUpState - currentNominalState );

//                 numericalPartialsHistory.at( currentTime ).block( 0, j, 7, 1 ) = currentPartial;
//                 modifiedAnalyticalStmHistory.at( currentTime ).block( 0, j, 7, 1 ) = analyticalStmHistory.at( currentTime ) * appliedStateDifferenceUp;
//             }
//         }
//     }


//     // PARTIALS WRT ANGULAR VELOCITY VECTOR
//     for( unsigned int j = 4; j < 7 ; j++ )
//     {
//         Eigen::Matrix< double, 7, 1 > appliedStateDifferenceUp, appliedStateDifferenceDown;

//         std::map< double, Eigen::VectorXd > upPerturbedStateHistory, downPerturbedStateHistory;
//         Eigen::VectorXd upPerturbedState, downPerturbedState;
//         perturbedState.setZero( );
//         perturbedState( j ) += statePerturbation( j );
//         upPerturbedStateHistory = executeRotationOnlySimulation< double, double >(
//                                    perturbedState, appliedStateDifferenceUp, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 ).second.first; //.rbegin()->second;
//                                 //    .second.at( 0 );
//         upPerturbedState = upPerturbedStateHistory.rbegin()->second;

//         std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > upPerturbedStateInterpolator =
//                 std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( upPerturbedStateHistory, 6 );

//         perturbedState.setZero( );
//         perturbedState( j ) -= statePerturbation( j );
//         downPerturbedStateHistory = executeRotationOnlySimulation< double, double >(
//                                      perturbedState, appliedStateDifferenceDown, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 ).second.first; //.rbegin()->second;
//                                     //  .second.at( 0 );
//         downPerturbedState = downPerturbedStateHistory.rbegin()->second;

//         std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > downPerturbedStateInterpolator =
//                 std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( downPerturbedStateHistory, 6 );

//         manualPartial.block( 0, j, 7, 1 ) =
//                 ( upPerturbedState.segment( 0, 7 ) - downPerturbedState.segment( 0, 7 ) ) / ( 2.0 * statePerturbation( j ) );

//         for ( unsigned int k = 0 ; k < refEpochs.size() ; k++ )
//         {
//             double currentTime = refEpochs[k];
//             Eigen::VectorXd currentUpState = upPerturbedStateInterpolator->interpolate(currentTime).segment( 0, 7 );
//             Eigen::VectorXd currentDownState = downPerturbedStateInterpolator->interpolate(currentTime).segment( 0, 7 );
//             Eigen::MatrixXd currentPartial = ( currentUpState - currentDownState ) / ( 2.0 * statePerturbation( j ) );

//             numericalPartialsHistory.at( currentTime ).block( 0, j, 7, 1 ) = currentPartial;
//         }
//     }

//     tudat::input_output::writeDataMapToTextFile(
//         numericalPartialsHistory,
//         "manualStateTransitionMatrixHistoryRotation.dat", "/Users/sam.fayolle/Downloads/" );

//     tudat::input_output::writeDataMapToTextFile(
//         modifiedAnalyticalStmHistory,
//         "stateTransitionMatrixHistoryRotation.dat", "/Users/sam.fayolle/Downloads/" );
//     tudat::input_output::writeDataMapToTextFile(
//         currentOutput.second.first,
//         "refStateHistoryRotation.dat", "/Users/sam.fayolle/Downloads/" );
//     tudat::input_output::writeDataMapToTextFile(
//         currentOutput.second.second,
//         "refDepVarHistoryRotation.dat", "/Users/sam.fayolle/Downloads/" );

//     std::cout << " ----------------------------------------------- " << std::endl;
//     std::cout << " PARTIALS WRT ANGULAR VELOCITY VECTOR " << std::endl;
//     std::cout << "manualPartial " << std::endl;
//     std::cout << manualPartial.block(0, 4, 7, 3) << std::endl;
//     std::cout << "state transition matrix " << std::endl;
//     std::cout << stateTransitionAndSensitivityMatrixAtEpoch.block(0, 4, 7, 3) << std::endl;

//     Eigen::MatrixXd differences = manualPartial.block(0, 4, 7, 3) - stateTransitionAndSensitivityMatrixAtEpoch.block(0, 4, 7, 3);
//     for ( unsigned int j = 0 ; j < 7 ; j++ )
//     {
//         for ( unsigned int k = 4 ; k < 4+3 ; k++ )
//         {
//             if ( stateTransitionAndSensitivityMatrixAtEpoch( j,k ) != 0.0 )
//             {
//                 differences( j,k-4 ) /= stateTransitionAndSensitivityMatrixAtEpoch( j,k );
//             }
//         }
//     }
//     std::cout << "relative differences" << std::endl;
//     std::cout << differences << std::endl;

// }

// BOOST_AUTO_TEST_CASE( testMassRateVariationalEquations )
// {
//     // Load spice kernels.
//     spice_interface::loadStandardSpiceKernels( );

//     // Set simulation time settings.
//     const double simulationStartEpoch = 0.0;
//     const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

//     // Define body settings for simulation.
//     std::vector< std::string > bodiesToCreate;
//     bodiesToCreate.push_back( "Earth" );

//     // Create body objects.
//     BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, "Earth", "J2000" );
//     SystemOfBodies bodies = createSystemOfBodies( bodySettings );
//     bodies.createEmptyBody( "Asterix" );
//     double initialBodyMass = 2000.0;
//     bodies.getBody( "Asterix" )->setConstantBodyMass( initialBodyMass );
//     bodies.getBody( "Asterix" )
//             ->setRotationalEphemeris( createRotationModel(
//                     std::make_shared< OrbitalStateBasedRotationSettings >( "Earth", true, false, "J2000", "BodyFixed" ),
//                     "Asterix",
//                     bodies ) );

//     Eigen::MatrixXd finalStateTransitionTranslationalOnly;
//     Eigen::MatrixXd finalStateTransitionCoupled;

//     for( int test = 0; test < 2; test++ )
//     {
//         // Define propagator settings variables.
//         SelectedAccelerationMap accelerationMap;
//         std::vector< std::string > bodiesToPropagate;
//         std::vector< std::string > centralBodies;

//         // Define propagation settings.
//         std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
//         accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
//         addEngineModel( "Asterix", "MainEngine", std::make_shared< ConstantThrustMagnitudeSettings >( 1.0E-4, 300.0 ), bodies );

//         accelerationsOfAsterix[ "Asterix" ].push_back( std::make_shared< ThrustAccelerationSettings >( "MainEngine" ) );

//         accelerationMap[ "Asterix" ] = accelerationsOfAsterix;
//         bodiesToPropagate.push_back( "Asterix" );
//         centralBodies.push_back( "Earth" );

//         basic_astrodynamics::AccelerationMap accelerationModelMap =
//                 createAccelerationModelsMap( bodies, accelerationMap, bodiesToPropagate, centralBodies );

//         ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//         ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
//         ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//         // Set Keplerian elements for Asterix.
//         Eigen::Vector6d asterixInitialStateInKeplerianElements = ( Eigen::Vector6d( ) << 7500.0E3,
//                                                                    0.1,
//                                                                    unit_conversions::convertDegreesToRadians( 85.3 ),
//                                                                    unit_conversions::convertDegreesToRadians( 235.7 ),
//                                                                    unit_conversions::convertDegreesToRadians( 23.4 ),
//                                                                    unit_conversions::convertDegreesToRadians( 139.87 ) )
//                                                                          .finished( );

//         double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
//         const Eigen::Vector6d asterixInitialState =
//                 convertKeplerianToCartesianElements( asterixInitialStateInKeplerianElements, earthGravitationalParameter );

//         std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings;
//         std::shared_ptr< SingleArcPropagatorSettings< double > > translationalPropagatorSettings =
//                 std::make_shared< TranslationalStatePropagatorSettings< double > >(
//                         centralBodies, accelerationModelMap, bodiesToPropagate, asterixInitialState, simulationEndEpoch );
//         std::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings;
//         if( test == 0 )
//         {
//             propagatorSettings = translationalPropagatorSettings;
//         }
//         else
//         {
//             std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
//             massRateModels[ "Asterix" ] =
//                     createMassRateModel( "Asterix", std::make_shared< FromThrustMassRateSettings >( 1 ), bodies, accelerationModelMap );
//             massPropagatorSettings = std::make_shared< MassPropagatorSettings< double > >(
//                     std::vector< std::string >{ "Asterix" },
//                     massRateModels,
//                     ( Eigen::VectorXd( 1 ) << initialBodyMass ).finished( ),
//                     std::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch ) );

//             std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
//             propagatorSettingsList.push_back( translationalPropagatorSettings );
//             propagatorSettingsList.push_back( massPropagatorSettings );
//             propagatorSettings = std::make_shared< MultiTypePropagatorSettings< double > >(
//                     propagatorSettingsList, std::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch ) );
//         }

//         const double fixedStepSize = 5.0;
//         std::shared_ptr< IntegratorSettings<> > integratorSettings =
//                 std::make_shared< IntegratorSettings<> >( rungeKutta4, 0.0, fixedStepSize );

//         ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//         ///////////////////////    DEFINE PARAMETERS FOR WHICH SENSITIVITY IS TO BE COMPUTED   ////////////////////////////////
//         ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//         // Define list of parameters to estimate.
//         std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
//         parameterNames = getInitialStateParameterSettings< double >( propagatorSettings, bodies );

//         // Create parameters
//         std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
//                 createParametersToEstimate( parameterNames, bodies );

//         ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//         ///////////////////////             PROPAGATE ORBIT AND VARIATIONAL EQUATIONS         /////////////////////////////////
//         ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//         // Create simulation object and propagate dynamics.
//         SingleArcVariationalEquationsSolver<> variationalEquationsSimulator(
//                 bodies,
//                 integratorSettings,
//                 propagatorSettings,
//                 parametersToEstimate,
//                 true,
//                 std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
//                 false,
//                 true );

//         std::map< double, Eigen::MatrixXd > stateTransitionResult = variationalEquationsSimulator.getStateTransitionMatrixSolution( );
//         std::map< double, Eigen::MatrixXd > sensitivityResult = variationalEquationsSimulator.getSensitivityMatrixSolution( );
//         std::map< double, Eigen::VectorXd > integrationResult =
//                 variationalEquationsSimulator.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );

//         if( test == 0 )
//         {
//             finalStateTransitionTranslationalOnly = stateTransitionResult.rbegin( )->second;
//         }
//         else
//         {
//             finalStateTransitionCoupled = stateTransitionResult.rbegin( )->second;
//             for( int i = 0; i < 6; i++ )
//             {
//                 BOOST_CHECK_EQUAL( finalStateTransitionCoupled( 6, i ), 0.0 );
//             }
//             BOOST_CHECK_EQUAL( finalStateTransitionCoupled( 6, 6 ), 1.0 );
//         }

//         Eigen::MatrixXd initialStateTransition = stateTransitionResult.begin( )->second;
//         int numberOfStateEntries = ( test == 0 ) ? 6 : 7;

//         for( int i = 0; i < numberOfStateEntries; i++ )
//         {
//             for( int j = 0; j < numberOfStateEntries; j++ )
//             {
//                 if( i == j )
//                 {
//                     BOOST_CHECK_EQUAL( initialStateTransition( i, j ), 1.0 );
//                 }
//                 else
//                 {
//                     BOOST_CHECK_EQUAL( initialStateTransition( i, j ), 0.0 );
//                 }
//             }
//         }

//         if( test == 1 )
//         {
//             Eigen::VectorXd upPerturbedInitialState, downPerturbedInitialState;
//             double massPerturbation = 1.0;
//             Eigen::VectorXd perturbedInitialMass = Eigen::VectorXd::Zero( 1 );
//             {
//                 perturbedInitialMass( 0 ) = initialBodyMass + massPerturbation;
//                 bodies.getBody( "Asterix" )->setConstantBodyMass( perturbedInitialMass( 0 ) );
//                 massPropagatorSettings->resetInitialStates( perturbedInitialMass );
//                 std::dynamic_pointer_cast< MultiTypePropagatorSettings< double > >( propagatorSettings )->recomputeInitialStates( );
//                 SingleArcDynamicsSimulator<> dynamicsSimulator( bodies, integratorSettings, propagatorSettings );
//                 upPerturbedInitialState = dynamicsSimulator.getEquationsOfMotionNumericalSolution( ).rbegin( )->second;
//             }

//             {
//                 perturbedInitialMass( 0 ) = initialBodyMass - massPerturbation;
//                 bodies.getBody( "Asterix" )->setConstantBodyMass( perturbedInitialMass( 0 ) );
//                 massPropagatorSettings->resetInitialStates( perturbedInitialMass );
//                 std::dynamic_pointer_cast< MultiTypePropagatorSettings< double > >( propagatorSettings )->recomputeInitialStates( );
//                 SingleArcDynamicsSimulator<> dynamicsSimulator( bodies, integratorSettings, propagatorSettings );
//                 downPerturbedInitialState = dynamicsSimulator.getEquationsOfMotionNumericalSolution( ).rbegin( )->second;
//             }
//             Eigen::VectorXd numericalStatePartialWrtMass =
//                     ( ( upPerturbedInitialState - downPerturbedInitialState ) / ( 2.0 * massPerturbation ) ).block( 0, 0, 6, 1 );
//             Eigen::VectorXd analyticalStatePartialWrtMass = finalStateTransitionCoupled.block( 0, 6, 6, 1 );
//             TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalStatePartialWrtMass, analyticalStatePartialWrtMass, 1.0E-4 );
//         }
//     }

//     TUDAT_CHECK_MATRIX_CLOSE_FRACTION( finalStateTransitionCoupled.block( 0, 0, 6, 6 ), finalStateTransitionTranslationalOnly, 1.0E-6 );
// }


// template< typename TimeType = double, typename StateScalarType = double >
// std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
//            std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >, std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >
// executePlanetRotationSimulation( const Eigen::Matrix< StateScalarType, 12, 1 > initialStateDifference,
//                            Eigen::Matrix< StateScalarType, 12, 1 >& appliedStateDifference,
//                            const bool propagateVariationalEquations = 1 )
// {

//     double initialTime = 0.0;
//     double finalTime = 20.0 * 86400.0; // 10.0 * 3600.0; //86400.0;

//     std::string globalFrameOrigin = "Io";
//     std::string globalFrameOrientation = "J2000";

//     std::vector< std::string > bodiesToCreate = { "Jupiter", "Io" }; 

//     double muJupiter = getBodyGravitationalParameter( "Jupiter" );
//     double muIo = getBodyGravitationalParameter( "Io" );
//     double muEffective = muJupiter + muIo;

//     // Get body settings.
//     BodyListSettings bodySettings =
//             getDefaultBodySettings( bodiesToCreate, initialTime - 86400.0, finalTime + 86400.0, globalFrameOrigin, globalFrameOrientation );

//     // // Set Io's rotation to synchronous        
//     // bodySettings.get( "Io" )->rotationModelSettings = simulation_setup::synchronousRotationModelSettings( "Jupiter", "J2000", "IAU_Io" );

//     // Create system of bodies
//     SystemOfBodies bodies = createSystemOfBodies( bodySettings );

//     // Set Io's mean moment of inertia
//     std::shared_ptr< SphericalHarmonicsGravityField > ioGravityFieldModel = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) );
//     std::shared_ptr< SphericalHarmonicsGravityField > jupiterGravityFieldModel = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Jupiter" )->getGravityFieldModel( ) );
//     double scaledMeanMomentOfInertia = 0.26;
//     jupiterGravityFieldModel->setScaledMeanMomentOfInertia( scaledMeanMomentOfInertia );

//     // Eigen::MatrixXd ioCosineCoefs = ioGravityFieldModel->getCosineCoefficients( );
//     // Eigen::MatrixXd ioSineCoefs = ioGravityFieldModel->getSineCoefficients( );
//     // ioCosineCoefs(2,1) = 1.0e-3;
//     // ioSineCoefs(2,1) = 1.0e-3;
//     // ioSineCoefs(2,2) = 2.0e-3;
//     // ioGravityFieldModel->setCosineCoefficients( ioCosineCoefs );
//     // ioGravityFieldModel->setSineCoefficients( ioSineCoefs );
//     // std::cout << "ioCosineCoefs" << std::endl;
//     // std::cout << ioCosineCoefs << std::endl;
//     // std::cout << "ioSineCoefs" << std::endl;
//     // std::cout << ioSineCoefs << std::endl;

//     // std::shared_ptr< SphericalHarmonicsGravityField > newIoGravityFieldModel = std::make_shared< SphericalHarmonicsGravityField >( 
//         // ioGravityFieldModel->getGravitationalParameter( ), ioGravityFieldModel->getReferenceRadius( ), ioCosineCoefs, ioSineCoefs, "IAU_Io", scaledMeanMomentOfInertia );
//     // bodies.at( "Io" )->setGravityFieldModel( newIoGravityFieldModel );

//     // Retrieve static gravity coefficients
//     Eigen::MatrixXd jupiterCosineCoefs = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Jupiter" )->getGravityFieldModel( ) )->getCosineCoefficients( );
//     Eigen::MatrixXd jupiterSineCoefs = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Jupiter" )->getGravityFieldModel( ) )->getSineCoefficients( );
//     Eigen::Vector5d normalisedStaticGravity = Eigen::Vector5d::Zero( );
//     normalisedStaticGravity[0] = jupiterCosineCoefs(2,0);
//     // normalisedStaticGravity[1] = 1.0e-4; // ioCosineCoefs(2,1);
//     normalisedStaticGravity[2] = jupiterCosineCoefs(2,2);
//     // normalisedStaticGravity[3] = 1.0e-4; // ioSineCoefs(2,1);
//     // normalisedStaticGravity[4] = 1.0e-4; // ioSineCoefs(2,2);
//     std::cout << "normalisedStaticGravity " << normalisedStaticGravity.transpose() << std::endl;
//     bodies.at( "Jupiter" )->setStaticDegreeTwoCoefficients( normalisedStaticGravity );

//     double ioSemiMajorAxis = 4.2e8;
//     // bodies.at( "Io" )
//     //         ->setEphemeris( tudat::ephemerides::getTabulatedEphemeris(
//     //                 std::make_shared< ephemerides::KeplerEphemeris >( ioKeplerElements, 0.0, muEffective, "Jupiter", "J2000" ),
//     //                 initialTime - 3600.0,
//     //                 finalTime + 3600,
//     //                 60.0 ) );

//     // // Set Io's parameters - TO BE MODIFIED
//     // double muIo = 5959924010272.5136719;
//     // double muJupiter = 126686534196012800.0;
//     // double muEffective = 126692494120023072.0;


//     double rotationRateJupiter = 1.76e-4; // std::sqrt( muEffective / ( ioSemiMajorAxis * ioSemiMajorAxis * ioSemiMajorAxis ) );

//     // Retrieve initial rotational state from ephemeris
//     Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRotationState = getInitialRotationalStateOfBody( "Jupiter", "J2000",  bodies, initialTime );
//     // initialRotationState[6] = rotationRateJupiter;

//     // ///// ALTERNATIVE DEFINITION OF IO'S ROTATIONAL STATE
//     // Eigen::Quaterniond noRotationQuaternion = Eigen::Quaterniond( Eigen::AngleAxisd( 1.0E-0, Eigen::Vector3d::UnitZ( ) ) *
//     //                                                               Eigen::AngleAxisd( 2.0E-0, Eigen::Vector3d::UnitX( ) ) *
//     //                                                               Eigen::AngleAxisd( -0.5E-0, Eigen::Vector3d::UnitZ( ) ) );

//     // initialRotationState( 0 ) = noRotationQuaternion.w( );
//     // initialRotationState( 1 ) = noRotationQuaternion.x( );
//     // initialRotationState( 2 ) = noRotationQuaternion.y( );
//     // initialRotationState( 3 ) = noRotationQuaternion.z( );
//     // initialRotationState( 4 ) = 1.0e-5;
//     // initialRotationState( 5 ) = -1.0e-5;
//     // initialRotationState( 6 ) = rotationRateIo;
//     // ////

//     // // initialRotationState = manualInitialState.segment( 0, 7 );

//     Eigen::Matrix< double, 7, 1 > originalRotationState = initialRotationState;
//     Eigen::Matrix< double, 7, 1 > stateDifferenceToAdd = initialStateDifference.segment( 0, 7 );
//     std::cout << "stateDifferenceToAdd " << stateDifferenceToAdd.transpose() << std::endl;

//     initialRotationState += stateDifferenceToAdd;
//     initialRotationState( 0 ) = originalRotationState( 0 ) / std::fabs( originalRotationState( 0 ) ) *
//             std::sqrt( 1.0 - std::pow( initialRotationState.segment( 1, 3 ).norm( ), 2.0 ) );

//     appliedStateDifference.segment( 0, 7 ) = initialRotationState - originalRotationState;
//     appliedStateDifference.segment( 7, 5 ) = initialStateDifference.segment( 7, 5 );
//     std::cout << "appliedStateDifference" << std::endl;
//     std::cout << appliedStateDifference.transpose() << std::endl;
    
    
//     std::map< double, Eigen::Matrix< double, 7, 1 > > dummyRotationMap;
//     dummyRotationMap[ -1.0E100 ] = initialRotationState;
//     dummyRotationMap[ 1.0E100 ] = initialRotationState;

//     std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > dummyInterpolator =
//             std::make_shared< interpolators::LinearInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( dummyRotationMap );
//     bodies.at( "Jupiter" )->setRotationalEphemeris(
//         std::make_shared< TabulatedRotationalEphemeris< double, double > >( dummyInterpolator, "J2000", "IAU_Jupiter" ) );

//     // Retrieve body objects for Io and Jupiter
//     std::shared_ptr< Body > io = bodies.at( "Io" );
//     std::shared_ptr< Body > jupiter = bodies.at( "Jupiter" );

//     // Update Jupiter and Io to current state
//     io->setStateFromEphemeris( initialTime );
//     jupiter->setStateFromEphemeris( initialTime );
//     io->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );
//     jupiter->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );

//     // Check translational intial state
//     Eigen::Vector6d initialTranslationalState = io->getState( ) - jupiter->getState( );
//     std::cout << "initialTranslationalState " << initialTranslationalState.transpose( ) << std::endl;
//     std::cout << "initialRotationState " << initialRotationState.transpose( ) << std::endl;    

//     SelectedTorqueMap torqueMap;
//     torqueMap[ "Jupiter" ][ "Io" ].push_back( std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );


//     // Create integrator settings
//     double timeStep = 100.0 / 2.0;
//     std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > > ( 
//         initialTime, timeStep, rungeKutta87DormandPrince, 1.0e-5, 1000.0, 1.0e-12, 1.0e-12 ); 

//     // Define propagator settings.
//     std::vector< std::string > bodiesToIntegrate;
//     bodiesToIntegrate.push_back( "Jupiter" );

//     // Create torque models
//     basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap( bodies, torqueMap, bodiesToIntegrate );

//     std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
//             std::make_shared< RotationalStatePropagatorSettings< double > >(
//                     torqueModelMap,
//                     bodiesToIntegrate,
//                     initialRotationState,
//                     initialTime,
//                     integratorSettings,
//                     std::make_shared< PropagationTimeTerminationSettings >( finalTime, true ) );

//     // Create gravity deformation model
//     double maxwellRelaxationTime = 179103.0;
//     double globalRelaxationTime = 24688.0;
//     double fluidLoveNumber = 1.5;
//     double radiusJupiter = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Jupiter" )->getGravityFieldModel( ) )->getReferenceRadius( ); 

//     std::vector< std::string > perturbingBody = {"Io"};
//     std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::make_shared< MaxwellDeformationSettings >( 
//         maxwellRelaxationTime, globalRelaxationTime, fluidLoveNumber, 2, 2, perturbingBody, normalisedStaticGravity, true, false );

//     std::map< std::string, std::vector< std::shared_ptr< GravityDeformationSettings > > > gravityDeformationModelMap;   
//     gravityDeformationModelMap[ "Jupiter" ] = { maxwellDeformationSettings };

//     basic_astrodynamics::GravityDeformationModelMap gravityDeformationModels = createGravityDeformationModelsMap(
//         bodies, gravityDeformationModelMap ); 

//     // Compute equilibrium coefficients to initialise gravity state   
//     Eigen::Vector5d equilibriumCoefficients = Eigen::Vector5d::Zero();

//     // Eigen::Vector6d initialTranslationalState = orbital_element_conversions::convertKeplerianToCartesianElements( ioKeplerElements, muEffective );
//     double distanceIo = initialTranslationalState.segment(0, 3).norm();
//     double radiusRatioPowerThree = ( radiusJupiter / distanceIo ) * ( radiusJupiter / distanceIo ) * ( radiusJupiter / distanceIo );
//     double gravitationalParametersRatio = muIo / muJupiter;

//     // std::cout << "translational " << initialTranslationalState.transpose() << std::endl;
//     std::cout << "rotational " << initialRotationState.transpose() << std::endl;
    
//     equilibriumCoefficients[ 0 ] = fluidLoveNumber / 2.0 * gravitationalParametersRatio * radiusRatioPowerThree 
//                 * ( 3.0 * std::sin( 0.0 ) * std::sin( 0.0 ) - 1.0 );
//                 // - fluidLoveNumber * rotationRateIo * rotationRateIo * radiusIo * radiusIo * radiusIo / ( 3.0 * muIo );
//     equilibriumCoefficients[ 2 ] = fluidLoveNumber / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
//                 ( 1.0 - std::sin( 0.0 ) * std::sin( 0.0 ) ) * std::cos( 2.0 * 0.0 );
//     equilibriumCoefficients[ 4 ] = fluidLoveNumber / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
//                 ( 1.0 - std::sin( 0.0 ) * std::sin( 0.0 ) ) * std::sin( 2.0 * 0.0 );
        
//     equilibriumCoefficients[ 1 ] = - fluidLoveNumber * gravitationalParametersRatio * radiusRatioPowerThree 
//                 * ( - std::cos( 0.0 ) * std::sin( 0.0 ) ) * std::cos( 0.0 );
//     equilibriumCoefficients[ 3 ] = - fluidLoveNumber * gravitationalParametersRatio * radiusRatioPowerThree 
//                 * ( - std::cos( 0.0 ) * std::sin( 0.0 ) ) * std::sin( 0.0 );

//     std::cout << "equilibriumCoefficients " << equilibriumCoefficients.transpose() << std::endl;

//     Eigen::Vector5d initialGravityState = equilibriumCoefficients + initialStateDifference.segment(7,5); 
//     // initialGravityState = manualInitialState.segment( 7, 5 ) + initialStateDifference.segment(7,5);    
//     std::cout << "initialGravityState " << initialGravityState.transpose() << std::endl;    

//     std::shared_ptr< GravityDeformationPropagatorSettings< > > gravityPropagatorSettings = 
//         std::make_shared< GravityDeformationPropagatorSettings< > >( bodiesToIntegrate, gravityDeformationModels, initialGravityState, integratorSettings,
//         std::make_shared< PropagationTimeTerminationSettings >( finalTime, true ) );  


//     std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
//     propagatorSettingsList.push_back( rotationalPropagatorSettings );
//     propagatorSettingsList.push_back( gravityPropagatorSettings );

//     // Define dependent variables
//     std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
//     dependentVariablesList.push_back(
//             std::make_shared< SingleDependentVariableSaveSettings >( body_fixed_relative_spherical_position, "Io", "Jupiter" ) );

//     std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings = std::make_shared< MultiTypePropagatorSettings< double > >(
//             propagatorSettingsList, std::make_shared< PropagationTimeTerminationSettings >( finalTime, true ), dependentVariablesList );   
    

//     // Define parameters.
//     std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
//     {
//         parameterNames = getInitialStateParameterSettings< StateScalarType >( propagatorSettings, bodies );
//     }

//     // Create parameters
//     std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
//             createParametersToEstimate( parameterNames, bodies );
//     // printEstimatableParameterEntries( parametersToEstimate );


//     // Perturb parameters.
//     Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > parameterVector =
//             parametersToEstimate->template getFullParameterValues< StateScalarType >( );
//     // parameterVector.block( 12, 0, numberOfParametersToEstimate, 1 ) += parameterPerturbation;
//     //    std::cout<<"Parameter perturbation "<<
//     //               ( parametersToEstimate->template getFullParameterValues< StateScalarType >( ) -
//     //                 parameterVector ).transpose( )<<std::endl<<
//     //               parameterPerturbation.transpose( )<<std::endl;;
//     parametersToEstimate->resetParameterValues( parameterVector );

//     std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
//                std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >, 
//                         std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > > results;

//     {
//         // Create dynamics simulator
//         propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType > dynamicsSimulator =
//                 propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType >(
//                         bodies,
//                         integratorSettings,
//                         propagatorSettings,
//                         parametersToEstimate,
//                         1,
//                         std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
//                         0,
//                         0 );

//         // Propagate requested equations.
//         if( propagateVariationalEquations )
//         {
//             dynamicsSimulator.integrateVariationalAndDynamicalEquations( propagatorSettings->getInitialStates( ), 1 );
//         }
//         else
//         {
//             dynamicsSimulator.integrateDynamicalEquationsOfMotionOnly( propagatorSettings->getInitialStates( ) );
//         }

//             //    tudat::input_output::writeDataMapToTextFile(
//             //                dynamicsSimulator.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( ),
//             //                "rotPropTest.dat" );

//         // Retrieve test data
//         double testEpoch = finalTime; // initialTime + 1.0 * 3600.0;
//         Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > testStates = Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( 12 );
//         testStates.segment( 0, 7 ) = dynamicsSimulator.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( ).rbegin()->second.segment( 0, 7 ); //bodies.at( "Io" )->getRotationalEphemeris( )->getRotationStateVector( testEpoch );
//         testStates.segment( 7, 5 ) = dynamicsSimulator.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( ).rbegin()->second.segment( 7, 5 );
//         std::cout << "final state " << testStates.transpose() << std::endl;\

//         if( propagateVariationalEquations )
//         {
//             // results.first.push_back(
//                     // dynamicsSimulator.getStateTransitionMatrixInterface( )->getCombinedStateTransitionAndSensitivityMatrix( testEpoch ) );
//             results.first = dynamicsSimulator.getStateTransitionMatrixSolution( );
//             Eigen::MatrixXd testMatrixDirect =
//                     dynamicsSimulator.getStateTransitionMatrixInterface( )->getCombinedStateTransitionAndSensitivityMatrix( testEpoch );
//             Eigen::MatrixXd testMatrixFull =
//                     dynamicsSimulator.getStateTransitionMatrixInterface( )->getFullCombinedStateTransitionAndSensitivityMatrix( testEpoch );
//             TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testMatrixDirect, testMatrixFull, std::numeric_limits< double >::epsilon( ) );
//         }
//         // results.second.push_back( testStates );
//         results.second.first = dynamicsSimulator.getEquationsOfMotionSolution( );

//         // Retrieve dependent variables history
//         std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > simulator = dynamicsSimulator.getDynamicsSimulator( );
//         std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > dependentVariablesHistory = simulator->getDependentVariableHistory();
//         results.second.second = dependentVariablesHistory;

//     }
//     return results;
// }

// BOOST_AUTO_TEST_CASE( testDeformationAndRotationalVariationalEquations )
// {
//     // Load spice kernels.
//     spice_interface::loadStandardSpiceKernels( );

//     std::pair< std::map< double, Eigen::MatrixXd >, std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > > currentOutput;

//     // Define variables for numerical differentiation
//     Eigen::Matrix< double, 12, 1 > perturbedState;
//     Eigen::Matrix< double, 12, 1 > statePerturbation;

//     // Define parameter perturbation
//     // int numberOfParametersToEstimate = 8;
//     // double sphericalHarmonicsPerturbation = 1.0E-4;
//     // Eigen::Matrix< double, 8, 1 > perturbedParameter;
//     // Eigen::Matrix< double, 8, 1 > parameterPerturbation;
//     // parameterPerturbation = Eigen::Matrix< double, 8, 1 >::Constant( sphericalHarmonicsPerturbation );
//     // parameterPerturbation( 0 ) = 1.0E-4;

//     // Compute state transition and sensitivity matrices
//     Eigen::Matrix< double, 12, 1 > appliedStateDifference;
//     currentOutput = executePlanetRotationSimulation< double, double >( Eigen::Matrix< double, 12, 1 >::Zero( ), appliedStateDifference );
//     Eigen::MatrixXd stateTransitionAndSensitivityMatrixAtEpoch = currentOutput.first.rbegin()->second; // currentOutput.first.at( 0 );
//     Eigen::VectorXd nominalState = currentOutput.second.first.rbegin()->second; // currentOutput.second.at( 0 );

//     std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > nominalStateInterpolator =
//                 std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( currentOutput.second.first, 6 );

//     std::vector< double > refEpochs;
//     for ( auto it : currentOutput.first )
//     {
//         refEpochs.push_back( it.first );
//     }

//     std::map< double, Eigen::MatrixXd > analyticalStmHistory = currentOutput.first;
//     std::map< double, Eigen::MatrixXd > modifiedAnalyticalStmHistory = analyticalStmHistory;


//     //    std::cout<<"Nominal "<<std::endl<<std::endl<<
//     //               stateTransitionAndSensitivityMatrixAtEpoch<<std::endl;
//     // Define state perturbation
//     statePerturbation = ( Eigen::Matrix< double, 12, 1 >( ) << 
//                           1.0E-8,
//                           1.0E-8,
//                           1.0E-8,
//                           1.0E-8,
//                           1.0E-10,
//                           1.0E-10,
//                           1.0E-9, 
//                           1.0E-6, 
//                           1.0E-6,
//                           1.0E-6,
//                           1.0E-6,
//                           1.0E-6 )
//                                 .finished( );

//     Eigen::MatrixXd manualPartial = Eigen::MatrixXd::Zero( 12, 12 );

//     // Numerically compute state transition matrix
//     std::map< double, Eigen::MatrixXd > numericalPartialsHistory;
//     for ( unsigned int k = 0 ; k < refEpochs.size() ; k++ )
//     {
//         numericalPartialsHistory[ refEpochs[k] ] = Eigen::MatrixXd::Zero( 12, 12 );
//     }

//     for( unsigned int test = 0; test < 1; test++ )
//     {
//         double perturbationMultiplier = ( test == 0 ? 1.0 : 1.0E-3 );
//         for( unsigned int j = 0; j < 4; j++ )
//         {
//             Eigen::Matrix< double, 12, 1 > appliedStateDifferenceUp, appliedStateDifferenceDown;

//             std::map< double, Eigen::VectorXd > upPerturbedStateHistory, downPerturbedStateHistory;
//             Eigen::VectorXd upPerturbedState, downPerturbedState;
//             perturbedState.setZero( );
//             perturbedState( j ) += perturbationMultiplier * statePerturbation( j );
//             upPerturbedStateHistory = executePlanetRotationSimulation< double, double >(
//                                        perturbedState, appliedStateDifferenceUp, 0 ).second.first;
//                                     //    .second.at( 0 );
//             upPerturbedState = upPerturbedStateHistory.rbegin()->second;

//              std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > upPerturbedsStateInterpolator =
//                 std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( upPerturbedStateHistory, 6 );


//             // std::cout << "upPerturbedState " << upPerturbedState.transpose() << std::endl;

//             Eigen::VectorXd stateDifferenceUp = upPerturbedState - nominalState;
//             Eigen::VectorXd stateDifferenceStm = stateTransitionAndSensitivityMatrixAtEpoch * appliedStateDifferenceUp;

//                        std::cout<<"Test output "<<test<<" "<<j<<"stateDifferenceUp"<<std::endl<<
//                                   stateDifferenceUp<<std::endl<<std::endl<<
//                                 //   "stateTransitionAndSensitivityMatrixAtEpoch"<<std::endl<<
//                                                         //  stateTransitionAndSensitivityMatrixAtEpoch<<std::endl<<std::endl<<
//                                 //   "appliedStateDifferenceUp"<<std::endl<<
//                                                         //  appliedStateDifferenceUp<<std::endl<<std::endl<<
//                                   "( stateTransitionAndSensitivityMatrixAtEpoch * appliedStateDifferenceUp )"<<std::endl<<
//                                                          stateDifferenceStm <<std::endl<<std::endl;

//             Eigen::VectorXd relativeErrors = Eigen::VectorXd::Zero(12);
//             for ( unsigned int k = 0 ; k < 12 ; k++ )
//             {
//                 relativeErrors[k] = ( stateDifferenceUp[k] - stateDifferenceStm[k] ) / stateDifferenceUp[k];
//             }

//             // std::cout << "nominalState " << nominalState.transpose() << std::endl;
//             std::cout << "relativeErrors" << std::endl;
//             std::cout << relativeErrors << std::endl;


//             for ( unsigned int k = 0 ; k < refEpochs.size() ; k++ )
//             {
//                 double currentTime = refEpochs[k];
//                 Eigen::VectorXd currentUpState = upPerturbedsStateInterpolator->interpolate(currentTime).segment( 0, 12 );
//                 Eigen::VectorXd currentNominalState = nominalStateInterpolator->interpolate(currentTime).segment( 0, 12 );
//                 Eigen::MatrixXd currentPartial = ( currentUpState - currentNominalState );

//                 numericalPartialsHistory.at( currentTime ).block( 0, j, 12, 1 ) = currentPartial;
//                 modifiedAnalyticalStmHistory.at( currentTime ).block( 0, j, 12, 1 ) = analyticalStmHistory.at( currentTime ) * appliedStateDifferenceUp;
//             }

//         }
//     }

    
//     // PARTIALS WRT GRAVITY DEFORMATION
//     for( unsigned int j = 7; j < 12; j++ )
//     {
//         // std::map< double, Eigen::MatrixXd > partialsWrtGravityDeformationHistory, interpolatedUpPerturbedStateHistory, interpolatedDownPerturbedStateHistory;

//         Eigen::Matrix< double, 12, 1 > appliedStateDifferenceUp, appliedStateDifferenceDown;

//         std::map< double, Eigen::VectorXd > upPerturbedStateHistory, downPerturbedStateHistory;
//         Eigen::VectorXd upPerturbedState, downPerturbedState;
//         perturbedState.setZero( );
//         perturbedState( j ) += statePerturbation( j );
//         upPerturbedStateHistory = executePlanetRotationSimulation< double, double >(
//                                    perturbedState, appliedStateDifferenceUp, 0 ).second.first; //.rbegin()->second;
//                                 //    .second.at( 0 );
//         upPerturbedState = upPerturbedStateHistory.rbegin()->second;

//         std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > upPerturbedStateInterpolator =
//                 std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( upPerturbedStateHistory, 6 );

//         perturbedState.setZero( );
//         perturbedState( j ) -= statePerturbation( j );
//         downPerturbedStateHistory = executePlanetRotationSimulation< double, double >(
//                                      perturbedState, appliedStateDifferenceDown, 0 ).second.first; //.rbegin()->second;
//                                     //  .second.at( 0 );
//         downPerturbedState = downPerturbedStateHistory.rbegin()->second;

//         std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > downPerturbedStateInterpolator =
//                 std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( downPerturbedStateHistory, 6 );

//         for ( unsigned int k = 0 ; k < refEpochs.size() ; k++ )
//         {
//             double currentTime = refEpochs[k];
//             Eigen::VectorXd currentUpState = upPerturbedStateInterpolator->interpolate(currentTime).segment( 0, 12 );
//             Eigen::VectorXd currentDownState = downPerturbedStateInterpolator->interpolate(currentTime).segment( 0, 12 );
//             Eigen::MatrixXd currentPartial = ( currentUpState - currentDownState ) / ( 2.0 * statePerturbation( j ) );

//             numericalPartialsHistory.at( currentTime ).block( 0, j, 12, 1 ) = currentPartial;
//         }

//         manualPartial.block( 0, j, 12, 1 ) =
//                 ( upPerturbedState.segment( 0, 12 ) - downPerturbedState.segment( 0, 12 ) ) / ( 2.0 * statePerturbation( j ) );
//     }


//     // PARTIALS WRT ANGULAR VELOCITY VECTOR
//     for( unsigned int j = 4; j < 7 ; j++ )
//     {
//         Eigen::Matrix< double, 12, 1 > appliedStateDifferenceUp, appliedStateDifferenceDown;

//         std::map< double, Eigen::VectorXd > upPerturbedStateHistory, downPerturbedStateHistory;
//         Eigen::VectorXd upPerturbedState, downPerturbedState;
//         perturbedState.setZero( );
//         perturbedState( j ) += statePerturbation( j );
//         upPerturbedStateHistory = executePlanetRotationSimulation< double, double >(
//                                    perturbedState, appliedStateDifferenceUp, 0 ).second.first; //.rbegin()->second;
//                                 //    .second.at( 0 );
//         upPerturbedState = upPerturbedStateHistory.rbegin()->second;

//         std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > upPerturbedStateInterpolator =
//                 std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( upPerturbedStateHistory, 6 );

//         perturbedState.setZero( );
//         perturbedState( j ) -= statePerturbation( j );
//         downPerturbedStateHistory = executePlanetRotationSimulation< double, double >(
//                                      perturbedState, appliedStateDifferenceDown, 0 ).second.first; //.rbegin()->second;
//                                     //  .second.at( 0 );
//         downPerturbedState = downPerturbedStateHistory.rbegin()->second;

//         std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > downPerturbedStateInterpolator =
//                 std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( downPerturbedStateHistory, 6 );

//         manualPartial.block( 0, j, 12, 1 ) =
//                 ( upPerturbedState.segment( 0, 12 ) - downPerturbedState.segment( 0, 12 ) ) / ( 2.0 * statePerturbation( j ) );

//         for ( unsigned int k = 0 ; k < refEpochs.size() ; k++ )
//         {
//             double currentTime = refEpochs[k];
//             Eigen::VectorXd currentUpState = upPerturbedStateInterpolator->interpolate(currentTime).segment( 0, 12 );
//             Eigen::VectorXd currentDownState = downPerturbedStateInterpolator->interpolate(currentTime).segment( 0, 12 );
//             Eigen::MatrixXd currentPartial = ( currentUpState - currentDownState ) / ( 2.0 * statePerturbation( j ) );

//             numericalPartialsHistory.at( currentTime ).block( 0, j, 12, 1 ) = currentPartial;

//         }
//     }

//     // tudat::input_output::writeDataMapToTextFile(
//     //     numericalPartialsHistory,
//     //     "manualStateTransitionMatrixHistory.dat", "/Users/sam.fayolle/Downloads/" );

//     // tudat::input_output::writeDataMapToTextFile(
//     //     modifiedAnalyticalStmHistory,
//     //     "stateTransitionMatrixHistory.dat", "/Users/sam.fayolle/Downloads/" );
//     // tudat::input_output::writeDataMapToTextFile(
//     //     currentOutput.second.first,
//     //     "refStateHistory.dat", "/Users/sam.fayolle/Downloads/" );
//     // tudat::input_output::writeDataMapToTextFile(
//     //     currentOutput.second.second,
//     //     "refDepVarHistory.dat", "/Users/sam.fayolle/Downloads/" );

//     std::cout << " ----------------------------------------------- " << std::endl;
//     std::cout << " PARTIALS WRT ANGULAR VELOCITY VECTOR " << std::endl;
//     std::cout << "manualPartial " << std::endl;
//     std::cout << manualPartial.block(0, 4, 12, 3) << std::endl;
//     std::cout << "state transition matrix " << std::endl;
//     std::cout << stateTransitionAndSensitivityMatrixAtEpoch.block(0, 4, 12, 3) << std::endl;

//     Eigen::MatrixXd differences = manualPartial.block(0, 4, 12, 3) - stateTransitionAndSensitivityMatrixAtEpoch.block(0, 4, 12, 3);
//     for ( unsigned int j = 0 ; j < 12 ; j++ )
//     {
//         for ( unsigned int k = 4 ; k < 4+3 ; k++ )
//         {
//             if ( stateTransitionAndSensitivityMatrixAtEpoch( j,k ) != 0.0 )
//             {
//                 differences( j,k-4 ) /= stateTransitionAndSensitivityMatrixAtEpoch( j,k );
//             }
//         }
//     }
//     std::cout << "relative differences" << std::endl;
//     std::cout << differences << std::endl;

//     std::cout << " ----------------------------------------------- " << std::endl;
//     std::cout << " PARTIALS WRT GRAVITY STATE DEFORMATION " << std::endl;
//     std::cout << "manualPartial " << std::endl;
//     std::cout << manualPartial.block(0, 7, 12, 5) << std::endl;
//     std::cout << "state transition matrix " << std::endl;
//     std::cout << stateTransitionAndSensitivityMatrixAtEpoch.block(0, 7, 12, 5) << std::endl;

//     Eigen::MatrixXd differences2 = manualPartial.block(0, 7, 12, 5) - stateTransitionAndSensitivityMatrixAtEpoch.block(0, 7, 12, 5);
//     for ( unsigned int j = 0 ; j < 12 ; j++ )
//     {
//         for ( unsigned int k = 7 ; k < 7+5 ; k++ )
//         {
//             if ( stateTransitionAndSensitivityMatrixAtEpoch( j,k ) != 0.0 )
//             {
//                 differences2( j,k-7 ) /= stateTransitionAndSensitivityMatrixAtEpoch( j,k );
//             }
//         }
//     }
//     std::cout << "relative differences" << std::endl;
//     std::cout << differences2 << std::endl;

// }


BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
