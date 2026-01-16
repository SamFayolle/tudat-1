/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <limits>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/ephemerides/keplerEphemeris.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/simulation/estimation.h"
#include "tudat/astro/orbit_determination/gravity_deformation_partials/maxwellDeformationPartial.h"
#include "tudat/math/basic/coordinateConversions.h"

namespace tudat
{
namespace unit_tests
{

// Using declarations.
using namespace tudat::ephemerides;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::propagators;
using namespace tudat::coordinate_conversions;
using namespace tudat::reference_frames;
using namespace tudat::observation_models;
using namespace tudat::estimatable_parameters;
using namespace tudat::orbit_determination;

BOOST_AUTO_TEST_SUITE( test_gravity_deformation_estimation )

Eigen::Matrix3d skew(const Eigen::Vector3d& w)
{
    Eigen::Matrix3d W;
    W <<     0, -w(2),  w(1),
          w(2),     0, -w(0),
         -w(1),  w(0),     0;
    return W;
}

Eigen::VectorXd computeEquilibriumCoefficients( 
        const double muRatio,
        const double kf,
        const double radius,
        const double distance,
        const double longitude,
        const double latitude ) 
 {
        double radiusRatioPowerThree = ( radius / distance ) *  ( radius / distance ) *  ( radius / distance ); 
        Eigen::VectorXd equilibriumCoefficients = Eigen::VectorXd::Zero( 5 );

        equilibriumCoefficients[ 0 ] = kf / 2.0 * muRatio * radiusRatioPowerThree 
            * ( 3.0 * std::sin( latitude ) * std::sin( latitude ) - 1.0 ); 
        equilibriumCoefficients[ 2 ] = kf / 4.0 * muRatio * radiusRatioPowerThree * 
            ( 1.0 - std::sin( latitude ) * std::sin( latitude ) ) * std::cos( 2.0 * longitude );
        equilibriumCoefficients[ 4 ] = kf / 4.0 * muRatio * radiusRatioPowerThree * 
            ( 1.0 - std::sin( latitude ) * std::sin( latitude ) ) * std::sin( 2.0 * longitude );
        equilibriumCoefficients[ 1 ] = - kf * muRatio * radiusRatioPowerThree 
        * ( - std::cos( latitude ) * std::sin( latitude ) ) * std::cos( longitude );
        equilibriumCoefficients[ 3 ] = - kf * muRatio * radiusRatioPowerThree 
        * ( - std::cos( latitude ) * std::sin( latitude ) ) * std::sin( longitude );
        
        return equilibriumCoefficients;
 }

 Eigen::MatrixXd equilibriumCoefficientsPartials( 
        const double kf, 
        const double muRatio,
        const double radius,
        const double currentRadialDistance,
        const double currentLongitude,
        const double currentLatitude )
 {
        Eigen::MatrixXd partial = Eigen::MatrixXd::Zero( 5, 3 );

    double radiusOverDistance = radius / currentRadialDistance;
    double radiusOverDistancePowerThree = radiusOverDistance * radiusOverDistance * radiusOverDistance;

    double c20eq_wrt_r = - 3.0 * kf / 2.0 * muRatio * radiusOverDistancePowerThree / currentRadialDistance 
        * ( 3.0 * std::sin( currentLatitude ) * std::sin( currentLatitude ) - 1.0 );
    double c21eq_wrt_r = - 3.0 * kf * muRatio * radiusOverDistancePowerThree / currentRadialDistance * std::cos( currentLongitude ) 
        * std::cos( currentLatitude ) * std::sin( currentLatitude );
    double c22eq_wrt_r = - 3.0 * kf / 4.0 * muRatio * radiusOverDistancePowerThree / currentRadialDistance * std::cos( 2.0 * currentLongitude )
        * ( 1.0 - std::sin( currentLatitude ) * std::sin( currentLatitude ) );
    double s21eq_wrt_r = - 3.0 * kf * muRatio * radiusOverDistancePowerThree / currentRadialDistance * std::sin( currentLongitude ) 
        * std::cos( currentLatitude ) * std::sin( currentLatitude );
    double s22eq_wrt_r = - 3.0 * kf / 4.0 * muRatio * radiusOverDistancePowerThree / currentRadialDistance * std::sin( 2.0 * currentLongitude )
        * ( 1.0 - std::sin( currentLatitude ) * std::sin( currentLatitude ) );

    double c20eq_wrt_longitude = 0.0;
    double c21eq_wrt_longitude = kf * muRatio * radiusOverDistancePowerThree * ( - std::sin( currentLongitude ) )
        * std::cos( currentLatitude ) * std::sin( currentLatitude );
    double c22eq_wrt_longitude = kf / 4.0 * muRatio * radiusOverDistancePowerThree * ( - 2.0 * std::sin( 2.0 * currentLongitude ) )
        * ( 1.0 - std::sin( currentLatitude ) * std::sin( currentLatitude ) );
    double s21eq_wrt_longitude = kf * muRatio * radiusOverDistancePowerThree * ( std::cos( currentLongitude ) ) 
        * std::cos( currentLatitude ) * std::sin( currentLatitude );
    double s22eq_wrt_longitude = kf / 4.0 * muRatio * radiusOverDistancePowerThree * ( 2.0 * std::cos( 2.0 * currentLongitude ) )
        * ( 1.0 - std::sin( currentLatitude ) * std::sin( currentLatitude ) );

    double c20eq_wrt_latitude = kf / 2.0 * muRatio * radiusOverDistancePowerThree 
        * ( 3.0 * 2.0 * std::sin( currentLatitude ) * std::cos( currentLatitude ) );
    double c21eq_wrt_latitude = kf * muRatio * radiusOverDistancePowerThree 
        * ( std::cos( currentLatitude ) * std::cos( currentLatitude ) - std::sin( currentLatitude ) * std::sin( currentLatitude ) ) * std::cos( currentLongitude );
    double c22eq_wrt_latitude = kf / 4.0 * muRatio * radiusOverDistancePowerThree 
        * ( - 2.0 * std::sin( currentLatitude ) * std::cos( currentLatitude ) ) * std::cos( 2.0 * currentLongitude );
    double s21eq_wrt_latitude = kf * muRatio * radiusOverDistancePowerThree 
        * ( std::cos( currentLatitude ) * std::cos( currentLatitude ) - std::sin( currentLatitude ) * std::sin( currentLatitude ) ) * std::sin( currentLongitude );
    double s22eq_wrt_latitude = kf / 4.0 * muRatio * radiusOverDistancePowerThree 
        * ( - 2.0 * std::sin( currentLatitude ) * std::cos( currentLatitude ) ) * std::sin( 2.0 * currentLongitude );

        partial(0, 0) = c20eq_wrt_r;
        partial(1, 0) = c21eq_wrt_r;
        partial(2, 0) = c22eq_wrt_r;
        partial(3, 0) = s21eq_wrt_r;
        partial(4, 0) = s22eq_wrt_r;

        partial(0, 1) = c20eq_wrt_longitude;
        partial(1, 1) = c21eq_wrt_longitude;
        partial(2, 1) = c22eq_wrt_longitude;
        partial(3, 1) = s21eq_wrt_longitude;
        partial(4, 1) = s22eq_wrt_longitude;

        partial(0, 2) = c20eq_wrt_latitude;
        partial(1, 2) = c21eq_wrt_latitude;
        partial(2, 2) = c22eq_wrt_latitude;
        partial(3, 2) = s21eq_wrt_latitude;
        partial(4, 2) = s22eq_wrt_latitude;

        return partial;
 }

 Eigen::Vector3d getSphericalVelocityComponents( 
    Eigen::Vector3d bodyFixedPosition,
    Eigen::Vector3d bodyFixedVelocity )
 {
    double distance = bodyFixedPosition.segment( 0, 3 ).norm( );
    double distanceDerivative = 
                ( bodyFixedPosition[ 0 ] * bodyFixedVelocity[ 0 ] 
                + bodyFixedPosition[ 1 ] * bodyFixedVelocity[ 1 ] 
                + bodyFixedPosition[ 2 ] * bodyFixedVelocity[ 2 ] ) / distance;

    double longitudeDerivative = ( bodyFixedVelocity[ 1 ] * bodyFixedPosition[ 0 ] - bodyFixedVelocity[ 0 ] * bodyFixedPosition[ 1 ] ) 
                / ( bodyFixedPosition[ 0 ] * bodyFixedPosition[ 0 ] + bodyFixedPosition[ 1 ] * bodyFixedPosition[ 1 ] );
        
    double latitudeDerivative =
        ( bodyFixedVelocity[2] * distance - bodyFixedPosition[2] * distanceDerivative ) /
        ( distance * std::sqrt( bodyFixedPosition[ 0 ] * bodyFixedPosition[ 0 ] + bodyFixedPosition[ 1 ] * bodyFixedPosition[ 1 ] ) );

    Eigen::Vector3d sphericalVelocity = ( Eigen::Vector3d( ) << distanceDerivative, - latitudeDerivative, longitudeDerivative ).finished( );

    return sphericalVelocity;
 }

// BOOST_AUTO_TEST_CASE( test_GravityDeformationEstimation )
// {
//      std::cout.precision( 20 );

//    // Load spice kernels.
//     spice_interface::loadStandardSpiceKernels( );

//     // Specify initial time
//     double initialTime = 0.0;
//     double finalTime = 1.0 * physical_constants::JULIAN_DAY / 10.0; 

//     std::string globalFrameOrigin = "Jupiter";
//     std::string globalFrameOrientation = "J2000";

//     std::vector< std::string > bodiesToCreate = { "Jupiter", "Io" }; 

//     // Get body settings.
//     BodyListSettings bodySettings =
//             getDefaultBodySettings( bodiesToCreate, initialTime - 86400.0, finalTime + 86400.0, globalFrameOrigin, globalFrameOrientation );
            
//         bodySettings.at( "Jupiter" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ), "SSB", globalFrameOrientation ); 
//         // bodySettings.at( "Jupiter" )->gravityFieldSettings = get_gravitational_field( planet, 'IAU_Jupiter' )

//         // Set Jupiter's rotation as constant (no precession)
//         double rightAscensionPole = ( 358.054324066462 - 90.0 ) * mathematical_constants::PI / 180.0;
//         double declinationPole = ( 90.0 - 25.5034135739821 ) * mathematical_constants::PI / 180.0;
//         double primeMeridian = ( 284.95 ) * mathematical_constants::PI / 180.0;
//         double rotationRateJupiter = ( 870.536 * mathematical_constants::PI / 180.0 ) / 86400.0;

//         // bodySettings.get( "Io" )->rotationModelSettings = simulation_setup::synchronousRotationModelSettings( "Jupiter", "J2000", "IAU_Io" );

//         double muIo = 5959924010272.5136719;
//         double muJupiter = 126686534196012800.0;
//         double muEffective = 126692494120023072.0;

//         double orbitalPeriodIo = 2.0 * mathematical_constants::PI * std::sqrt( 4.2e8 * 4.2e8 * 4.2e8 / muIo );
//         double rotationRateIo = std::sqrt( muEffective / ( 4.2e8 * 4.2e8 * 4.2e8 ) );

        
//         Eigen::Vector6d initialKeplerianState = ( Eigen::Vector6d( ) << 4.2e8, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( );
//         // std::shared_ptr< KeplerEphemerisSettings > keplerEphemerisSettings = std::make_shared< KeplerEphemerisSettings >( initialKeplerianState, 0.0, muEffective, "Jupiter", "J2000" );
//         // bodySettings.at( "Io" )->ephemerisSettings = keplerEphemerisSettings;

//         bodySettings.at( "Io" )->ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
//             ( Eigen::Vector6d( ) << 1.0 * 421.8E6, 1.0 * 0.004, 0.0, 0.0, 0.0, 0.0 ).finished( ),
//             0.0,
//             getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Io" ),
//             "Jupiter",
//             "J2000" );

//         // Create bodies needed in simulation
//         SystemOfBodies bodies = createSystemOfBodies( bodySettings );

//         bodies.at("Jupiter")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//                 rightAscensionPole, declinationPole, primeMeridian, rotationRateJupiter, initialTime, globalFrameOrientation, "IAU_Jupiter" ) );

//         bodies.at("Io")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//                 0.0, mathematical_constants::PI / 180.0, 0.0, rotationRateIo, initialTime, globalFrameOrientation, "IAU_Io" ) );

//         double scaledMeanMomentOfInertia = 0.37685;
//         std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )->setScaledMeanMomentOfInertia( scaledMeanMomentOfInertia );

//         double maxwellRelaxationTime = 179103.0;
//         double globalRelaxationTime = 24688.0;
//         double fluidLoveNumber = 1.5;
//         std::vector< std::string > perturbingBody = {"Jupiter"};
//         std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::make_shared< MaxwellDeformationSettings >( 
//                 maxwellRelaxationTime, globalRelaxationTime, fluidLoveNumber, 2, 2, perturbingBody );

    

// //     // orbital period
// //     double orbitalPeriodIo = 2.0 * mathematical_constants::PI * std::sqrt( 4.2e8 * 4.2e8 * 4.2e8 / muIo );
// //     double rotationRateIo = std::sqrt( muEffective / ( 4.2e8 * 4.2e8 * 4.2e8 ) );
// //     std::cout << "rotationRateIo " << rotationRateIo << std::endl;

// // //     Eigen::Matrix3d initialOrientation = Eigen::Matrix3d::Identity( );
// // //     initialOrientation( 0, 0 ) = - 1.0;
// // //     initialOrientation( 1, 1 ) = - 1.0;
// //     bodySettings.get( "Io" )->rotationModelSettings = simulation_setup::synchronousRotationModelSettings( "Jupiter", "J2000", "IAU_Io" );
// // //     bodySettings.at( "Io" )->rotationModelSettings = std::make_shared< simulation_setup::SimpleRotationModelSettings >( 
// //         // "J2000", "IAU_Io", Eigen::Quaterniond( initialOrientation ), initialTime, rotationRateIo );

// //     // Create bodies needed in simulation
// //     SystemOfBodies bodies = createSystemOfBodies( bodySettings );

//     bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );
//     bodies.at( "Io" )->setStateFromEphemeris<>( initialTime );
//     bodies.at( "Jupiter" )->setStateFromEphemeris<>( initialTime );
//     // bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );

// //     double scaledMeanMomentOfInertiaIo = 0.37685;
// //     std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( 
// //         bodies.at( "Io" )->getGravityFieldModel( ) )->setScaledMeanMomentOfInertia( scaledMeanMomentOfInertiaIo );
// //     std::cout << "initial inertia tensor " << std::endl;
// //     std::cout << bodies.at( "Io" )->getGravityFieldModel( )->getInertiaTensor( ) << std::endl;

// //     const double maxwellRelaxationTime = 1000.0;
// //     const double globalRelaxationTime = 2000.0;
// //     const double loveNumber = 0.4;
// //     // const double rotationRate = ( 2.0 * mathematical_constants::PI ) / ( 1.77 * 86400.0 );
// //     const int maximumDegree = 2;
// //     const int maximumOrder = 2;

// //     // std::cout << "initial rotation rate " << rotationRate << std::endl;

// //     std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::make_shared< MaxwellDeformationSettings >( 
// //         maxwellRelaxationTime, globalRelaxationTime, loveNumber, /*rotationRateIo,*/ maximumDegree, maximumOrder, "Jupiter" );

// //     // std::shared_ptr< Body > deformingBody = bodies.at( "Io" );
// //     // std::shared_ptr< Body > perturbingBody = bodies.at( "Jupiter" );
// //     // std::shared_ptr< basic_astrodynamics::MaxwellGravityDeformationModel > maxwellDeformationModel = 
// //     //     createMaxwellGravityFieldDeformationModel( deformingBody, perturbingBody, "Io", "Jupiter", maxwellDeformationSettings );

// //     // std::map< std::string, std::shared_ptr< basic_astrodynamics::GravityDeformationModel > > gravityDeformationModels;
// //     // gravityDeformationModels[ "Io" ] = maxwellDeformationModel;


//     std::map< std::string, std::vector< std::shared_ptr< GravityDeformationSettings > > > gravityDeformationModelMap;   
//     gravityDeformationModelMap[ "Io" ] = { maxwellDeformationSettings };

//     basic_astrodynamics::GravityDeformationModelMap gravityDeformationModels = createGravityDeformationModelsMap(
//         bodies, gravityDeformationModelMap );

// //     std::map< std::string, std::shared_ptr< basic_astrodynamics::GravityDeformationModel > > gravityDeformationModels;
// //     gravityDeformationModels[ "Io" ] = deformationModels.at( "Io" )[ 0 ];



//     double timeStep = 10.0;
//     std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > > ( 
//         initialTime, timeStep, rungeKutta87DormandPrince, timeStep, timeStep );  

//     std::vector< std::string > bodiesToPropagate = { "Io" };
// //     // Eigen::Matrix< double, Eigen::Dynamic, 1 > initialBodyGravity = Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( 3, 1 );

// //     Eigen::Vector3d computedEquilibriumCoefficients = Eigen::Vector3d::Zero( );
// //     std::shared_ptr< SphericalHarmonicsGravityField > shModel = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) );
// //     double distance = 4.2e8;
// //     double radius = shModel->getReferenceRadius( );
// //     double muJup = bodies.at( "Jupiter" )->getGravitationalParameter( );  
// //     double ratioDistanceRadiusPowerThree = radius * radius * radius / ( distance * distance * distance );
// //     double muRatio = muJup / muIo;
// //     computedEquilibriumCoefficients[ 0 ] = - loveNumber * ( rotationRateIo * rotationRateIo * radius * radius * radius / ( 3.0 * muIo ) + 0.5 * muRatio * ratioDistanceRadiusPowerThree );
// //     computedEquilibriumCoefficients[ 1 ] = loveNumber / 4 * muRatio * ratioDistanceRadiusPowerThree;
// //     computedEquilibriumCoefficients[ 2 ] = 0.0;
// //     std::cout << "computedEquilibriumCoefficients " << computedEquilibriumCoefficients.transpose( ) << std::endl;

// //     Eigen::MatrixXd originalCosineMatrix = shModel->getCosineCoefficients();
// //     Eigen::MatrixXd originalSineMatrix = shModel->getSineCoefficients();
// //     Eigen::Vector3d computedInitialCoefficients = Eigen::Vector3d::Zero( );
// //     computedInitialCoefficients[ 0 ] = originalCosineMatrix( 2, 0 );
// //     computedInitialCoefficients[ 1 ] = originalCosineMatrix( 2, 2 );
// //     computedInitialCoefficients[ 2 ] = originalSineMatrix( 2, 2 );
// //     std::cout << "computedInitialCoefficients " << computedInitialCoefficients.transpose( ) << std::endl;

// //     Eigen::Vector3d computedInitialTransientCoefficients = Eigen::Vector3d::Zero( );
// //     computedInitialTransientCoefficients = ( 1.0 / ( globalRelaxationTime - maxwellRelaxationTime ) ) * 
// //         ( globalRelaxationTime * computedInitialCoefficients - maxwellRelaxationTime * computedEquilibriumCoefficients );
// //     // std::cout << "computedInitialTransientCoefficients: " << computedInitialTransientCoefficients.transpose( ) << std::endl;

//         Eigen::VectorXd computedInitialCoefficients = Eigen::VectorXd::Zero( 5 );
//         Eigen::VectorXd perturbedInitialCoefficients = Eigen::VectorXd::Zero( 5 );
//         perturbedInitialCoefficients[0] = 1.0e-2;
//         perturbedInitialCoefficients[2] = 1.0e-2;
//         perturbedInitialCoefficients[4] = 1.0e-2;

//     std::shared_ptr< GravityDeformationPropagatorSettings< > > gravityPropagatorSettings = 
//         std::make_shared< GravityDeformationPropagatorSettings< > >( bodiesToPropagate, gravityDeformationModels, computedInitialCoefficients, integratorSettings,
//         std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

//     // Translational dynamics propagator
//     std::vector< std::string > centralBodies = { "Jupiter" };
//     SelectedAccelerationMap accelerationSettingsMap;
//     accelerationSettingsMap[ "Io" ][ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

//     std::shared_ptr< SingleArcPropagatorProcessingSettings > outputSettings =
//             std::make_shared< SingleArcPropagatorProcessingSettings >( );
//     outputSettings->setIntegratedResult( false );

//     AccelerationMap accelerationsMap = createAccelerationModelsMap( bodies, accelerationSettingsMap, bodiesToPropagate, centralBodies );
//     Eigen::Vector6d initialState = orbital_element_conversions::convertKeplerianToCartesianElements( initialKeplerianState, muEffective );
//     std::shared_ptr< TranslationalStatePropagatorSettings<  > > translationalPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >( 
//         centralBodies, accelerationsMap, bodiesToPropagate, initialState, initialTime, integratorSettings, 
//         std::make_shared< PropagationTimeTerminationSettings >( finalTime ), cowell );

// //     // Create torque models
// //     SelectedTorqueMap torqueSettings;
// //     torqueSettings[ "Io" ][ "Jupiter" ].push_back(
// //                             std::make_shared< SphericalHarmonicTorqueSettings >( 2, 2 ) );
// //     basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap( bodies, torqueSettings, bodiesToPropagate );

// //     Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRotationState = getInitialRotationalStateOfBody(
// //         "Io", "J2000",  bodies, initialTime );
// //     // std::cout << "initialRotationState " << initialRotationState << std::endl;

// //     // Create propagator settings for rotational dynamics
// //     std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
// //             std::make_shared< RotationalStatePropagatorSettings< double > >
// //             ( torqueModelMap, bodiesToPropagate, initialRotationState, initialTime, integratorSettings, 
// //             std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

//     std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > >  propagatorSettingsList;
//     propagatorSettingsList.push_back( translationalPropagatorSettings );
// //     propagatorSettingsList.push_back( rotationalPropagatorSettings );
//     propagatorSettingsList.push_back( gravityPropagatorSettings );
//     std::shared_ptr< MultiTypePropagatorSettings< > > fullPropagatorSettings = std::make_shared< MultiTypePropagatorSettings< > >(
//             propagatorSettingsList, integratorSettings, initialTime, std::make_shared< PropagationTimeTerminationSettings >( finalTime ),
//             std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ), outputSettings );
    

//     std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
//             getInitialStateParameterSettings< double, double >( fullPropagatorSettings, bodies );
//     // Create parameters
//     std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
//             createParametersToEstimate< double, double >( parameterNames, bodies );
// printEstimatableParameterEntries( parametersToEstimate );

// std::shared_ptr< SingleArcVariationalEquationsSolver< double, double > > variationalEquationsSolver =
//             std::make_shared< SingleArcVariationalEquationsSolver< double, double > >(
//                     bodies, fullPropagatorSettings, parametersToEstimate, true );

// std::map< double, Eigen::VectorXd > results = variationalEquationsSolver->getEquationsOfMotionSolution();
// std::cout << "initial state " << results.begin()->second.transpose() << std::endl;
//     std::cout << "final state " << results.rbegin()->second.transpose() << std::endl;
// // for ( auto it : results )
// // {
// // std::cout << it.second.transpose( ) << std::endl;
// // }

// std::map< double, Eigen::MatrixXd > stateTransitionMatrixHistory = variationalEquationsSolver->getStateTransitionMatrixSolution( );
// std::cout << "stateTransitionMatrixHistory size " << stateTransitionMatrixHistory.size( ) << std::endl;

// std::cout << "initial STM " << std::endl;
// std::cout << stateTransitionMatrixHistory.begin( )->second << std::endl;
// std::cout << "-----------------------" << std::endl;
// std::cout << "final STM " << std::endl;
// std::cout << stateTransitionMatrixHistory.rbegin( )->second << std::endl;

// std::cout << "test " << std::endl;
// Eigen::VectorXd testPerturbation = Eigen::VectorXd::Zero(11);
// testPerturbation[0] = 0.001 * initialState[0];
// std::cout << "testPerturbation " << testPerturbation.transpose() << std::endl;
// std::cout << ( stateTransitionMatrixHistory.rbegin( )->second * testPerturbation ).transpose( ) << std::endl;

// Eigen::Vector6d perturbedInitialState = initialState;
// perturbedInitialState[0] *= 1.001;

// Eigen::VectorXd fullPerturbedInitialState = Eigen::VectorXd::Zero( 11 );
// fullPerturbedInitialState.segment( 0, 6 ) = perturbedInitialState;
// fullPerturbedInitialState.segment( 6, 5 ) = computedInitialCoefficients;

// std::cout << "fullPerturbedInitialState " << fullPerturbedInitialState.transpose( ) << std::endl;

// fullPropagatorSettings->resetInitialStates( fullPerturbedInitialState );
// SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, fullPropagatorSettings ); 
// std::map< double, Eigen::VectorXd > results2 = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

// std::cout << "initial state " << results2.begin()->second.transpose() << std::endl;
// std::cout << "final state " << results2.rbegin()->second.transpose() << std::endl;
// // for ( auto it : results2 )
// // {
// // std::cout << it.second.transpose( ) << std::endl;
// // }

// std::map< double, Eigen::VectorXd > stateVariation;
// for ( auto it : results )
// {
//         stateVariation[ it.first ] = ( results2.at( it.first ) - it.second );
// }

// std::cout << "final state variation " << stateVariation.rbegin()->second.transpose( ) << std::endl;
// // for ( auto it : stateVariation )
// // {
// //         std::cout << it.second.transpose( ) << std::endl;
// // }

// // TEST COMPUTATION EQUILIBRIUM COEFFICIENTS

// double radius = 1800.0e3;
// double distance = 4.2e8;
// double latitude = mathematical_constants::PI / 180.0 * 2.0;
// double longitude = mathematical_constants::PI / 180.0 * 15.0;

// Eigen::VectorXd equilibriumCoefficients = computeEquilibriumCoefficients( 
//         muJupiter / muIo, fluidLoveNumber, radius, distance, longitude, latitude );
// std::cout << "equilibriumCoefficients " << equilibriumCoefficients.transpose() << std::endl;

// double perturbation_r = 0.001 * distance;
// double perturbation_longitude = 0.01 * longitude;
// double perturbation_latitude = 0.01 * latitude;

// Eigen::VectorXd perturbedEquilibriumCoefficients_r = computeEquilibriumCoefficients( 
//         muJupiter / muIo, fluidLoveNumber, radius, distance + perturbation_r, longitude, latitude );
// Eigen::VectorXd perturbedEquilibriumCoefficients_longitude = computeEquilibriumCoefficients( 
//         muJupiter / muIo, fluidLoveNumber, radius, distance, longitude + perturbation_longitude, latitude );
// Eigen::VectorXd perturbedEquilibriumCoefficients_latitude = computeEquilibriumCoefficients( 
//         muJupiter / muIo, fluidLoveNumber, radius, distance, longitude, latitude + perturbation_latitude );

// Eigen::VectorXd variation_r = ( perturbedEquilibriumCoefficients_r - equilibriumCoefficients ) / perturbation_r;
// std::cout << "variation - r" << variation_r.transpose( ) << std::endl;
// Eigen::VectorXd variation_longitude = ( perturbedEquilibriumCoefficients_longitude - equilibriumCoefficients ) / perturbation_longitude;
// std::cout << "variation - longitude" << variation_longitude.transpose( ) << std::endl;
// Eigen::VectorXd variation_latitude = ( perturbedEquilibriumCoefficients_latitude - equilibriumCoefficients ) / perturbation_latitude;
// std::cout << "variation - latitude" << variation_latitude.transpose( ) << std::endl;

//  Eigen::MatrixXd partial = equilibriumCoefficientsPartials( 
//         fluidLoveNumber, muJupiter / muIo, radius, distance, longitude, latitude );
// std::cout << "partial " << partial.transpose( ) << std::endl;


// Eigen::Vector3d initialPosition = initialState.segment( 0, 3 );
// Eigen::Vector3d sphericalPosition = 
//                 coordinate_conversions::convertCartesianToSpherical( initialPosition );

// Eigen::Vector3d perturbed_x = initialPosition;
// perturbed_x[0] += 1000.0;
// Eigen::Vector3d perturbedSphericalPosition_x = 
//                 coordinate_conversions::convertCartesianToSpherical( perturbed_x );
// Eigen::Vector3d variationSphericalPosition_x = (perturbedSphericalPosition_x - sphericalPosition) / 1000.0;
// std::cout << "variationSphericalPosition_x : " << variationSphericalPosition_x.transpose( ) << std::endl;

// Eigen::Vector3d perturbed_y = initialPosition;
// perturbed_y[1] += 1000.0;
// Eigen::Vector3d perturbedSphericalPosition_y = 
//                 coordinate_conversions::convertCartesianToSpherical( perturbed_y );
// Eigen::Vector3d variationSphericalPosition_y = (perturbedSphericalPosition_y - sphericalPosition) / 1000.0;
// std::cout << "variationSphericalPosition_y : " << variationSphericalPosition_y.transpose( ) << std::endl;

// Eigen::Vector3d perturbed_z = initialPosition;
// perturbed_z[2] += 1000.0;
// Eigen::Vector3d perturbedSphericalPosition_z = 
//                 coordinate_conversions::convertCartesianToSpherical( perturbed_z );
// Eigen::Vector3d variationSphericalPosition_z = (perturbedSphericalPosition_z - sphericalPosition) / 1000.0;
// std::cout << "variationSphericalPosition_z : " << variationSphericalPosition_z.transpose( ) << std::endl;

// Eigen::Matrix3d partialSphericalPosition = acceleration_partials::computeSphericalJacobian(initialPosition);
// std::cout << "partialSphericalPosition " << partialSphericalPosition.transpose( ) << std::endl;

// // Test full partial of body-fixed spherical state wrt body-fixed cartesian state
// bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );
// bodies.at( "Io" )->setStateFromEphemeris<>( initialTime );
// bodies.at( "Jupiter" )->setStateFromEphemeris<>( initialTime );

// Eigen::Vector6d currentStateIo;
// bodies.at( "Io" )->getStateByReference( currentStateIo );
// std::cout << "currentStateIo " << currentStateIo.transpose() << std::endl;

// Eigen::Vector3d bodyFixedAngularVelocity = bodies.at( "Io" )->getCurrentAngularVelocityVectorInLocalFrame( );
// // std::cout << "bodyFixedAngularVelocity " << bodyFixedAngularVelocity.transpose( ) << std::endl;
// Eigen::Matrix3d rotationMatrixToLocalFrame = ( bodies.at( "Io" )->getCurrentRotationToGlobalFrame( ) ).toRotationMatrix( ).transpose( );
// // std::cout << "rotationMatrixToLocalFrame " << rotationMatrixToLocalFrame << std::endl;
// Eigen::Matrix3d rotationMatrixDerivativeToLocalFrame = bodies.at( "Io" )->getCurrentRotationMatrixDerivativeToLocalFrame( );
// // std::cout << "rotationMatrixDerivativeToLocalFrame " << rotationMatrixDerivativeToLocalFrame << std::endl;

// Eigen::Vector3d bodyFixedCartesianPosition = rotationMatrixToLocalFrame * currentStateIo.segment( 0, 3 );
// Eigen::Vector3d bodyFixedCartesianVelocity = rotationMatrixToLocalFrame * currentStateIo.segment( 3, 3 )
//                 + rotationMatrixDerivativeToLocalFrame * currentStateIo.segment( 0, 3 );

// Eigen::Matrix6d fullSphericalWrtCartesianStatePartials;
// acceleration_partials::computeFullSphericalStatePartials( bodyFixedCartesianPosition, 
//                                    bodyFixedCartesianVelocity,
//                                    bodyFixedAngularVelocity,
//                                    fullSphericalWrtCartesianStatePartials );

// std::cout << "computeFullSphericalStatePartials " << std::endl;
// std::cout << fullSphericalWrtCartesianStatePartials << std::endl;

// // Perturb body-fixed cartesian position
// double deltaPosition = 1000.0;
// Eigen::Vector3d perturbedBodyFixedCartesianPositionX = bodyFixedCartesianPosition;
// perturbedBodyFixedCartesianPositionX[0] += deltaPosition;
// Eigen::Vector3d perturbedBodyFixedCartesianPositionY = bodyFixedCartesianPosition;
// perturbedBodyFixedCartesianPositionY[1] += deltaPosition;
// Eigen::Vector3d perturbedBodyFixedCartesianPositionZ = bodyFixedCartesianPosition;
// perturbedBodyFixedCartesianPositionZ[2] += deltaPosition;

// // Perturb body-fixed cartesian velocity
// double deltaVelocity = 1.0;
// Eigen::Vector3d perturbedBodyFixedCartesianVelocityX = bodyFixedCartesianVelocity;
// perturbedBodyFixedCartesianVelocityX[0] += deltaVelocity;
// Eigen::Vector3d perturbedBodyFixedCartesianVelocityY = bodyFixedCartesianVelocity;
// perturbedBodyFixedCartesianVelocityY[1] += deltaVelocity;
// Eigen::Vector3d perturbedBodyFixedCartesianVelocityZ = bodyFixedCartesianVelocity;
// perturbedBodyFixedCartesianVelocityZ[2] += deltaVelocity;

// // Compute nominal body-fixed spherical position and velocity
// Eigen::Vector3d nominalSphericalPosition = coordinate_conversions::convertCartesianToSpherical( bodyFixedCartesianPosition );
// Eigen::Vector3d nominalSphericalVelocity = getSphericalVelocityComponents( bodyFixedCartesianPosition, bodyFixedCartesianVelocity );

// Eigen::Matrix6d numericalPartial = Eigen::Matrix6d::Zero();

// // Perturb body-fixed cartesian position
// Eigen::Vector3d perturbedSphericalPositionX = coordinate_conversions::convertCartesianToSpherical( perturbedBodyFixedCartesianPositionX );
// Eigen::Vector3d variationSphericalPositionX = (perturbedSphericalPositionX - nominalSphericalPosition) / deltaPosition;

// Eigen::Vector3d perturbedSphericalPositionY = coordinate_conversions::convertCartesianToSpherical( perturbedBodyFixedCartesianPositionY );
// Eigen::Vector3d variationSphericalPositionY = (perturbedSphericalPositionY - nominalSphericalPosition) / deltaPosition;

// Eigen::Vector3d perturbedSphericalPositionZ = coordinate_conversions::convertCartesianToSpherical( perturbedBodyFixedCartesianPositionZ );
// Eigen::Vector3d variationSphericalPositionZ = (perturbedSphericalPositionZ - nominalSphericalPosition) / deltaPosition;

// Eigen::Vector3d perturbedSphericalVelocityX = getSphericalVelocityComponents( perturbedBodyFixedCartesianPositionX, bodyFixedCartesianVelocity );
// Eigen::Vector3d variationSphericalVelocityX = (perturbedSphericalVelocityX - nominalSphericalVelocity) / deltaPosition;

// Eigen::Vector3d perturbedSphericalVelocityY = getSphericalVelocityComponents( perturbedBodyFixedCartesianPositionY, bodyFixedCartesianVelocity );
// Eigen::Vector3d variationSphericalVelocityY = (perturbedSphericalVelocityY - nominalSphericalVelocity) / deltaPosition;

// Eigen::Vector3d perturbedSphericalVelocityZ = getSphericalVelocityComponents( perturbedBodyFixedCartesianPositionZ, bodyFixedCartesianVelocity );
// Eigen::Vector3d variationSphericalVelocityZ = (perturbedSphericalVelocityZ - nominalSphericalVelocity) / deltaPosition;

// numericalPartial.block(0,0,3,1) = variationSphericalPositionX;
// numericalPartial.block(0,1,3,1) = variationSphericalPositionY;
// numericalPartial.block(0,2,3,1) = variationSphericalPositionZ;

// numericalPartial.block(3,0,3,1) = variationSphericalVelocityX;
// numericalPartial.block(3,1,3,1) = variationSphericalVelocityY;
// numericalPartial.block(3,2,3,1) = variationSphericalVelocityZ;

// // Perturb body-fixed cartesian velocity
// perturbedSphericalVelocityX = getSphericalVelocityComponents( bodyFixedCartesianPosition, perturbedBodyFixedCartesianVelocityX );
// variationSphericalVelocityX = (perturbedSphericalVelocityX - nominalSphericalVelocity) / deltaVelocity;

// perturbedSphericalVelocityY = getSphericalVelocityComponents( bodyFixedCartesianPosition, perturbedBodyFixedCartesianVelocityY );
// variationSphericalVelocityY = (perturbedSphericalVelocityY - nominalSphericalVelocity) / deltaVelocity;

// perturbedSphericalVelocityZ = getSphericalVelocityComponents( bodyFixedCartesianPosition, perturbedBodyFixedCartesianVelocityZ );
// variationSphericalVelocityZ = (perturbedSphericalVelocityZ - nominalSphericalVelocity) / deltaVelocity;

// numericalPartial.block(3,3,3,1) = variationSphericalVelocityX;
// numericalPartial.block(3,4,3,1) = variationSphericalVelocityY;
// numericalPartial.block(3,5,3,1) = variationSphericalVelocityZ;
// std::cout << "numericalPartial " << std::endl;
// std::cout << numericalPartial << std::endl;

// //  Eigen::Vector3d getSphericalVelocityComponents( 
// //     Eigen::Vector3d bodyFixedPosition,
// //     Eigen::Vector3d bodyFixedVelocity )

// }

// BOOST_AUTO_TEST_CASE( test_GravityDeformationEstimation_SynchronousRotation )
// {
//      std::cout.precision( 20 );

//    // Load spice kernels.
//     spice_interface::loadStandardSpiceKernels( );

//     // Specify initial time
//     double initialTime = 0.0;
//     double finalTime = 1.0 * physical_constants::JULIAN_DAY / 10.0; 

//     std::string globalFrameOrigin = "Jupiter";
//     std::string globalFrameOrientation = "J2000";

//     std::vector< std::string > bodiesToCreate = { "Jupiter", "Io" }; 

//     // Get body settings.
//     BodyListSettings bodySettings =
//             getDefaultBodySettings( bodiesToCreate, initialTime - 86400.0, finalTime + 86400.0, globalFrameOrigin, globalFrameOrientation );
            
//         bodySettings.at( "Jupiter" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ), "SSB", globalFrameOrientation ); 
//         // bodySettings.at( "Jupiter" )->gravityFieldSettings = get_gravitational_field( planet, 'IAU_Jupiter' )

//         // Set Jupiter's rotation as constant (no precession)
//         double rightAscensionPole = ( 358.054324066462 - 90.0 ) * mathematical_constants::PI / 180.0;
//         double declinationPole = ( 90.0 - 25.5034135739821 ) * mathematical_constants::PI / 180.0;
//         double primeMeridian = ( 284.95 ) * mathematical_constants::PI / 180.0;
//         double rotationRateJupiter = ( 870.536 * mathematical_constants::PI / 180.0 ) / 86400.0;

//         bodySettings.get( "Io" )->rotationModelSettings = simulation_setup::synchronousRotationModelSettings( "Jupiter", "J2000", "IAU_Io" );

//         double muIo = 5959924010272.5136719;
//         double muJupiter = 126686534196012800.0;
//         double muEffective = 126692494120023072.0;

//         double orbitalPeriodIo = 2.0 * mathematical_constants::PI * std::sqrt( 4.2e8 * 4.2e8 * 4.2e8 / muIo );
//         double rotationRateIo = std::sqrt( muEffective / ( 4.2e8 * 4.2e8 * 4.2e8 ) );

        
//         Eigen::Vector6d initialKeplerianState = ( Eigen::Vector6d( ) << 4.2e8, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( );
//         // std::shared_ptr< KeplerEphemerisSettings > keplerEphemerisSettings = std::make_shared< KeplerEphemerisSettings >( initialKeplerianState, 0.0, muEffective, "Jupiter", "J2000" );
//         // bodySettings.at( "Io" )->ephemerisSettings = keplerEphemerisSettings;

//         bodySettings.at( "Io" )->ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
//             ( Eigen::Vector6d( ) << 1.0 * 421.8E6, 1.0 * 0.004, 0.0, 0.0, 0.0, 0.0 ).finished( ),
//             0.0,
//             getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Io" ),
//             "Jupiter",
//             "J2000" );

//         // Create bodies needed in simulation
//         SystemOfBodies bodies = createSystemOfBodies( bodySettings );

//         bodies.at("Jupiter")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//                 rightAscensionPole, declinationPole, primeMeridian, rotationRateJupiter, initialTime, globalFrameOrientation, "IAU_Jupiter" ) );

//         // bodies.at("Io")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//                 // 0.0, mathematical_constants::PI / 180.0, 0.0, rotationRateIo, initialTime, globalFrameOrientation, "IAU_Io" ) );

//         double scaledMeanMomentOfInertia = 0.37685;
//         std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )->setScaledMeanMomentOfInertia( scaledMeanMomentOfInertia );

//         double maxwellRelaxationTime = 179103.0;
//         double globalRelaxationTime = 24688.0;
//         double fluidLoveNumber = 1.5;
//         std::vector< std::string > perturbingBody = {"Jupiter"};
//         std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::make_shared< MaxwellDeformationSettings >( 
//                 maxwellRelaxationTime, globalRelaxationTime, fluidLoveNumber, 2, 2, perturbingBody );

//     bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );
//     bodies.at( "Io" )->setStateFromEphemeris<>( initialTime );
//     bodies.at( "Jupiter" )->setStateFromEphemeris<>( initialTime );
//     // bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );


//     std::map< std::string, std::vector< std::shared_ptr< GravityDeformationSettings > > > gravityDeformationModelMap;   
//     gravityDeformationModelMap[ "Io" ] = { maxwellDeformationSettings };

//     basic_astrodynamics::GravityDeformationModelMap gravityDeformationModels = createGravityDeformationModelsMap(
//         bodies, gravityDeformationModelMap );



//     double timeStep = 10.0;
//     std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > > ( 
//         initialTime, timeStep, rungeKutta87DormandPrince, timeStep, timeStep );  

//     std::vector< std::string > bodiesToPropagate = { "Io" };

//         Eigen::VectorXd computedInitialCoefficients = Eigen::VectorXd::Zero( 5 );
//         Eigen::VectorXd perturbedInitialCoefficients = Eigen::VectorXd::Zero( 5 );
//         perturbedInitialCoefficients[0] = 1.0e-2;
//         perturbedInitialCoefficients[2] = 1.0e-2;
//         perturbedInitialCoefficients[4] = 1.0e-2;

//     std::shared_ptr< GravityDeformationPropagatorSettings< > > gravityPropagatorSettings = 
//         std::make_shared< GravityDeformationPropagatorSettings< > >( bodiesToPropagate, gravityDeformationModels, computedInitialCoefficients, integratorSettings,
//         std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

//     // Translational dynamics propagator
//     std::vector< std::string > centralBodies = { "Jupiter" };
//     SelectedAccelerationMap accelerationSettingsMap;
//     accelerationSettingsMap[ "Io" ][ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

//     std::shared_ptr< SingleArcPropagatorProcessingSettings > outputSettings =
//             std::make_shared< SingleArcPropagatorProcessingSettings >( );
//     outputSettings->setIntegratedResult( false );

//     AccelerationMap accelerationsMap = createAccelerationModelsMap( bodies, accelerationSettingsMap, bodiesToPropagate, centralBodies );
//     Eigen::Vector6d initialState = orbital_element_conversions::convertKeplerianToCartesianElements( initialKeplerianState, muEffective );
//     std::shared_ptr< TranslationalStatePropagatorSettings<  > > translationalPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >( 
//         centralBodies, accelerationsMap, bodiesToPropagate, initialState, initialTime, integratorSettings, 
//         std::make_shared< PropagationTimeTerminationSettings >( finalTime ), cowell );

//     std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > >  propagatorSettingsList;
//     propagatorSettingsList.push_back( translationalPropagatorSettings );
//     propagatorSettingsList.push_back( gravityPropagatorSettings );
//     std::shared_ptr< MultiTypePropagatorSettings< > > fullPropagatorSettings = std::make_shared< MultiTypePropagatorSettings< > >(
//             propagatorSettingsList, integratorSettings, initialTime, std::make_shared< PropagationTimeTerminationSettings >( finalTime ),
//             std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ), outputSettings );
    

//     std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
//             getInitialStateParameterSettings< double, double >( fullPropagatorSettings, bodies );
//     // Create parameters
//     std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
//             createParametersToEstimate< double, double >( parameterNames, bodies );
// printEstimatableParameterEntries( parametersToEstimate );

// std::shared_ptr< SingleArcVariationalEquationsSolver< double, double > > variationalEquationsSolver =
//             std::make_shared< SingleArcVariationalEquationsSolver< double, double > >(
//                     bodies, fullPropagatorSettings, parametersToEstimate, true );

// std::map< double, Eigen::VectorXd > results = variationalEquationsSolver->getEquationsOfMotionSolution();
// std::cout << "initial state " << results.begin()->second.transpose() << std::endl;
//     std::cout << "final state " << results.rbegin()->second.transpose() << std::endl;
// // for ( auto it : results )
// // {
// // std::cout << it.second.transpose( ) << std::endl;
// // }

// std::map< double, Eigen::MatrixXd > stateTransitionMatrixHistory = variationalEquationsSolver->getStateTransitionMatrixSolution( );
// std::cout << "stateTransitionMatrixHistory size " << stateTransitionMatrixHistory.size( ) << std::endl;

// std::cout << "initial STM " << std::endl;
// std::cout << stateTransitionMatrixHistory.begin( )->second << std::endl;
// std::cout << "-----------------------" << std::endl;
// std::cout << "final STM " << std::endl;
// std::cout << stateTransitionMatrixHistory.rbegin( )->second << std::endl;

// std::cout << "test " << std::endl;
// Eigen::VectorXd testPerturbation = Eigen::VectorXd::Zero(11);
// testPerturbation[0] = 0.001 * initialState[0];
// testPerturbation[4] = 0.001 * initialState[4];
// std::cout << "testPerturbation " << testPerturbation.transpose() << std::endl;
// std::cout << ( stateTransitionMatrixHistory.rbegin( )->second * testPerturbation ).transpose( ) << std::endl;

// Eigen::Vector6d perturbedInitialState = initialState;
// perturbedInitialState[0] *= 1.001;
// perturbedInitialState[4] *= 1.001;

// Eigen::VectorXd fullPerturbedInitialState = Eigen::VectorXd::Zero( 11 );
// fullPerturbedInitialState.segment( 0, 6 ) = perturbedInitialState;
// fullPerturbedInitialState.segment( 6, 5 ) = computedInitialCoefficients;

// std::cout << "fullPerturbedInitialState " << fullPerturbedInitialState.transpose( ) << std::endl;

// fullPropagatorSettings->resetInitialStates( fullPerturbedInitialState );
// SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, fullPropagatorSettings ); 
// std::map< double, Eigen::VectorXd > results2 = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

// std::cout << "initial state " << results2.begin()->second.transpose() << std::endl;
// std::cout << "final state " << results2.rbegin()->second.transpose() << std::endl;
// // for ( auto it : results2 )
// // {
// // std::cout << it.second.transpose( ) << std::endl;
// // }

// std::map< double, Eigen::VectorXd > stateVariation;
// for ( auto it : results )
// {
//         stateVariation[ it.first ] = ( results2.at( it.first ) - it.second );
// }

// std::cout << "final state variation " << stateVariation.rbegin()->second.transpose( ) << std::endl;
// // for ( auto it : stateVariation )
// // {
// //         std::cout << it.second.transpose( ) << std::endl;
// // }

// // TEST COMPUTATION EQUILIBRIUM COEFFICIENTS

// double radius = 1800.0e3;
// double distance = 4.2e8;
// double latitude = mathematical_constants::PI / 180.0 * 2.0;
// double longitude = mathematical_constants::PI / 180.0 * 15.0;

// Eigen::VectorXd equilibriumCoefficients = computeEquilibriumCoefficients( 
//         muJupiter / muIo, fluidLoveNumber, radius, distance, longitude, latitude );
// std::cout << "equilibriumCoefficients " << equilibriumCoefficients.transpose() << std::endl;

// double perturbation_r = 0.001 * distance;
// double perturbation_longitude = 0.01 * longitude;
// double perturbation_latitude = 0.01 * latitude;

// Eigen::VectorXd perturbedEquilibriumCoefficients_r = computeEquilibriumCoefficients( 
//         muJupiter / muIo, fluidLoveNumber, radius, distance + perturbation_r, longitude, latitude );
// Eigen::VectorXd perturbedEquilibriumCoefficients_longitude = computeEquilibriumCoefficients( 
//         muJupiter / muIo, fluidLoveNumber, radius, distance, longitude + perturbation_longitude, latitude );
// Eigen::VectorXd perturbedEquilibriumCoefficients_latitude = computeEquilibriumCoefficients( 
//         muJupiter / muIo, fluidLoveNumber, radius, distance, longitude, latitude + perturbation_latitude );

// Eigen::VectorXd variation_r = ( perturbedEquilibriumCoefficients_r - equilibriumCoefficients ) / perturbation_r;
// std::cout << "variation - r" << variation_r.transpose( ) << std::endl;
// Eigen::VectorXd variation_longitude = ( perturbedEquilibriumCoefficients_longitude - equilibriumCoefficients ) / perturbation_longitude;
// std::cout << "variation - longitude" << variation_longitude.transpose( ) << std::endl;
// Eigen::VectorXd variation_latitude = ( perturbedEquilibriumCoefficients_latitude - equilibriumCoefficients ) / perturbation_latitude;
// std::cout << "variation - latitude" << variation_latitude.transpose( ) << std::endl;

//  Eigen::MatrixXd partial = equilibriumCoefficientsPartials( 
//         fluidLoveNumber, muJupiter / muIo, radius, distance, longitude, latitude );
// std::cout << "partial " << partial.transpose( ) << std::endl;


// Eigen::Vector3d initialPosition = initialState.segment( 0, 3 );
// Eigen::Vector3d sphericalPosition = 
//                 coordinate_conversions::convertCartesianToSpherical( initialPosition );

// Eigen::Vector3d perturbed_x = initialPosition;
// perturbed_x[0] += 1000.0;
// Eigen::Vector3d perturbedSphericalPosition_x = 
//                 coordinate_conversions::convertCartesianToSpherical( perturbed_x );
// Eigen::Vector3d variationSphericalPosition_x = (perturbedSphericalPosition_x - sphericalPosition) / 1000.0;
// std::cout << "variationSphericalPosition_x : " << variationSphericalPosition_x.transpose( ) << std::endl;

// Eigen::Vector3d perturbed_y = initialPosition;
// perturbed_y[1] += 1000.0;
// Eigen::Vector3d perturbedSphericalPosition_y = 
//                 coordinate_conversions::convertCartesianToSpherical( perturbed_y );
// Eigen::Vector3d variationSphericalPosition_y = (perturbedSphericalPosition_y - sphericalPosition) / 1000.0;
// std::cout << "variationSphericalPosition_y : " << variationSphericalPosition_y.transpose( ) << std::endl;

// Eigen::Vector3d perturbed_z = initialPosition;
// perturbed_z[2] += 1000.0;
// Eigen::Vector3d perturbedSphericalPosition_z = 
//                 coordinate_conversions::convertCartesianToSpherical( perturbed_z );
// Eigen::Vector3d variationSphericalPosition_z = (perturbedSphericalPosition_z - sphericalPosition) / 1000.0;
// std::cout << "variationSphericalPosition_z : " << variationSphericalPosition_z.transpose( ) << std::endl;

// Eigen::Matrix3d partialSphericalPosition = acceleration_partials::computeSphericalJacobian(initialPosition);
// std::cout << "partialSphericalPosition " << partialSphericalPosition.transpose( ) << std::endl;

// // Test full partial of body-fixed spherical state wrt body-fixed cartesian state
// bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );
// bodies.at( "Io" )->setStateFromEphemeris<>( initialTime );
// bodies.at( "Jupiter" )->setStateFromEphemeris<>( initialTime );

// Eigen::Vector6d currentStateIo;
// bodies.at( "Io" )->getStateByReference( currentStateIo );
// std::cout << "currentStateIo " << currentStateIo.transpose() << std::endl;

// Eigen::Vector3d bodyFixedAngularVelocity = bodies.at( "Io" )->getCurrentAngularVelocityVectorInLocalFrame( );
// // std::cout << "bodyFixedAngularVelocity " << bodyFixedAngularVelocity.transpose( ) << std::endl;
// Eigen::Matrix3d rotationMatrixToLocalFrame = ( bodies.at( "Io" )->getCurrentRotationToGlobalFrame( ) ).toRotationMatrix( ).transpose( );
// // std::cout << "rotationMatrixToLocalFrame " << rotationMatrixToLocalFrame << std::endl;
// Eigen::Matrix3d rotationMatrixDerivativeToLocalFrame = bodies.at( "Io" )->getCurrentRotationMatrixDerivativeToLocalFrame( );
// // std::cout << "rotationMatrixDerivativeToLocalFrame " << rotationMatrixDerivativeToLocalFrame << std::endl;

// Eigen::Vector3d bodyFixedCartesianPosition = rotationMatrixToLocalFrame * currentStateIo.segment( 0, 3 );
// Eigen::Vector3d bodyFixedCartesianVelocity = rotationMatrixToLocalFrame * currentStateIo.segment( 3, 3 )
//                 + rotationMatrixDerivativeToLocalFrame * currentStateIo.segment( 0, 3 );

// Eigen::Matrix6d fullSphericalWrtCartesianStatePartials;
// acceleration_partials::computeFullSphericalStatePartials( bodyFixedCartesianPosition, 
//                                    bodyFixedCartesianVelocity,
//                                    bodyFixedAngularVelocity,
//                                    fullSphericalWrtCartesianStatePartials );

// std::cout << "computeFullSphericalStatePartials " << std::endl;
// std::cout << fullSphericalWrtCartesianStatePartials << std::endl;

// // Perturb body-fixed cartesian position
// double deltaPosition = 1000.0;
// Eigen::Vector3d perturbedBodyFixedCartesianPositionX = bodyFixedCartesianPosition;
// perturbedBodyFixedCartesianPositionX[0] += deltaPosition;
// Eigen::Vector3d perturbedBodyFixedCartesianPositionY = bodyFixedCartesianPosition;
// perturbedBodyFixedCartesianPositionY[1] += deltaPosition;
// Eigen::Vector3d perturbedBodyFixedCartesianPositionZ = bodyFixedCartesianPosition;
// perturbedBodyFixedCartesianPositionZ[2] += deltaPosition;

// // Perturb body-fixed cartesian velocity
// double deltaVelocity = 1.0;
// Eigen::Vector3d perturbedBodyFixedCartesianVelocityX = bodyFixedCartesianVelocity;
// perturbedBodyFixedCartesianVelocityX[0] += deltaVelocity;
// Eigen::Vector3d perturbedBodyFixedCartesianVelocityY = bodyFixedCartesianVelocity;
// perturbedBodyFixedCartesianVelocityY[1] += deltaVelocity;
// Eigen::Vector3d perturbedBodyFixedCartesianVelocityZ = bodyFixedCartesianVelocity;
// perturbedBodyFixedCartesianVelocityZ[2] += deltaVelocity;

// // Compute nominal body-fixed spherical position and velocity
// Eigen::Vector3d nominalSphericalPosition = coordinate_conversions::convertCartesianToSpherical( bodyFixedCartesianPosition );
// Eigen::Vector3d nominalSphericalVelocity = getSphericalVelocityComponents( bodyFixedCartesianPosition, bodyFixedCartesianVelocity );

// Eigen::Matrix6d numericalPartial = Eigen::Matrix6d::Zero();

// // Perturb body-fixed cartesian position
// Eigen::Vector3d perturbedSphericalPositionX = coordinate_conversions::convertCartesianToSpherical( perturbedBodyFixedCartesianPositionX );
// Eigen::Vector3d variationSphericalPositionX = (perturbedSphericalPositionX - nominalSphericalPosition) / deltaPosition;

// Eigen::Vector3d perturbedSphericalPositionY = coordinate_conversions::convertCartesianToSpherical( perturbedBodyFixedCartesianPositionY );
// Eigen::Vector3d variationSphericalPositionY = (perturbedSphericalPositionY - nominalSphericalPosition) / deltaPosition;

// Eigen::Vector3d perturbedSphericalPositionZ = coordinate_conversions::convertCartesianToSpherical( perturbedBodyFixedCartesianPositionZ );
// Eigen::Vector3d variationSphericalPositionZ = (perturbedSphericalPositionZ - nominalSphericalPosition) / deltaPosition;

// Eigen::Vector3d perturbedSphericalVelocityX = getSphericalVelocityComponents( perturbedBodyFixedCartesianPositionX, bodyFixedCartesianVelocity );
// Eigen::Vector3d variationSphericalVelocityX = (perturbedSphericalVelocityX - nominalSphericalVelocity) / deltaPosition;

// Eigen::Vector3d perturbedSphericalVelocityY = getSphericalVelocityComponents( perturbedBodyFixedCartesianPositionY, bodyFixedCartesianVelocity );
// Eigen::Vector3d variationSphericalVelocityY = (perturbedSphericalVelocityY - nominalSphericalVelocity) / deltaPosition;

// Eigen::Vector3d perturbedSphericalVelocityZ = getSphericalVelocityComponents( perturbedBodyFixedCartesianPositionZ, bodyFixedCartesianVelocity );
// Eigen::Vector3d variationSphericalVelocityZ = (perturbedSphericalVelocityZ - nominalSphericalVelocity) / deltaPosition;

// numericalPartial.block(0,0,3,1) = variationSphericalPositionX;
// numericalPartial.block(0,1,3,1) = variationSphericalPositionY;
// numericalPartial.block(0,2,3,1) = variationSphericalPositionZ;

// numericalPartial.block(3,0,3,1) = variationSphericalVelocityX;
// numericalPartial.block(3,1,3,1) = variationSphericalVelocityY;
// numericalPartial.block(3,2,3,1) = variationSphericalVelocityZ;

// // Perturb body-fixed cartesian velocity
// perturbedSphericalVelocityX = getSphericalVelocityComponents( bodyFixedCartesianPosition, perturbedBodyFixedCartesianVelocityX );
// variationSphericalVelocityX = (perturbedSphericalVelocityX - nominalSphericalVelocity) / deltaVelocity;

// perturbedSphericalVelocityY = getSphericalVelocityComponents( bodyFixedCartesianPosition, perturbedBodyFixedCartesianVelocityY );
// variationSphericalVelocityY = (perturbedSphericalVelocityY - nominalSphericalVelocity) / deltaVelocity;

// perturbedSphericalVelocityZ = getSphericalVelocityComponents( bodyFixedCartesianPosition, perturbedBodyFixedCartesianVelocityZ );
// variationSphericalVelocityZ = (perturbedSphericalVelocityZ - nominalSphericalVelocity) / deltaVelocity;

// numericalPartial.block(3,3,3,1) = variationSphericalVelocityX;
// numericalPartial.block(3,4,3,1) = variationSphericalVelocityY;
// numericalPartial.block(3,5,3,1) = variationSphericalVelocityZ;
// std::cout << "numericalPartial " << std::endl;
// std::cout << numericalPartial << std::endl;

// //  Eigen::Vector3d getSphericalVelocityComponents( 
// //     Eigen::Vector3d bodyFixedPosition,
// //     Eigen::Vector3d bodyFixedVelocity )

// }

// BOOST_AUTO_TEST_CASE( test_GravityDeformationEstimation_RotationalState )
// {
//      std::cout.precision( 20 );

//    // Load spice kernels.
//     spice_interface::loadStandardSpiceKernels( );

//     // Specify initial time
//     double initialTime = 0.0;
//     double finalTime = 1.0 * physical_constants::JULIAN_DAY / 10.0; 

//     std::string globalFrameOrigin = "Jupiter";
//     std::string globalFrameOrientation = "J2000";

//     std::vector< std::string > bodiesToCreate = { "Jupiter", "Io" }; 

//     // Get body settings.
//     BodyListSettings bodySettings =
//             getDefaultBodySettings( bodiesToCreate, initialTime - 86400.0, finalTime + 86400.0, globalFrameOrigin, globalFrameOrientation );
            
//         bodySettings.at( "Jupiter" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ), "SSB", globalFrameOrientation ); 
//         // bodySettings.at( "Jupiter" )->gravityFieldSettings = get_gravitational_field( planet, 'IAU_Jupiter' )

//         // Set Jupiter's rotation as constant (no precession)
//         double rightAscensionPole = ( 358.054324066462 - 90.0 ) * mathematical_constants::PI / 180.0;
//         double declinationPole = ( 90.0 - 25.5034135739821 ) * mathematical_constants::PI / 180.0;
//         double primeMeridian = ( 284.95 ) * mathematical_constants::PI / 180.0;
//         double rotationRateJupiter = ( 870.536 * mathematical_constants::PI / 180.0 ) / 86400.0;

//         bodySettings.get( "Io" )->rotationModelSettings = simulation_setup::synchronousRotationModelSettings( "Jupiter", "J2000", "IAU_Io" );

//         double muIo = 5959924010272.5136719;
//         double muJupiter = 126686534196012800.0;
//         double muEffective = 126692494120023072.0;

//         double orbitalPeriodIo = 2.0 * mathematical_constants::PI * std::sqrt( 4.2e8 * 4.2e8 * 4.2e8 / muIo );
//         double rotationRateIo = std::sqrt( muEffective / ( 4.2e8 * 4.2e8 * 4.2e8 ) );

        
//         Eigen::Vector6d initialKeplerianState = ( Eigen::Vector6d( ) << 4.2e8, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( );
//         // std::shared_ptr< KeplerEphemerisSettings > keplerEphemerisSettings = std::make_shared< KeplerEphemerisSettings >( initialKeplerianState, 0.0, muEffective, "Jupiter", "J2000" );
//         // bodySettings.at( "Io" )->ephemerisSettings = keplerEphemerisSettings;

//         bodySettings.at( "Io" )->ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
//             ( Eigen::Vector6d( ) << 1.0 * 421.8E6, 1.0 * 0.004, 0.0, 0.0, 0.0, 0.0 ).finished( ),
//             0.0,
//             getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Io" ),
//             "Jupiter",
//             "J2000" );

//         // Create bodies needed in simulation
//         SystemOfBodies bodies = createSystemOfBodies( bodySettings );

//         bodies.at("Jupiter")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//                 rightAscensionPole, declinationPole, primeMeridian, rotationRateJupiter, initialTime, globalFrameOrientation, "IAU_Jupiter" ) );

//         // bodies.at("Io")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//                 // 0.0, mathematical_constants::PI / 180.0, 0.0, rotationRateIo, initialTime, globalFrameOrientation, "IAU_Io" ) );

//         double scaledMeanMomentOfInertia = 0.37685;
//         std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )->setScaledMeanMomentOfInertia( scaledMeanMomentOfInertia );

//         double maxwellRelaxationTime = 179103.0;
//         double globalRelaxationTime = 24688.0;
//         double fluidLoveNumber = 1.5;
//         std::vector< std::string > perturbingBody = {"Jupiter"};
//         std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::make_shared< MaxwellDeformationSettings >( 
//                 maxwellRelaxationTime, globalRelaxationTime, fluidLoveNumber, 2, 2, perturbingBody );

//     bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );
//     bodies.at( "Io" )->setStateFromEphemeris<>( initialTime );
//     bodies.at( "Jupiter" )->setStateFromEphemeris<>( initialTime );
//     // bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );


//     std::map< std::string, std::vector< std::shared_ptr< GravityDeformationSettings > > > gravityDeformationModelMap;   
//     gravityDeformationModelMap[ "Io" ] = { maxwellDeformationSettings };

//     basic_astrodynamics::GravityDeformationModelMap gravityDeformationModels = createGravityDeformationModelsMap(
//         bodies, gravityDeformationModelMap );



//     double timeStep = 10.0;
//     std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > > ( 
//         initialTime, timeStep, rungeKutta87DormandPrince, timeStep, timeStep );  

//     std::vector< std::string > bodiesToPropagate = { "Io" };

//         Eigen::VectorXd computedInitialCoefficients = Eigen::VectorXd::Zero( 5 );
//         Eigen::VectorXd perturbedInitialCoefficients = Eigen::VectorXd::Zero( 5 );
//         perturbedInitialCoefficients[0] = 1.0e-2;
//         perturbedInitialCoefficients[2] = 1.0e-2;
//         perturbedInitialCoefficients[4] = 1.0e-2;

//     std::shared_ptr< GravityDeformationPropagatorSettings< > > gravityPropagatorSettings = 
//         std::make_shared< GravityDeformationPropagatorSettings< > >( bodiesToPropagate, gravityDeformationModels, computedInitialCoefficients, integratorSettings,
//         std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

//     // // Translational dynamics propagator
//     // std::vector< std::string > centralBodies = { "Jupiter" };
//     // SelectedAccelerationMap accelerationSettingsMap;
//     // accelerationSettingsMap[ "Io" ][ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

//     std::shared_ptr< SingleArcPropagatorProcessingSettings > outputSettings =
//             std::make_shared< SingleArcPropagatorProcessingSettings >( );
//     outputSettings->setIntegratedResult( false );

//     // AccelerationMap accelerationsMap = createAccelerationModelsMap( bodies, accelerationSettingsMap, bodiesToPropagate, centralBodies );
//     Eigen::Vector6d initialTranslationalState = orbital_element_conversions::convertKeplerianToCartesianElements( initialKeplerianState, muEffective );
//     // std::shared_ptr< TranslationalStatePropagatorSettings<  > > translationalPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >( 
//     //     centralBodies, accelerationsMap, bodiesToPropagate, initialState, initialTime, integratorSettings, 
//     //     std::make_shared< PropagationTimeTerminationSettings >( finalTime ), cowell );

//     // Create torque models
//     SelectedTorqueMap torqueSettings;
//     torqueSettings[ "Io" ][ "Jupiter" ].push_back(
//                             std::make_shared< SphericalHarmonicTorqueSettings >( 1, 1 ) );
//     basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap( bodies, torqueSettings, bodiesToPropagate );

//     Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRotationState = getInitialRotationalStateOfBody( "Io", "J2000",  bodies, initialTime );
//     initialRotationState[6] = rotationRateIo;
//     // std::cout << "initialRotationState " << initialRotationState << std::endl;

//     // Create propagator settings for rotational dynamics
//     std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
//             std::make_shared< RotationalStatePropagatorSettings< double > >
//             ( torqueModelMap, bodiesToPropagate, initialRotationState, initialTime, integratorSettings, 
//             std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

//     std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > >  propagatorSettingsList;
//     propagatorSettingsList.push_back( rotationalPropagatorSettings );
//     propagatorSettingsList.push_back( gravityPropagatorSettings );
//     std::shared_ptr< MultiTypePropagatorSettings< > > fullPropagatorSettings = std::make_shared< MultiTypePropagatorSettings< > >(
//             propagatorSettingsList, integratorSettings, initialTime, std::make_shared< PropagationTimeTerminationSettings >( finalTime ),
//             std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ), outputSettings );
    

//     std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
//             getInitialStateParameterSettings< double, double >( fullPropagatorSettings, bodies );
//     // Create parameters
//     std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
//             createParametersToEstimate< double, double >( parameterNames, bodies );
// printEstimatableParameterEntries( parametersToEstimate );

// std::shared_ptr< SingleArcVariationalEquationsSolver< double, double > > variationalEquationsSolver =
//             std::make_shared< SingleArcVariationalEquationsSolver< double, double > >(
//                     bodies, fullPropagatorSettings, parametersToEstimate, true );

// std::map< double, Eigen::VectorXd > results = variationalEquationsSolver->getEquationsOfMotionSolution();
// std::cout << "initial state " << results.begin()->second.transpose() << std::endl;
//     std::cout << "final state " << results.rbegin()->second.transpose() << std::endl;
// // for ( auto it : results )
// // {
// // std::cout << it.second.transpose( ) << std::endl;
// // }

// std::map< double, Eigen::MatrixXd > stateTransitionMatrixHistory = variationalEquationsSolver->getStateTransitionMatrixSolution( );
// std::cout << "stateTransitionMatrixHistory size " << stateTransitionMatrixHistory.size( ) << std::endl;

// std::cout << "initial STM " << std::endl;
// std::cout << stateTransitionMatrixHistory.begin( )->second << std::endl;
// std::cout << "-----------------------" << std::endl;
// std::cout << "final STM " << std::endl;
// std::cout << stateTransitionMatrixHistory.rbegin( )->second << std::endl;

// std::cout << "test " << std::endl;
// Eigen::VectorXd testPerturbation = Eigen::VectorXd::Zero(7+5);
// // testPerturbation[0] = 0.001 * initialRotationState[0];
// testPerturbation[6] = 0.001 * initialRotationState[6];
// // std::cout << "testPerturbation " << testPerturbation.transpose() << std::endl;
// std::cout << ( stateTransitionMatrixHistory.rbegin( )->second * testPerturbation ).transpose( ) << std::endl;

// Eigen::VectorXd perturbedInitialState = initialRotationState;
// // perturbedInitialState[0] *= 1.001;
// perturbedInitialState[6] *= 1.001;

// Eigen::VectorXd fullPerturbedInitialState = Eigen::VectorXd::Zero( 7+5 );
// fullPerturbedInitialState.segment( 0, 7 ) = perturbedInitialState;
// fullPerturbedInitialState.segment( 7, 5 ) = computedInitialCoefficients;

// std::cout << "fullPerturbedInitialState " << fullPerturbedInitialState.transpose( ) << std::endl;

// fullPropagatorSettings->resetInitialStates( fullPerturbedInitialState );
// SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, fullPropagatorSettings ); 
// std::map< double, Eigen::VectorXd > results2 = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

// std::cout << "initial state " << results2.begin()->second.transpose() << std::endl;
// std::cout << "final state " << results2.rbegin()->second.transpose() << std::endl;
// // for ( auto it : results2 )
// // {
// // std::cout << it.second.transpose( ) << std::endl;
// // }

// std::map< double, Eigen::VectorXd > stateVariation;
// for ( auto it : results )
// {
//         stateVariation[ it.first ] = ( results2.at( it.first ) - it.second );
// }

// std::cout << "final state variation " << stateVariation.rbegin()->second.transpose( ) << std::endl;

// // MANUAL TEST

// bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );
// bodies.at( "Io" )->setStateFromEphemeris<>( initialTime );
// bodies.at( "Jupiter" )->setStateFromEphemeris<>( initialTime );

// // Compute analytical partials
//  std::shared_ptr< observation_partials::RotationMatrixPartialWrtQuaternion > rotationMatrixPartialObject = 
//     std::make_shared< observation_partials::RotationMatrixPartialWrtQuaternion >( 
//         std::bind( &Body::getCurrentRotationToGlobalFrame, bodies.at( "Io" ) ) );

// std::vector< Eigen::Matrix3d > rotationMatrixPartial = 
//     rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter( initialTime );

// Eigen::Vector6d currentInertialState = - initialTranslationalState; 
// Eigen::Vector3d currentInertialPosition = currentInertialState.segment( 0, 3 );
// Eigen::Vector3d currentInertialVelocity = currentInertialState.segment( 3, 3 );
// Eigen::Vector3d angularVelocityBodyFixedFrame = initialRotationState.segment(4, 3);  
// Eigen::Matrix3d rotationLocalToGlobal = Eigen::Quaterniond( 
//     initialRotationState[0], initialRotationState[1], initialRotationState[2], initialRotationState[3] ).toRotationMatrix();  
// Eigen::Matrix3d rotationGlobalToLocal = rotationLocalToGlobal.transpose();
// Eigen::Vector3d currentBodyFixedPosition = rotationGlobalToLocal * currentInertialPosition;  

// Eigen::MatrixXd partialBodyFixedTranslationalWrtRotationalState = Eigen::MatrixXd::Zero( 6, 7 );
    
// for( unsigned int i = 0; i < 4; i++ )
// {
//     // d R / d q0, q1, q2, q3
//     Eigen::Matrix3d currentRotationMatrixPartial = rotationMatrixPartial[ i ].transpose();

//     Eigen::Vector3d temp = currentRotationMatrixPartial * currentInertialPosition;
//     // accelerationPartial.block( 0, i, 3, 1) = temp;
//     partialBodyFixedTranslationalWrtRotationalState.block( 0, i, 3, 1) = temp;

//     Eigen::Vector3d term = currentRotationMatrixPartial * currentInertialVelocity;
//     Eigen::Vector3d crossTerm = angularVelocityBodyFixedFrame.cross( currentRotationMatrixPartial * currentInertialPosition );
//     // drdotb_dq.col(i) = term - crossTerm;

//     partialBodyFixedTranslationalWrtRotationalState.block(3, i, 3, 1) = term - crossTerm; 
// }

//     Eigen::Matrix3d drb_domega = Eigen::Matrix3d::Zero();
//     Eigen::Matrix3d drdotb_domega = linear_algebra::getCrossProductMatrix(currentBodyFixedPosition); // -[r_b]_x
//     partialBodyFixedTranslationalWrtRotationalState.block(3, 4, 3, 3) = drdotb_domega;

// std::cout << "analytical partials" << std::endl;
// std::cout << partialBodyFixedTranslationalWrtRotationalState << std::endl;

// // Compute numerical partials
// Eigen::VectorXd nominalRotationalState = initialRotationState;
// Eigen::Quaterniond nominalQuaterniond = Eigen::Quaterniond( 
//     nominalRotationalState[0], nominalRotationalState[1], nominalRotationalState[2], nominalRotationalState[3] );
// Eigen::Vector3d nominalAngularVelocity = nominalRotationalState.segment( 4, 3 );

// Eigen::Matrix3d nominalRotationMatrix = nominalQuaterniond.toRotationMatrix().transpose();
// Eigen::Matrix3d nominalRotationMatrixDerivative = - linear_algebra::getCrossProductMatrix( nominalAngularVelocity ) *
//       nominalRotationMatrix;

// Eigen::Vector3d nominalBodyFixedPosition = nominalRotationMatrix * currentInertialPosition;
// Eigen::Vector3d nominalBodyFixedVelocity = 
//     nominalRotationMatrix * currentInertialVelocity + nominalRotationMatrixDerivative * currentInertialPosition;

// double perturbation = 1.0e-8;
// Eigen::MatrixXd numericalPartials = Eigen::MatrixXd::Zero( 6, 7 );

// for ( unsigned int k = 0 ; k < 7 ; k++ )
// {
//     // Perturb initial rotational state
//     Eigen::Matrix<double, 7, 1> perturbedRotationalState = nominalRotationalState;

//     if ( k < 4 )
//     {
        // // Perturb quaternion slightly
        // Eigen::Vector4d perturbedQuaternion = nominalRotationalState.segment<4>(0);
        // perturbedQuaternion(k) += perturbation;  // or whichever component you want to perturb
        // perturbedQuaternion.normalize();    // renormalize to preserve unit quaternion constraint
        // std::cout << "perturbed q " << perturbedQuaternion.transpose() << std::endl;
        //  perturbedRotationalState.segment<4>(0) = perturbedQuaternion;
//     }
//     else
//     {
//         // Optionally perturb angular velocity too
//         perturbedRotationalState(k) += 1.0e-8;  // perturb wx, for example
//     }

//     Eigen::Quaterniond perturbedQuaterniond = Eigen::Quaterniond( 
//         perturbedRotationalState[0], perturbedRotationalState[1], perturbedRotationalState[2], perturbedRotationalState[3] );
//     Eigen::Vector3d perturbedAngularVelocity = perturbedRotationalState.segment( 4, 3 );

//     Eigen::Matrix3d perturbedRotationMatrix = perturbedQuaterniond.toRotationMatrix().transpose();
//     Eigen::Matrix3d perturbedRotationMatrixDerivative = - linear_algebra::getCrossProductMatrix( perturbedAngularVelocity ) *
//         perturbedRotationMatrix;

//     Eigen::Vector3d perturbedBodyFixedPosition = perturbedRotationMatrix * currentInertialPosition;
//     Eigen::Vector3d perturbedBodyFixedVelocity = 
//         perturbedRotationMatrix * currentInertialVelocity + perturbedRotationMatrixDerivative * currentInertialPosition;

//     Eigen::Vector3d diffPosition = ( perturbedBodyFixedPosition - nominalBodyFixedPosition ) / perturbation;
//     Eigen::Vector3d diffVelocity = ( perturbedBodyFixedVelocity - nominalBodyFixedVelocity ) / perturbation;
//     numericalPartials.block( 0, k, 3, 1 ) = diffPosition;
//     numericalPartials.block( 3, k, 3, 1 ) = diffVelocity;
// }

// std::cout << "numerical partials" << std::endl;
// std::cout << numericalPartials << std::endl;










// // std::cout << "diffPositionWrtQuaternion " << diffPositionWrtQuaternion << std::endl;
// // for ( auto it : stateVariation )
// // {
// //         std::cout << it.second.transpose( ) << std::endl;
// // }

// // // TEST COMPUTATION EQUILIBRIUM COEFFICIENTS

// // double radius = 1800.0e3;
// // double distance = 4.2e8;
// // double latitude = mathematical_constants::PI / 180.0 * 2.0;
// // double longitude = mathematical_constants::PI / 180.0 * 15.0;

// // Eigen::VectorXd equilibriumCoefficients = computeEquilibriumCoefficients( 
// //         muJupiter / muIo, fluidLoveNumber, radius, distance, longitude, latitude );
// // std::cout << "equilibriumCoefficients " << equilibriumCoefficients.transpose() << std::endl;

// // double perturbation_r = 0.001 * distance;
// // double perturbation_longitude = 0.01 * longitude;
// // double perturbation_latitude = 0.01 * latitude;

// // Eigen::VectorXd perturbedEquilibriumCoefficients_r = computeEquilibriumCoefficients( 
// //         muJupiter / muIo, fluidLoveNumber, radius, distance + perturbation_r, longitude, latitude );
// // Eigen::VectorXd perturbedEquilibriumCoefficients_longitude = computeEquilibriumCoefficients( 
// //         muJupiter / muIo, fluidLoveNumber, radius, distance, longitude + perturbation_longitude, latitude );
// // Eigen::VectorXd perturbedEquilibriumCoefficients_latitude = computeEquilibriumCoefficients( 
// //         muJupiter / muIo, fluidLoveNumber, radius, distance, longitude, latitude + perturbation_latitude );

// // Eigen::VectorXd variation_r = ( perturbedEquilibriumCoefficients_r - equilibriumCoefficients ) / perturbation_r;
// // std::cout << "variation - r" << variation_r.transpose( ) << std::endl;
// // Eigen::VectorXd variation_longitude = ( perturbedEquilibriumCoefficients_longitude - equilibriumCoefficients ) / perturbation_longitude;
// // std::cout << "variation - longitude" << variation_longitude.transpose( ) << std::endl;
// // Eigen::VectorXd variation_latitude = ( perturbedEquilibriumCoefficients_latitude - equilibriumCoefficients ) / perturbation_latitude;
// // std::cout << "variation - latitude" << variation_latitude.transpose( ) << std::endl;

// //  Eigen::MatrixXd partial = equilibriumCoefficientsPartials( 
// //         fluidLoveNumber, muJupiter / muIo, radius, distance, longitude, latitude );
// // std::cout << "partial " << partial.transpose( ) << std::endl;


// // Eigen::Vector3d initialPosition = initialState.segment( 0, 3 );
// // Eigen::Vector3d sphericalPosition = 
// //                 coordinate_conversions::convertCartesianToSpherical( initialPosition );

// // Eigen::Vector3d perturbed_x = initialPosition;
// // perturbed_x[0] += 1000.0;
// // Eigen::Vector3d perturbedSphericalPosition_x = 
// //                 coordinate_conversions::convertCartesianToSpherical( perturbed_x );
// // Eigen::Vector3d variationSphericalPosition_x = (perturbedSphericalPosition_x - sphericalPosition) / 1000.0;
// // std::cout << "variationSphericalPosition_x : " << variationSphericalPosition_x.transpose( ) << std::endl;

// // Eigen::Vector3d perturbed_y = initialPosition;
// // perturbed_y[1] += 1000.0;
// // Eigen::Vector3d perturbedSphericalPosition_y = 
// //                 coordinate_conversions::convertCartesianToSpherical( perturbed_y );
// // Eigen::Vector3d variationSphericalPosition_y = (perturbedSphericalPosition_y - sphericalPosition) / 1000.0;
// // std::cout << "variationSphericalPosition_y : " << variationSphericalPosition_y.transpose( ) << std::endl;

// // Eigen::Vector3d perturbed_z = initialPosition;
// // perturbed_z[2] += 1000.0;
// // Eigen::Vector3d perturbedSphericalPosition_z = 
// //                 coordinate_conversions::convertCartesianToSpherical( perturbed_z );
// // Eigen::Vector3d variationSphericalPosition_z = (perturbedSphericalPosition_z - sphericalPosition) / 1000.0;
// // std::cout << "variationSphericalPosition_z : " << variationSphericalPosition_z.transpose( ) << std::endl;

// // Eigen::Matrix3d partialSphericalPosition = acceleration_partials::computeSphericalJacobian(initialPosition);
// // std::cout << "partialSphericalPosition " << partialSphericalPosition.transpose( ) << std::endl;

// // // Test full partial of body-fixed spherical state wrt body-fixed cartesian state
// // bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );
// // bodies.at( "Io" )->setStateFromEphemeris<>( initialTime );
// // bodies.at( "Jupiter" )->setStateFromEphemeris<>( initialTime );

// // Eigen::Vector6d currentStateIo;
// // bodies.at( "Io" )->getStateByReference( currentStateIo );
// // std::cout << "currentStateIo " << currentStateIo.transpose() << std::endl;

// // Eigen::Vector3d bodyFixedAngularVelocity = bodies.at( "Io" )->getCurrentAngularVelocityVectorInLocalFrame( );
// // // std::cout << "bodyFixedAngularVelocity " << bodyFixedAngularVelocity.transpose( ) << std::endl;
// // Eigen::Matrix3d rotationMatrixToLocalFrame = ( bodies.at( "Io" )->getCurrentRotationToGlobalFrame( ) ).toRotationMatrix( ).transpose( );
// // // std::cout << "rotationMatrixToLocalFrame " << rotationMatrixToLocalFrame << std::endl;
// // Eigen::Matrix3d rotationMatrixDerivativeToLocalFrame = bodies.at( "Io" )->getCurrentRotationMatrixDerivativeToLocalFrame( );
// // // std::cout << "rotationMatrixDerivativeToLocalFrame " << rotationMatrixDerivativeToLocalFrame << std::endl;

// // Eigen::Vector3d bodyFixedCartesianPosition = rotationMatrixToLocalFrame * currentStateIo.segment( 0, 3 );
// // Eigen::Vector3d bodyFixedCartesianVelocity = rotationMatrixToLocalFrame * currentStateIo.segment( 3, 3 )
// //                 + rotationMatrixDerivativeToLocalFrame * currentStateIo.segment( 0, 3 );

// // Eigen::Matrix6d fullSphericalWrtCartesianStatePartials;
// // acceleration_partials::computeFullSphericalStatePartials( bodyFixedCartesianPosition, 
// //                                    bodyFixedCartesianVelocity,
// //                                    bodyFixedAngularVelocity,
// //                                    fullSphericalWrtCartesianStatePartials );

// // std::cout << "computeFullSphericalStatePartials " << std::endl;
// // std::cout << fullSphericalWrtCartesianStatePartials << std::endl;

// // // Perturb body-fixed cartesian position
// // double deltaPosition = 1000.0;
// // Eigen::Vector3d perturbedBodyFixedCartesianPositionX = bodyFixedCartesianPosition;
// // perturbedBodyFixedCartesianPositionX[0] += deltaPosition;
// // Eigen::Vector3d perturbedBodyFixedCartesianPositionY = bodyFixedCartesianPosition;
// // perturbedBodyFixedCartesianPositionY[1] += deltaPosition;
// // Eigen::Vector3d perturbedBodyFixedCartesianPositionZ = bodyFixedCartesianPosition;
// // perturbedBodyFixedCartesianPositionZ[2] += deltaPosition;

// // // Perturb body-fixed cartesian velocity
// // double deltaVelocity = 1.0;
// // Eigen::Vector3d perturbedBodyFixedCartesianVelocityX = bodyFixedCartesianVelocity;
// // perturbedBodyFixedCartesianVelocityX[0] += deltaVelocity;
// // Eigen::Vector3d perturbedBodyFixedCartesianVelocityY = bodyFixedCartesianVelocity;
// // perturbedBodyFixedCartesianVelocityY[1] += deltaVelocity;
// // Eigen::Vector3d perturbedBodyFixedCartesianVelocityZ = bodyFixedCartesianVelocity;
// // perturbedBodyFixedCartesianVelocityZ[2] += deltaVelocity;

// // // Compute nominal body-fixed spherical position and velocity
// // Eigen::Vector3d nominalSphericalPosition = coordinate_conversions::convertCartesianToSpherical( bodyFixedCartesianPosition );
// // Eigen::Vector3d nominalSphericalVelocity = getSphericalVelocityComponents( bodyFixedCartesianPosition, bodyFixedCartesianVelocity );

// // Eigen::Matrix6d numericalPartial = Eigen::Matrix6d::Zero();

// // // Perturb body-fixed cartesian position
// // Eigen::Vector3d perturbedSphericalPositionX = coordinate_conversions::convertCartesianToSpherical( perturbedBodyFixedCartesianPositionX );
// // Eigen::Vector3d variationSphericalPositionX = (perturbedSphericalPositionX - nominalSphericalPosition) / deltaPosition;

// // Eigen::Vector3d perturbedSphericalPositionY = coordinate_conversions::convertCartesianToSpherical( perturbedBodyFixedCartesianPositionY );
// // Eigen::Vector3d variationSphericalPositionY = (perturbedSphericalPositionY - nominalSphericalPosition) / deltaPosition;

// // Eigen::Vector3d perturbedSphericalPositionZ = coordinate_conversions::convertCartesianToSpherical( perturbedBodyFixedCartesianPositionZ );
// // Eigen::Vector3d variationSphericalPositionZ = (perturbedSphericalPositionZ - nominalSphericalPosition) / deltaPosition;

// // Eigen::Vector3d perturbedSphericalVelocityX = getSphericalVelocityComponents( perturbedBodyFixedCartesianPositionX, bodyFixedCartesianVelocity );
// // Eigen::Vector3d variationSphericalVelocityX = (perturbedSphericalVelocityX - nominalSphericalVelocity) / deltaPosition;

// // Eigen::Vector3d perturbedSphericalVelocityY = getSphericalVelocityComponents( perturbedBodyFixedCartesianPositionY, bodyFixedCartesianVelocity );
// // Eigen::Vector3d variationSphericalVelocityY = (perturbedSphericalVelocityY - nominalSphericalVelocity) / deltaPosition;

// // Eigen::Vector3d perturbedSphericalVelocityZ = getSphericalVelocityComponents( perturbedBodyFixedCartesianPositionZ, bodyFixedCartesianVelocity );
// // Eigen::Vector3d variationSphericalVelocityZ = (perturbedSphericalVelocityZ - nominalSphericalVelocity) / deltaPosition;

// // numericalPartial.block(0,0,3,1) = variationSphericalPositionX;
// // numericalPartial.block(0,1,3,1) = variationSphericalPositionY;
// // numericalPartial.block(0,2,3,1) = variationSphericalPositionZ;

// // numericalPartial.block(3,0,3,1) = variationSphericalVelocityX;
// // numericalPartial.block(3,1,3,1) = variationSphericalVelocityY;
// // numericalPartial.block(3,2,3,1) = variationSphericalVelocityZ;

// // // Perturb body-fixed cartesian velocity
// // perturbedSphericalVelocityX = getSphericalVelocityComponents( bodyFixedCartesianPosition, perturbedBodyFixedCartesianVelocityX );
// // variationSphericalVelocityX = (perturbedSphericalVelocityX - nominalSphericalVelocity) / deltaVelocity;

// // perturbedSphericalVelocityY = getSphericalVelocityComponents( bodyFixedCartesianPosition, perturbedBodyFixedCartesianVelocityY );
// // variationSphericalVelocityY = (perturbedSphericalVelocityY - nominalSphericalVelocity) / deltaVelocity;

// // perturbedSphericalVelocityZ = getSphericalVelocityComponents( bodyFixedCartesianPosition, perturbedBodyFixedCartesianVelocityZ );
// // variationSphericalVelocityZ = (perturbedSphericalVelocityZ - nominalSphericalVelocity) / deltaVelocity;

// // numericalPartial.block(3,3,3,1) = variationSphericalVelocityX;
// // numericalPartial.block(3,4,3,1) = variationSphericalVelocityY;
// // numericalPartial.block(3,5,3,1) = variationSphericalVelocityZ;
// // std::cout << "numericalPartial " << std::endl;
// // std::cout << numericalPartial << std::endl;

// // //  Eigen::Vector3d getSphericalVelocityComponents( 
// // //     Eigen::Vector3d bodyFixedPosition,
// // //     Eigen::Vector3d bodyFixedVelocity )

// }


BOOST_AUTO_TEST_CASE( test_GravityDeformationEstimation_RotationalState ) // IMPORTANT TEST FOR ROTATIONAL STATE
{
     std::cout.precision( 20 );

   // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Specify initial time
    double initialTime = 0.0;
    double finalTime = 10.0 * physical_constants::JULIAN_DAY; 

    std::string globalFrameOrigin = "Jupiter";
    std::string globalFrameOrientation = "J2000";

    std::vector< std::string > bodiesToCreate = { "Jupiter", "Io" }; 

    // Get body settings.
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodiesToCreate, initialTime - 86400.0, finalTime + 86400.0, globalFrameOrigin, globalFrameOrientation );
            
        // bodySettings.at( "Jupiter" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ), "SSB", globalFrameOrientation ); 
        // bodySettings.at( "Jupiter" )->gravityFieldSettings = get_gravitational_field( planet, 'IAU_Jupiter' )

        // Set Jupiter's rotation as constant (no precession)
        double rightAscensionPole = ( 358.054324066462 - 90.0 ) * mathematical_constants::PI / 180.0;
        double declinationPole = ( 90.0 - 25.5034135739821 ) * mathematical_constants::PI / 180.0;
        double primeMeridian = ( 284.95 ) * mathematical_constants::PI / 180.0;
        double rotationRateJupiter = ( 870.536 * mathematical_constants::PI / 180.0 ) / 86400.0;

        bodySettings.get( "Io" )->rotationModelSettings = simulation_setup::synchronousRotationModelSettings( "Jupiter", "J2000", "IAU_Io" );

        double muIo = 5959924010272.5136719;
        double muJupiter = 126686534196012800.0;
        double muEffective = 126692494120023072.0;

        double orbitalPeriodIo = 2.0 * mathematical_constants::PI * std::sqrt( 4.2e8 * 4.2e8 * 4.2e8 / muIo );
        double rotationRateIo = std::sqrt( muEffective / ( 4.2e8 * 4.2e8 * 4.2e8 ) );

        
        Eigen::Vector6d initialKeplerianState = ( Eigen::Vector6d( ) << 4.2e8, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( );
        // std::shared_ptr< KeplerEphemerisSettings > keplerEphemerisSettings = std::make_shared< KeplerEphemerisSettings >( initialKeplerianState, 0.0, muEffective, "Jupiter", "J2000" );
        // bodySettings.at( "Io" )->ephemerisSettings = keplerEphemerisSettings;

        // bodySettings.at( "Io" )->ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
        //     ( Eigen::Vector6d( ) << 1.0 * 421.8E6, 1.0 * 0.004, 0.0, 0.0, 0.0, 0.0 ).finished( ),
        //     0.0,
        //     getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Io" ),
        //     "Jupiter",
        //     "J2000" );

        // Create bodies needed in simulation
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        bodies.at("Jupiter")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
                rightAscensionPole, declinationPole, primeMeridian, rotationRateJupiter, initialTime, globalFrameOrientation, "IAU_Jupiter" ) );

        bodies.at("Io")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
                0.0, mathematical_constants::PI / 180.0, 0.0, rotationRateIo, initialTime, globalFrameOrientation, "IAU_Io" ) );

        double scaledMeanMomentOfInertia = 0.37685;
        std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )->setScaledMeanMomentOfInertia( scaledMeanMomentOfInertia );
        double radiusIo = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )->getReferenceRadius( );    

        Eigen::Vector6d initialTranslationalState = orbital_element_conversions::convertKeplerianToCartesianElements( initialKeplerianState, muEffective );

        double maxwellRelaxationTime = 179103.0;
        double globalRelaxationTime = 24688.0;
        double fluidLoveNumber = 1.5;
        std::vector< std::string > perturbingBody = {"Jupiter"};
        std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::make_shared< MaxwellDeformationSettings >( 
                maxwellRelaxationTime, globalRelaxationTime, fluidLoveNumber, 2, 2, perturbingBody, Eigen::Vector5d::Zero( ), true, false );

    bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );
    bodies.at( "Io" )->setStateFromEphemeris<>( initialTime );
    bodies.at( "Jupiter" )->setStateFromEphemeris<>( initialTime );
    // bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );


    std::map< std::string, std::vector< std::shared_ptr< GravityDeformationSettings > > > gravityDeformationModelMap;   
    gravityDeformationModelMap[ "Io" ] = { maxwellDeformationSettings };

    basic_astrodynamics::GravityDeformationModelMap gravityDeformationModels = createGravityDeformationModelsMap(
        bodies, gravityDeformationModelMap );



    double timeStep = 100.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > > ( 
        initialTime, timeStep, rungeKutta87DormandPrince, 1.0e-5, 1000.0 );  

    // std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< IntegratorSettings< > >( rungeKutta4, 0.0, timeStep );

    std::vector< std::string > bodiesToPropagate = { "Io" };

    Eigen::VectorXd equilibriumCoefficients = Eigen::VectorXd::Zero(5);

    double distanceIo = initialTranslationalState.segment(0, 3).norm();
    double radiusRatioPowerThree = ( radiusIo / distanceIo ) * ( radiusIo / distanceIo ) * ( radiusIo / distanceIo );
    double gravitationalParametersRatio = muJupiter / muIo;
    
    equilibriumCoefficients[ 0 ] = fluidLoveNumber / 2.0 * gravitationalParametersRatio * radiusRatioPowerThree 
                * ( 3.0 * std::sin( 0.0 ) * std::sin( 0.0 ) - 1.0 );
                // - fluidLoveNumber * rotationRateIo * rotationRateIo * radiusIo * radiusIo * radiusIo / ( 3.0 * muIo ); 
    equilibriumCoefficients[ 2 ] = fluidLoveNumber / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
                ( 1.0 - std::sin( 0.0 ) * std::sin( 0.0 ) ) * std::cos( 2.0 * 0.0 );
    equilibriumCoefficients[ 4 ] = fluidLoveNumber / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
                ( 1.0 - std::sin( 0.0 ) * std::sin( 0.0 ) ) * std::sin( 2.0 * 0.0 );

    std::cout << "c20 - 1 " << fluidLoveNumber / 2.0 * gravitationalParametersRatio * radiusRatioPowerThree 
                * ( 3.0 * std::sin( 0.0 ) * std::sin( 0.0 ) - 1.0 ) << std::endl;
    std::cout << "c20 - 2 " << - fluidLoveNumber * rotationRateIo * rotationRateIo * radiusIo * radiusIo * radiusIo / ( 3.0 * muIo ) << std::endl;
        
    equilibriumCoefficients[ 1 ] = - fluidLoveNumber * gravitationalParametersRatio * radiusRatioPowerThree 
                * ( - std::cos( 0.0 ) * std::sin( 0.0 ) ) * std::cos( 0.0 );
    equilibriumCoefficients[ 3 ] = - fluidLoveNumber * gravitationalParametersRatio * radiusRatioPowerThree 
                * ( - std::cos( 0.0 ) * std::sin( 0.0 ) ) * std::sin( 0.0 );

        Eigen::VectorXd computedInitialCoefficients = equilibriumCoefficients; // Eigen::VectorXd::Zero( 5 );
        // Eigen::VectorXd perturbedInitialCoefficients = Eigen::VectorXd::Zero( 5 );
        // perturbedInitialCoefficients[0] = 1.0e-2;
        // perturbedInitialCoefficients[2] = 1.0e-2;
        // perturbedInitialCoefficients[4] = 1.0e-2;



    std::shared_ptr< GravityDeformationPropagatorSettings< > > gravityPropagatorSettings = 
        std::make_shared< GravityDeformationPropagatorSettings< > >( bodiesToPropagate, gravityDeformationModels, equilibriumCoefficients, integratorSettings,
        std::make_shared< PropagationTimeTerminationSettings >( finalTime, true ) );

    // // Translational dynamics propagator
    // std::vector< std::string > centralBodies = { "Jupiter" };
    // SelectedAccelerationMap accelerationSettingsMap;
    // accelerationSettingsMap[ "Io" ][ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

    std::shared_ptr< SingleArcPropagatorProcessingSettings > outputSettings =
            std::make_shared< SingleArcPropagatorProcessingSettings >( );
    outputSettings->setIntegratedResult( false );

    // AccelerationMap accelerationsMap = createAccelerationModelsMap( bodies, accelerationSettingsMap, bodiesToPropagate, centralBodies );
    // std::shared_ptr< TranslationalStatePropagatorSettings<  > > translationalPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >( 
    //     centralBodies, accelerationsMap, bodiesToPropagate, initialState, initialTime, integratorSettings, 
    //     std::make_shared< PropagationTimeTerminationSettings >( finalTime ), cowell );

    // Create torque models
    SelectedTorqueMap torqueSettings;
    torqueSettings[ "Io" ][ "Jupiter" ].push_back( std::make_shared< TorqueSettings >( basic_astrodynamics::second_order_gravitational_torque ) );
    // torqueSettings[ "Io" ][ "Jupiter" ].push_back( std::make_shared<SphericalHarmonicTorqueSettings>(2,2) );
    basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap( bodies, torqueSettings, bodiesToPropagate );

    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRotationState = getInitialRotationalStateOfBody( "Io", "J2000",  bodies, initialTime );
    initialRotationState[6] = rotationRateIo;
    // std::cout << "initialRotationState " << initialRotationState << std::endl;

    // Create propagator settings for rotational dynamics
    std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
            std::make_shared< RotationalStatePropagatorSettings< double > >
            ( torqueModelMap, bodiesToPropagate, initialRotationState, initialTime, integratorSettings, 
            std::make_shared< PropagationTimeTerminationSettings >( finalTime, true ) );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > >  propagatorSettingsList;
    propagatorSettingsList.push_back( rotationalPropagatorSettings );
    propagatorSettingsList.push_back( gravityPropagatorSettings );
    std::shared_ptr< MultiTypePropagatorSettings< > > fullPropagatorSettings = std::make_shared< MultiTypePropagatorSettings< > >(
            propagatorSettingsList, integratorSettings, initialTime, std::make_shared< PropagationTimeTerminationSettings >( finalTime, true ),
            std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ), outputSettings );
    

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialStateParameterSettings< double, double >( fullPropagatorSettings, bodies );
    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies );
printEstimatableParameterEntries( parametersToEstimate );

// bodies.at( "Io" )->getMassProperties( )->updateInertiaTensorDerivative( 
        // ( Eigen::Vector5d( ) << 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ) ); 

Eigen::VectorXd initialTranslationalState2 =
        getInitialStatesOfBodies( bodiesToPropagate, std::vector< std::string >({"Jupiter"}), bodies, initialTime );
std::cout << "initialTranslationalState2 " << initialTranslationalState2.transpose() << std::endl;

std::shared_ptr< SingleArcVariationalEquationsSolver< double, double > > variationalEquationsSolver =
            std::make_shared< SingleArcVariationalEquationsSolver< double, double > >(
                    bodies, fullPropagatorSettings, parametersToEstimate, true );

std::map< double, Eigen::VectorXd > results = variationalEquationsSolver->getEquationsOfMotionSolution();
std::cout << "initial state " << results.begin()->second.transpose() << std::endl;
    std::cout << "final state " << results.rbegin()->second.transpose() << std::endl;
// for ( auto it : results )
// {
// std::cout << it.second.transpose( ) << std::endl;
// }

std::map< double, Eigen::MatrixXd > stateTransitionMatrixHistory = variationalEquationsSolver->getStateTransitionMatrixSolution( );
std::cout << "stateTransitionMatrixHistory size " << stateTransitionMatrixHistory.size( ) << std::endl;

std::cout << "initial STM " << std::endl;
std::cout << stateTransitionMatrixHistory.begin( )->second << std::endl;
std::cout << "-----------------------" << std::endl;
std::cout << "final STM " << std::endl;
std::cout << stateTransitionMatrixHistory.rbegin( )->second << std::endl;

Eigen::VectorXd fullInitialState = Eigen::VectorXd::Zero( 7+5 );
fullInitialState.segment( 0, 7 ) = initialRotationState;
fullInitialState.segment( 7, 5 ) = computedInitialCoefficients;

// Eigen::VectorXd testPerturbation = Eigen::VectorXd::Zero(7+5);
// testPerturbation[0] = 1.0e-8;
// testPerturbation[6] = 0.001 * initialRotationState[6];
// std::cout << "testPerturbation " << testPerturbation.transpose() << std::endl;


// // Perturb quaternion slightly
// double perturbation = 1.0e-8;
Eigen::VectorXd perturbedRotationalState = initialRotationState;

// Eigen::Vector4d perturbedQuaternion = perturbedRotationalState.segment<4>(0);
// perturbedQuaternion[0] += perturbation;  
// perturbedQuaternion.normalize();    // re-normalise
// perturbedRotationalState.segment<4>(0) = perturbedQuaternion;

// Perturb rotation

// // Define a very small rotation perturbation (in radians) about a chosen axis, e.g. body X
// Eigen::Vector3d deltaTheta = Eigen::Vector3d::UnitZ() * 1.0e-3;

// // Build a small rotation quaternion: q_delta ≈ [1, 0.5*deltaTheta]
// Eigen::Quaterniond q_delta(1.0, 0.5 * deltaTheta.x(), 0.5 * deltaTheta.y(), 0.5 * deltaTheta.z());
// q_delta.normalize();  // ensure unit length, though nearly unnecessary for small angles

// // Apply perturbation (post-multiply for body-frame rotation, pre-multiply for inertial-frame)
// Eigen::Quaterniond rotationToIntegrationFrame = Eigen::Quaterniond( initialRotationState[0], initialRotationState[1], initialRotationState[2], initialRotationState[3] );
// Eigen::Quaterniond perturbedRotation = ( rotationToIntegrationFrame * q_delta ).normalized();

// perturbedRotationalState.segment<4>(0) = linear_algebra::convertQuaternionToVectorFormat( perturbedRotation );

// Eigen::VectorXd rotationalStatePerturbation = perturbedRotationalState - initialRotationState;
// std::cout << "rotationalStatePerturbation " << rotationalStatePerturbation.transpose() << std::endl;


////////////////////////

//  Eigen::Matrix< double, 7, 1 > unitRotationState = Eigen::Matrix< double, 7, 1 >::Zero( );
//     unitRotationState( 0 ) = noRotationQuaternion.w( );
//     unitRotationState( 1 ) = noRotationQuaternion.x( );
//     unitRotationState( 2 ) = noRotationQuaternion.y( );
//     unitRotationState( 3 ) = noRotationQuaternion.z( );
//     unitRotationState( 4 ) = 0.0;
//     unitRotationState( 5 ) = 0.0;
//     unitRotationState( 6 ) = rotationRateIo + 1.0E-6;

//     Eigen::Matrix< double, 7, 1 > originalRotationState = unitRotationState;
//     Eigen::Matrix< double, 7, 1 > stateDifferenceToAdd = initialStateDifference.segment( 0, 7 );


    perturbedRotationalState[1] += 1.0e-5; // += stateDifferenceToAdd;
    perturbedRotationalState[0] = ( initialRotationState( 0 ) > 0 ? 1.0 : -1.0 ) * 
            std::sqrt( 1.0 - std::pow( perturbedRotationalState.segment( 1, 3 ).norm( ), 2.0 ) );
    perturbedRotationalState.segment( 0, 4 ) = perturbedRotationalState.segment( 0, 4 ).normalized( );
    
    Eigen::VectorXd rotationalStatePerturbation = perturbedRotationalState - initialRotationState;
    std::cout << "rotationalStatePerturbation " << rotationalStatePerturbation.transpose() << std::endl;
    // appliedStateDifference.segment( 7, 5 ) = initialStateDifference.segment( 7, 5 );

    ///////////////

// Reconstruct full initial state (rotation + deformation)
Eigen::VectorXd fullPerturbedInitialState = Eigen::VectorXd::Zero( 7+5 );
fullPerturbedInitialState.segment( 0, 7 ) = perturbedRotationalState;
fullPerturbedInitialState.segment( 7, 5 ) = computedInitialCoefficients;
std::cout << "fullPerturbedInitialState " << fullPerturbedInitialState.transpose( ) << std::endl;

Eigen::VectorXd fullStatePerturbation = fullPerturbedInitialState - fullInitialState;
std::cout << "fullStatePerturbation " << fullStatePerturbation.transpose( ) << std::endl;

std::cout << "from state transition matrix" << std::endl;
std::cout << ( stateTransitionMatrixHistory.rbegin( )->second * fullStatePerturbation ).transpose( ) << std::endl;

fullPropagatorSettings->resetInitialStates( fullPerturbedInitialState );
// bodies.at( "Io" )->getMassProperties( )->updateInertiaTensorDerivative( 
        // ( Eigen::Vector5d( ) << 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ) ); 
// bodies.at( "Io" )->setCurrentPropagatedGravityField( fullPerturbedInitialState.segment(7,5) ); 
// SingleArcDynamicsSimulator< > dynamicsSimulator_test( bodies, fullPropagatorSettings ); 
// std::map< double, Eigen::VectorXd > results = dynamicsSimulator_test.getEquationsOfMotionNumericalSolution( );

SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, fullPropagatorSettings ); 
std::map< double, Eigen::VectorXd > results2 = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

// std::cout << "initial state " << results2.begin()->second.transpose() << std::endl;
// std::cout << "final state " << results2.rbegin()->second.transpose() << std::endl;
// for ( auto it : results2 )
// {
// std::cout << it.second.transpose( ) << std::endl;
// }

std::map< double, Eigen::VectorXd > stateVariation;
std::cout << "numerical variation " << std::endl;
// for ( auto it : results )
// {
//         stateVariation[ it.first ] = ( results2.at( it.first ) - it.second );
//         // std::cout << stateVariation[ it.first ].transpose( ) << std::endl; 
// }

std::cout << ( results2.rbegin()->second - results.rbegin()->second ).transpose( ) << std::endl;    

}

// BOOST_AUTO_TEST_CASE( test_GravityDeformationEstimation_TranslationalState ) // IMPORTANT?
// {
//      std::cout.precision( 20 );

//    // Load spice kernels.
//     spice_interface::loadStandardSpiceKernels( );

//     // Specify initial time
//     double initialTime = 0.0;
//     double finalTime = 1.0 * physical_constants::JULIAN_DAY; 

//     std::string globalFrameOrigin = "Jupiter";
//     std::string globalFrameOrientation = "J2000";

//     std::vector< std::string > bodiesToCreate = { "Jupiter", "Io" }; 

//     // Get body settings.
//     BodyListSettings bodySettings =
//             getDefaultBodySettings( bodiesToCreate, initialTime - 86400.0, finalTime + 86400.0, globalFrameOrigin, globalFrameOrientation );
            
//         bodySettings.at( "Jupiter" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ), "SSB", globalFrameOrientation ); 
//         // bodySettings.at( "Jupiter" )->gravityFieldSettings = get_gravitational_field( planet, 'IAU_Jupiter' )

//         // Set Jupiter's rotation as constant (no precession)
//         double rightAscensionPole = ( 358.054324066462 - 90.0 ) * mathematical_constants::PI / 180.0;
//         double declinationPole = ( 90.0 - 25.5034135739821 ) * mathematical_constants::PI / 180.0;
//         double primeMeridian = ( 284.95 ) * mathematical_constants::PI / 180.0;
//         double rotationRateJupiter = ( 870.536 * mathematical_constants::PI / 180.0 ) / 86400.0;

//         bodySettings.get( "Io" )->rotationModelSettings = simulation_setup::synchronousRotationModelSettings( "Jupiter", "J2000", "IAU_Io" );

//         double muIo = 5959924010272.5136719;
//         double muJupiter = 126686534196012800.0;
//         double muEffective = 126692494120023072.0;

//         double orbitalPeriodIo = 2.0 * mathematical_constants::PI * std::sqrt( 4.2e8 * 4.2e8 * 4.2e8 / muIo );
//         double rotationRateIo = std::sqrt( muEffective / ( 4.2e8 * 4.2e8 * 4.2e8 ) );

        
//         Eigen::Vector6d initialKeplerianState = ( Eigen::Vector6d( ) << 4.2e8, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( );
//         // std::shared_ptr< KeplerEphemerisSettings > keplerEphemerisSettings = std::make_shared< KeplerEphemerisSettings >( initialKeplerianState, 0.0, muEffective, "Jupiter", "J2000" );
//         // bodySettings.at( "Io" )->ephemerisSettings = keplerEphemerisSettings;

//         bodySettings.at( "Io" )->ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
//             ( Eigen::Vector6d( ) << 1.0 * 421.8E6, 1.0 * 0.004, 0.0, 0.0, 0.0, 0.0 ).finished( ),
//             0.0,
//             getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Io" ),
//             "Jupiter",
//             "J2000" );

//         // Create bodies needed in simulation
//         SystemOfBodies bodies = createSystemOfBodies( bodySettings );

//         bodies.at("Jupiter")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//                 rightAscensionPole, declinationPole, primeMeridian, rotationRateJupiter, initialTime, globalFrameOrientation, "IAU_Jupiter" ) );

//         bodies.at("Io")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//                 0.0, mathematical_constants::PI / 180.0, 0.0, rotationRateIo, initialTime, globalFrameOrientation, "IAU_Io" ) );

//         double scaledMeanMomentOfInertia = 0.37685;
//         std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )->setScaledMeanMomentOfInertia( scaledMeanMomentOfInertia );
//         double radiusIo = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )->getReferenceRadius( );    

//         Eigen::Vector6d initialTranslationalState = orbital_element_conversions::convertKeplerianToCartesianElements( initialKeplerianState, muEffective );

//         double maxwellRelaxationTime = 179103.0;
//         double globalRelaxationTime = 24688.0;
//         double fluidLoveNumber = 1.5;
//         std::vector< std::string > perturbingBody = {"Jupiter"};
//         std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::make_shared< MaxwellDeformationSettings >( 
//                 maxwellRelaxationTime, globalRelaxationTime, fluidLoveNumber, 2, 2, perturbingBody );

//     bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );
//     bodies.at( "Io" )->setStateFromEphemeris<>( initialTime );
//     bodies.at( "Jupiter" )->setStateFromEphemeris<>( initialTime );
//     // bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );


//     std::map< std::string, std::vector< std::shared_ptr< GravityDeformationSettings > > > gravityDeformationModelMap;   
//     gravityDeformationModelMap[ "Io" ] = { maxwellDeformationSettings };

//     basic_astrodynamics::GravityDeformationModelMap gravityDeformationModels = createGravityDeformationModelsMap(
//         bodies, gravityDeformationModelMap );



//     double timeStep = 10.0 / 2.0;
//     std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > > ( 
//         initialTime, timeStep, rungeKutta87DormandPrince, timeStep, timeStep );  

//     // std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< IntegratorSettings< > >( rungeKutta4, 0.0, timeStep );

//     std::vector< std::string > bodiesToPropagate = { "Io" };

//     Eigen::VectorXd equilibriumCoefficients = Eigen::VectorXd::Zero(5);

//     double distanceIo = initialTranslationalState.segment(0, 3).norm();
//     double radiusRatioPowerThree = ( radiusIo / distanceIo ) * ( radiusIo / distanceIo ) * ( radiusIo / distanceIo );
//     double gravitationalParametersRatio = muJupiter / muIo;
    
//     equilibriumCoefficients[ 0 ] = fluidLoveNumber / 2.0 * gravitationalParametersRatio * radiusRatioPowerThree 
//                 * ( 3.0 * std::sin( 0.0 ) * std::sin( 0.0 ) - 1.0 ); 
//     equilibriumCoefficients[ 2 ] = fluidLoveNumber / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
//                 ( 1.0 - std::sin( 0.0 ) * std::sin( 0.0 ) ) * std::cos( 2.0 * 0.0 );
//     equilibriumCoefficients[ 4 ] = fluidLoveNumber / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
//                 ( 1.0 - std::sin( 0.0 ) * std::sin( 0.0 ) ) * std::sin( 2.0 * 0.0 );
        
//     equilibriumCoefficients[ 1 ] = - fluidLoveNumber * gravitationalParametersRatio * radiusRatioPowerThree 
//                 * ( - std::cos( 0.0 ) * std::sin( 0.0 ) ) * std::cos( 0.0 );
//     equilibriumCoefficients[ 3 ] = - fluidLoveNumber * gravitationalParametersRatio * radiusRatioPowerThree 
//                 * ( - std::cos( 0.0 ) * std::sin( 0.0 ) ) * std::sin( 0.0 );

//         Eigen::VectorXd computedInitialCoefficients = equilibriumCoefficients; // Eigen::VectorXd::Zero( 5 );
//         // Eigen::VectorXd perturbedInitialCoefficients = Eigen::VectorXd::Zero( 5 );
//         // perturbedInitialCoefficients[0] = 1.0e-2;
//         // perturbedInitialCoefficients[2] = 1.0e-2;
//         // perturbedInitialCoefficients[4] = 1.0e-2;



//     std::shared_ptr< GravityDeformationPropagatorSettings< > > gravityPropagatorSettings = 
//         std::make_shared< GravityDeformationPropagatorSettings< > >( bodiesToPropagate, gravityDeformationModels, equilibriumCoefficients, integratorSettings,
//         std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

//     // Translational dynamics propagator
//     std::vector< std::string > centralBodies = { "Jupiter" };
//     SelectedAccelerationMap accelerationSettingsMap;
//     accelerationSettingsMap[ "Io" ][ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

//     std::shared_ptr< SingleArcPropagatorProcessingSettings > outputSettings =
//             std::make_shared< SingleArcPropagatorProcessingSettings >( );
//     outputSettings->setIntegratedResult( false );

//     Eigen::Vector6d initialState = orbital_element_conversions::convertKeplerianToCartesianElements( initialKeplerianState, muEffective );

//     AccelerationMap accelerationsMap = createAccelerationModelsMap( bodies, accelerationSettingsMap, bodiesToPropagate, centralBodies );
//     std::shared_ptr< TranslationalStatePropagatorSettings<  > > translationalPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >( 
//         centralBodies, accelerationsMap, bodiesToPropagate, initialState, initialTime, integratorSettings, 
//         std::make_shared< PropagationTimeTerminationSettings >( finalTime ), cowell );

//     std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > >  propagatorSettingsList;
//     propagatorSettingsList.push_back( translationalPropagatorSettings );
//     propagatorSettingsList.push_back( gravityPropagatorSettings );
//     std::shared_ptr< MultiTypePropagatorSettings< > > fullPropagatorSettings = std::make_shared< MultiTypePropagatorSettings< > >(
//             propagatorSettingsList, integratorSettings, initialTime, std::make_shared< PropagationTimeTerminationSettings >( finalTime ),
//             std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ), outputSettings );
    

//     std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
//             getInitialStateParameterSettings< double, double >( fullPropagatorSettings, bodies );
//     // Create parameters
//     std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
//             createParametersToEstimate< double, double >( parameterNames, bodies );
// printEstimatableParameterEntries( parametersToEstimate );

// bodies.at( "Io" )->getMassProperties( )->updateInertiaTensorDerivative( 
//         ( Eigen::Vector5d( ) << 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ) ); 

// std::shared_ptr< SingleArcVariationalEquationsSolver< double, double > > variationalEquationsSolver =
//             std::make_shared< SingleArcVariationalEquationsSolver< double, double > >(
//                     bodies, fullPropagatorSettings, parametersToEstimate, true );

// std::map< double, Eigen::VectorXd > results = variationalEquationsSolver->getEquationsOfMotionSolution();
// std::cout << "initial state " << results.begin()->second.transpose() << std::endl;
//     std::cout << "final state " << results.rbegin()->second.transpose() << std::endl;
// // for ( auto it : results )
// // {
// // std::cout << it.second.transpose( ) << std::endl;
// // }

// std::map< double, Eigen::MatrixXd > stateTransitionMatrixHistory = variationalEquationsSolver->getStateTransitionMatrixSolution( );
// std::cout << "stateTransitionMatrixHistory size " << stateTransitionMatrixHistory.size( ) << std::endl;

// std::cout << "initial STM " << std::endl;
// std::cout << stateTransitionMatrixHistory.begin( )->second << std::endl;
// std::cout << "-----------------------" << std::endl;
// std::cout << "final STM " << std::endl;
// std::cout << stateTransitionMatrixHistory.rbegin( )->second << std::endl;

// Eigen::VectorXd fullInitialState = Eigen::VectorXd::Zero( 6+5 );
// fullInitialState.segment( 0, 6 ) = initialState;
// fullInitialState.segment( 6, 5 ) = computedInitialCoefficients;

// // Eigen::VectorXd testPerturbation = Eigen::VectorXd::Zero(6+5);
// // testPerturbation[0] = 0.0001 * fullInitialState[0];
// // testPerturbation[6] = 0.001 * initialRotationState[6];
// // std::cout << "testPerturbation " << testPerturbation.transpose() << std::endl;


// // // Perturb quaternion slightly
// // double perturbation = 1.0e-8;
// // Eigen::VectorXd perturbedRotationalState = initialRotationState;

// // Eigen::Vector4d perturbedQuaternion = perturbedRotationalState.segment<4>(0);
// // perturbedQuaternion[0] += perturbation;  
// // perturbedQuaternion.normalize();    // re-normalise
// // perturbedRotationalState.segment<4>(0) = perturbedQuaternion;

// // Perturb rotation

// // // Define a very small rotation perturbation (in radians) about a chosen axis, e.g. body X
// // Eigen::Vector3d deltaTheta = Eigen::Vector3d::UnitZ() * 1.0e-3;

// // // Build a small rotation quaternion: q_delta ≈ [1, 0.5*deltaTheta]
// // Eigen::Quaterniond q_delta(1.0, 0.5 * deltaTheta.x(), 0.5 * deltaTheta.y(), 0.5 * deltaTheta.z());
// // q_delta.normalize();  // ensure unit length, though nearly unnecessary for small angles

// // // Apply perturbation (post-multiply for body-frame rotation, pre-multiply for inertial-frame)
// // Eigen::Quaterniond rotationToIntegrationFrame = Eigen::Quaterniond( initialRotationState[0], initialRotationState[1], initialRotationState[2], initialRotationState[3] );
// // Eigen::Quaterniond perturbedRotation = ( rotationToIntegrationFrame * q_delta ).normalized();

// // perturbedRotationalState.segment<4>(0) = linear_algebra::convertQuaternionToVectorFormat( perturbedRotation );

// // Eigen::VectorXd rotationalStatePerturbation = perturbedRotationalState - initialRotationState;
// // std::cout << "rotationalStatePerturbation " << rotationalStatePerturbation.transpose() << std::endl;


// ////////////////////////

// //  Eigen::Matrix< double, 7, 1 > unitRotationState = Eigen::Matrix< double, 7, 1 >::Zero( );
// //     unitRotationState( 0 ) = noRotationQuaternion.w( );
// //     unitRotationState( 1 ) = noRotationQuaternion.x( );
// //     unitRotationState( 2 ) = noRotationQuaternion.y( );
// //     unitRotationState( 3 ) = noRotationQuaternion.z( );
// //     unitRotationState( 4 ) = 0.0;
// //     unitRotationState( 5 ) = 0.0;
// //     unitRotationState( 6 ) = rotationRateIo + 1.0E-6;

// //     Eigen::Matrix< double, 7, 1 > originalRotationState = unitRotationState;
// //     Eigen::Matrix< double, 7, 1 > stateDifferenceToAdd = initialStateDifference.segment( 0, 7 );


//     // perturbedRotationalState[1] += 1.0e-8; // += stateDifferenceToAdd;
//     // perturbedRotationalState[0] = initialRotationState[0] / std::fabs( initialRotationState[0] ) *
//     //         std::sqrt( 1.0 - std::pow( perturbedRotationalState.segment( 1, 3 ).norm( ), 2.0 ) );

//     // Eigen::VectorXd rotationalStatePerturbation = perturbedRotationalState - initialRotationState;
//     // std::cout << "rotationalStatePerturbation " << rotationalStatePerturbation.transpose() << std::endl;
//     // // appliedStateDifference.segment( 7, 5 ) = initialStateDifference.segment( 7, 5 );


//     ///////////////

// // Reconstruct full initial state (rotation + deformation)
// Eigen::VectorXd fullPerturbedInitialState = fullInitialState;
// fullPerturbedInitialState[0] += 0.0001 * fullInitialState[0];
// // fullPerturbedInitialState.segment( 0, 7 ) = perturbedRotationalState;
// // fullPerturbedInitialState.segment( 7, 5 ) = computedInitialCoefficients;
// // std::cout << "fullPerturbedInitialState " << fullPerturbedInitialState.transpose( ) << std::endl;

// Eigen::VectorXd fullStatePerturbation = fullPerturbedInitialState - fullInitialState;
// std::cout << "fullStatePerturbation " << fullStatePerturbation.transpose( ) << std::endl;

// std::cout << "from state transition matrix" << std::endl;
// std::cout << ( stateTransitionMatrixHistory.rbegin( )->second * fullStatePerturbation ).transpose( ) << std::endl;

// fullPropagatorSettings->resetInitialStates( fullPerturbedInitialState );
// // bodies.at( "Io" )->getMassProperties( )->updateInertiaTensorDerivative( 
//         // ( Eigen::Vector5d( ) << 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ) ); 
// // bodies.at( "Io" )->setCurrentPropagatedGravityField( fullPerturbedInitialState.segment(7,5) ); 
// // SingleArcDynamicsSimulator< > dynamicsSimulator_test( bodies, fullPropagatorSettings ); 
// // std::map< double, Eigen::VectorXd > results = dynamicsSimulator_test.getEquationsOfMotionNumericalSolution( );

// SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, fullPropagatorSettings ); 
// std::map< double, Eigen::VectorXd > results2 = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

// // std::cout << "initial state " << results2.begin()->second.transpose() << std::endl;
// // std::cout << "final state " << results2.rbegin()->second.transpose() << std::endl;
// // for ( auto it : results2 )
// // {
// // std::cout << it.second.transpose( ) << std::endl;
// // }

// std::map< double, Eigen::VectorXd > stateVariation;
// std::cout << "numerical variation " << std::endl;
// for ( auto it : results )
// {
//         stateVariation[ it.first ] = ( results2.at( it.first ) - it.second );
//         // std::cout << stateVariation[ it.first ].transpose( ) << std::endl; 
// }

// std::cout << stateVariation.rbegin()->second.transpose( ) << std::endl;    

// }

// SystemOfBodies createIoDynamicalEnvironment( const double initialTime, const double finalTime )
// {
//     std::string globalFrameOrigin = "Jupiter";
//     std::string globalFrameOrientation = "J2000";

//     // Load spice kernels.
//     spice_interface::loadStandardSpiceKernels( );

//     // Create bodies
//     SystemOfBodies bodies = SystemOfBodies( globalFrameOrigin, globalFrameOrientation );

//     // Create Jupiter
//     bodies.createEmptyBody( "Jupiter", false );
    
//     // Set Jupiter's ephemeris
//     bodies.at( "Jupiter" )->setEphemeris(
//             std::make_shared< ephemerides::ConstantEphemeris >( [ = ]( ) { return Eigen::Vector6d::Zero( ); }, "SSB", globalFrameOrientation ) );

//     // Set Jupiter's gravity field
//     double muJupiter = spice_interface::getBodyGravitationalParameter( "Jupiter" );
//     bodies.at( "Jupiter" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( muJupiter ) );

//     // Set Jupiter's rotation as constant (no precession)
//     double rightAscensionPole = ( 358.054324066462 - 90.0 ) * mathematical_constants::PI / 180.0;
//     double declinationPole = ( 90.0 - 25.5034135739821 ) * mathematical_constants::PI / 180.0;
//     double primeMeridian = ( 284.95 ) * mathematical_constants::PI / 180.0;
//     double rotationRateJupiter = ( 870.536 * mathematical_constants::PI / 180.0 ) / 86400.0;   
    
//     bodies.at("Jupiter")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//                 rightAscensionPole, declinationPole, primeMeridian, rotationRateJupiter, initialTime, globalFrameOrientation, "IAU_Jupiter" ) );
    
//     // Create Io
//     bodies.createEmptyBody( "Io" );

//     // Set Io's gravity field
//     double muIo = spice_interface::getBodyGravitationalParameter( "Io" ); 
//     double radiusIo = 1821.6e3;
//     double scaledMeanMomentOfInertia = 0.37685;
//     Eigen::MatrixXd ioCosineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );
//     Eigen::MatrixXd ioSineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );
//     ioCosineCoefficients( 0, 0 ) = 1.0;
//     ioCosineCoefficients( 2, 0 ) = -1845.9E-6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
//     ioCosineCoefficients( 2, 2 ) = 553.7E-6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
          
//     bodies.at( "Io" )->setGravityFieldModel( std::make_shared< gravitation::SphericalHarmonicsGravityField >( 
//         muIo, radiusIo, ioCosineCoefficients, ioSineCoefficients, "IAU_Io", scaledMeanMomentOfInertia ) );

//     // Retrieve Io's default inertia tensor
//     Eigen::Matrix3d ioInertiaTensor = bodies.at("Io")->getBodyInertiaTensor();
//     std::cout << "ioInertiaTensor" << std::endl;
//     std::cout << ioInertiaTensor << std::endl;

//     // Set Io's ephemeris
//     double muEffective = muIo + muJupiter;
//     Eigen::Vector6d ioKeplerElements = Eigen::Vector6d::Zero( );
//     double ioSemiMajorAxis = 4.2e8;
//     ioKeplerElements( 0 ) = ioSemiMajorAxis;
//     ioKeplerElements( 1 ) = 0.004;
//     bodies.at( "Io" )->setEphemeris( std::make_shared< ephemerides::KeplerEphemeris >( ioKeplerElements, 0.0, muEffective, globalFrameOrigin, globalFrameOrientation ) );

//     // Set Io's rotation model 
//     double rotationRateIo = std::sqrt( muEffective / ( 4.2e8 * 4.2e8 * 4.2e8 ) );
// //     bodies.at( "Io" )->setRotationalEphemeris( std::make_shared< SynchronousRotationalEphemeris >(
//         // createRelativeStateFunction( bodies, "Io", "Jupiter" ), "Jupiter", "J2000", "IAU_Io" ) );
//     bodies.at("Io")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//         0.0, mathematical_constants::PI / 2.0, 0.0, rotationRateIo, initialTime, "J2000", "IAU_Io" ) );       

//     return bodies;
 
// }

// template< typename TimeType = double, typename StateScalarType = double >
// std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
//            std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
// executeIoRotationDeformationSimulation( 
//     const Eigen::Matrix< StateScalarType, 12, 1 > initialStateDifference,
//     Eigen::Matrix< StateScalarType, 12, 1 >& appliedStateDifference,
//     const Eigen::VectorXd parameterPerturbation = Eigen::VectorXd::Zero( 8 ),
//     const bool propagateVariationalEquations = 1 )
// {
//     double initialTime = 0.0;
//     double finalTime = 86400.0 / 10.0;
//     int numberOfParametersToEstimate = 8;

//     // Create Io-Jupiter dynamical environment
//     SystemOfBodies bodies = createIoDynamicalEnvironment( initialTime, finalTime );

//     // Retrieve Io and Jupiter's Body objects
//     std::shared_ptr< Body > io = bodies.at( "Io" );
//     std::shared_ptr< Body > jupiter = bodies.at( "Jupiter" );

//     // Update Io and Jupiter to initial time
//     io->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );
//     io->setStateFromEphemeris<>( initialTime );
//     jupiter->setStateFromEphemeris<>( initialTime );


//     // Retrieve relevant dynamical parameters
//     double radiusIo = std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( io->getGravityFieldModel( ) )->getReferenceRadius( );
//     double muIo = io->getGravitationalParameter();
//     double muJupiter = jupiter->getGravitationalParameter();
//     double muEffective = muIo + muJupiter;
//     double distanceIo = io->getState( ).segment( 0, 3 ).norm( );
//     double rotationRateIo = std::sqrt( muEffective / ( 4.2e8 * 4.2e8 * 4.2e8 ) );
//     std::cout << "distanceIo " << distanceIo << std::endl;


//     Eigen::Quaterniond noRotationQuaternion = Eigen::Quaterniond( Eigen::AngleAxisd( 1.0E-0, Eigen::Vector3d::UnitZ( ) ) *
//                                                                   Eigen::AngleAxisd( 2.0E-0, Eigen::Vector3d::UnitX( ) ) *
//                                                                   Eigen::AngleAxisd( -0.5E-0, Eigen::Vector3d::UnitZ( ) ) );
    // Eigen::Matrix< double, 7, 1 > unitRotationState = Eigen::Matrix< double, 7, 1 >::Zero( );
    // unitRotationState( 0 ) = noRotationQuaternion.w( );
    // unitRotationState( 1 ) = noRotationQuaternion.x( );
    // unitRotationState( 2 ) = noRotationQuaternion.y( );
    // unitRotationState( 3 ) = noRotationQuaternion.z( );
    // unitRotationState( 4 ) = 0.0;
    // unitRotationState( 5 ) = 0.0;
    // unitRotationState( 6 ) = rotationRateIo + 1.0E-6;

    // Eigen::Matrix< double, 7, 1 > originalRotationState = unitRotationState;
    // Eigen::Matrix< double, 7, 1 > stateDifferenceToAdd = initialStateDifference.segment( 0, 7 );

    // unitRotationState += stateDifferenceToAdd;
    // unitRotationState( 0 ) = originalRotationState( 0 ) / std::fabs( originalRotationState( 0 ) ) *
    //         std::sqrt( 1.0 - std::pow( unitRotationState.segment( 1, 3 ).norm( ), 2.0 ) );

    // appliedStateDifference.segment( 0, 7 ) = unitRotationState - originalRotationState;
    // appliedStateDifference.segment( 7, 5 ) = initialStateDifference.segment( 7, 5 );

//     std::map< double, Eigen::Matrix< double, 7, 1 > > dummyRotationMap;
//     dummyRotationMap[ -1.0E100 ] = unitRotationState;
//     dummyRotationMap[ 1.0E100 ] = unitRotationState;

//     std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > dummyInterpolator =
//             std::make_shared< interpolators::LinearInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( dummyRotationMap );
//     bodies.at( "Io" )
//             ->setRotationalEphemeris(
//                     std::make_shared< TabulatedRotationalEphemeris< double, double > >( dummyInterpolator, "J2000", "IAU_Io" ) );

    
    

//     // SelectedAccelerationMap accelerationMap;
//     // std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
//     // //    accelerationMap[ "Phobos" ][ "Mars" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
//     // accelerationMap[ "Phobos" ][ "Mars" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 2, 2 ) );

//     // std::vector< std::string > translationalBodiesToIntegrate;
//     // std::vector< std::string > translationalCentralBodies;

//     // translationalBodiesToIntegrate.push_back( "Phobos" );
//     // translationalCentralBodies.push_back( "Mars" );

//     // AccelerationMap accelerationModelMap =
//     //         createAccelerationModelsMap( bodies, accelerationMap, translationalBodiesToIntegrate, translationalCentralBodies );

//     // Define integrator settings
//     double timeStep = 10.0;
//     // std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > > ( 
//         // initialTime, timeStep, rungeKutta87DormandPrince, timeStep, timeStep );  
//     std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< IntegratorSettings< > >( rungeKutta4, 0.0, timeStep );

//      // Create torque models
//     std::vector< std::string > bodiesToPropagate = { "Io" }; 

//     SelectedTorqueMap torqueSettings;
//     // torqueSettings[ "Io" ][ "Jupiter" ].push_back( std::make_shared< TorqueSettings >( basic_astrodynamics::second_order_gravitational_torque ) );
//     // torqueSettings[ "Io" ][ "Jupiter" ].push_back( std::make_shared<SphericalHarmonicTorqueSettings>(2,2) );
//     basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap( bodies, torqueSettings, bodiesToPropagate );

//     Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRotationState = getInitialRotationalStateOfBody( "Io", "J2000",  bodies, initialTime );
//     initialRotationState[6] = rotationRateIo;

//     // Create propagator settings for rotational dynamics
//     std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
//             std::make_shared< RotationalStatePropagatorSettings< double > >
//             ( torqueModelMap, bodiesToPropagate, initialRotationState, initialTime, integratorSettings, 
//             std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );


//     // Create deformation model
//     double maxwellRelaxationTime = 179103.0;
//     double globalRelaxationTime = 24688.0;
//     double fluidLoveNumber = 1.5;
//     std::vector< std::string > perturbingBody = {"Jupiter"};
//     std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::make_shared< MaxwellDeformationSettings >( 
//             maxwellRelaxationTime, globalRelaxationTime, fluidLoveNumber, 2, 2, perturbingBody );

//     std::map< std::string, std::vector< std::shared_ptr< GravityDeformationSettings > > > gravityDeformationModelMap;   
//     gravityDeformationModelMap[ "Io" ] = { maxwellDeformationSettings };

//     basic_astrodynamics::GravityDeformationModelMap gravityDeformationModels = createGravityDeformationModelsMap(
//         bodies, gravityDeformationModelMap );

//     // Initialise gravity state    
//     Eigen::VectorXd equilibriumCoefficients = Eigen::VectorXd::Zero(5);
//     // double distanceIo = initialTranslationalState.segment(0, 3).norm();
//     double radiusRatioPowerThree = ( radiusIo / distanceIo ) * ( radiusIo / distanceIo ) * ( radiusIo / distanceIo );
//     double gravitationalParametersRatio = muJupiter / muIo;
    
//     equilibriumCoefficients[ 0 ] = fluidLoveNumber / 2.0 * gravitationalParametersRatio * radiusRatioPowerThree 
//                 * ( 3.0 * std::sin( 0.0 ) * std::sin( 0.0 ) - 1.0 ); 
//     equilibriumCoefficients[ 2 ] = fluidLoveNumber / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
//                 ( 1.0 - std::sin( 0.0 ) * std::sin( 0.0 ) ) * std::cos( 2.0 * 0.0 );
//     equilibriumCoefficients[ 4 ] = fluidLoveNumber / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
//                 ( 1.0 - std::sin( 0.0 ) * std::sin( 0.0 ) ) * std::sin( 2.0 * 0.0 );
        
//     equilibriumCoefficients[ 1 ] = - fluidLoveNumber * gravitationalParametersRatio * radiusRatioPowerThree 
//                 * ( - std::cos( 0.0 ) * std::sin( 0.0 ) ) * std::cos( 0.0 );
//     equilibriumCoefficients[ 3 ] = - fluidLoveNumber * gravitationalParametersRatio * radiusRatioPowerThree 
//                 * ( - std::cos( 0.0 ) * std::sin( 0.0 ) ) * std::sin( 0.0 );

//     Eigen::VectorXd initialGravityCoefficients = equilibriumCoefficients; // Eigen::VectorXd::Zero( 5 );
//         // Eigen::VectorXd perturbedInitialCoefficients = Eigen::VectorXd::Zero( 5 );
//         // perturbedInitialCoefficients[0] = 1.0e-2;
//         // perturbedInitialCoefficients[2] = 1.0e-2;
//         // perturbedInitialCoefficients[4] = 1.0e-2;

//     // Create propagator settings for gravity deformation    
//     std::shared_ptr< GravityDeformationPropagatorSettings< > > gravityPropagatorSettings = 
//         std::make_shared< GravityDeformationPropagatorSettings< > >( bodiesToPropagate, gravityDeformationModels, initialGravityCoefficients, integratorSettings,
//         std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

    

//     // // Define propagator settings.
//     // std::vector< std::string > bodiesToIntegrate;
//     // bodiesToIntegrate.push_back( "Phobos" );

//     // // Create torque models
//     // basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap( bodies, torqueMap, bodiesToIntegrate );

//     // std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
//     //         std::make_shared< RotationalStatePropagatorSettings< double > >(
//     //                 torqueModelMap,
//     //                 bodiesToIntegrate,
//     //                 unitRotationState,
//     //                 std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime ) );

//     // Eigen::VectorXd initialTranslationalState;
//     // initialTranslationalState =
//     //         getInitialStatesOfBodies( translationalBodiesToIntegrate, translationalCentralBodies, bodies, initialEphemerisTime );

//     // initialTranslationalState += initialStateDifference.segment( 0, 6 );
//     // std::shared_ptr< TranslationalStatePropagatorSettings<> > translationalPropagatorSettings =
//     //         std::make_shared< TranslationalStatePropagatorSettings<> >( translationalCentralBodies,
//     //                                                                     accelerationModelMap,
//     //                                                                     translationalBodiesToIntegrate,
//     //                                                                     initialTranslationalState,
//     //                                                                     finalEphemerisTime,
//     //                                                                     cowell );

//     std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
//     propagatorSettingsList.push_back( rotationalPropagatorSettings );
//     propagatorSettingsList.push_back( gravityPropagatorSettings );

//     std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings = std::make_shared< MultiTypePropagatorSettings< double > >(
//             propagatorSettingsList, std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

//     // // Create integrator settings
//     // std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
//     //         std::make_shared< IntegratorSettings< TimeType > >( rungeKutta4, TimeType( initialEphemerisTime ), 15.0 );

//     // Define parameters.
//     std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
//     {
//         parameterNames = getInitialStateParameterSettings< double, double >( propagatorSettings, bodies );

//         // parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Io", mean_moment_of_inertia ) );
//         // parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
//         //         1, 0, 2, 2, "Io", spherical_harmonics_cosine_coefficient_block ) );
//         // parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
//         //         2, 1, 2, 2, "Io", spherical_harmonics_sine_coefficient_block ) );
//     }

//     // Create parameters
//     std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
//             createParametersToEstimate( parameterNames, bodies );
//     // std::cout << "********************************************* " << std::endl;
//     printEstimatableParameterEntries( parametersToEstimate );

//     // Eigen::MatrixXd constraintStateMultiplier;
//     // Eigen::VectorXd constraintRightHandSide;
//     // parametersToEstimate->getConstraints( constraintStateMultiplier, constraintRightHandSide );
//     // //    std::cout<<"Unit rotation: "<<std::endl<<unitRotationState.transpose( )<<std::endl;
//     // //    std::cout<<"Constraints: "<<std::endl<<constraintStateMultiplier.transpose( )<<std::endl;

//     // TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( constraintStateMultiplier.block( 0, 0, 1, 4 ) ),
//     //                                    ( unitRotationState.segment( 0, 4 ) ).transpose( ),
//     //                                    std::numeric_limits< double >::epsilon( ) );
//     // TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( constraintStateMultiplier.block( 0, 4, 1, 9 ) ),
//     //                                    ( Eigen::MatrixXd::Zero( 1, 9 ) ),
//     //                                    std::numeric_limits< double >::epsilon( ) );
//     // TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//     //         ( constraintRightHandSide.block( 0, 0, 1, 1 ) ), ( Eigen::MatrixXd::Zero( 1, 1 ) ), std::numeric_limits< double >::epsilon( ) );

//     // // Perturb parameters.
//     // Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > parameterVector =
//     //         parametersToEstimate->template getFullParameterValues< StateScalarType >( );
//     // parameterVector.block( 13, 0, numberOfParametersToEstimate, 1 ) += parameterPerturbation;
//     // //    std::cout<<"Parameter perturbation "<<
//     // //               ( parametersToEstimate->template getFullParameterValues< StateScalarType >( ) -
//     // //                 parameterVector ).transpose( )<<std::endl<<
//     // //               parameterPerturbation.transpose( )<<std::endl;;
//     // parametersToEstimate->resetParameterValues( parameterVector );

//     std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
//                std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > results;

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

//             //             const simulation_setup::SystemOfBodies& bodies,
//             // const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
//             // const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
//             // const bool integrateDynamicalAndVariationalEquationsConcurrently = true,
//             // const bool integrateEquationsOnCreation = true ):

//         // Propagate requested equations.
//         if( propagateVariationalEquations )
//         {
//             dynamicsSimulator.integrateVariationalAndDynamicalEquations( propagatorSettings->getInitialStates( ), 1 );
//         }
//         else
//         {
//             dynamicsSimulator.integrateDynamicalEquationsOfMotionOnly( propagatorSettings->getInitialStates( ) );
//         }

//         //        tudat::input_output::writeDataMapToTextFile(
//         //                    dynamicsSimulator.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( ),
//         //                    "rotPropTest.dat" );

//         // Retrieve test data
//         double testEpoch = initialTime + 2.0 * 3600.0;
//         Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > testStates = Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( 12 );
//         testStates.segment( 0, 7 ) = bodies.at( "Io" )->getRotationalEphemeris( )->getRotationStateVector( testEpoch );
//         testStates.segment( 7, 5 ) = bodies.at( "Io" )->getPropagatedDegreeTwoCoefficients( );
//         std::cout << "testStates " << testStates.transpose() << std::endl;

//         if( propagateVariationalEquations )
//         {
//             results.first.push_back(
//                     dynamicsSimulator.getStateTransitionMatrixInterface( )->getCombinedStateTransitionAndSensitivityMatrix( testEpoch ) );
//             Eigen::MatrixXd testMatrixDirect =
//                     dynamicsSimulator.getStateTransitionMatrixInterface( )->getCombinedStateTransitionAndSensitivityMatrix( testEpoch );
//             Eigen::MatrixXd testMatrixFull =
//                     dynamicsSimulator.getStateTransitionMatrixInterface( )->getFullCombinedStateTransitionAndSensitivityMatrix( testEpoch );
//             TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testMatrixDirect, testMatrixFull, std::numeric_limits< double >::epsilon( ) );
//         }
//         results.second.push_back( testStates );
//     }
//     return results;
// }

// BOOST_AUTO_TEST_CASE( testIoRotationDeformationVariationalEquationCalculation )
// {
//     // Load spice kernels.
//     spice_interface::loadStandardSpiceKernels( );

//     std::pair< std::vector< Eigen::MatrixXd >, std::vector< Eigen::VectorXd > > currentOutput;

//     // Define variables for numerical differentiation
//     Eigen::Matrix< double, 12, 1 > perturbedState;
//     Eigen::Matrix< double, 12, 1 > statePerturbation;

//     // Define parameter perturbation
//     int numberOfParametersToEstimate = 8;
//     double sphericalHarmonicsPerturbation = 1.0E-4;
//     Eigen::Matrix< double, 8, 1 > perturbedParameter;
//     Eigen::Matrix< double, 8, 1 > parameterPerturbation;
//     parameterPerturbation = Eigen::Matrix< double, 8, 1 >::Constant( sphericalHarmonicsPerturbation );
//     parameterPerturbation( 0 ) = 1.0E-4;

//     // Compute state transition and sensitivity matrices
//     Eigen::Matrix< double, 12, 1 > appliedStateDifference;
//     currentOutput = executeIoRotationDeformationSimulation< double, double >(  Eigen::Matrix< double, 12, 1 >::Zero( ), appliedStateDifference );
//     Eigen::MatrixXd stateTransitionAndSensitivityMatrixAtEpoch = currentOutput.first.at( 0 );
//     Eigen::VectorXd nominalState = currentOutput.second.at( 0 );

//     //    std::cout<<"Nominal "<<std::endl<<std::endl<<
//     //               stateTransitionAndSensitivityMatrixAtEpoch<<std::endl;
//     // Define state perturbation
//     statePerturbation <<
//         1.0E-5,
//         1.0E-5,
//         1.0E-5,
//         1.0E-5,
//         1.0E-7,
//         1.0E-7,
//         1.0E-7,
//         0.0,
//         0.0,
//         0.0,
//         0.0,
//         0.0;

//     Eigen::MatrixXd manualPartial = Eigen::MatrixXd::Zero( 12, 12 + numberOfParametersToEstimate );

//     // Numerically compute state transition matrix
//     for( unsigned int test = 0; test <= 1; test++ )
//     {
//         double perturbationMultiplier = ( test == 0 ? 1.0 : 1.0E-3 );
//         for( unsigned int j = 0; j < 4; j++ )
//         {
//             Eigen::Matrix< double, 12, 1 > appliedStateDifferenceUp, appliedStateDifferenceDown;

//             Eigen::VectorXd upPerturbedState, downPerturbedState;
//             perturbedState.setZero( );
//             perturbedState( j ) += perturbationMultiplier * statePerturbation( j );
//             upPerturbedState = executeIoRotationDeformationSimulation< double, double >(
//                                        perturbedState, appliedStateDifferenceUp, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 )
//                                        .second.at( 0 );

//             std::cout << "upPerturbedState " <<    upPerturbedState.transpose() << std::endl; 
//             std::cout << "nominalState " <<    nominalState.transpose() << std::endl;                        
//             Eigen::VectorXd stateDifferenceUp = upPerturbedState - nominalState;

//             //            std::cout<<"Test output "<<test<<" "<<j<<"stateDifferenceUp"<<std::endl<<
//             //                       stateDifferenceUp<<std::endl<<std::endl<<
//             //                       "stateTransitionAndSensitivityMatrixAtEpoch"<<std::endl<<
//             //                                              stateTransitionAndSensitivityMatrixAtEpoch<<std::endl<<std::endl<<
//             //                       "appliedStateDifferenceUp"<<std::endl<<
//             //                                              appliedStateDifferenceUp<<std::endl<<std::endl<<
//             //                       "( stateTransitionAndSensitivityMatrixAtEpoch * appliedStateDifferenceUp )"<<std::endl<<
//             //                                              ( stateTransitionAndSensitivityMatrixAtEpoch * appliedStateDifferenceUp
//             //                                              )<<std::endl<<std::endl;
//             if( test == 0 )
//             {
//                 Eigen::VectorXd testMatrix =
//                         ( stateTransitionAndSensitivityMatrixAtEpoch.block( 0, 0, 12, 12 ) * appliedStateDifferenceUp );
//                 std::cout << "testMatrix.segment( 7, 5 )" << std::endl;
//                 std::cout << testMatrix.segment( 7, 5 ) << std::endl;
//                 TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( testMatrix.segment( 7, 5 ) ), ( stateDifferenceUp.segment( 7, 5 ) ), 1.5E-3 );
//                 std::cout << "( stateDifferenceUp.segment( 7, 5 ) )" << std::endl;
//                 std::cout << ( stateDifferenceUp.segment( 7, 5 ) ) << std::endl;
//             }
//             else
//             {
//                 Eigen::VectorXd testMatrix =
//                         ( stateTransitionAndSensitivityMatrixAtEpoch.block( 0, 0, 12, 12 ) * appliedStateDifferenceUp );
//                 TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( testMatrix.segment( 0, 7 ) ), ( stateDifferenceUp.segment( 0, 7 ) ), 1.0E-5 );
//             }
//         }
//     }

    // for( unsigned int j = 0; j < 6; j++ )
    // {
    //     Eigen::Matrix< double, 12, 1 > appliedStateDifferenceUp, appliedStateDifferenceDown;

    //     Eigen::VectorXd upPerturbedState, downPerturbedState;
    //     perturbedState.setZero( );
    //     perturbedState( j ) += statePerturbation( j );
    //     upPerturbedState = executePhobosRotationSimulation< double, double >(
    //                                perturbedState, appliedStateDifferenceUp, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 )
    //                                .second.at( 0 );

    //     perturbedState.setZero( );
    //     perturbedState( j ) -= statePerturbation( j );
    //     downPerturbedState = executePhobosRotationSimulation< double, double >(
    //                                  perturbedState, appliedStateDifferenceDown, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 )
    //                                  .second.at( 0 );

    //     manualPartial.block( 0, j, 13, 1 ) =
    //             ( upPerturbedState.segment( 0, 13 ) - downPerturbedState.segment( 0, 13 ) ) / ( 2.0 * statePerturbation( j ) );
    // }

    // TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
    //         ( manualPartial.block( 0, 0, 13, 6 ) ), ( stateTransitionAndSensitivityMatrixAtEpoch.block( 0, 0, 13, 6 ) ), 1.0E-4 );

    // for( unsigned int j = 10; j < 13; j++ )
    // {
    //     Eigen::Matrix< double, 13, 1 > appliedStateDifferenceUp, appliedStateDifferenceDown;

    //     Eigen::VectorXd upPerturbedState, downPerturbedState;
    //     perturbedState.setZero( );
    //     perturbedState( j ) += statePerturbation( j );
    //     upPerturbedState = executePhobosRotationSimulation< double, double >(
    //                                perturbedState, appliedStateDifferenceUp, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 )
    //                                .second.at( 0 );

    //     perturbedState.setZero( );
    //     perturbedState( j ) -= statePerturbation( j );
    //     downPerturbedState = executePhobosRotationSimulation< double, double >(
    //                                  perturbedState, appliedStateDifferenceDown, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 )
    //                                  .second.at( 0 );

    //     manualPartial.block( 0, j, 13, 1 ) =
    //             ( upPerturbedState.segment( 0, 13 ) - downPerturbedState.segment( 0, 13 ) ) / ( 2.0 * statePerturbation( j ) );
    // }

    // // Check element separately
    // BOOST_CHECK_SMALL( std::fabs( manualPartial( 3, 1 + 10 ) - stateTransitionAndSensitivityMatrixAtEpoch( 3, 1 + 10 ) ), 1.0E-2 );
    // manualPartial( 3, 1 + 10 ) = stateTransitionAndSensitivityMatrixAtEpoch( 3, 1 + 10 );

    // TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
    //         ( manualPartial.block( 0, 10, 6, 3 ) ), ( stateTransitionAndSensitivityMatrixAtEpoch.block( 0, 10, 6, 3 ) ), 1.0E-3 );
    // TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
    //         ( manualPartial.block( 6, 10, 7, 3 ) ), ( stateTransitionAndSensitivityMatrixAtEpoch.block( 6, 10, 7, 3 ) ), 1.0E-5 );

    // // Numerically compute sensitivity matrix
    // for( int j = 0; j < numberOfParametersToEstimate; j++ )
    // {
    //     Eigen::Matrix< double, 13, 1 > appliedStateDifference;

    //     Eigen::VectorXd upPerturbedState, downPerturbedState;
    //     perturbedState.setZero( );
    //     perturbedParameter.setZero( );
    //     perturbedParameter( j ) += parameterPerturbation( j );

    //     //        std::cout<<"Test "<<j<<" "<<perturbedParameter.transpose( )<<std::endl;

    //     upPerturbedState =
    //             executePhobosRotationSimulation< double, double >( perturbedState, appliedStateDifference, perturbedParameter, 0 )
    //                     .second.at( 0 );

    //     perturbedParameter.setZero( );
    //     perturbedParameter( j ) -= parameterPerturbation( j );
    //     downPerturbedState =
    //             executePhobosRotationSimulation< double, double >( perturbedState, appliedStateDifference, perturbedParameter, 0 )
    //                     .second.at( 0 );

    //     manualPartial.block( 0, j + 13, 13, 1 ) =
    //             ( upPerturbedState.segment( 0, 13 ) - downPerturbedState.segment( 0, 13 ) ) / ( 2.0 * parameterPerturbation( j ) );
    // }
    // //    std::cout<<manualPartial<<std::endl<<std::endl
    // //            <<stateTransitionAndSensitivityMatrixAtEpoch<<std::endl<<std::endl<<
    // //              ( manualPartial - stateTransitionAndSensitivityMatrixAtEpoch ).cwiseQuotient(
    // //                  stateTransitionAndSensitivityMatrixAtEpoch )<<std::endl;

    // // Check three values separately: could not find perturbations for which all partials are sufficiently within the linear regime

    // BOOST_CHECK_SMALL( std::fabs( manualPartial( 4, 5 + 13 ) - stateTransitionAndSensitivityMatrixAtEpoch( 4, 5 + 13 ) ), 1.0E-4 );
    // BOOST_CHECK_SMALL( std::fabs( manualPartial( 11, 4 + 13 ) - stateTransitionAndSensitivityMatrixAtEpoch( 11, 4 + 13 ) ), 1.0E-5 );
    // BOOST_CHECK_SMALL( std::fabs( manualPartial( 12, 3 + 13 ) - stateTransitionAndSensitivityMatrixAtEpoch( 12, 3 + 13 ) ), 1.0E-2 );

    // manualPartial( 4, 5 + 13 ) = stateTransitionAndSensitivityMatrixAtEpoch( 4, 5 + 13 );
    // manualPartial( 11, 4 + 13 ) = stateTransitionAndSensitivityMatrixAtEpoch( 11, 4 + 13 );
    // manualPartial( 12, 3 + 13 ) = stateTransitionAndSensitivityMatrixAtEpoch( 12, 3 + 13 );
    // std::cout << manualPartial << std::endl
    //           << std::endl
    //           << ( manualPartial - stateTransitionAndSensitivityMatrixAtEpoch ).cwiseQuotient( stateTransitionAndSensitivityMatrixAtEpoch )
    //           << std::endl;
    // TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
    //         ( manualPartial.block( 0, 13, 13, 8 ) ), ( stateTransitionAndSensitivityMatrixAtEpoch.block( 0, 13, 13, 8 ) ), 7.5E-3 );
// }


// BOOST_AUTO_TEST_CASE( test_TranslationalDeformation )
// {
//      std::cout.precision( 20 );

//    // Load spice kernels.
//     spice_interface::loadStandardSpiceKernels( );

//     // Specify initial time
//     double initialTime = 0.0;
//     double finalTime = 1.0 * physical_constants::JULIAN_DAY / 10.0; 

//     std::string globalFrameOrigin = "Jupiter";
//     std::string globalFrameOrientation = "J2000";

//     std::vector< std::string > bodiesToCreate = { "Jupiter", "Io" }; 

//     // Get body settings.
//     BodyListSettings bodySettings =
//             getDefaultBodySettings( bodiesToCreate, initialTime - 86400.0, finalTime + 86400.0, globalFrameOrigin, globalFrameOrientation );
            
//         bodySettings.at( "Jupiter" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ), "SSB", globalFrameOrientation ); 
//         // bodySettings.at( "Jupiter" )->gravityFieldSettings = get_gravitational_field( planet, 'IAU_Jupiter' )

//         // Set Jupiter's rotation as constant (no precession)
//         double rightAscensionPole = ( 358.054324066462 - 90.0 ) * mathematical_constants::PI / 180.0;
//         double declinationPole = ( 90.0 - 25.5034135739821 ) * mathematical_constants::PI / 180.0;
//         double primeMeridian = ( 284.95 ) * mathematical_constants::PI / 180.0;
//         double rotationRateJupiter = ( 870.536 * mathematical_constants::PI / 180.0 ) / 86400.0;

//         // bodySettings.get( "Io" )->rotationModelSettings = simulation_setup::synchronousRotationModelSettings( "Jupiter", "J2000", "IAU_Io" );

//         double muIo = 5959924010272.5136719;
//         double muJupiter = 126686534196012800.0;
//         double muEffective = 126692494120023072.0;

//         double orbitalPeriodIo = 2.0 * mathematical_constants::PI * std::sqrt( 4.2e8 * 4.2e8 * 4.2e8 / muIo );
//         double rotationRateIo = std::sqrt( muEffective / ( 4.2e8 * 4.2e8 * 4.2e8 ) );

        
//         Eigen::Vector6d initialKeplerianState = ( Eigen::Vector6d( ) << 4.2e8, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( );
//         // std::shared_ptr< KeplerEphemerisSettings > keplerEphemerisSettings = std::make_shared< KeplerEphemerisSettings >( initialKeplerianState, 0.0, muEffective, "Jupiter", "J2000" );
//         // bodySettings.at( "Io" )->ephemerisSettings = keplerEphemerisSettings;

//         bodySettings.at( "Io" )->ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
//             ( Eigen::Vector6d( ) << 1.0 * 421.8E6, 1.0 * 0.004, 0.0, 0.0, 0.0, 0.0 ).finished( ),
//             0.0,
//             getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Io" ),
//             "Jupiter",
//             "J2000" );

//         // Create bodies needed in simulation
//         SystemOfBodies bodies = createSystemOfBodies( bodySettings );

//         bodies.at("Jupiter")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//                 rightAscensionPole, declinationPole, primeMeridian, rotationRateJupiter, initialTime, globalFrameOrientation, "IAU_Jupiter" ) );

//         bodies.at("Io")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//                 0.0, mathematical_constants::PI / 180.0, 0.0, rotationRateIo, initialTime, globalFrameOrientation, "IAU_Io" ) );

//         double scaledMeanMomentOfInertia = 0.37685;
//         std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )->setScaledMeanMomentOfInertia( scaledMeanMomentOfInertia );

//         double maxwellRelaxationTime = 179103.0;
//         double globalRelaxationTime = 24688.0;
//         double fluidLoveNumber = 1.5;
//         std::vector< std::string > perturbingBody = {"Jupiter"};
//         std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::make_shared< MaxwellDeformationSettings >( 
//                 maxwellRelaxationTime, globalRelaxationTime, fluidLoveNumber, 2, 2, perturbingBody );

    

// //     // orbital period
// //     double orbitalPeriodIo = 2.0 * mathematical_constants::PI * std::sqrt( 4.2e8 * 4.2e8 * 4.2e8 / muIo );
// //     double rotationRateIo = std::sqrt( muEffective / ( 4.2e8 * 4.2e8 * 4.2e8 ) );
// //     std::cout << "rotationRateIo " << rotationRateIo << std::endl;

// // //     Eigen::Matrix3d initialOrientation = Eigen::Matrix3d::Identity( );
// // //     initialOrientation( 0, 0 ) = - 1.0;
// // //     initialOrientation( 1, 1 ) = - 1.0;
// //     bodySettings.get( "Io" )->rotationModelSettings = simulation_setup::synchronousRotationModelSettings( "Jupiter", "J2000", "IAU_Io" );
// // //     bodySettings.at( "Io" )->rotationModelSettings = std::make_shared< simulation_setup::SimpleRotationModelSettings >( 
// //         // "J2000", "IAU_Io", Eigen::Quaterniond( initialOrientation ), initialTime, rotationRateIo );

// //     // Create bodies needed in simulation
// //     SystemOfBodies bodies = createSystemOfBodies( bodySettings );

//     bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );
//     bodies.at( "Io" )->setStateFromEphemeris<>( initialTime );
//     bodies.at( "Jupiter" )->setStateFromEphemeris<>( initialTime );
//     // bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );

// //     double scaledMeanMomentOfInertiaIo = 0.37685;
// //     std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( 
// //         bodies.at( "Io" )->getGravityFieldModel( ) )->setScaledMeanMomentOfInertia( scaledMeanMomentOfInertiaIo );
// //     std::cout << "initial inertia tensor " << std::endl;
// //     std::cout << bodies.at( "Io" )->getGravityFieldModel( )->getInertiaTensor( ) << std::endl;

// //     const double maxwellRelaxationTime = 1000.0;
// //     const double globalRelaxationTime = 2000.0;
// //     const double loveNumber = 0.4;
// //     // const double rotationRate = ( 2.0 * mathematical_constants::PI ) / ( 1.77 * 86400.0 );
// //     const int maximumDegree = 2;
// //     const int maximumOrder = 2;

// //     // std::cout << "initial rotation rate " << rotationRate << std::endl;

// //     std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::make_shared< MaxwellDeformationSettings >( 
// //         maxwellRelaxationTime, globalRelaxationTime, loveNumber, /*rotationRateIo,*/ maximumDegree, maximumOrder, "Jupiter" );

// //     // std::shared_ptr< Body > deformingBody = bodies.at( "Io" );
// //     // std::shared_ptr< Body > perturbingBody = bodies.at( "Jupiter" );
// //     // std::shared_ptr< basic_astrodynamics::MaxwellGravityDeformationModel > maxwellDeformationModel = 
// //     //     createMaxwellGravityFieldDeformationModel( deformingBody, perturbingBody, "Io", "Jupiter", maxwellDeformationSettings );

// //     // std::map< std::string, std::shared_ptr< basic_astrodynamics::GravityDeformationModel > > gravityDeformationModels;
// //     // gravityDeformationModels[ "Io" ] = maxwellDeformationModel;


//     std::map< std::string, std::vector< std::shared_ptr< GravityDeformationSettings > > > gravityDeformationModelMap;   
//     gravityDeformationModelMap[ "Io" ] = { maxwellDeformationSettings };

//     basic_astrodynamics::GravityDeformationModelMap gravityDeformationModels = createGravityDeformationModelsMap(
//         bodies, gravityDeformationModelMap );

// //     std::map< std::string, std::shared_ptr< basic_astrodynamics::GravityDeformationModel > > gravityDeformationModels;
// //     gravityDeformationModels[ "Io" ] = deformationModels.at( "Io" )[ 0 ];



//     double timeStep = 10.0;
//     std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > > ( 
//         initialTime, timeStep, rungeKutta87DormandPrince, timeStep, timeStep );  

//     std::vector< std::string > bodiesToPropagate = { "Io" };
// //     // Eigen::Matrix< double, Eigen::Dynamic, 1 > initialBodyGravity = Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( 3, 1 );

// //     Eigen::Vector3d computedEquilibriumCoefficients = Eigen::Vector3d::Zero( );
// //     std::shared_ptr< SphericalHarmonicsGravityField > shModel = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) );
// //     double distance = 4.2e8;
// //     double radius = shModel->getReferenceRadius( );
// //     double muJup = bodies.at( "Jupiter" )->getGravitationalParameter( );  
// //     double ratioDistanceRadiusPowerThree = radius * radius * radius / ( distance * distance * distance );
// //     double muRatio = muJup / muIo;
// //     computedEquilibriumCoefficients[ 0 ] = - loveNumber * ( rotationRateIo * rotationRateIo * radius * radius * radius / ( 3.0 * muIo ) + 0.5 * muRatio * ratioDistanceRadiusPowerThree );
// //     computedEquilibriumCoefficients[ 1 ] = loveNumber / 4 * muRatio * ratioDistanceRadiusPowerThree;
// //     computedEquilibriumCoefficients[ 2 ] = 0.0;
// //     std::cout << "computedEquilibriumCoefficients " << computedEquilibriumCoefficients.transpose( ) << std::endl;

// //     Eigen::MatrixXd originalCosineMatrix = shModel->getCosineCoefficients();
// //     Eigen::MatrixXd originalSineMatrix = shModel->getSineCoefficients();
// //     Eigen::Vector3d computedInitialCoefficients = Eigen::Vector3d::Zero( );
// //     computedInitialCoefficients[ 0 ] = originalCosineMatrix( 2, 0 );
// //     computedInitialCoefficients[ 1 ] = originalCosineMatrix( 2, 2 );
// //     computedInitialCoefficients[ 2 ] = originalSineMatrix( 2, 2 );
// //     std::cout << "computedInitialCoefficients " << computedInitialCoefficients.transpose( ) << std::endl;

// //     Eigen::Vector3d computedInitialTransientCoefficients = Eigen::Vector3d::Zero( );
// //     computedInitialTransientCoefficients = ( 1.0 / ( globalRelaxationTime - maxwellRelaxationTime ) ) * 
// //         ( globalRelaxationTime * computedInitialCoefficients - maxwellRelaxationTime * computedEquilibriumCoefficients );
// //     // std::cout << "computedInitialTransientCoefficients: " << computedInitialTransientCoefficients.transpose( ) << std::endl;

//         Eigen::VectorXd computedInitialCoefficients = Eigen::VectorXd::Zero( 5 );
//         Eigen::VectorXd perturbedInitialCoefficients = Eigen::VectorXd::Zero( 5 );
//         perturbedInitialCoefficients[0] = 1.0e-2;
//         perturbedInitialCoefficients[2] = 1.0e-2;
//         perturbedInitialCoefficients[4] = 1.0e-2;

//     std::shared_ptr< GravityDeformationPropagatorSettings< > > gravityPropagatorSettings = 
//         std::make_shared< GravityDeformationPropagatorSettings< > >( bodiesToPropagate, gravityDeformationModels, computedInitialCoefficients, integratorSettings,
//         std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

//     // Translational dynamics propagator
//     std::vector< std::string > centralBodies = { "Jupiter" };
//     SelectedAccelerationMap accelerationSettingsMap;
//     accelerationSettingsMap[ "Io" ][ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
//     // accelerationSettingsMap["Io"][ "Jupiter" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 0, 0, 2, 2 ) );

//     std::shared_ptr< SingleArcPropagatorProcessingSettings > outputSettings =
//             std::make_shared< SingleArcPropagatorProcessingSettings >( );
//     outputSettings->setIntegratedResult( false );

//     AccelerationMap accelerationsMap = createAccelerationModelsMap( bodies, accelerationSettingsMap, bodiesToPropagate, centralBodies );
//     Eigen::Vector6d initialState = orbital_element_conversions::convertKeplerianToCartesianElements( initialKeplerianState, muEffective );
//     std::shared_ptr< TranslationalStatePropagatorSettings<  > > translationalPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >( 
//         centralBodies, accelerationsMap, bodiesToPropagate, initialState, initialTime, integratorSettings, 
//         std::make_shared< PropagationTimeTerminationSettings >( finalTime ), cowell );

// //     // Create torque models
// //     SelectedTorqueMap torqueSettings;
// //     torqueSettings[ "Io" ][ "Jupiter" ].push_back(
// //                             std::make_shared< SphericalHarmonicTorqueSettings >( 2, 2 ) );
// //     basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap( bodies, torqueSettings, bodiesToPropagate );

// //     Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRotationState = getInitialRotationalStateOfBody(
// //         "Io", "J2000",  bodies, initialTime );
// //     // std::cout << "initialRotationState " << initialRotationState << std::endl;

// //     // Create propagator settings for rotational dynamics
// //     std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
// //             std::make_shared< RotationalStatePropagatorSettings< double > >
// //             ( torqueModelMap, bodiesToPropagate, initialRotationState, initialTime, integratorSettings, 
// //             std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

//     std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > >  propagatorSettingsList;
//     propagatorSettingsList.push_back( translationalPropagatorSettings );
// //     propagatorSettingsList.push_back( rotationalPropagatorSettings );
//     propagatorSettingsList.push_back( gravityPropagatorSettings );
//     std::shared_ptr< MultiTypePropagatorSettings< > > fullPropagatorSettings = std::make_shared< MultiTypePropagatorSettings< > >(
//             propagatorSettingsList, integratorSettings, initialTime, std::make_shared< PropagationTimeTerminationSettings >( finalTime ),
//             std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ), outputSettings );
    

//     std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
//             getInitialStateParameterSettings< double, double >( fullPropagatorSettings, bodies );
//     // Create parameters
//     std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
//             createParametersToEstimate< double, double >( parameterNames, bodies );
// printEstimatableParameterEntries( parametersToEstimate );

// std::shared_ptr< SingleArcVariationalEquationsSolver< double, double > > variationalEquationsSolver =
//             std::make_shared< SingleArcVariationalEquationsSolver< double, double > >(
//                     bodies, fullPropagatorSettings, parametersToEstimate, true );

// std::map< double, Eigen::VectorXd > results = variationalEquationsSolver->getEquationsOfMotionSolution();
// std::cout << "initial state " << results.begin()->second.transpose() << std::endl;
//     std::cout << "final state " << results.rbegin()->second.transpose() << std::endl;
// // for ( auto it : results )
// // {
// // std::cout << it.second.transpose( ) << std::endl;
// // }

// std::map< double, Eigen::MatrixXd > stateTransitionMatrixHistory = variationalEquationsSolver->getStateTransitionMatrixSolution( );
// std::cout << "stateTransitionMatrixHistory size " << stateTransitionMatrixHistory.size( ) << std::endl;

// std::cout << "initial STM " << std::endl;
// std::cout << stateTransitionMatrixHistory.begin( )->second << std::endl;
// std::cout << "-----------------------" << std::endl;
// std::cout << "final STM " << std::endl;
// std::cout << stateTransitionMatrixHistory.rbegin( )->second << std::endl;

// std::cout << "test " << std::endl;
// Eigen::VectorXd testPerturbation = Eigen::VectorXd::Zero(11);
// testPerturbation[0] = 0.001 * initialState[0];
// // std::cout << "testPerturbation " << testPerturbation.transpose() << std::endl;
// std::cout << ( stateTransitionMatrixHistory.rbegin( )->second * testPerturbation ).transpose( ) << std::endl;

// Eigen::Vector6d perturbedInitialState = initialState;
// perturbedInitialState[0] *= 1.001;

// Eigen::VectorXd fullPerturbedInitialState = Eigen::VectorXd::Zero( 11 );
// fullPerturbedInitialState.segment( 0, 6 ) = perturbedInitialState;
// fullPerturbedInitialState.segment( 6, 5 ) = computedInitialCoefficients;

// // std::cout << "fullPerturbedInitialState " << fullPerturbedInitialState.transpose( ) << std::endl;

// fullPropagatorSettings->resetInitialStates( fullPerturbedInitialState );
// SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, fullPropagatorSettings ); 
// std::map< double, Eigen::VectorXd > results2 = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

// // std::cout << "initial state " << results2.begin()->second.transpose() << std::endl;
// // std::cout << "final state " << results2.rbegin()->second.transpose() << std::endl;
// // for ( auto it : results2 )
// // {
// // std::cout << it.second.transpose( ) << std::endl;
// // }

// std::map< double, Eigen::VectorXd > stateVariation;
// for ( auto it : results )
// {
//         stateVariation[ it.first ] = ( results2.at( it.first ) - it.second );
// }

// std::cout << "final state variation " << std::endl;
// std::cout << stateVariation.rbegin()->second.transpose( ) << std::endl;
// // for ( auto it : stateVariation )
// // {
// //         std::cout << it.second.transpose( ) << std::endl;
// // }

// }


// BOOST_AUTO_TEST_CASE( test_GravityDeformationEstimation_RotationalState )
// {
//      std::cout.precision( 20 );

//    // Load spice kernels.
//     spice_interface::loadStandardSpiceKernels( );

//     // Specify initial time
//     double initialTime = 0.0;
//     double finalTime = 1.0 * physical_constants::JULIAN_DAY / 10.0; 

//     std::string globalFrameOrigin = "Jupiter";
//     std::string globalFrameOrientation = "J2000";

//     std::vector< std::string > bodiesToCreate = { "Jupiter", "Io" }; 

//     // Get body settings.
//     BodyListSettings bodySettings =
//             getDefaultBodySettings( bodiesToCreate, initialTime - 86400.0, finalTime + 86400.0, globalFrameOrigin, globalFrameOrientation );
            
//         bodySettings.at( "Jupiter" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ), "SSB", globalFrameOrientation ); 
//         // bodySettings.at( "Jupiter" )->gravityFieldSettings = get_gravitational_field( planet, 'IAU_Jupiter' )

//         // Set Jupiter's rotation as constant (no precession)
//         double rightAscensionPole = ( 358.054324066462 - 90.0 ) * mathematical_constants::PI / 180.0;
//         double declinationPole = ( 90.0 - 25.5034135739821 ) * mathematical_constants::PI / 180.0;
//         double primeMeridian = ( 284.95 ) * mathematical_constants::PI / 180.0;
//         double rotationRateJupiter = ( 870.536 * mathematical_constants::PI / 180.0 ) / 86400.0;

//         bodySettings.get( "Io" )->rotationModelSettings = simulation_setup::synchronousRotationModelSettings( "Jupiter", "J2000", "IAU_Io" );

//         double muIo = 5959924010272.5136719;
//         double muJupiter = 126686534196012800.0;
//         double muEffective = 126692494120023072.0;

//         double orbitalPeriodIo = 2.0 * mathematical_constants::PI * std::sqrt( 4.2e8 * 4.2e8 * 4.2e8 / muIo );
//         double rotationRateIo = std::sqrt( muEffective / ( 4.2e8 * 4.2e8 * 4.2e8 ) );

        
//         Eigen::Vector6d initialKeplerianState = ( Eigen::Vector6d( ) << 4.2e8, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( );
//         // std::shared_ptr< KeplerEphemerisSettings > keplerEphemerisSettings = std::make_shared< KeplerEphemerisSettings >( initialKeplerianState, 0.0, muEffective, "Jupiter", "J2000" );
//         // bodySettings.at( "Io" )->ephemerisSettings = keplerEphemerisSettings;

//         bodySettings.at( "Io" )->ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
//             ( Eigen::Vector6d( ) << 1.0 * 421.8E6, 1.0 * 0.004, 0.0, 0.0, 0.0, 0.0 ).finished( ),
//             0.0,
//             getBodyGravitationalParameter( "Jupiter" ) + getBodyGravitationalParameter( "Io" ),
//             "Jupiter",
//             "J2000" );

//         // Create bodies needed in simulation
//         SystemOfBodies bodies = createSystemOfBodies( bodySettings );

//         bodies.at("Jupiter")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//                 rightAscensionPole, declinationPole, primeMeridian, rotationRateJupiter, initialTime, globalFrameOrientation, "IAU_Jupiter" ) );

//         bodies.at("Io")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//                 0.0, mathematical_constants::PI / 180.0, 0.0, rotationRateIo, initialTime, globalFrameOrientation, "IAU_Io" ) );

//         double scaledMeanMomentOfInertia = 0.37685;
//         std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )->setScaledMeanMomentOfInertia( scaledMeanMomentOfInertia );
//         double radiusIo = std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Io" )->getGravityFieldModel( ) )->getReferenceRadius( );    

//         Eigen::Vector6d initialTranslationalState = orbital_element_conversions::convertKeplerianToCartesianElements( initialKeplerianState, muEffective );

//         double maxwellRelaxationTime = 179103.0;
//         double globalRelaxationTime = 24688.0;
//         double fluidLoveNumber = 1.5;
//         std::vector< std::string > perturbingBody = {"Jupiter"};
//         std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::make_shared< MaxwellDeformationSettings >( 
//                 maxwellRelaxationTime, globalRelaxationTime, fluidLoveNumber, 2, 2, perturbingBody );

//     bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );
//     bodies.at( "Io" )->setStateFromEphemeris<>( initialTime );
//     bodies.at( "Jupiter" )->setStateFromEphemeris<>( initialTime );
//     // bodies.at( "Io" )->setCurrentRotationalStateToLocalFrameFromEphemeris( initialTime );

//     bodies.at( "Io" )->getMassProperties( )->updateInertiaTensorDerivative( 
//         ( Eigen::Vector5d( ) << 1.0e-8, 2.0e-8, 1.0e-8, -1.5e-8, 2.0e-8 ).finished( ) ); 


//     std::map< std::string, std::vector< std::shared_ptr< GravityDeformationSettings > > > gravityDeformationModelMap;   
//     gravityDeformationModelMap[ "Io" ] = { maxwellDeformationSettings };

//     basic_astrodynamics::GravityDeformationModelMap gravityDeformationModels = createGravityDeformationModelsMap(
//         bodies, gravityDeformationModelMap );



//     double timeStep = 10.0;
//     std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< > > ( 
//         initialTime, timeStep, rungeKutta87DormandPrince, timeStep, timeStep );  

//     // std::shared_ptr< IntegratorSettings< > > integratorSettings = std::make_shared< IntegratorSettings< > >( rungeKutta4, 0.0, timeStep );

//     std::vector< std::string > bodiesToPropagate = { "Io" };

//     Eigen::VectorXd equilibriumCoefficients = Eigen::VectorXd::Zero(5);

//     double distanceIo = initialTranslationalState.segment(0, 3).norm();
//     double radiusRatioPowerThree = ( radiusIo / distanceIo ) * ( radiusIo / distanceIo ) * ( radiusIo / distanceIo );
//     double gravitationalParametersRatio = muJupiter / muIo;
    
//     equilibriumCoefficients[ 0 ] = fluidLoveNumber / 2.0 * gravitationalParametersRatio * radiusRatioPowerThree 
//                 * ( 3.0 * std::sin( 0.0 ) * std::sin( 0.0 ) - 1.0 ); 
//     equilibriumCoefficients[ 2 ] = fluidLoveNumber / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
//                 ( 1.0 - std::sin( 0.0 ) * std::sin( 0.0 ) ) * std::cos( 2.0 * 0.0 );
//     equilibriumCoefficients[ 4 ] = fluidLoveNumber / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
//                 ( 1.0 - std::sin( 0.0 ) * std::sin( 0.0 ) ) * std::sin( 2.0 * 0.0 );
        
//     equilibriumCoefficients[ 1 ] = - fluidLoveNumber * gravitationalParametersRatio * radiusRatioPowerThree 
//                 * ( - std::cos( 0.0 ) * std::sin( 0.0 ) ) * std::cos( 0.0 );
//     equilibriumCoefficients[ 3 ] = - fluidLoveNumber * gravitationalParametersRatio * radiusRatioPowerThree 
//                 * ( - std::cos( 0.0 ) * std::sin( 0.0 ) ) * std::sin( 0.0 );

//         Eigen::VectorXd computedInitialCoefficients = equilibriumCoefficients; // Eigen::VectorXd::Zero( 5 );
//         // Eigen::VectorXd perturbedInitialCoefficients = Eigen::VectorXd::Zero( 5 );
//         // perturbedInitialCoefficients[0] = 1.0e-2;
//         // perturbedInitialCoefficients[2] = 1.0e-2;
//         // perturbedInitialCoefficients[4] = 1.0e-2;



//     // std::shared_ptr< GravityDeformationPropagatorSettings< > > gravityPropagatorSettings = 
//         // std::make_shared< GravityDeformationPropagatorSettings< > >( bodiesToPropagate, gravityDeformationModels, equilibriumCoefficients, integratorSettings,
//         // std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

//     // // Translational dynamics propagator
//     // std::vector< std::string > centralBodies = { "Jupiter" };
//     // SelectedAccelerationMap accelerationSettingsMap;
//     // accelerationSettingsMap[ "Io" ][ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

//     std::shared_ptr< SingleArcPropagatorProcessingSettings > outputSettings =
//             std::make_shared< SingleArcPropagatorProcessingSettings >( );
//     outputSettings->setIntegratedResult( false );

//     // AccelerationMap accelerationsMap = createAccelerationModelsMap( bodies, accelerationSettingsMap, bodiesToPropagate, centralBodies );
//     // std::shared_ptr< TranslationalStatePropagatorSettings<  > > translationalPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >( 
//     //     centralBodies, accelerationsMap, bodiesToPropagate, initialState, initialTime, integratorSettings, 
//     //     std::make_shared< PropagationTimeTerminationSettings >( finalTime ), cowell );

//     // Create torque models
//     SelectedTorqueMap torqueSettings;
//     // torqueSettings[ "Io" ][ "Jupiter" ].push_back( std::make_shared< TorqueSettings >( basic_astrodynamics::second_order_gravitational_torque ) );
//     torqueSettings[ "Io" ][ "Jupiter" ].push_back( std::make_shared<SphericalHarmonicTorqueSettings>(2,2) );
//     basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap( bodies, torqueSettings, bodiesToPropagate );

//     Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRotationState = getInitialRotationalStateOfBody( "Io", "J2000",  bodies, initialTime );
//     initialRotationState[6] = rotationRateIo;
//     // std::cout << "initialRotationState " << initialRotationState << std::endl;

//     // Create propagator settings for rotational dynamics
//     std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
//             std::make_shared< RotationalStatePropagatorSettings< double > >
//             ( torqueModelMap, bodiesToPropagate, initialRotationState, initialTime, integratorSettings, 
//             std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

//     std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > >  propagatorSettingsList;
//     propagatorSettingsList.push_back( rotationalPropagatorSettings );
//     // propagatorSettingsList.push_back( gravityPropagatorSettings );
//     std::shared_ptr< MultiTypePropagatorSettings< > > fullPropagatorSettings = std::make_shared< MultiTypePropagatorSettings< > >(
//             propagatorSettingsList, integratorSettings, initialTime, std::make_shared< PropagationTimeTerminationSettings >( finalTime ),
//             std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ), outputSettings );
    

//     std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
//             getInitialStateParameterSettings< double, double >( fullPropagatorSettings, bodies );
//     // Create parameters
//     std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
//             createParametersToEstimate< double, double >( parameterNames, bodies );
// printEstimatableParameterEntries( parametersToEstimate );

// std::shared_ptr< SingleArcVariationalEquationsSolver< double, double > > variationalEquationsSolver =
//             std::make_shared< SingleArcVariationalEquationsSolver< double, double > >(
//                     bodies, fullPropagatorSettings, parametersToEstimate, true );

// std::map< double, Eigen::VectorXd > results = variationalEquationsSolver->getEquationsOfMotionSolution();
// std::cout << "initial state " << results.begin()->second.transpose() << std::endl;
//     std::cout << "final state " << results.rbegin()->second.transpose() << std::endl;
// // for ( auto it : results )
// // {
// // std::cout << it.second.transpose( ) << std::endl;
// // }

// std::map< double, Eigen::MatrixXd > stateTransitionMatrixHistory = variationalEquationsSolver->getStateTransitionMatrixSolution( );
// std::cout << "stateTransitionMatrixHistory size " << stateTransitionMatrixHistory.size( ) << std::endl;

// std::cout << "initial STM " << std::endl;
// std::cout << stateTransitionMatrixHistory.begin( )->second << std::endl;
// std::cout << "-----------------------" << std::endl;
// std::cout << "final STM " << std::endl;
// std::cout << stateTransitionMatrixHistory.rbegin( )->second << std::endl;

// Eigen::VectorXd fullInitialState = Eigen::VectorXd::Zero( 7 );
// fullInitialState.segment( 0, 7 ) = initialRotationState;
// // fullInitialState.segment( 7, 5 ) = computedInitialCoefficients;

// Eigen::VectorXd testPerturbation = Eigen::VectorXd::Zero(7);
// testPerturbation[0] = 1.0e-8;
// // testPerturbation[6] = 0.001 * initialRotationState[6];
// // std::cout << "testPerturbation " << testPerturbation.transpose() << std::endl;


// // // Perturb quaternion slightly
// // double perturbation = 1.0e-8;
// Eigen::VectorXd perturbedRotationalState = initialRotationState;

// // Eigen::Vector4d perturbedQuaternion = perturbedRotationalState.segment<4>(0);
// // perturbedQuaternion[0] += perturbation;  
// // perturbedQuaternion.normalize();    // re-normalise
// // perturbedRotationalState.segment<4>(0) = perturbedQuaternion;

// // Perturb rotation

// // Define a very small rotation perturbation (in radians) about a chosen axis, e.g. body X
// Eigen::Vector3d deltaTheta = Eigen::Vector3d::UnitZ() * 1.0e-3;

// // Build a small rotation quaternion: q_delta ≈ [1, 0.5*deltaTheta]
// Eigen::Quaterniond q_delta(1.0, 0.5 * deltaTheta.x(), 0.5 * deltaTheta.y(), 0.5 * deltaTheta.z());
// q_delta.normalize();  // ensure unit length, though nearly unnecessary for small angles

// // Apply perturbation (post-multiply for body-frame rotation, pre-multiply for inertial-frame)
// Eigen::Quaterniond rotationToIntegrationFrame = Eigen::Quaterniond( initialRotationState[0], initialRotationState[1], initialRotationState[2], initialRotationState[3] );
// Eigen::Quaterniond perturbedRotation = ( rotationToIntegrationFrame * q_delta ).normalized();

// perturbedRotationalState.segment<4>(0) = linear_algebra::convertQuaternionToVectorFormat( perturbedRotation );

// Eigen::VectorXd rotationalStatePerturbation = perturbedRotationalState - initialRotationState;
// std::cout << "rotationalStatePerturbation " << rotationalStatePerturbation.transpose() << std::endl;

// // Reconstruct full initial state (rotation + deformation)
// Eigen::VectorXd fullPerturbedInitialState = Eigen::VectorXd::Zero( 7 );
// fullPerturbedInitialState.segment( 0, 7 ) = perturbedRotationalState;
// // fullPerturbedInitialState.segment( 7, 5 ) = computedInitialCoefficients;
// std::cout << "fullPerturbedInitialState " << fullPerturbedInitialState.transpose( ) << std::endl;

// Eigen::VectorXd fullStatePerturbation = fullPerturbedInitialState - fullInitialState;

// std::cout << "from state transition matrix" << std::endl;
// std::cout << ( stateTransitionMatrixHistory.rbegin( )->second * fullStatePerturbation ).transpose( ) << std::endl;

// fullPropagatorSettings->resetInitialStates( fullPerturbedInitialState );
// SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, fullPropagatorSettings ); 
// std::map< double, Eigen::VectorXd > results2 = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

// // std::cout << "initial state " << results2.begin()->second.transpose() << std::endl;
// // std::cout << "final state " << results2.rbegin()->second.transpose() << std::endl;
// // for ( auto it : results2 )
// // {
// // std::cout << it.second.transpose( ) << std::endl;
// // }

// std::map< double, Eigen::VectorXd > stateVariation;
// for ( auto it : results )
// {
//         stateVariation[ it.first ] = ( results2.at( it.first ) - it.second );
// }

// std::cout << "numerical variation " << std::endl;
// std::cout << stateVariation.rbegin()->second.transpose( ) << std::endl;
// }



BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
