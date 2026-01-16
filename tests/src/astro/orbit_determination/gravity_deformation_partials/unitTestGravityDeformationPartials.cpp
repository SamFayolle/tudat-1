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
#include <string>
#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"

#include <boost/test/unit_test.hpp>

#include <boost/lambda/lambda.hpp>

#include "tudat/astro/basic_astro/sphericalStateConversions.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/ephemerides/keplerEphemeris.h"
#include "tudat/astro/relativity/metric.h"
#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravitationalParameter.h"
#include "tudat/simulation/estimation_setup/createGravityDeformationPartials.h"
#include "tudat/simulation/estimation_setup/createTorquePartials.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/propagation_setup/createGravityDeformationModels.h"
#include "tudat/simulation/propagation_setup/createTorqueModel.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::relativity;
using namespace tudat::gravitation;
using namespace tudat::aerodynamics;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::unit_conversions;
using namespace tudat::orbit_determination;
using namespace tudat::acceleration_partials;
using namespace tudat::spice_interface;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::electromagnetism;
using namespace tudat::basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_deformation_partials )



Eigen::Vector5d computeEquilibriumCoefficientsFromSphericalCoordinates(
        const double k2,
        const double gravitationalParametersRatio,
        const double gravitationalParameter,
        const double radius,
        const Eigen::Vector6d& sphericalBodyFixedState,
        const Eigen::Vector3d angularVelocityVector = Eigen::Vector3d::Zero( ),
        const Eigen::Vector3d angularVelocityVectorDerivative = Eigen::Vector3d::Zero( ) )
{
        double distance = sphericalBodyFixedState[0];
        double latitude = mathematical_constants::PI / 2.0 - sphericalBodyFixedState[1];
        double longitude = sphericalBodyFixedState[2];

        double radiusRatioPowerThree = ( radius * radius * radius ) / ( distance * distance * distance );

        double angularVelocity = angularVelocityVector.norm();
        
        Eigen::Vector5d equilibriumCoefficients = Eigen::Vector5d::Zero();
        
        equilibriumCoefficients[ 0 ] = k2 * ( 
                - angularVelocity * angularVelocity * radius * radius * radius / ( 3.0 * gravitationalParameter ) 
                + 0.5 * gravitationalParametersRatio * radiusRatioPowerThree 
                * ( 3.0 * std::sin( latitude ) * std::sin( latitude ) - 1.0 ) ); 
        equilibriumCoefficients[ 1 ] = - k2 * gravitationalParametersRatio * radiusRatioPowerThree 
                * ( - std::cos( latitude ) * std::sin( latitude ) ) * std::cos( longitude );
        equilibriumCoefficients[ 2 ] = k2 / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
                ( 1.0 - std::sin( latitude ) * std::sin( latitude ) ) * std::cos( 2.0 * longitude );
        
        
        equilibriumCoefficients[ 3 ] = - k2 * gravitationalParametersRatio * radiusRatioPowerThree 
                * ( - std::cos( latitude ) * std::sin( latitude ) ) * std::sin( longitude );
        equilibriumCoefficients[ 4 ] = k2 / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
                ( 1.0 - std::sin( latitude ) * std::sin( latitude ) ) * std::sin( 2.0 * longitude );
                

        return equilibriumCoefficients;
}

Eigen::Vector5d computeEquilibriumCoefficientsDerivativeFromSphericalCoordinates(
        const double k2,
        const double gravitationalParametersRatio,
        const double gravitationalParameter,
        const double radius,
        const Eigen::Vector6d& sphericalBodyFixedState,
        const Eigen::Vector3d angularVelocityVector = Eigen::Vector3d::Zero( ),
        const Eigen::Vector3d angularVelocityVectorDerivative = Eigen::Vector3d::Zero( ) )
{
        double distance = sphericalBodyFixedState[0];
        double latitude = mathematical_constants::PI / 2.0 - sphericalBodyFixedState[1];
        double longitude = sphericalBodyFixedState[2];

        double distanceDerivative = sphericalBodyFixedState[3];
        double latitudeDerivative = - sphericalBodyFixedState[4];
        double longitudeDerivative = sphericalBodyFixedState[5];

        double radiusRatioPowerThree = ( radius * radius * radius ) / ( distance * distance * distance );

        double angularVelocity = angularVelocityVector.norm();
        double angularVelocityDerivative = 0.0;
        if ( angularVelocity > 0.0 )
        {
                angularVelocityDerivative = ( angularVelocityVector[0] * angularVelocityVectorDerivative[0]
                        + angularVelocityVector[1] * angularVelocityVectorDerivative[1] + angularVelocityVector[2] * angularVelocityVectorDerivative[2] ) / angularVelocity;
        } 

        Eigen::Vector5d equilibriumCoefficientsDerivative = Eigen::Vector5d::Zero();

        equilibriumCoefficientsDerivative[ 0 ] = - k2 * ( 
            2.0 * angularVelocity * radius * radius * radius / ( 3.0 * gravitationalParameter ) * angularVelocityDerivative
            + 1.0 / 2.0 * gravitationalParametersRatio * radiusRatioPowerThree 
            * 3.0 * distanceDerivative / distance * ( 3.0 * std::sin( latitude ) * std::sin( latitude ) - 1.0 )
            - 1.0 / 2.0 * gravitationalParametersRatio * radiusRatioPowerThree 
            * ( 6.0 * latitudeDerivative * std::sin( latitude ) * std::cos( latitude ) )  );

        equilibriumCoefficientsDerivative[ 1 ] = - k2 * gravitationalParametersRatio * radiusRatioPowerThree * (
            3.0 * distanceDerivative / distance 
            * std::sin( latitude ) * std::cos( latitude ) * std::cos( longitude )
            + longitudeDerivative * std::sin( latitude ) * std::cos( latitude ) * std::sin( longitude )
            - latitudeDerivative * 
            ( std::cos( latitude ) * std::cos( latitude ) - std::sin( latitude ) * std::sin( latitude ) ) 
            * std::cos( longitude ) );

        equilibriumCoefficientsDerivative[ 2 ] = - k2 / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * (
            3.0 * distanceDerivative / distance 
            * ( 1.0 - std::sin( latitude ) * std::sin( latitude ) ) * std::cos( 2.0 * longitude )
            + 2.0 * longitudeDerivative * std::sin( 2.0 * longitude ) * ( 1.0 - std::sin( latitude ) * std::sin( latitude ) ) 
            + 2.0 * std::cos( 2.0 * longitude ) * latitudeDerivative * std::sin( latitude ) * std::cos( latitude ) );

        equilibriumCoefficientsDerivative[ 3 ] = - k2 * gravitationalParametersRatio * radiusRatioPowerThree * (
            3.0 * distanceDerivative / distance 
            * std::sin( latitude ) * std::cos( latitude ) * std::sin( longitude )
            - longitudeDerivative * std::sin( latitude ) * std::cos( latitude ) * std::cos( longitude )
            - latitudeDerivative * 
            ( std::cos( latitude ) * std::cos( latitude ) - std::sin( latitude ) * std::sin( latitude ) ) 
            * std::sin( longitude ) );

        equilibriumCoefficientsDerivative[ 4 ] = - k2 / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * (
            3.0 * distanceDerivative / distance
            * ( 1.0 - std::sin( latitude ) * std::sin( latitude ) ) * std::sin( 2.0 * longitude )
            - 2.0 * longitudeDerivative * std::cos( 2.0 * longitude ) * ( 1.0 - std::sin( latitude ) * std::sin( latitude ) ) 
            + 2.0 * std::sin( 2.0 * longitude ) * latitudeDerivative * std::sin( latitude ) * std::cos( latitude ) );

        return equilibriumCoefficientsDerivative;
}

Eigen::Vector5d computeDeformationFromSphericalCoordinates(
        const double k2,
        const double gravitationalParametersRatio,
        const double gravitationalParameter,
        const double radius,
        const double globalRelaxationTime,
        const double maxwellRelaxationTime,
        const Eigen::Vector6d& sphericalBodyFixedState,
        const Eigen::Vector3d angularVelocityVector = Eigen::Vector3d::Zero( ),
        const Eigen::Vector3d angularVelocityVectorDerivative = Eigen::Vector3d::Zero( ) )
{
        Eigen::Vector5d currentCoefficients = Eigen::Vector5d::Zero(); // can be set to zero for the partials testing

        Eigen::Vector5d equilibriumCoefficients = computeEquilibriumCoefficientsFromSphericalCoordinates(
                k2, gravitationalParametersRatio, gravitationalParameter, radius, sphericalBodyFixedState, angularVelocityVector, angularVelocityVectorDerivative );
        Eigen::Vector5d equilibriumCoefficientsDerivative = computeEquilibriumCoefficientsDerivativeFromSphericalCoordinates(
                k2, gravitationalParametersRatio, gravitationalParameter, radius, sphericalBodyFixedState, angularVelocityVector, angularVelocityVectorDerivative );

        Eigen::Vector5d deformation = ( 1.0 / globalRelaxationTime ) * ( equilibriumCoefficients - currentCoefficients + maxwellRelaxationTime * equilibriumCoefficientsDerivative );

        return deformation;
}

Eigen::Vector6d computeSphericalCoordinatesFromBodyFixedState( 
        const Eigen::Vector6d& bodyFixedState )
{
        Eigen::Vector6d sphericalState = Eigen::Vector6d::Zero( );
        Eigen::Vector3d position = bodyFixedState.segment( 0, 3 );
        Eigen::Vector3d velocity = bodyFixedState.segment( 3, 3 );
        
        // Compute spherical position (r, pi/2 - latitude, longitude)
        sphericalState.segment( 0, 3 ) = coordinate_conversions::convertCartesianToSpherical( position );
        
        // Compute spherical "velocity" 
        double distance = position.segment( 0, 3 ).norm( );
        double distanceDerivative = ( position[ 0 ] * velocity[ 0 ] + position[ 1 ] * velocity[ 1 ] + position[ 2 ] * velocity[ 2 ] ) / distance;

        double longitudeDerivative = ( velocity[ 1 ] * position[ 0 ]  - velocity[ 0 ] * position[ 1 ] ) / 
                ( position[ 0 ] * position[ 0 ] + position[ 1 ] * position[ 1 ] ) ;

        double latitudeDerivative = ( velocity[2] * distance - position[2] * distanceDerivative ) /
                ( distance * std::sqrt( position[ 0 ] * position[ 0 ] + position[ 1 ] * position[ 1 ] ) );

        sphericalState[3] = distanceDerivative;
        sphericalState[4] = - latitudeDerivative;
        sphericalState[5] = longitudeDerivative;

        return sphericalState;
}

Eigen::Vector6d computeBodyFixedFromInertialState( 
        const Eigen::Vector6d& inertialState,
        const Eigen::Matrix3d& rotationMatrixToBodyFixed,
        const Eigen::Matrix3d& rotationMatrixToBodyFixedDerivative )
{
        Eigen::Vector6d bodyFixedState = Eigen::Vector6d::Zero();
        bodyFixedState.segment( 0, 3 ) = rotationMatrixToBodyFixed * inertialState.segment( 0, 3 );
        bodyFixedState.segment( 3, 3 ) = rotationMatrixToBodyFixed * inertialState.segment( 3, 3 )
                + rotationMatrixToBodyFixedDerivative * inertialState.segment( 0, 3 );

        return bodyFixedState;
}


// BOOST_AUTO_TEST_CASE( testMaxwellEquilibriumCoefficientsPartials )
// {
//     std::cout.precision(20);

//     double initialTime = 0.0;
        
//     // Load spice kernels.
//     spice_interface::loadStandardSpiceKernels( );

//     // Create bodies
//     SystemOfBodies bodies = SystemOfBodies( "Jupiter", "J2000" );

//     // Create Jupiter
//     bodies.createEmptyBody( "Jupiter", false );

//     // Set Jupiter's ephemeris to
//     bodies.at( "Jupiter" )->setEphemeris(
//             std::make_shared< ephemerides::ConstantEphemeris >( [ = ]( ) { return Eigen::Vector6d::Zero( ); }, "SSB", "J2000" ) );

//     double muJupiter = spice_interface::getBodyGravitationalParameter( "Jupiter" );
//     bodies.at( "Jupiter" )->setGravityFieldModel(
//             std::make_shared< gravitation::GravityFieldModel >( muJupiter ) );

//     // Set Jupiter's rotation model
//     double rightAscensionPole = ( 358.054324066462 - 90.0 ) * mathematical_constants::PI / 180.0;
//     double declinationPole = ( 90.0 - 25.5034135739821 ) * mathematical_constants::PI / 180.0;
//     double primeMeridian = ( 284.95 ) * mathematical_constants::PI / 180.0;
//     double rotationRateJupiter = ( 870.536 * mathematical_constants::PI / 180.0 ) / 86400.0;   
    
//     bodies.at("Jupiter")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//                 rightAscensionPole, declinationPole, primeMeridian, rotationRateJupiter, initialTime, "J2000", "IAU_Jupiter" ) );
    
//     // Create Io
//     bodies.createEmptyBody( "Io" );
//     std::shared_ptr< Body > io = bodies.at( "Io" );
//     std::shared_ptr< Body > jupiter = bodies.at( "Jupiter" );

//     // Set Io gravity field
//     double muIo = spice_interface::getBodyGravitationalParameter( "Io" ); 
//     double radiusIo = 1821.6E3;
//     double scaledMeanMomentOfInertia = 0.37685;
//     Eigen::MatrixXd ioCosineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );
//     Eigen::MatrixXd ioSineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );
//     ioCosineCoefficients( 0, 0 ) = 1.0;
//     ioCosineCoefficients( 2, 0 ) = -1845.9E-6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
//     ioCosineCoefficients( 2, 2 ) = 553.7E-6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
          
//     bodies.at( "Io" )->setGravityFieldModel( std::make_shared< gravitation::SphericalHarmonicsGravityField >( 
//         muIo, radiusIo, ioCosineCoefficients, ioSineCoefficients, "IAU_Io", scaledMeanMomentOfInertia ) );


//     // Set Io's ephemeris
//     double muEffective = muIo + muJupiter;
//     Eigen::Vector6d ioKeplerElements = Eigen::Vector6d::Zero( );
//     double ioSemiMajorAxis = 4.2e8;
//     ioKeplerElements( 0 ) = ioSemiMajorAxis;
//     bodies.at( "Io" )->setEphemeris( std::make_shared< ephemerides::KeplerEphemeris >( ioKeplerElements, 0.0, muEffective, "Jupiter", "J2000" ) );
    
//     double rotationRateIo = std::sqrt( muEffective / ( 4.2e8 * 4.2e8 * 4.2e8 ) );

//     // Two test cases:
//     // test case 0: state-independent rotation model (i.e., simple rotation model)
//     // test case 1: state-dependent rotation model (i.e., synchronous rotation model)

//     for ( unsigned int testCase = 0 ; testCase < 1 ; testCase++ )
//     {
//         if ( testCase == 0 ) // state-independent rotation model
//         {
//                 // Set Io's rotation model to constant
//                 bodies.at("Io")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//                         rightAscensionPole, declinationPole, primeMeridian, rotationRateIo, initialTime, "J2000", "IAU_Io" ) );
//         }
//         else if ( testCase == 1 ) // state dependent rotation model
//         {
//                 // Set Io's rotation model to synchronous
//                 bodies.at( "Io" )->setRotationalEphemeris( std::make_shared< SynchronousRotationalEphemeris >(
//                         createRelativeStateFunction( bodies, "Io", "Jupiter" ), "Jupiter", "J2000", "IAU_Io" ) );
//         }

//         // Update Jupiter and Io to testTime
//         double testTime = 1000.0;
//         io->setStateFromEphemeris( testTime );
//         jupiter->setStateFromEphemeris( testTime );
//         io->setCurrentRotationalStateToLocalFrameFromEphemeris( testTime );

//         // Create Maxwell deformation model 
//         double maxwellRelaxationTime = 179103.0;
//         double globalRelaxationTime = 24688.0;
//         double fluidLoveNumber = 1.5;
//         std::vector< std::string > perturbingBody = { "Jupiter" };
//         std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::make_shared< MaxwellDeformationSettings >( 
//                 maxwellRelaxationTime, globalRelaxationTime, fluidLoveNumber, 2, 2, perturbingBody, Eigen::Vector5d::Zero( ), true, true );

//         std::vector< std::shared_ptr< simulation_setup::Body > > perturbingBodies = { jupiter };
//         std::shared_ptr< basic_astrodynamics::MaxwellGravityDeformationModel > maxwellDeformationModel = createMaxwellGravityFieldDeformationModel(
//                 bodies.at("Io"), perturbingBodies, "Io", std::vector< std::string >( {"Jupiter"} ), maxwellDeformationSettings );

//         // Create parameter objects (necessary to the automatic computation of the rotation matrix partials)
//         Eigen::Vector6d initialTranslationalState = io->getState( );
//         Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRotationState = propagators::getInitialRotationalStateOfBody( "Io", "J2000",  bodies, initialTime );
//         initialRotationState[6] = rotationRateIo;

//         std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames; 
//         parameterNames.push_back( std::make_shared< estimatable_parameters::InitialTranslationalStateEstimatableParameterSettings< double > >( 
//                 "Io", initialTranslationalState, "Jupiter" ) ); 
//         parameterNames.push_back( std::make_shared< estimatable_parameters::InitialRotationalStateEstimatableParameterSettings< double > >( "Io", initialRotationState, "J2000" ) ); 
//         if ( testCase == 0 ) // if constant rotation model 
//         {
//                 parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Io", constant_rotation_rate ) );
//         }
//         std::shared_ptr< EstimatableParameterSet< double > > parameterSet = createParametersToEstimate( parameterNames, bodies );
//         printEstimatableParameterEntries( parameterSet );

//         // Create deformation partial.
//         std::shared_ptr< MaxwellDeformationPartial > deformationPartial = std::dynamic_pointer_cast< MaxwellDeformationPartial >(
//                 createAnalyticalGravityDeformationPartial( maxwellDeformationModel, std::make_pair( "Io", io ), bodies, parameterSet ) );   

//         // Get Io's and Jupiter's state
//         Eigen::Vector6d stateIo = io->getStateInBaseFrameFromEphemeris( testTime );
//         Eigen::Vector6d stateJupiter = jupiter->getStateInBaseFrameFromEphemeris( testTime );
//         Eigen::Vector6d inertialState = stateJupiter - stateIo;        

//         // Update deformation partial        
//         deformationPartial->update( testTime );


//         // Get rotation matrix to body-fixed frame (and its derivative)
//         Eigen::Quaterniond originalQuaternion = io->getCurrentRotationToGlobalFrame( );
//         Eigen::Matrix3d rotationMatrixToBodyFixed = io->getCurrentRotationToGlobalFrame( ).toRotationMatrix( ).inverse( );
//         Eigen::Matrix3d rotationMatrixToBodyFixedDerivative = io->getCurrentRotationMatrixDerivativeToLocalFrame( );
//         Eigen::Vector3d angularVelocityBodyFixed = io->getCurrentAngularVelocityVectorInLocalFrame( );
//         Eigen::Vector3d angularVelocityDerivative = io->getCurrentAngularVelocityDerivativeVectorInLocalFrame( );

//         // std::cout << "rotationMatrixToBodyFixedDerivative" << std::endl;
//         // std::cout << rotationMatrixToBodyFixedDerivative << std::endl;
//         // std::cout << "angularVelocityBodyFixed " << angularVelocityBodyFixed.transpose() << std::endl;
//         // std::cout << "rotationMatrixToBodyFixed" << std::endl;
//         // std::cout << rotationMatrixToBodyFixed << std::endl;

//         // Check partials of body-fixed state wrt inertial state
//         Eigen::Matrix6d partialsBodyFixedWrtInertialState = Eigen::Matrix6d::Zero();    

//         Eigen::Vector6d statePerturbation = ( Eigen::Vector6d( ) << 100.0, 100.0, 100.0, 0.1, 0.1, 0.1 ).finished( );
//         for ( unsigned int i = 0 ; i < 6 ; i++ )
//         {
//                 Eigen::Matrix3d perturbedRotationMatrixToBodyFixed = rotationMatrixToBodyFixed;
//                 Eigen::Matrix3d perturbedRotationMatrixToBodyFixedDerivative = rotationMatrixToBodyFixedDerivative;

//                 // Up-perturb inertial state
//                 Eigen::Vector6d upPerturbedInertialState = inertialState;
//                 upPerturbedInertialState[ i ] += statePerturbation[ i ];

//                 // Compute indirect effect on rotation matrix
//                 if ( testCase == 1 ) // synchronous rotation model test case
//                 {
//                         std::shared_ptr< SynchronousRotationalEphemeris > synchronousRotation = std::dynamic_pointer_cast< SynchronousRotationalEphemeris >(
//                                 io->getRotationalEphemeris( ) );

//                         perturbedRotationMatrixToBodyFixed = ( synchronousRotation->getFullyLockedRotationToBaseFrame( - upPerturbedInertialState, testTime ) 
//                                 * synchronousRotation->getLibrationRotation( - upPerturbedInertialState, testTime ) ).inverse();
//                 }

//                 Eigen::Vector6d upPerturbedBodyFixedState = computeBodyFixedFromInertialState( 
//                         upPerturbedInertialState, perturbedRotationMatrixToBodyFixed, perturbedRotationMatrixToBodyFixedDerivative );


//                 // Down-perturb inertial state        
//                 Eigen::Vector6d downPerturbedInertialState = inertialState;
//                 downPerturbedInertialState[ i ] -= statePerturbation[ i ];

//                 // Compute indirect effect on rotation matrix
//                 if ( testCase == 1 ) // synchronous rotation model test case
//                 {
//                         std::shared_ptr< SynchronousRotationalEphemeris > synchronousRotation = std::dynamic_pointer_cast< SynchronousRotationalEphemeris >(
//                                 io->getRotationalEphemeris( ) );

//                         perturbedRotationMatrixToBodyFixed = ( synchronousRotation->getFullyLockedRotationToBaseFrame( - downPerturbedInertialState, testTime ) 
//                                 * synchronousRotation->getLibrationRotation( - downPerturbedInertialState, testTime ) ).inverse();
//                 }

//                 Eigen::Vector6d downPerturbedBodyFixedState = computeBodyFixedFromInertialState( 
//                         downPerturbedInertialState, perturbedRotationMatrixToBodyFixed, perturbedRotationMatrixToBodyFixedDerivative );

//                 partialsBodyFixedWrtInertialState.block( 0, i, 6, 1 ) = ( upPerturbedBodyFixedState - downPerturbedBodyFixedState ) / ( 2.0 * statePerturbation[i] );
//         }

//         Eigen::Matrix6d analyticalPartialBodyFixedWrtInertial = deformationPartial->bodyFixedWrtGlobalState( testTime, true );
//         std::cout << "partialsBodyFixedWrtInertialState" << std::endl;
//         std::cout << partialsBodyFixedWrtInertialState << std::endl;
//         std::cout << "analyticalPartialBodyFixedWrtInertial" << std::endl;
//         std::cout << analyticalPartialBodyFixedWrtInertial << std::endl;

//         for ( unsigned int i = 0 ; i < 6 ; i++ )
//         {
//                 for ( unsigned int j = 0 ; j < 6 ; j++ )
//                 {
//                         BOOST_CHECK_SMALL( partialsBodyFixedWrtInertialState( i, j ) - analyticalPartialBodyFixedWrtInertial( i, j ), 1.0E-9 );
//                 }
//         }


//         // Check partials of body-fixed state wrt quaternions (body-fixed to inertial)
//         Eigen::Vector6d bodyFixedState = computeBodyFixedFromInertialState( 
//                         inertialState, rotationMatrixToBodyFixed, rotationMatrixToBodyFixedDerivative );

//         Eigen::Vector4d orientationPerturbation = ( Eigen::Vector4d( ) << 1.0E-9, 1.0E-9, 1.0E-9, 1.0E-9 ).finished();
//         std::vector< Eigen::Vector4d > appliedQuaternionPerturbation( 4 );
        
//         Eigen::MatrixXd bodyFixedStateDeviations = Eigen::MatrixXd::Zero( 6, 3 );

//         Eigen::Vector4d originalQuaternionVector = linear_algebra::convertQuaternionToVectorFormat( originalQuaternion );
//         Eigen::Vector4d perturbedQuaternionVector = originalQuaternionVector;

//         // Calculate perturbed body-fixed state
//         for( int i = 1; i < 4; i++ )
//         {
//                 perturbedQuaternionVector( i ) += orientationPerturbation( i );
//                 perturbedQuaternionVector( 0 ) = ( originalQuaternionVector( 0 ) > 0 ? 1.0 : -1.0 ) *
//                         std::sqrt( 1.0 - std::pow( perturbedQuaternionVector.segment( 1, 3 ).norm( ), 2 ) );

//                 appliedQuaternionPerturbation[ i ] = perturbedQuaternionVector.segment( 0, 4 ).normalized( ) - originalQuaternionVector.segment( 0, 4 );

//                 Eigen::Quaterniond perturbedQuaternion = linear_algebra::convertVectorToQuaternionFormat( perturbedQuaternionVector );

//                 Eigen::Matrix3d perturbedRotationMatrixToInertial = perturbedQuaternion.toRotationMatrix( );
//                 Eigen::Matrix3d perturbedRotationMatrixToBodyFixed = perturbedRotationMatrixToInertial.inverse( );

//                 Eigen::Matrix3d perturbedRotationMatrixToBodyFixedDerivative = 
//                         - linear_algebra::getCrossProductMatrix( angularVelocityBodyFixed ) * perturbedRotationMatrixToBodyFixed;

//                 Eigen::Vector6d perturbedBodyFixedState = computeBodyFixedFromInertialState( 
//                         inertialState, perturbedRotationMatrixToBodyFixed, perturbedRotationMatrixToBodyFixedDerivative );

//                 bodyFixedStateDeviations.block( 0, i - 1, 6, 1 ) = perturbedBodyFixedState - bodyFixedState;

//                 perturbedQuaternionVector = originalQuaternionVector;
                        
//         }

//         // std::cout << "appliedQuaternionPerturbation " << std::endl;
//         // std::cout << appliedQuaternionPerturbation[0].transpose() << std::endl;
//         // std::cout << appliedQuaternionPerturbation[1].transpose() << std::endl;
//         // std::cout << appliedQuaternionPerturbation[2].transpose() << std::endl;
//         // std::cout << appliedQuaternionPerturbation[3].transpose() << std::endl;

//         // Compute analytical partials
//         Eigen::MatrixXd analyticalPartialBodyFixedWrtRotationalState = deformationPartial->bodyFixedWrtRotational();
//         Eigen::MatrixXd analyticalPartialBodyFixedWrtQuaternion = analyticalPartialBodyFixedWrtRotationalState.block( 0, 0, 6, 4 );


//         // Compare numerical and analytical results.
//         for ( int index = 1 ; index < 4 ; index++ )
//         {
//                 Eigen::Vector6d numericalChangeInBodyFixedState = bodyFixedStateDeviations.block( 0, index - 1, 6, 1 );
//                 Eigen::Vector6d analyticalChangeInBodyFixedState =
//                         analyticalPartialBodyFixedWrtQuaternion.block( 0, 0, 6, 1 ) * appliedQuaternionPerturbation[ index ]( 0 ) +
//                         analyticalPartialBodyFixedWrtQuaternion.block( 0, index, 6, 1 ) * appliedQuaternionPerturbation[ index ]( index );

//                 std::cout << "numericalChangeInBodyFixedState" << std::endl;
//                 std::cout << numericalChangeInBodyFixedState.transpose() << std::endl;
//                 std::cout << "analyticalChangeInBodyFixedState" << std::endl;
//                 std::cout << analyticalChangeInBodyFixedState.transpose() << std::endl;

//                 for ( unsigned int j = 0 ; j < 6 ; j++ )
//                 {
//                         BOOST_CHECK_SMALL( numericalChangeInBodyFixedState[j] - analyticalChangeInBodyFixedState[j], 1.0E-6 );
//                 }
//         }


//         // Check partials of body-fixed state wrt body-fixed angular velocity and deformation wrt body-fixed angular velocity (due to centrifugal potential)
//         Eigen::MatrixXd partialsBodyFixedWrtAngularVelocity = Eigen::MatrixXd::Zero( 6, 3 );    

//         Eigen::Vector3d angularVelocityPerturbation = ( Eigen::Vector3d( ) << 1.0e-6, 1.0e-6, 1.0e-6 ).finished( );
//         for ( unsigned int i = 0 ; i < 3 ; i++ )
//         {
//                 // Up-perturbation
//                 Eigen::Vector3d upPerturbedAngularVelocity = angularVelocityBodyFixed;
//                 upPerturbedAngularVelocity[ i ] += angularVelocityPerturbation[ i ];

//                 Eigen::Matrix3d upPerturbedRotationMatrixToBodyFixedDerivative = 
//                         - linear_algebra::getCrossProductMatrix( upPerturbedAngularVelocity ) * rotationMatrixToBodyFixed;

//                 Eigen::Vector6d upPerturbedBodyFixedState = computeBodyFixedFromInertialState( 
//                         inertialState, rotationMatrixToBodyFixed, upPerturbedRotationMatrixToBodyFixedDerivative );     

//                 // Down-perturbation        
//                 Eigen::Vector3d downPerturbedAngularVelocity = angularVelocityBodyFixed;
//                 downPerturbedAngularVelocity[ i ] -= angularVelocityPerturbation[ i ];
                
//                 Eigen::Matrix3d downPerturbedRotationMatrixToBodyFixedDerivative = 
//                         - linear_algebra::getCrossProductMatrix( downPerturbedAngularVelocity ) * rotationMatrixToBodyFixed;

//                 Eigen::Vector6d downPerturbedBodyFixedState = computeBodyFixedFromInertialState( 
//                         inertialState, rotationMatrixToBodyFixed, downPerturbedRotationMatrixToBodyFixedDerivative );

//                 // Compute numerical partials        
//                 partialsBodyFixedWrtAngularVelocity.block( 0, i, 6, 1 ) = 
//                         ( upPerturbedBodyFixedState - downPerturbedBodyFixedState ) / ( 2.0 * angularVelocityPerturbation[i] );
//         }

//         Eigen::MatrixXd analyticalPartialBodyFixedWrtAngularVelocity = 
//                 analyticalPartialBodyFixedWrtRotationalState.block( 0, 4, 6, 3 );

//         std::cout << "partialsBodyFixedWrtAngularVelocity" << std::endl;
//         std::cout << partialsBodyFixedWrtAngularVelocity << std::endl;
//         std::cout << "analyticalPartialBodyFixedWrtAngularVelocity" << std::endl;
//         std::cout << analyticalPartialBodyFixedWrtAngularVelocity << std::endl; 

//         for ( unsigned int i = 0 ; i < 6 ; i++ )
//         {
//                 for ( unsigned int j = 0 ; j < 3 ; j++ )
//                 {
//                         BOOST_CHECK_SMALL( partialsBodyFixedWrtAngularVelocity(i,j) - analyticalPartialBodyFixedWrtAngularVelocity(i,j), 1.0E-6 );
//                 }
//         }


//         // Check partials of body-fixed state wrt constant angular rate (for constant rotation model only)
//         if ( testCase == 0 )
//         {
//                 std::shared_ptr< SimpleRotationalEphemeris > constantRotationModel = 
//                         std::dynamic_pointer_cast< SimpleRotationalEphemeris >( io->getRotationalEphemeris() );
//                 double angularRatePerturbation = 1.0e-7;

//                 // Up-perturbation
//                 constantRotationModel->resetRotationRate( rotationRateIo + angularRatePerturbation );
//                 Eigen::Matrix3d upPerturbedRotationMatrixToBodyFixed = 
//                         constantRotationModel->getRotationToBaseFrame( testTime ).toRotationMatrix( ).transpose();
//                 Eigen::Matrix3d upPerturbedRotationMatrixToBodyFixedDerivative = 
//                         constantRotationModel->getDerivativeOfRotationToTargetFrame( testTime );

//                 // std::cout << "upPerturbedRotationMatrixToBodyFixedDerivative" << std::endl;
//                 // std::cout << upPerturbedRotationMatrixToBodyFixedDerivative << std::endl;

//                 Eigen::Vector6d upPerturbedBodyFixedState = computeBodyFixedFromInertialState( 
//                         inertialState, upPerturbedRotationMatrixToBodyFixed, upPerturbedRotationMatrixToBodyFixedDerivative );

//                 // Down-perturbation        
//                 constantRotationModel->resetRotationRate( rotationRateIo - angularRatePerturbation );
//                 Eigen::Matrix3d downPerturbedRotationMatrixToBodyFixed = 
//                         constantRotationModel->getRotationToBaseFrame( testTime ).toRotationMatrix( ).transpose();
//                 Eigen::Matrix3d downPerturbedRotationMatrixToBodyFixedDerivative = 
//                         constantRotationModel->getDerivativeOfRotationToTargetFrame( testTime );

//                 // std::cout << "downPerturbedRotationMatrixToBodyFixedDerivative" << std::endl;
//                 // std::cout << downPerturbedRotationMatrixToBodyFixedDerivative << std::endl;

//                 Eigen::Vector6d downPerturbedBodyFixedState = computeBodyFixedFromInertialState( 
//                         inertialState, downPerturbedRotationMatrixToBodyFixed, downPerturbedRotationMatrixToBodyFixedDerivative );

//                 Eigen::MatrixXd partialsBodyFixedWrtAngularRate = 
//                         ( upPerturbedBodyFixedState - downPerturbedBodyFixedState ) / ( 2.0 * angularRatePerturbation );
                
//                 constantRotationModel->resetRotationRate( rotationRateIo );

//                 Eigen::MatrixXd analyticalPartialBodyFixedWrtAngularRate = deformationPartial->bodyFixedWrtRotationParameter( constant_rotation_rate, "" );        


//                 std::cout << "partialsBodyFixedWrtAngularRate" << std::endl;
//                 std::cout << partialsBodyFixedWrtAngularRate << std::endl;
//                 std::cout << "analyticalPartialBodyFixedWrtAngularRate" << std::endl;
//                 std::cout << analyticalPartialBodyFixedWrtAngularRate << std::endl; 

//                 for ( unsigned int i = 0 ; i < 6 ; i++ )
//                 {
//                         if ( i == 2 || i == 5 ) // check z and vz separately (small values)
//                         {
//                                 BOOST_CHECK_SMALL( partialsBodyFixedWrtAngularRate(i) - analyticalPartialBodyFixedWrtAngularRate(i), 1.0E-5 );
//                         }
//                         else
//                         {
//                                 BOOST_CHECK_CLOSE_FRACTION( partialsBodyFixedWrtAngularRate(i), analyticalPartialBodyFixedWrtAngularRate(i), 1.0E-8 );
//                         }
                        
//                 }
//         }


//         // Check partials of spherical wrt cartesian state coordinates 
//         Eigen::Matrix6d partialsSphericalWrtCartesianState = Eigen::Matrix6d::Zero();    

//         for ( unsigned int i = 0 ; i < 6 ; i++ )
//         {
//                 Eigen::Vector6d upPerturbedBodyFixedState = bodyFixedState;
//                 upPerturbedBodyFixedState[ i ] += statePerturbation[ i ];
//                 Eigen::Vector6d upPerturbedSphericalState = computeSphericalCoordinatesFromBodyFixedState( upPerturbedBodyFixedState );

//                 Eigen::Vector6d downPerturbedBodyFixedState = bodyFixedState;
//                 downPerturbedBodyFixedState[ i ] -= statePerturbation[ i ];
//                 Eigen::Vector6d downPerturbedSphericalState = computeSphericalCoordinatesFromBodyFixedState( downPerturbedBodyFixedState );

//                 partialsSphericalWrtCartesianState.block( 0, i, 6, 1 ) = ( upPerturbedSphericalState - downPerturbedSphericalState ) / ( 2.0 * statePerturbation[i] );
//         }

//         Eigen::Matrix6d analyticalPartialSphericalWrtCartesianState = deformationPartial->sphericalWrtCartesianBodyFixedState();
//         std::cout << "partialsSphericalWrtCartesianState" << std::endl;
//         std::cout << partialsSphericalWrtCartesianState << std::endl;
//         std::cout << "analyticalPartialSphericalWrtCartesianState" << std::endl;
//         std::cout << analyticalPartialSphericalWrtCartesianState << std::endl;

//         for ( unsigned int i = 0 ; i < 6 ; i++ )
//         {
//                 for ( unsigned int j = 0 ; j < 6 ; j++ )
//                 {
//                         BOOST_CHECK_SMALL( partialsSphericalWrtCartesianState(i,j) - analyticalPartialSphericalWrtCartesianState(i,j), 1.0E-9 );
//                 }
//         }


//         // Check partials of Maxwell equilibrium coefficients, derivatives, and full deformation wrt (body-fixed) spherical state coordinates 
//         Eigen::Vector6d sphericalState = computeSphericalCoordinatesFromBodyFixedState( bodyFixedState );

//         double gravitationalParametersRatio = jupiter->getGravitationalParameter() / io->getGravitationalParameter();
//         double invGlobalRelaxationTime = 1.0 / globalRelaxationTime;

//         Eigen::MatrixXd partialsEquilibriumCoefficientsWrtSphericalState = Eigen::MatrixXd::Zero( 5, 6 );   
//         Eigen::MatrixXd partialsEquilibriumCoefficientsDerivativeWrtSphericalState = Eigen::MatrixXd::Zero( 5, 6 ); 
//         Eigen::MatrixXd partialsDeformationWrtSphericalState = Eigen::MatrixXd::Zero( 5, 6 );  

//         for ( unsigned int i = 0 ; i < 6 ; i++ )
//         {
//                 double perturbation = sphericalState[i] * 0.0001;
//                 if ( perturbation == 0.0 )
//                 {
//                         perturbation = 1.0e-10;
//                 }
//                 if ( testCase == 1 && i == 2 ) // specific latitude perturbation if rotation is synchronous (otherwise too small to be detectable in the numerical partials)
//                 {
//                         perturbation = 1.0e-12;
//                 }        

//                 Eigen::Vector6d upPerturbedSphericalState = sphericalState;
//                 upPerturbedSphericalState[ i ] += perturbation;
                
//                 Eigen::Vector5d upPerturbedEquilibriumCoefficients = computeEquilibriumCoefficientsFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, upPerturbedSphericalState );
//                 Eigen::Vector5d upPerturbedEquilibriumCoefficientsDerivative = computeEquilibriumCoefficientsDerivativeFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, upPerturbedSphericalState );
//                 Eigen::Vector5d upPerturbedDeformation = computeDeformationFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, globalRelaxationTime, maxwellRelaxationTime, upPerturbedSphericalState );

//                 Eigen::Vector6d downPerturbedSphericalState = sphericalState;
//                 downPerturbedSphericalState[ i ] -= perturbation;
                
//                 Eigen::Vector5d downPerturbedEquilibriumCoefficients = computeEquilibriumCoefficientsFromSphericalCoordinates(
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, downPerturbedSphericalState );
//                 Eigen::Vector5d downPerturbedEquilibriumCoefficientsDerivative = computeEquilibriumCoefficientsDerivativeFromSphericalCoordinates(
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, downPerturbedSphericalState );
//                 Eigen::Vector5d downPerturbedDeformation = computeDeformationFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, globalRelaxationTime, maxwellRelaxationTime, downPerturbedSphericalState );

//                 partialsEquilibriumCoefficientsWrtSphericalState.block( 0, i, 5, 1 ) = 
//                         ( upPerturbedEquilibriumCoefficients - downPerturbedEquilibriumCoefficients ) / ( 2.0 * perturbation );
//                 partialsEquilibriumCoefficientsDerivativeWrtSphericalState.block( 0, i, 5, 1 ) = 
//                         ( upPerturbedEquilibriumCoefficientsDerivative - downPerturbedEquilibriumCoefficientsDerivative ) / ( 2.0 * perturbation );
//                 partialsDeformationWrtSphericalState.block( 0, i, 5, 1 ) = 
//                         ( upPerturbedDeformation - downPerturbedDeformation ) / ( 2.0 * perturbation );
//         }

//         // Compute analytical partials
//         Eigen::MatrixXd analyticalPartialEquilibriumCoefficientsWrtSphericalState = deformationPartial->equilibriumCoefficientsWrtSphericalBodyFixedState();
//         Eigen::MatrixXd analyticalPartialEquilibriumCoefficientsDerivativeWrtSphericalState = deformationPartial->equilibriumCoefficientsDerivativeWrtSphericalBodyFixedState();
//         Eigen::MatrixXd analyticalPartialDeformationWrtSphericalState = deformationPartial->deformationWrtSphericalBodyFixedState();

//         std::cout << "--------------------------" << std::endl;
//         std::cout << "partials eq. coefficients wrt spherical body-fixed state" << std::endl;
//         std::cout << "partialsEquilibriumCoefficientsWrtSphericalState" << std::endl;
//         std::cout << partialsEquilibriumCoefficientsWrtSphericalState << std::endl;
//         std::cout << "analyticalPartialEquilibriumCoefficientsWrtSphericalState" << std::endl;
//         std::cout << analyticalPartialEquilibriumCoefficientsWrtSphericalState << std::endl;

//         std::cout << "--------------------------" << std::endl;
//         std::cout << "partials eq. coefficients **derivative** wrt spherical body-fixed state" << std::endl;
//         std::cout << "partialsEquilibriumCoefficientsDerivativeWrtSphericalState" << std::endl;
//         std::cout << partialsEquilibriumCoefficientsDerivativeWrtSphericalState << std::endl;
//         std::cout << "analyticalPartialEquilibriumCoefficientsDerivativeWrtSphericalState" << std::endl;
//         std::cout << analyticalPartialEquilibriumCoefficientsDerivativeWrtSphericalState << std::endl;

//         std::cout << "--------------------------" << std::endl;
//         std::cout << "partials Maxwell deformation wrt spherical body-fixed state" << std::endl;
//         std::cout << "partialsDeformationWrtSphericalState" << std::endl;
//         std::cout << partialsDeformationWrtSphericalState << std::endl;
//         std::cout << "analyticalPartialDeformationWrtSphericalState" << std::endl;
//         std::cout << analyticalPartialDeformationWrtSphericalState << std::endl;

//         for ( unsigned int i = 0 ; i < 5 ; i++ )
//         {
//                 for ( unsigned int j = 0 ; j < 6 ; j++ )
//                 {
//                         BOOST_CHECK_SMALL( partialsEquilibriumCoefficientsWrtSphericalState(i,j) - analyticalPartialEquilibriumCoefficientsWrtSphericalState(i,j), 1.0E-10 );
//                         BOOST_CHECK_SMALL( partialsEquilibriumCoefficientsDerivativeWrtSphericalState(i,j) - analyticalPartialEquilibriumCoefficientsDerivativeWrtSphericalState(i,j), 1.0E-10 );
//                         BOOST_CHECK_SMALL( partialsDeformationWrtSphericalState(i,j) - analyticalPartialDeformationWrtSphericalState(i,j), 1.0E-10 );
//                 }
//         }


//         // Check partials of Maxwell equilibrium coefficients, derivatives, and full deformation wrt (body-fixed) angular velocity (due to centrifugal potential)
//         Eigen::MatrixXd partialsEquilibriumCoefficientsWrtAngularVelocity = Eigen::MatrixXd::Zero( 5, 3 );   
//         Eigen::MatrixXd partialsEquilibriumCoefficientsDerivativeWrtAngularVelocity = Eigen::MatrixXd::Zero( 5, 3 );
//         Eigen::MatrixXd partialsDeformationWrtAngularVelocity = Eigen::MatrixXd::Zero( 5, 3 );

//         for ( unsigned int i = 0 ; i < 3 ; i++ )
//         {
//                 // Up-perturbation
//                 Eigen::Vector3d upPerturbedAngularVelocity = angularVelocityBodyFixed;
//                 upPerturbedAngularVelocity[ i ] += angularVelocityPerturbation[ i ];

//                 Eigen::Vector5d upPerturbedEquilibriumCoefficients = computeEquilibriumCoefficientsFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, sphericalState, upPerturbedAngularVelocity );
//                 Eigen::Vector5d upPerturbedEquilibriumCoefficientsDerivative = computeEquilibriumCoefficientsDerivativeFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, sphericalState, upPerturbedAngularVelocity );
//                 Eigen::Vector5d upPerturbedDeformation = computeDeformationFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, globalRelaxationTime, maxwellRelaxationTime, sphericalState, upPerturbedAngularVelocity );       

//                 // Down-perturbation        
//                 Eigen::Vector3d downPerturbedAngularVelocity = angularVelocityBodyFixed;
//                 downPerturbedAngularVelocity[ i ] -= angularVelocityPerturbation[ i ];
                
//                 Eigen::Vector5d downPerturbedEquilibriumCoefficients = computeEquilibriumCoefficientsFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, sphericalState, downPerturbedAngularVelocity );
//                 Eigen::Vector5d downPerturbedEquilibriumCoefficientsDerivative = computeEquilibriumCoefficientsDerivativeFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, sphericalState, downPerturbedAngularVelocity );
//                 Eigen::Vector5d downPerturbedDeformation = computeDeformationFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, globalRelaxationTime, maxwellRelaxationTime, sphericalState, downPerturbedAngularVelocity ); 

//                 // Compute numerical partials        
//                 partialsEquilibriumCoefficientsWrtAngularVelocity.block( 0, i, 5, 1 ) = 
//                         ( upPerturbedEquilibriumCoefficients - downPerturbedEquilibriumCoefficients ) / ( 2.0 * angularVelocityPerturbation[i] );
//                 partialsEquilibriumCoefficientsDerivativeWrtAngularVelocity.block( 0, i, 5, 1 ) = 
//                         ( upPerturbedEquilibriumCoefficientsDerivative - downPerturbedEquilibriumCoefficientsDerivative ) / ( 2.0 * angularVelocityPerturbation[i] );
//                 partialsDeformationWrtAngularVelocity.block( 0, i, 5, 1 ) = 
//                         ( upPerturbedDeformation - downPerturbedDeformation ) / ( 2.0 * angularVelocityPerturbation[i] );
//         }

//         Eigen::MatrixXd analyticalPartialEquilibriumCoefficientsWrtAngularVelocity = deformationPartial->equilibriumCoefficientsWrtAngularVelocityVector();
//         Eigen::MatrixXd analyticalPartialEquilibriumCoefficientsDerivativeWrtAngularVelocity = deformationPartial->equilibriumCoefficientsDerivativeWrtAngularVelocityVector();
//         Eigen::MatrixXd analyticalPartialDeformationWrtRotationalState = deformationPartial->deformationWrtRotationalState();
//         Eigen::MatrixXd analyticalPartialDeformationWrtAngularVelocity = analyticalPartialDeformationWrtRotationalState.block( 0, 4, 5, 3 );


//         std::cout << "--------------------------" << std::endl;
//         std::cout << "partials eq. coefficients wrt angular velocity vector" << std::endl;
//         std::cout << "partialsEquilibriumCoefficientsWrtAngularVelocity" << std::endl;
//         std::cout << partialsEquilibriumCoefficientsWrtAngularVelocity << std::endl;
//         std::cout << "analyticalPartialEquilibriumCoefficientsWrtAngularVelocity" << std::endl;
//         std::cout << analyticalPartialEquilibriumCoefficientsWrtAngularVelocity << std::endl;

//         std::cout << "--------------------------" << std::endl;
//         std::cout << "partials eq. coefficients **derivative** wrt angular velocity vector" << std::endl;
//         std::cout << "partialsEquilibriumCoefficientsDerivativeWrtAngularVelocity" << std::endl;
//         std::cout << partialsEquilibriumCoefficientsDerivativeWrtAngularVelocity << std::endl;
//         std::cout << "analyticalPartialEquilibriumCoefficientsDerivativeWrtAngularVelocity" << std::endl;
//         std::cout << analyticalPartialEquilibriumCoefficientsDerivativeWrtAngularVelocity << std::endl;

//         std::cout << "--------------------------" << std::endl;
//         std::cout << "partials Maxwell deformation wrt angular velocity vector" << std::endl;
//         std::cout << "partialsDeformationWrtAngularVelocity" << std::endl;
//         std::cout << partialsDeformationWrtAngularVelocity << std::endl;
//         std::cout << "analyticalPartialDeformationWrtAngularVelocity" << std::endl;
//         std::cout << analyticalPartialDeformationWrtAngularVelocity << std::endl;

//         for ( unsigned int i = 0 ; i < 5 ; i++ )
//         {
//                 for ( unsigned int j = 0 ; j < 3 ; j++ )
//                 {
//                         if ( partialsEquilibriumCoefficientsWrtAngularVelocity(i,j) > 1.0e-14 && analyticalPartialEquilibriumCoefficientsWrtAngularVelocity(i,j) > 1.0e-14 )
//                         {
//                                 BOOST_CHECK_SMALL( partialsEquilibriumCoefficientsWrtAngularVelocity(i,j) - analyticalPartialEquilibriumCoefficientsWrtAngularVelocity(i,j), 1.0E-10 );
//                         }
//                         if ( partialsEquilibriumCoefficientsDerivativeWrtAngularVelocity(i,j) > 1.0e-14 && analyticalPartialEquilibriumCoefficientsDerivativeWrtAngularVelocity(i,j) > 1.0e-14 )
//                         {
//                                 BOOST_CHECK_SMALL( partialsEquilibriumCoefficientsDerivativeWrtAngularVelocity(i,j) - analyticalPartialEquilibriumCoefficientsDerivativeWrtAngularVelocity(i,j), 1.0E-10 );
//                         }
//                         if ( partialsDeformationWrtAngularVelocity(i,j) > 1.0e-14 && analyticalPartialDeformationWrtAngularVelocity(i,j) > 1.0e-14 )
//                         {
//                                 BOOST_CHECK_SMALL( partialsDeformationWrtAngularVelocity(i,j) - analyticalPartialDeformationWrtAngularVelocity(i,j), 1.0E-10 );
//                         }
//                 }
//         }

//         TUDAT_CHECK_MATRIX_CLOSE_FRACTION( analyticalPartialDeformationWrtRotationalState.block( 0, 0, 5, 4 ), Eigen::MatrixXd::Zero( 5, 4 ), std::numeric_limits< double >::epsilon( ) );


//         // Check partials of Maxwell equilibrium coefficients and full deformation wrt constant angular rate (due to centrifugal potential) (for constant rotation model only)
//         if ( testCase == 0 )
//         {
//                 Eigen::MatrixXd partialsEquilibriumCoefficientsWrtRotationRate = Eigen::MatrixXd::Zero( 5, 1 );   
//                 Eigen::MatrixXd partialsDeformationWrtRotationRate = Eigen::MatrixXd::Zero( 5, 1 );

//                 std::shared_ptr< SimpleRotationalEphemeris > constantRotationModel = 
//                         std::dynamic_pointer_cast< SimpleRotationalEphemeris >( io->getRotationalEphemeris() );
//                 double angularRatePerturbation = 1.0e-7;

//                 // Up-perturbation
//                 // constantRotationModel->resetRotationRate( rotationRateIo + angularRatePerturbation );
//                 Eigen::Vector3d upPerturbedAngularVelocity = ( Eigen::Vector3d( ) << 0.0, 0.0, rotationRateIo + angularRatePerturbation ).finished( ); 

//                 Eigen::Vector5d upPerturbedEquilibriumCoefficients = computeEquilibriumCoefficientsFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, sphericalState, upPerturbedAngularVelocity );
//                 Eigen::Vector5d upPerturbedDeformation = computeDeformationFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, globalRelaxationTime, maxwellRelaxationTime, sphericalState, upPerturbedAngularVelocity ); 

//                 // Down-perturbation        
//                 Eigen::Vector3d downPerturbedAngularVelocity = ( Eigen::Vector3d( ) << 0.0, 0.0, rotationRateIo - angularRatePerturbation ).finished( );

//                 Eigen::Vector5d downPerturbedEquilibriumCoefficients = computeEquilibriumCoefficientsFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, sphericalState, downPerturbedAngularVelocity );
//                 Eigen::Vector5d downPerturbedDeformation = computeDeformationFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, globalRelaxationTime, maxwellRelaxationTime, sphericalState, downPerturbedAngularVelocity ); 

//                 // Compute numerical partials        
//                 partialsEquilibriumCoefficientsWrtRotationRate.block( 0, 0, 5, 1 ) = 
//                         ( upPerturbedEquilibriumCoefficients - downPerturbedEquilibriumCoefficients ) / ( 2.0 * angularRatePerturbation );
//                 partialsDeformationWrtRotationRate.block( 0, 0, 5, 1 ) = 
//                         ( upPerturbedDeformation - downPerturbedDeformation ) / ( 2.0 * angularRatePerturbation );        

                
//                 Eigen::MatrixXd analyticalPartialEquilibriumCoefficientsWrtRotationRate = deformationPartial->equilibriumCoefficientsWrtRotationRate();
//                 Eigen::MatrixXd analyticalPartialDeformationWrtRotationRate = deformationPartial->deformationWrtRotationRate();        


//                 std::cout << "--------------------------" << std::endl;
//                 std::cout << "partials eq. coefficients wrt rotation rate " << std::endl;
//                 std::cout << "partialsEquilibriumCoefficientsWrtRotationRate" << std::endl;
//                 std::cout << partialsEquilibriumCoefficientsWrtRotationRate << std::endl;
//                 std::cout << "analyticalPartialEquilibriumCoefficientsWrtRotationRate" << std::endl;
//                 std::cout << analyticalPartialEquilibriumCoefficientsWrtRotationRate << std::endl;

//                 std::cout << "--------------------------" << std::endl;
//                 std::cout << "partials Maxwell deformation wrt rotation rate" << std::endl;
//                 std::cout << "partialsDeformationWrtRotationRate" << std::endl;
//                 std::cout << partialsDeformationWrtRotationRate << std::endl;
//                 std::cout << "analyticalPartialDeformationWrtRotationRate" << std::endl;
//                 std::cout << analyticalPartialDeformationWrtRotationRate << std::endl;

//                 for ( unsigned int i = 0 ; i < 5 ; i++ )
//                 {
//                         BOOST_CHECK_SMALL( partialsEquilibriumCoefficientsWrtRotationRate(i,0) - analyticalPartialEquilibriumCoefficientsWrtRotationRate(i,0), 1.0E-10 );
//                         BOOST_CHECK_SMALL( partialsDeformationWrtRotationRate(i,0) - analyticalPartialDeformationWrtRotationRate(i,0), 1.0E-10 );
//                 }
//         }


//         // Check partials of Maxwell equilibrium coefficients derivatives and full deformation wrt
//         // (body-fixed) angular velocity derivative (due to centrifugal potential)
//         Eigen::MatrixXd partialsEquilibriumCoefficientsDerivativeWrtAngularVelocityDerivative = Eigen::MatrixXd::Zero( 5, 3 );
//         Eigen::MatrixXd partialsDeformationWrtAngularVelocityDerivative = Eigen::MatrixXd::Zero( 5, 3 );

//         for ( unsigned int i = 0 ; i < 3 ; i++ )
//         {
//                 // Up-perturbation
//                 Eigen::Vector3d upPerturbedAngularVelocityDerivative = angularVelocityDerivative;
//                 upPerturbedAngularVelocityDerivative[ i ] += angularVelocityPerturbation[ i ];

//                 Eigen::Vector5d upPerturbedEquilibriumCoefficientsDerivative = computeEquilibriumCoefficientsDerivativeFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, sphericalState, 
//                         angularVelocityBodyFixed, upPerturbedAngularVelocityDerivative );
//                 Eigen::Vector5d upPerturbedDeformation = computeDeformationFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, globalRelaxationTime, maxwellRelaxationTime, sphericalState, 
//                         angularVelocityBodyFixed, upPerturbedAngularVelocityDerivative );       

//                 // Down-perturbation        
//                 Eigen::Vector3d downPerturbedAngularVelocityDerivative = angularVelocityDerivative;
//                 downPerturbedAngularVelocityDerivative[ i ] -= angularVelocityPerturbation[ i ];
                
//                 Eigen::Vector5d downPerturbedEquilibriumCoefficientsDerivative = computeEquilibriumCoefficientsDerivativeFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, sphericalState,
//                         angularVelocityBodyFixed, downPerturbedAngularVelocityDerivative );
//                 Eigen::Vector5d downPerturbedDeformation = computeDeformationFromSphericalCoordinates( 
//                         fluidLoveNumber, gravitationalParametersRatio, io->getGravitationalParameter(), radiusIo, globalRelaxationTime, maxwellRelaxationTime, sphericalState,
//                         angularVelocityBodyFixed, downPerturbedAngularVelocityDerivative ); 

//                 // Compute numerical partials        
//                 partialsEquilibriumCoefficientsDerivativeWrtAngularVelocityDerivative.block( 0, i, 5, 1 ) = 
//                         ( upPerturbedEquilibriumCoefficientsDerivative - downPerturbedEquilibriumCoefficientsDerivative ) / ( 2.0 * angularVelocityPerturbation[i] );
//                 partialsDeformationWrtAngularVelocityDerivative.block( 0, i, 5, 1 ) = 
//                         ( upPerturbedDeformation - downPerturbedDeformation ) / ( 2.0 * angularVelocityPerturbation[i] );
//         }

//         Eigen::MatrixXd analyticalPartialEquilibriumCoefficientsDerivativeWrtAngularVelocityDerivative = 
//                 deformationPartial->equilibriumCoefficientsDerivativeWrtAngularVelocityVectorDerivative();
//         Eigen::MatrixXd analyticalPartialDeformationWrtAngularVelocityDerivative = deformationPartial->deformationWrtAngularVelocityVectorDerivative();

//         std::cout << "--------------------------" << std::endl;
//         std::cout << "partials eq. coefficients **derivative** wrt angular velocity vector derivative" << std::endl;
//         std::cout << "partialsEquilibriumCoefficientsDerivativeWrtAngularVelocityDerivative" << std::endl;
//         std::cout << partialsEquilibriumCoefficientsDerivativeWrtAngularVelocityDerivative << std::endl;
//         std::cout << "analyticalPartialEquilibriumCoefficientsDerivativeWrtAngularVelocityDerivative" << std::endl;
//         std::cout << analyticalPartialEquilibriumCoefficientsDerivativeWrtAngularVelocityDerivative << std::endl;

//         std::cout << "--------------------------" << std::endl;
//         std::cout << "partials Maxwell deformation wrt angular velocity vector derivative" << std::endl;
//         std::cout << "partialsDeformationWrtAngularVelocityDerivative" << std::endl;
//         std::cout << partialsDeformationWrtAngularVelocityDerivative << std::endl;
//         std::cout << "analyticalPartialDeformationWrtAngularVelocityDerivative" << std::endl;
//         std::cout << analyticalPartialDeformationWrtAngularVelocityDerivative << std::endl;

//         for ( unsigned int i = 0 ; i < 5 ; i++ )
//         {
//                 for ( unsigned int j = 0 ; j < 3 ; j++ )
//                 {
//                         BOOST_CHECK_SMALL( partialsEquilibriumCoefficientsDerivativeWrtAngularVelocityDerivative(i,j) - analyticalPartialEquilibriumCoefficientsDerivativeWrtAngularVelocityDerivative(i,j), 1.0E-10 );
//                         BOOST_CHECK_SMALL( partialsDeformationWrtAngularVelocityDerivative(i,j) - analyticalPartialDeformationWrtAngularVelocityDerivative(i,j), 1.0E-10 );
//                 }
//         }


//     }

// }


// BOOST_AUTO_TEST_CASE( testMaxwellDeformationPartials )
// {
//     std::cout << "START testMaxwellDeformationPartials" << std::endl;
//     std::cout.precision(20);

//     double initialTime = 0.0;
        
//     // Load spice kernels.
//     spice_interface::loadStandardSpiceKernels( );

//     // Create bodies
//     SystemOfBodies bodies = SystemOfBodies( "Jupiter", "J2000" );

//     // Create Jupiter
//     bodies.createEmptyBody( "Jupiter", false );

//     // Set Jupiter's ephemeris to
//     bodies.at( "Jupiter" )->setEphemeris(
//             std::make_shared< ephemerides::ConstantEphemeris >( [ = ]( ) { return Eigen::Vector6d::Zero( ); }, "SSB", "J2000" ) );

//     double muJupiter = spice_interface::getBodyGravitationalParameter( "Jupiter" );
//     bodies.at( "Jupiter" )->setGravityFieldModel(
//             std::make_shared< gravitation::GravityFieldModel >( muJupiter ) );

//     // Set Jupiter's rotation model
//     double rightAscensionPole = ( 358.054324066462 - 90.0 ) * mathematical_constants::PI / 180.0;
//     double declinationPole = ( 90.0 - 25.5034135739821 ) * mathematical_constants::PI / 180.0;
//     double primeMeridian = ( 284.95 ) * mathematical_constants::PI / 180.0;
//     double rotationRateJupiter = ( 870.536 * mathematical_constants::PI / 180.0 ) / 86400.0;   
    
//     bodies.at("Jupiter")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//                 rightAscensionPole, declinationPole, primeMeridian, rotationRateJupiter, initialTime, "J2000", "IAU_Jupiter" ) );
    
//     // Create Io
//     bodies.createEmptyBody( "Io" );
//     std::shared_ptr< Body > io = bodies.at( "Io" );
//     std::shared_ptr< Body > jupiter = bodies.at( "Jupiter" );

//     // Set Io's gravity field
//     double muIo = spice_interface::getBodyGravitationalParameter( "Io" ); 
//     double radiusIo = 1821.6E3;
//     double scaledMeanMomentOfInertia = 0.37685;
//     Eigen::MatrixXd ioCosineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );
//     Eigen::MatrixXd ioSineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );
//     ioCosineCoefficients( 0, 0 ) = 1.0;
//     ioCosineCoefficients( 2, 0 ) = -1845.9E-6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
//     ioCosineCoefficients( 2, 2 ) = 553.7E-6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
          
//     bodies.at( "Io" )->setGravityFieldModel( std::make_shared< gravitation::SphericalHarmonicsGravityField >( 
//         muIo, radiusIo, ioCosineCoefficients, ioSineCoefficients, "IAU_Io", scaledMeanMomentOfInertia ) );

//     Eigen::Vector5d ioGravityDeformation = Eigen::Vector5d::Zero();
//     ioGravityDeformation.segment(0, 3) = ioCosineCoefficients.block(2,0,1,3);
//     ioGravityDeformation.segment(3, 2) = ioSineCoefficients.block(2,1,1,2);
//     std::cout << "ioGravityDeformation " << ioGravityDeformation.transpose() << std::endl;
//     io->setCurrentPropagatedGravityField( ioGravityDeformation );

//     // Set Io's ephemeris
//     double muEffective = muIo + muJupiter;
//     Eigen::Vector6d ioKeplerElements = Eigen::Vector6d::Zero( );
//     double ioSemiMajorAxis = 4.2e8;
//     ioKeplerElements( 0 ) = ioSemiMajorAxis;
// //     ioKeplerElements( 1 ) = 0.004;
//     bodies.at( "Io" )->setEphemeris( std::make_shared< ephemerides::KeplerEphemeris >( ioKeplerElements, 0.0, muEffective, "Jupiter", "J2000" ) );
    
//      double rotationRateIo = std::sqrt( muEffective / ( 4.2e8 * 4.2e8 * 4.2e8 ) );


//      // Define state parameters perturbations 
//      Eigen::Vector5d gravityDeformationPerturbation;
//      gravityDeformationPerturbation << 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6;
    
//      Eigen::Vector4d orientationPerturbation;
//      orientationPerturbation << 1.0E-9, 1.0E-9, 1.0E-9, 1.0E-9;
     
//      Eigen::Vector7d rotationalStatePerturbation;
//      rotationalStatePerturbation << 1.0E-9, 1.0E-9, 1.0E-9, 1.0E-9, 1.0E-6, 1.0E-6, 1.0E-6;
     
//      Eigen::Vector6d translationalStatePerturbation;
//      translationalStatePerturbation << 1.0, 1.0, 100.0, 1.0E-3, 1.0E-3, 1.0E-3;


//      // Two test cases:
//     // test case 0: state-independent rotation model (i.e., simple rotation model)
//     // test case 1: state-dependent rotation model (i.e., synchronous rotation model)

//     for ( unsigned int testCase = 0 ; testCase < 2 ; testCase++ )
//     {
//         if ( testCase == 0 ) // state-independent rotation model
//         {
//                 // Set Io's rotation model to constant
//                 bodies.at("Io")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
//                         rightAscensionPole, declinationPole, primeMeridian, rotationRateIo, initialTime, "J2000", "IAU_Io" ) );
//         }
//         else if ( testCase == 1 ) // state dependent rotation model
//         {
//                 // Set Io's rotation model to synchronous
//                 std::function< Eigen::Vector6d( const double, bool ) > relativeStateFunction = [ io ]( const double, bool ) { return io->getState(); };
//                 bodies.at( "Io" )->setRotationalEphemeris( std::make_shared< SynchronousRotationalEphemeris >(
//                         relativeStateFunction, /*createRelativeStateFunction( bodies, "Io", "Jupiter" ),*/ "Jupiter", "J2000", "IAU_Io" ) );
//         }
   

//         // Update Jupiter and Io to current state
//         double testTime = 1000.0;
//         io->setStateFromEphemeris( testTime );
//         jupiter->setStateFromEphemeris( testTime );
//         // io->setCurrentRotationToLocalFrameFromEphemeris( testTime );
//         io->setCurrentRotationalStateToLocalFrameFromEphemeris( testTime );


//         // Create Maxwell deformation model 
//         double maxwellRelaxationTime = 179103.0;
//         double globalRelaxationTime = 24688.0;
//         double fluidLoveNumber = 1.5;
//         std::vector< std::string > perturbingBody = { "Jupiter" };
//         std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::make_shared< MaxwellDeformationSettings >( 
//                 maxwellRelaxationTime, globalRelaxationTime, fluidLoveNumber, 2, 2, perturbingBody, Eigen::VectorXd::Zero( 5 ), 
//                 true, true );

//         std::vector< std::shared_ptr< simulation_setup::Body > > perturbingBodies = { jupiter };
//         std::shared_ptr< basic_astrodynamics::MaxwellGravityDeformationModel > maxwellDeformationModel = createMaxwellGravityFieldDeformationModel(
//                 bodies.at("Io"), perturbingBodies, "Io", std::vector< std::string >( {"Jupiter"} ), maxwellDeformationSettings );


//         // Create parameter objects (necessary to the automatic computation of the rotation matrix partials)
//         Eigen::Vector6d initialTranslationalState = io->getState( );
//         Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRotationState = propagators::getInitialRotationalStateOfBody( "Io", "J2000",  bodies, initialTime );
//         initialRotationState[6] = rotationRateIo;

//         std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames; 
//         parameterNames.push_back( std::make_shared< estimatable_parameters::InitialTranslationalStateEstimatableParameterSettings< double > >( 
//                 "Io", initialTranslationalState, "Jupiter" ) ); 
//         parameterNames.push_back( std::make_shared< estimatable_parameters::InitialRotationalStateEstimatableParameterSettings< double > >( "Io", initialRotationState, "J2000" ) ); 
//         parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Io", gravitational_parameter ) );
//         parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Jupiter", gravitational_parameter ) );
//         if ( testCase == 0 ) // for constant rotation model only
//         {
//                 parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Io", constant_rotation_rate ) );
//         }

//         parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
//                 1, 0, 2, 2, "Io", spherical_harmonics_cosine_coefficient_block ) );
//         parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
//                 1, 1, 2, 2, "Io", spherical_harmonics_sine_coefficient_block ) );

//         std::shared_ptr< EstimatableParameterSet< double > > parameterSet = createParametersToEstimate( parameterNames, bodies );
//         printEstimatableParameterEntries( parameterSet );

//         // Create deformation partial.
//         std::shared_ptr< DeformationPartial > deformationPartial =
//                 createAnalyticalGravityDeformationPartial( maxwellDeformationModel, std::make_pair( "Io", io ), bodies, parameterSet );

//         std::shared_ptr< EstimatableParameter< double > > ioGravitationalParameterParameter =
//                 parameterSet->getEstimatedDoubleParameters( ).at( 0 );
//         std::shared_ptr< EstimatableParameter< double > > jupiterGravitationalParameterParameter = 
//                 parameterSet->getEstimatedDoubleParameters( ).at( 1 );
//         std::shared_ptr< EstimatableParameter< double > > ioConstantRotationRateParameter;
//         if ( testCase == 0 )
//         {
//                 ioConstantRotationRateParameter = parameterSet->getEstimatedDoubleParameters( ).at( 2 );
//         } 
//         std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > ioCosineCoefficientsParameter =
//                 parameterSet->getEstimatedVectorParameters( ).at( 0 );
//         std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > ioSineCoefficientsParameter =
//                 parameterSet->getEstimatedVectorParameters( ).at( 1 );

//         // Calculate analytical partials.
//         deformationPartial->update( testTime );

//         // Compute analytical partials wrt Io's own gravity state
//         Eigen::MatrixXd partialWrtGravityDeformation = Eigen::MatrixXd::Zero( 5, 5 );
//         deformationPartial->wrtStateOfDeformingBody( partialWrtGravityDeformation.block( 0, 0, 5, 5 ) );

//         Eigen::Matrix< double, 5, 5 > testPartialWrtGravityDeformation = Eigen::Matrix< double, 5, 5 >::Zero( );
//         std::function< void( Eigen::VectorXd ) > ioGravityDeformationSetFunction = 
//                 std::bind( &Body::setCurrentPropagatedGravityField, io, std::placeholders::_1 ); 

//         // Compute numerical partials wrt gravity state   
//         testPartialWrtGravityDeformation = calculateDeformationWrtGravityDeformationStatePartials( 
//                 ioGravityDeformationSetFunction,
//                 maxwellDeformationModel,
//                 ioGravityDeformation,
//                 gravityDeformationPerturbation, 0 );

//         // Check partials wrt Io's own gravity deformation   
//         std::cout << "------ partials wrt gravity state ------" << std::endl;
//         std::cout << "testPartialWrtGravityDeformation " << std::endl;
//         std::cout << testPartialWrtGravityDeformation << std::endl;
//         std::cout << "partialWrtGravityDeformation " << std::endl;
//         std::cout << partialWrtGravityDeformation << std::endl;
//         TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtGravityDeformation, partialWrtGravityDeformation, 1.0E-10 );


//         if ( testCase == 0 )
//         {
//                 // Compute analytical partials wrt Io's rotational state
//                 Eigen::MatrixXd partialWrtIoRotationalState = Eigen::MatrixXd::Zero( 5, 7 );
//                 deformationPartial->wrtNonDeformationStateOfAdditionalBody( partialWrtIoRotationalState.block(0, 0, 5, 7), std::make_pair( "Io", "" ), propagators::rotational_state );

//                 Eigen::Matrix< double, 5, 7 > testPartialWrtRotationalState = Eigen::Matrix< double, 5, 7 >::Zero( );
//                 Eigen::Matrix< double, 5, 3 > testPartialWrtAngularVelocity = Eigen::Matrix< double, 5, 3 >::Zero( );
//                 std::function< void( Eigen::Vector7d ) > ioRotationalStateSetFunction =
//                         std::bind( &Body::setCurrentRotationalStateToLocalFrame, io, std::placeholders::_1 );

//                 // Compute numerical partials wrt Io's rotational state     
//                 testPartialWrtRotationalState = calculateDeformationWrtRotationalStatePartials( ioRotationalStateSetFunction,
//                                                                                                         maxwellDeformationModel,
//                                                                                                         io->getRotationalStateVector( ),
//                                                                                                         rotationalStatePerturbation,
//                                                                                                         0, 7 );
//                 testPartialWrtAngularVelocity = testPartialWrtRotationalState.block( 0, 4, 5, 3 );

//                 // Compute numerical partials wrt quaternion
//                 std::vector< Eigen::Vector4d > appliedQuaternionPerturbation;
//                 Eigen::MatrixXd deformationDeviations = calculateGravityDeformationDeviationDueToOrientationChange( 
//                         ioRotationalStateSetFunction,
//                         maxwellDeformationModel,
//                         io->getRotationalStateVector( ),
//                         orientationPerturbation,
//                         appliedQuaternionPerturbation );

//                 // Check partials wrt Io's orientation
//                 std::cout << "------ partials wrt quaternion orientation ------" << std::endl;
//                 for( int index = 1; index < 4; index++ )
//                 {
//                         Eigen::Vector5d numericalChangeInDeformation = deformationDeviations.block( 0, index - 1, 5, 1 );
//                         Eigen::Vector5d analyticalChangeInDeformation =
//                                 partialWrtIoRotationalState.block( 0, 0, 5, 1 ) * appliedQuaternionPerturbation[ index ]( 0 ) +
//                                 partialWrtIoRotationalState.block( 0, index, 5, 1 ) * appliedQuaternionPerturbation[ index ]( index );

//                         std::cout << "numericalChangeInDeformation" << std::endl;
//                         std::cout << numericalChangeInDeformation.transpose() << std::endl;
//                         std::cout << "analyticalChangeInDeformation" << std::endl;
//                         std::cout << analyticalChangeInDeformation.transpose() << std::endl;

//                         TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalChangeInDeformation, analyticalChangeInDeformation, 5.0e-5 );
//                 }

//                 // Check partials wrt Io's angular velocity vector
//                 std::cout << "------ partials wrt angular velocity vector ------" << std::endl;
//                 std::cout << "testPartialWrtAngularVelocityVector " << std::endl;
//                 std::cout << testPartialWrtAngularVelocity << std::endl;
//                 std::cout << "partialWrtIoAngularVelocityVector " << std::endl;
//                 std::cout << partialWrtIoRotationalState.block( 0, 4, 5, 3 ) << std::endl;
//                 for ( unsigned int i = 0 ; i < 5 ; i++ )
//                 {
//                         for ( unsigned int j = 0 ; j < 3 ; j++ )
//                         {
//                                 // if ( i == 0 && j == 2 ) // test partial of C20 wrt wz separately (very small)
//                                 // {
//                                 //         BOOST_CHECK_SMALL( std::fabs( partialWrtIoRotationalState( i, 4+j ) ), 1.0E-16 );
//                                 // }
//                                 // else
//                                 // {
//                                         BOOST_CHECK_CLOSE_FRACTION( testPartialWrtAngularVelocity( i, j ), partialWrtIoRotationalState( i, 4+j ), 1.0E-10 );
//                                 // }
//                         }
//                 }
//         }
        


//         // Compute analytical partials wrt Io's translational state
//         Eigen::MatrixXd partialWrtIoTranslationalState = Eigen::MatrixXd::Zero( 5, 6 );
//         deformationPartial->wrtNonDeformationStateOfAdditionalBody( partialWrtIoTranslationalState.block(0, 0, 5, 6), std::make_pair( "Io", "" ), propagators::translational_state );

//         Eigen::Matrix< double, 5, 6 > testPartialWrtTranslationalState = Eigen::Matrix< double, 5, 6 >::Zero( );
//         std::function< void( Eigen::Vector6d ) > ioTranslationalStateSetFunction = std::bind( &Body::setState, io, std::placeholders::_1 );

//         // Update rotational state function
//         std::function< void( ) > updateRotationalStateFunction =
//                 std::bind( &Body::setCurrentRotationalStateToLocalFrameFromEphemeris< double >, bodies.at( "Io" ), testTime );

//         // Compute numerical partials wrt Io's translational state
//         testPartialWrtTranslationalState = calculateDeformationWrtTranslationalStatePartials( ioTranslationalStateSetFunction,
//                                                                                                 maxwellDeformationModel,
//                                                                                                 io->getState( ),
//                                                                                                 translationalStatePerturbation,
//                                                                                                 0, 6, updateRotationalStateFunction );
        
//         // Check partials wrt Io's translational state
//         std::cout << "------ partials wrt Io's translational state ------" << std::endl;
//         std::cout << "testPartialWrtTranslationalState " << std::endl;
//         std::cout << testPartialWrtTranslationalState << std::endl;
//         std::cout << "partialWrtIoTranslationalState " << std::endl;
//         std::cout << partialWrtIoTranslationalState << std::endl;
//         for ( unsigned int i = 0 ; i < 5 ; i++ )
//         {
//                 for ( unsigned int j = 0 ; j < 6 ; j++ )
//                 {
//                         if ( testPartialWrtTranslationalState( i, j ) > 1.0e-30 && partialWrtIoTranslationalState( i, j ) > 1.0e-30 )
//                         {
//                                 BOOST_CHECK_CLOSE_FRACTION( testPartialWrtTranslationalState( i, j ), partialWrtIoTranslationalState( i, j ), 1.0E-6 );
//                         }
//                 }
//         }                                                                                        


//         // Compute analytical partials wrt Jupiter's translational state
//         Eigen::MatrixXd partialWrtJupiterTranslationalState = Eigen::MatrixXd::Zero( 5, 6 );
//         deformationPartial->wrtNonDeformationStateOfAdditionalBody( partialWrtJupiterTranslationalState.block(0, 0, 5, 6), std::make_pair( "Jupiter", "" ), propagators::translational_state );      

//         Eigen::Matrix< double, 5, 6 > testPartialWrtJupiterTranslationalState = Eigen::Matrix< double, 5, 6 >::Zero( );
//         std::function< void( Eigen::Vector6d ) > jupiterTranslationalStateSetFunction = std::bind( &Body::setState, jupiter, std::placeholders::_1 );

//         // Compute numerical partials wrt Jupiter's translational state
//         testPartialWrtJupiterTranslationalState = calculateDeformationWrtTranslationalStatePartials( jupiterTranslationalStateSetFunction,
//                                                                                                         maxwellDeformationModel,
//                                                                                                         jupiter->getState( ),
//                                                                                                         translationalStatePerturbation,
//                                                                                                         0, 6 );        
//         // Check partials wrt Jupiter's translational state
//         std::cout << "------ partials wrt Jupiter's translational state ------" << std::endl;
//         std::cout << "testPartialWrtJupiterTranslationalState " << std::endl;
//         std::cout << testPartialWrtJupiterTranslationalState << std::endl;
//         std::cout << "partialWrtJupiterTranslationalState " << std::endl;
//         std::cout << partialWrtJupiterTranslationalState << std::endl;
//         for ( unsigned int i = 0 ; i < 5 ; i++ )
//         {
//                 for ( unsigned int j = 0 ; j < 6 ; j++ )
//                 {
//                         BOOST_CHECK_CLOSE_FRACTION( testPartialWrtJupiterTranslationalState( i, j ), partialWrtJupiterTranslationalState( i, j ), 1.0E-5 );
//                 }
//         }        

//         // Compute analytical partials wrt gravitational parameters    
//         Eigen::Vector5d partialWrtIoGravitationalParameter = deformationPartial->wrtParameter( ioGravitationalParameterParameter );
//         Eigen::Vector5d partialWrtJupiterGravitationalParameter = deformationPartial->wrtParameter( jupiterGravitationalParameterParameter );

//         // Numerical partials wrt gravitational parameters
//         Eigen::Vector5d testPartialWrtIoGravitationalParameter =
//                 calculateDeformationWrtParameterPartials( ioGravitationalParameterParameter, maxwellDeformationModel, 1.0E8 );
//         Eigen::Vector5d testPartialWrtJupiterGravitationalParameter =
//                 calculateDeformationWrtParameterPartials( jupiterGravitationalParameterParameter, maxwellDeformationModel, 1.0E8 );
        
//         // Check partials wrt Io's gravitational parameter
//         std::cout << "------ partials wrt Io's gravitational parameter ------" << std::endl;
//         std::cout << "testPartialWrtIoGravitationalParameter " << std::endl;
//         std::cout << testPartialWrtIoGravitationalParameter.transpose() << std::endl;
//         std::cout << "partialWrtIoGravitationalParameter " << std::endl;
//         std::cout << partialWrtIoGravitationalParameter.transpose() << std::endl;

//         TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtIoGravitationalParameter, partialWrtIoGravitationalParameter, 1.0E-8 );

//         // Check partials wrt Jupiter's gravitational parameter
//         std::cout << "testPartialWrtJupiterGravitationalParameter " << std::endl;
//         std::cout << testPartialWrtJupiterGravitationalParameter.transpose() << std::endl;
//         std::cout << "partialWrtJupiterGravitationalParameter " << std::endl;
//         std::cout << partialWrtJupiterGravitationalParameter.transpose() << std::endl;

//         TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtJupiterGravitationalParameter, partialWrtJupiterGravitationalParameter, 1.0E-5 );


//         if ( testCase == 0 )
//         {
//                 // Compute analytical partials wrt Io's constant rotation rate (for constant rotation model only, 
//                 // to test partials wrt rotation parameter that is not the initial rotational state)
//                 Eigen::Vector5d partialWrtIoConstantRotationRate = deformationPartial->wrtParameter( ioConstantRotationRateParameter );
        

//                 // Numerical partials wrt Io's constant rotation rate (for constant rotation model only)
//                 Eigen::Vector5d testPartialWrtIoConstantRotationRate = Eigen::Vector5d::Zero();

//                 std::function< void( ) > updateFunctionRotationRate =
//                 std::bind( &Body::setCurrentRotationalStateToLocalFrameFromEphemeris< double >, bodies.at( "Io" ), testTime );
//                 testPartialWrtIoConstantRotationRate = calculateDeformationWrtParameterPartials( ioConstantRotationRateParameter, maxwellDeformationModel, 1.0E-8, updateFunctionRotationRate );

//                 // Check partials wrt Io's constant rotation rate  
//                 std::cout << "testPartialWrtIoConstantRotationRate" << std::endl;
//                 std::cout << testPartialWrtIoConstantRotationRate.transpose() << std::endl;
//                 std::cout << "partialWrtIoConstantRotationRate" << std::endl;
//                 std::cout << partialWrtIoConstantRotationRate.transpose() << std::endl;

//                 for ( unsigned int i = 0 ; i < 5 ; i++ )
//                 {
//                         if ( i == 0 ) // Check C20 deformation separately (very low value)
//                         {
//                                 BOOST_CHECK_SMALL( std::fabs( partialWrtIoTranslationalState( i ) ), 1.0E-14 );
//                         }
//                         else
//                         {
//                                 BOOST_CHECK_CLOSE_FRACTION( testPartialWrtIoConstantRotationRate( i ), partialWrtIoConstantRotationRate( i ), 1.0E-6 );  
//                         }
//                 }
//         }
        
//         // Compute analytical partials wrt Io's "static" gravity coefficients 
//         Eigen::MatrixXd partialWrtIoCosineCoefficients = deformationPartial->wrtParameter( ioCosineCoefficientsParameter );
//         Eigen::MatrixXd partialWrtIoSineCoefficients = deformationPartial->wrtParameter( ioSineCoefficientsParameter );

//         // Numerical partials wrt Io's cosine/sine gravity coefficients
//         std::function< void( ) > updateFunction =
//                 std::bind( &RigidBodyProperties::update, bodies.at( "Io" )->getMassProperties( ), testTime );
//         Eigen::MatrixXd testPartialWrtIoCosineCoefficients = calculateDeformationWrtParameterPartials(
//                 ioCosineCoefficientsParameter, maxwellDeformationModel,
//                 Eigen::VectorXd::Constant( ioCosineCoefficientsParameter->getParameterSize( ), 1.0e-6 ), updateFunction );
//         Eigen::MatrixXd testPartialWrtIoSineCoefficients = calculateDeformationWrtParameterPartials( 
//                 ioSineCoefficientsParameter, maxwellDeformationModel,
//                 Eigen::VectorXd::Constant( ioSineCoefficientsParameter->getParameterSize( ), 1.0e-6 ), updateFunction );

//         // Check partials wrt Io's cosine gravity coefficients
//         std::cout << "testPartialWrtIoCosineCoefficients" << std::endl;
//         std::cout << testPartialWrtIoCosineCoefficients << std::endl;
//         std::cout << "partialWrtIoCosineCoefficients" << std::endl;
//         std::cout << partialWrtIoCosineCoefficients << std::endl;

//         // Check partials wrt Io's sine gravity coefficients
//         std::cout << "testPartialWrtIoSineCoefficients" << std::endl;
//         std::cout << testPartialWrtIoSineCoefficients << std::endl;
//         std::cout << "partialWrtIoSineCoefficients" << std::endl;
//         std::cout << partialWrtIoSineCoefficients << std::endl;

//         TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtIoCosineCoefficients, partialWrtIoCosineCoefficients, 1.0E-9 );
//         TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtIoSineCoefficients, partialWrtIoSineCoefficients, 1.0E-9 );


//     }
    
// }

BOOST_AUTO_TEST_CASE( testStateDerivativeInterdependenciesPartials )
{
    std::cout.precision(20);

    double initialTime = 0.0;
        
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies
    SystemOfBodies bodies = SystemOfBodies( "Jupiter", "J2000" );

    // Create Jupiter
    bodies.createEmptyBody( "Jupiter", false );
    bodies.at( "Jupiter" )->setEphemeris(
            std::make_shared< ephemerides::ConstantEphemeris >( [ = ]( ) { return Eigen::Vector6d::Zero( ); }, "SSB", "J2000" ) );

    double muJupiter = spice_interface::getBodyGravitationalParameter( "Jupiter" );
    bodies.at( "Jupiter" )->setGravityFieldModel(
            std::make_shared< gravitation::GravityFieldModel >( muJupiter ) );

    // Set Jupiter's rotation as constant (no precession)
    double rightAscensionPole = ( 358.054324066462 - 90.0 ) * mathematical_constants::PI / 180.0;
    double declinationPole = ( 90.0 - 25.5034135739821 ) * mathematical_constants::PI / 180.0;
    double primeMeridian = ( 284.95 ) * mathematical_constants::PI / 180.0;
    double rotationRateJupiter = ( 870.536 * mathematical_constants::PI / 180.0 ) / 86400.0;   
    
    bodies.at("Jupiter")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
                rightAscensionPole, declinationPole, primeMeridian, rotationRateJupiter, initialTime, "J2000", "IAU_Jupiter" ) );
    
    // Create Io
    bodies.createEmptyBody( "Io" );
    std::shared_ptr< Body > io = bodies.at( "Io" );
    std::shared_ptr< Body > jupiter = bodies.at( "Jupiter" );

    // Set Io gravity field
    double muIo = spice_interface::getBodyGravitationalParameter( "Io" ); 
    double radiusIo = 1821.6E3;
    double scaledMeanMomentOfInertia = 0.37685;
    Eigen::MatrixXd ioCosineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );
    Eigen::MatrixXd ioSineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );
    ioCosineCoefficients( 0, 0 ) = 1.0;
    ioCosineCoefficients( 2, 0 ) = -1845.9E-6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
    ioCosineCoefficients( 2, 2 ) = 553.7E-6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
          
    bodies.at( "Io" )->setGravityFieldModel( std::make_shared< gravitation::SphericalHarmonicsGravityField >( 
        muIo, radiusIo, ioCosineCoefficients, ioSineCoefficients, "IAU_Io", scaledMeanMomentOfInertia ) );

//     // Retrieve Io's default inertia tensor
//     Eigen::Matrix3d ioInertiaTensor = bodies.at("Io")->getBodyInertiaTensor();
//     std::cout << "ioInertiaTensor" << std::endl;
//     std::cout << ioInertiaTensor << std::endl;

    // Set Io's ephemeris
    double muEffective = muIo + muJupiter;
    Eigen::Vector6d ioKeplerElements = Eigen::Vector6d::Zero( );
    double ioSemiMajorAxis = 4.2e8;
    ioKeplerElements( 0 ) = ioSemiMajorAxis;
    bodies.at( "Io" )->setEphemeris( std::make_shared< ephemerides::KeplerEphemeris >( ioKeplerElements, 0.0, muEffective, "Jupiter", "J2000" ) );
    
     double rotationRateIo = std::sqrt( muEffective / ( 4.2e8 * 4.2e8 * 4.2e8 ) );

    // Set Io's rotation model 
//     bodies.at( "Io" )->setRotationalEphemeris( std::make_shared< SynchronousRotationalEphemeris >(
        // createRelativeStateFunction( bodies, "Io", "Jupiter" ), "Jupiter", "J2000", "IAU_Io" ) );
    bodies.at("Io")->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >( 
            rightAscensionPole, declinationPole, primeMeridian, rotationRateIo, initialTime, "J2000", "IAU_Io" ) );

    // Update Jupiter and Io to current state
    double testTime = 1000.0;
    io->setStateFromEphemeris( testTime );
    jupiter->setStateFromEphemeris( testTime );
    io->setCurrentRotationalStateToLocalFrameFromEphemeris( testTime );

    // Define Io's initial gravity state        
    Eigen::Vector5d ioGravityDeformation = Eigen::Vector5d::Zero();
    ioGravityDeformation.segment(0, 3) = ioCosineCoefficients.block(2,0,1,3).transpose();
    ioGravityDeformation.segment(3, 2) = ioSineCoefficients.block(2,1,1,2).transpose();
    std::cout << "ioGravityDeformation " << ioGravityDeformation.transpose() << std::endl;
    io->setCurrentPropagatedGravityField( ioGravityDeformation );

    // Create Maxwell deformation model 
    double maxwellRelaxationTime = 179103.0;
    double globalRelaxationTime = 24688.0;
    double fluidLoveNumber = 1.5;
    std::vector< std::string > perturbingBody = { "Jupiter" };
    std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::make_shared< MaxwellDeformationSettings >( 
            maxwellRelaxationTime, globalRelaxationTime, fluidLoveNumber, 2, 2, perturbingBody );

    std::vector< std::shared_ptr< simulation_setup::Body > > perturbingBodies = { jupiter };
    std::shared_ptr< basic_astrodynamics::MaxwellGravityDeformationModel > maxwellDeformationModel = createMaxwellGravityFieldDeformationModel(
        bodies.at("Io"), perturbingBodies, "Io", std::vector< std::string >( {"Jupiter"} ), maxwellDeformationSettings );

    // Create inertial torque model for Io
    std::shared_ptr< InertialTorqueModel > inertialTorqueModel = createInertialTorqueModel( bodies.at( "Io" ), "Io" );
    inertialTorqueModel->updateMembers( 0.0 );  

    // Create parameter objects.
    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialRotationState = propagators::getInitialRotationalStateOfBody( "Io", "J2000",  bodies, initialTime );
    initialRotationState[6] = rotationRateIo;

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames; 
    parameterNames.push_back( std::make_shared< estimatable_parameters::InitialRotationalStateEstimatableParameterSettings< double > >( "Io", initialRotationState, "J2000" ) ); 
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Io", gravitational_parameter ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Io", mean_moment_of_inertia ) );

    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
            1, 0, 2, 2, "Io", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
            1, 1, 2, 2, "Io", spherical_harmonics_sine_coefficient_block ) );

    std::shared_ptr< EstimatableParameterSet< double > > parameterSet = createParametersToEstimate( parameterNames, bodies );

    // Create deformation partial.
    std::shared_ptr< DeformationPartial > deformationPartial =
            createAnalyticalGravityDeformationPartial( maxwellDeformationModel, std::make_pair( "Io", io ), bodies, parameterSet );

    // Create direct torque partial.
    std::shared_ptr< TorquePartial > torquePartial =
            createAnalyticalTorquePartial( inertialTorqueModel, std::make_pair( "Io", io ), std::make_pair( "Io", io ) );


    // Calculate analytical partials.
    Eigen::Vector5d nominalDeformation = basic_astrodynamics::updateAndGetDeformation( maxwellDeformationModel, testTime );
    std::cout << "nominalDeformation " << nominalDeformation.transpose() << std::endl;
    io->getMassProperties()->updateInertiaTensorDerivative( nominalDeformation );

    deformationPartial->update( testTime );
    torquePartial->update( testTime );

    // Retrieve Io's default inertia tensor
    Eigen::Matrix3d ioInertiaTensor = bodies.at("Io")->getBodyInertiaTensor();
    std::cout << "ioInertiaTensor" << std::endl;
    std::cout << ioInertiaTensor << std::endl;
    
    // Retrieve interdependency multiplying factor
    Eigen::MatrixXd multiplyingFactorTorque = ioInertiaTensor * std::dynamic_pointer_cast< InertialTorquePartial >( torquePartial )->wrtOtherStateDerivative();
    std::cout << "multiplyingFactorTorque" << std::endl;
    std::cout << multiplyingFactorTorque << std::endl; 
//     multiplyingFactorTorque = Eigen::MatrixXd::Zero( 3, 5 ); 

    Eigen::MatrixXd multiplyingFactorDeformation = std::dynamic_pointer_cast< MaxwellDeformationPartial >( deformationPartial )->wrtOtherStateDerivative();
    std::cout << "multiplyingFactorDeformation" << std::endl;
    std::cout << multiplyingFactorDeformation << std::endl; 
      

    // Check partials wrt Io's translational state    
    Eigen::MatrixXd deformationPartialWrtTranslationalState = Eigen::MatrixXd::Zero( 5, 6 );
    deformationPartial->wrtNonDeformationStateOfAdditionalBody( deformationPartialWrtTranslationalState.block(0, 0, 5, 6), std::make_pair( "Io", "" ), propagators::translational_state );
    
    Eigen::MatrixXd torquePartialWrtTranslationalState = Eigen::MatrixXd::Zero( 3, 6 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
            torquePartialWrtTranslationalState.block( 0, 0, 3, 6 ), std::make_pair( "Io", "" ), propagators::translational_state );

    torquePartialWrtTranslationalState += ( multiplyingFactorTorque * deformationPartialWrtTranslationalState );

    // Check partials with respect to Io's rotational state
    Eigen::MatrixXd deformationPartialWrtRotationalState = Eigen::MatrixXd::Zero( 5, 7 );
    deformationPartial->wrtNonDeformationStateOfAdditionalBody( deformationPartialWrtRotationalState.block(0, 0, 5, 7), std::make_pair( "Io", "" ), propagators::rotational_state );

    Eigen::MatrixXd directTorquePartialWrtRotationalState = Eigen::MatrixXd::Zero( 3, 7 );
    torquePartial->wrtRotationalStateOfAcceleratedBody( directTorquePartialWrtRotationalState.block( 0, 0, 3, 7 ) );

    Eigen::MatrixXd torquePartialWrtRotationalState = directTorquePartialWrtRotationalState + multiplyingFactorTorque * deformationPartialWrtRotationalState;  

    // Check partials with respect to Io's gravity deformation  
    Eigen::MatrixXd deformationPartialWrtGravityState = Eigen::MatrixXd::Zero( 5, 5 );
    deformationPartial->wrtStateOfDeformingBody( deformationPartialWrtGravityState.block(0, 0, 5, 5) );
    
    Eigen::MatrixXd directTorquePartialWrtGravityState = Eigen::MatrixXd::Zero( 3, 5 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
            directTorquePartialWrtGravityState.block( 0, 0, 3, 5 ), std::make_pair( "Io", "" ), propagators::gravity_deformation_state );

    Eigen::MatrixXd torquePartialWrtGravityState = directTorquePartialWrtGravityState + multiplyingFactorTorque * deformationPartialWrtGravityState;


   // Declare numerical partials.
   Eigen::Matrix< double, 3, 6 > testTorquePartialWrtTranslationalState = Eigen::Matrix< double, 3, 6 >::Zero( );
   Eigen::Matrix< double, 3, 4 > testTorquexPartialWrtOrientation = Eigen::Matrix< double, 3, 4 >::Zero( );
   Eigen::Matrix< double, 3, 3 > testTorquePartialWrtRotationalVelocity = Eigen::Matrix< double, 3, 3 >::Zero( );
   Eigen::Matrix< double, 3, 5 > testTorquePartialWrtGravityState = Eigen::Matrix< double, 3, 5 >::Zero( );

    // Declare perturbations for numerical partials
    Eigen::Vector4d orientationPerturbation;
    orientationPerturbation << 1.0E-9, 1.0E-9, 1.0E-9, 1.0E-9;
    Eigen::Vector3d rotationalVelocityPerturbation;
    rotationalVelocityPerturbation << 1.0E-6, 1.0E-6, 1.0E-6;

    Eigen::Vector6d translationalStatePerturbation;
    translationalStatePerturbation << 1.0, 1.0, 100.0, 1.0E-3, 1.0E-3, 1.0E-3;

    Eigen::Vector5d gravityDeformationPerturbation;
    gravityDeformationPerturbation << 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6;

    // Create state access/modification functions for bodies.
    std::function< void( Eigen::VectorXd ) > ioGravityDeformationSetFunction = 
        std::bind( &Body::setCurrentPropagatedGravityField, io, std::placeholders::_1 ); 
    std::function< void( Eigen::Vector7d ) > ioRotationalStateSetFunction =
            std::bind( &Body::setCurrentRotationalStateToLocalFrame, io, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > ioTranslationalStateSetFunction = std::bind( &Body::setState, io, std::placeholders::_1 );
       
   std::function< void( Eigen::VectorXd ) > updateAngularVelocityDerivative = 
        std::bind( &simulation_setup::Body::setCurrentAngularVelocityDerivativeVectorInLocalFrame, io, std::placeholders::_1 ); 

    // Calculate numerical partials.
    std::vector< Eigen::Vector4d > appliedQuaternionPerturbation;
    Eigen::MatrixXd torqueDeviationsDueToOrientation = calculateTorqueDeviationViaGravityDerivativeDependencyDueToOrientationChange( 
        ioRotationalStateSetFunction, maxwellDeformationModel, inertialTorqueModel, 
        io->getRotationalStateVector( ), nominalDeformation,
        orientationPerturbation, appliedQuaternionPerturbation, 
        std::bind( &RigidBodyProperties::updateInertiaTensorDerivative, io->getMassProperties( ), std::placeholders::_1 ) );

    testTorquePartialWrtTranslationalState = calculateTorqueDeviationDueToGravityDerivativeDependency( 
        ioTranslationalStateSetFunction,
        maxwellDeformationModel, inertialTorqueModel, io->getState( ), nominalDeformation, 
        translationalStatePerturbation, 6, 0,
        std::bind( &RigidBodyProperties::updateInertiaTensorDerivative, io->getMassProperties( ), std::placeholders::_1 ) );

    testTorquePartialWrtRotationalVelocity = calculateTorqueDeviationDueToGravityDerivativeDependency( 
        ioRotationalStateSetFunction,
        maxwellDeformationModel, inertialTorqueModel, io->getRotationalStateVector( ), nominalDeformation, 
        rotationalVelocityPerturbation, 3, 4,
        std::bind( &RigidBodyProperties::updateInertiaTensorDerivative, io->getMassProperties( ), std::placeholders::_1 ) );
 
    testTorquePartialWrtGravityState = calculateTorqueDeviationDueToGravityDerivativeDependency( 
        ioGravityDeformationSetFunction,
        maxwellDeformationModel, inertialTorqueModel, ioGravityDeformation, nominalDeformation, 
        gravityDeformationPerturbation, 5, 0,
        std::bind( &RigidBodyProperties::updateInertiaTensorDerivative, io->getMassProperties( ), std::placeholders::_1 ) );

    // Compare numerical and analytical results.
    for( int index = 1; index < 4; index++ )
    {
        Eigen::Vector3d numericalChangeInDeformation = torqueDeviationsDueToOrientation.block( 0, index - 1, 3, 1 );
        Eigen::Vector3d analyticalChangeInDeformation =
                torquePartialWrtRotationalState.block( 0, 0, 3, 1 ) * appliedQuaternionPerturbation[ index ]( 0 ) +
                torquePartialWrtRotationalState.block( 0, index, 3, 1 ) * appliedQuaternionPerturbation[ index ]( index );

        std::cout << "numericalChangeInDeformation" << std::endl;
        std::cout << numericalChangeInDeformation << std::endl;
        std::cout << "analyticalChangeInDeformation" << std::endl;
        std::cout << analyticalChangeInDeformation << std::endl;

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalChangeInDeformation, analyticalChangeInDeformation, 1.0e-6 );
    }

    std::cout << "testTorquePartialWrtTranslationalState" << std::endl;
    std::cout << testTorquePartialWrtTranslationalState << std::endl;

    std::cout << "deformationPartialWrtRotationalState" << std::endl;
    std::cout << deformationPartialWrtRotationalState << std::endl;
    std::cout << "direct torquePartialWrtRotationalState" << std::endl;
    std::cout << directTorquePartialWrtRotationalState << std::endl;
    std::cout << "full torquePartialWrtRotationalState" << std::endl;
    std::cout << torquePartialWrtRotationalState << std::endl;

    std::cout << "testTorquePartialWrtRotationalVelocity" << std::endl;
    std::cout << testTorquePartialWrtRotationalVelocity << std::endl;

    std::cout << "deformationPartialWrtGravityState" << std::endl;
    std::cout << deformationPartialWrtGravityState << std::endl;
    std::cout << "direct torquePartialWrtGravityState" << std::endl;
    std::cout << directTorquePartialWrtGravityState << std::endl;
    std::cout << "full torquePartialWrtGravityState" << std::endl;
    std::cout << torquePartialWrtGravityState << std::endl;

    std::cout << "testTorquePartialWrtGravityState" << std::endl;
    std::cout << testTorquePartialWrtGravityState << std::endl;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testTorquePartialWrtTranslationalState, torquePartialWrtTranslationalState, 1.0e-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testTorquePartialWrtRotationalVelocity, torquePartialWrtRotationalState.block(0, 4, 3, 3), 1.0e-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testTorquePartialWrtGravityState, torquePartialWrtGravityState, 1.0e-6 );


//     for( int i = 0; i < 3; i++ )
//     {
//         BOOST_CHECK_SMALL( std::fabs( partialWrtMeanMomentOfInertia( i, 0 ) - testPartialWrtMeanMomentOfInertia( i, 0 ) ), 1.0E1 );
//     }
}


BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
