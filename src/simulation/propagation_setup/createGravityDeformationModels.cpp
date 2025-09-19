/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <algorithm>
#include <functional>
#include <memory>

#include "tudat/simulation/propagation_setup/createGravityDeformationModels.h"


namespace tudat
{

namespace simulation_setup
{

using namespace aerodynamics;
using namespace gravitation;
using namespace basic_astrodynamics;
using namespace electromagnetism;
using namespace ephemerides;


std::shared_ptr< basic_astrodynamics::MaxwellGravityDeformationModel >
createMaxwellGravityFieldDeformationModel(
        const std::shared_ptr< simulation_setup::Body > deformingBody,
        const std::shared_ptr< simulation_setup::Body > perturbingBody,
        const std::string& nameOfDeformingBody,
        const std::string& nameOfPerturbingBody,
        const std::shared_ptr< GravityDeformationSettings > deformationSettings ) 
{ 
    // Declare pointer to return object
    std::shared_ptr< MaxwellGravityDeformationModel > deformationModel;

    // Dynamic cast deformation settings to required type and check consistency.
    std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings =
            std::dynamic_pointer_cast< MaxwellDeformationSettings >( deformationSettings );
    if( maxwellDeformationSettings == nullptr )
    {
        throw std::runtime_error( 
            std::string( "Error, deformation settings inconsistent ") + " making maxwell gravity deformation of " 
            + nameOfDeformingBody + " due to " + nameOfPerturbingBody );
    }
    else
    {
        // Get pointer to gravity field and rotational ephemeris of deforming body and cast to required type.
        std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityField =
            std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( deformingBody->getGravityFieldModel( ) );

        std::shared_ptr< RotationalEphemeris> rotationalEphemeris = deformingBody->getRotationalEphemeris( );
        if( sphericalHarmonicsGravityField == nullptr )
        {
            throw std::runtime_error(
                        std::string( "Error, spherical harmonic gravity field model not set when ")
                        + " creating Maxwell gravity deformation model of " + nameOfDeformingBody );
        }
        else
        {
            if( rotationalEphemeris == nullptr )
            {
                throw std::runtime_error( "Warning when creating Maxwell gravity deformation of " + nameOfDeformingBody 
                    + "no rotation model found for " + nameOfDeformingBody );  
            }

            // Create gravity deformation object.
            deformationModel = std::make_shared< MaxwellGravityDeformationModel >(
                    std::bind( &Body::getStateByReference, deformingBody, std::placeholders::_1 ),
                    nameOfPerturbingBody,
                    maxwellDeformationSettings->maxwellRelaxationTime_,
                    maxwellDeformationSettings->globalRelaxationTime_,
                    sphericalHarmonicsGravityField->getGravitationalParameter( ),
                    perturbingBody->getGravitationalParameter( ),
                    sphericalHarmonicsGravityField->getReferenceRadius( ),
                    std::bind( &Body::getCurrentAngularVelocityVectorInLocalFrame, deformingBody ),
                    std::bind( &Body::getCurrentAngularVelocityDerivativeVectorInLocalFrame, deformingBody ),
                    maxwellDeformationSettings->loveNumber_,
                    std::bind( &SphericalHarmonicsGravityField::getCosineCoefficientsBlock,
                                sphericalHarmonicsGravityField,
                                maxwellDeformationSettings->maximumDegree_,
                                maxwellDeformationSettings->maximumOrder_ ),
                    std::bind( &SphericalHarmonicsGravityField::getSineCoefficientsBlock,
                                sphericalHarmonicsGravityField,
                                maxwellDeformationSettings->maximumDegree_,
                                maxwellDeformationSettings->maximumOrder_ ),
                    std::bind( &Body::getStateByReference, perturbingBody, std::placeholders::_1 ),
                    std::bind( &Body::getCurrentRotationToGlobalFrame, deformingBody ),
                    std::bind( &Body::getCurrentRotationMatrixDerivativeToLocalFrame, deformingBody ),
                    maxwellDeformationSettings->staticCoefficients_,
                    maxwellDeformationSettings->includeOrder1_ );
        }
    }
    return deformationModel;
};

// //! Function to create acceleration model object.
// std::shared_ptr< AccelerationModel< Eigen::Vector3d > > createAccelerationModel(
//         const std::shared_ptr< Body > bodyUndergoingAcceleration,
//         const std::shared_ptr< Body > bodyExertingAcceleration,
//         const std::shared_ptr< AccelerationSettings > accelerationSettings,
//         const std::string& nameOfBodyUndergoingAcceleration,
//         const std::string& nameOfBodyExertingAcceleration,
//         const std::shared_ptr< Body > centralBody,
//         const std::string& nameOfCentralBody,
//         const SystemOfBodies& bodies )
// {
//     // Declare pointer to return object.
//     std::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModelPointer;

//     // Switch to call correct acceleration model type factory function.
//     switch( accelerationSettings->accelerationType_ )
//     {
//     case point_mass_gravity:
//     case spherical_harmonic_gravity:
//     case mutual_spherical_harmonic_gravity:
//     case polyhedron_gravity:
//     case ring_gravity:
//         accelerationModelPointer = createGravitationalAccelerationModel(
//                     bodyUndergoingAcceleration, bodyExertingAcceleration, accelerationSettings,
//                     nameOfBodyUndergoingAcceleration, nameOfBodyExertingAcceleration,
//                     centralBody, nameOfCentralBody );
//         break;
//     case aerodynamic:
//         accelerationModelPointer = createAerodynamicAcceleratioModel(
//                     bodyUndergoingAcceleration,
//                     bodyExertingAcceleration,
//                     nameOfBodyUndergoingAcceleration,
//                     nameOfBodyExertingAcceleration );
//         break;
//     case radiation_pressure:
//         accelerationModelPointer = createRadiationPressureAccelerationModel(
//                     bodyUndergoingAcceleration,
//                     bodyExertingAcceleration,
//                     nameOfBodyUndergoingAcceleration,
//                     nameOfBodyExertingAcceleration,
//                     bodies,
//                     accelerationSettings );
//         break;
//     case cannon_ball_radiation_pressure:
//         accelerationModelPointer = createCannonballRadiationPressureAcceleratioModel(
//             bodyUndergoingAcceleration,
//             bodyExertingAcceleration,
//             nameOfBodyUndergoingAcceleration,
//             nameOfBodyExertingAcceleration,
//             bodies);
//         break;
//     case thrust_acceleration:
//         accelerationModelPointer = createThrustAcceleratioModel(
//                     accelerationSettings, bodies,
//                     nameOfBodyUndergoingAcceleration );
//         break;
//     case relativistic_correction_acceleration:
//         accelerationModelPointer = createRelativisticCorrectionAcceleration(
//                     bodyUndergoingAcceleration,
//                     bodyExertingAcceleration,
//                     nameOfBodyUndergoingAcceleration,
//                     nameOfBodyExertingAcceleration,
//                     accelerationSettings, bodies );
//         break;
//     case empirical_acceleration:
//         accelerationModelPointer = createEmpiricalAcceleration(
//                     bodyUndergoingAcceleration,
//                     bodyExertingAcceleration,
//                     nameOfBodyUndergoingAcceleration,
//                     nameOfBodyExertingAcceleration,
//                     accelerationSettings );
//         break;
//     case direct_tidal_dissipation_in_central_body_acceleration:
//         accelerationModelPointer = createDirectTidalDissipationAcceleration(
//                     bodyUndergoingAcceleration,
//                     bodyExertingAcceleration,
//                     nameOfBodyUndergoingAcceleration,
//                     nameOfBodyExertingAcceleration,
//                     accelerationSettings );
//         break;
//     case direct_tidal_dissipation_in_orbiting_body_acceleration:
//         accelerationModelPointer = createDirectTidalDissipationAcceleration(
//                     bodyUndergoingAcceleration,
//                     bodyExertingAcceleration,
//                     nameOfBodyUndergoingAcceleration,
//                     nameOfBodyExertingAcceleration,
//                     accelerationSettings );
//         break;
//     case momentum_wheel_desaturation_acceleration:
//         accelerationModelPointer = createMomentumWheelDesaturationAcceleration(
//                     bodyUndergoingAcceleration,
//                     bodyExertingAcceleration,
//                     nameOfBodyUndergoingAcceleration,
//                     nameOfBodyExertingAcceleration,
//                     accelerationSettings );
//         break;
//     case yarkovsky_acceleration:
//         accelerationModelPointer = createYarkovskyAcceleration(
//                     bodyUndergoingAcceleration,
//                     bodyExertingAcceleration,
//                     nameOfBodyUndergoingAcceleration,
//                     nameOfBodyExertingAcceleration,
//                     accelerationSettings );
//         break;
//     case custom_acceleration:
//         accelerationModelPointer = createCustomAccelerationModel(
//                     accelerationSettings,
//                     nameOfBodyUndergoingAcceleration );
//         break;
//     default:
//         throw std::runtime_error(
//                     std::string( "Error, acceleration model ") +
//                     std::to_string( accelerationSettings->accelerationType_ ) +
//                     " not recognized when making acceleration model of" +
//                     nameOfBodyExertingAcceleration + " on " +
//                     nameOfBodyUndergoingAcceleration );
//         break;
//     }
//     return accelerationModelPointer;
// }

//! Function to create a list of mass rate models for a list of bodies.
basic_astrodynamics::GravityDeformationModelMap createGravityDeformationModelsMap(
        const SystemOfBodies& bodies,
        const SelectedGravityDeformationModelMap& gravityDeformationSettings )
{
    // Iterate over all bodies
    std::map< std::string, std::vector< std::shared_ptr< basic_astrodynamics::GravityDeformationModel > > > gravityDeformationModels;
    for( std::map< std::string, std::vector< std::shared_ptr< GravityDeformationSettings > > >::const_iterator settingsIterator =
         gravityDeformationSettings.begin( ); settingsIterator != gravityDeformationSettings.end( ); settingsIterator++)
    {
        // Iterate over all mass model settings for current body.
        for( unsigned int i = 0; i < settingsIterator->second.size( ); i++ )
        {
            switch ( settingsIterator->second.at( i )->deformationType_ )
            {
            case maxwell_deformation:
            {
                std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings = std::dynamic_pointer_cast< MaxwellDeformationSettings >( settingsIterator->second.at( i ) );
                std::string deformingBody = settingsIterator->first;
                std::string perturbingBody = maxwellDeformationSettings->perturbingBody_;
                gravityDeformationModels[ settingsIterator->first ].push_back( createMaxwellGravityFieldDeformationModel( bodies.at( deformingBody ), bodies.at( perturbingBody ),
                deformingBody, perturbingBody, settingsIterator->second.at( i ) ) );
                break;
            }
            default:
                break;
            }
        }
    }
    return gravityDeformationModels;

}


} // namespace simulation_setup

} // namespace tudat