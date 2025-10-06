/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEGRAVITYDEFORMATIONMODELS_H
#define TUDAT_CREATEGRAVITYDEFORMATIONMODELS_H

#include <vector>
#include <string>
#include <memory>
#include <functional>

#include "tudat/astro/basic_astro/gravityDeformationModel.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/propagation_setup/gravityDeformationSettings.h"

namespace tudat
{

namespace simulation_setup
{

using namespace basic_astrodynamics;
using namespace gravitation;
using namespace ephemerides;


//! Function to create Maxwell gravity deformation model.
/*!
 *  Function to create Maxwell gravity deformation model from perturbing and deforming bodies.
 *  \param deformingBody Pointer to object of deforming body.
 *  \param perturbingBody Pointer to object of perturbing body.
 *  \param nameOfDeformingBody Name of deforming body.
 *  \param nameOfPerturbingBody Name of perturbing body.
 *  \param deformationSettings Settings for gravity deformation model that is to be created.
 *  \return Maxwell gravity deformation model pointer.
 */
std::shared_ptr< basic_astrodynamics::MaxwellGravityDeformationModel >
createMaxwellGravityFieldDeformationModel(
        const std::shared_ptr< simulation_setup::Body > deformingBody,
        const std::vector< std::shared_ptr< simulation_setup::Body > > perturbingBody,
        const std::string& nameOfDeformingBody,
        const std::vector< std::string >& nameOfPerturbingBody,
        const std::shared_ptr< GravityDeformationSettings > deformationSettings ); 

basic_astrodynamics::GravityDeformationModelMap createGravityDeformationModelsMap(
        const SystemOfBodies& bodies,
        const SelectedGravityDeformationModelMap& gravityDeformationSettings );        
// { 
//     // Declare pointer to return object
//     std::shared_ptr< MaxwellGravityDeformationModel > deformationModel;

//     // Dynamic cast deformation settings to required type and check consistency.
//     std::shared_ptr< MaxwellDeformationSettings > maxwellDeformationSettings =
//             std::dynamic_pointer_cast< MaxwellDeformationSettings >( deformationSettings );
//     if( maxwellDeformationSettings == nullptr )
//     {
//         throw std::runtime_error( 
//             std::string( "Error, deformation settings inconsistent ") + " making maxwell gravity deformation of " 
//             + nameOfDeformingBody + " due to " + nameOfPerturbingBody );
//     }
//     else
//     {
//         // Get pointer to gravity field and rotational ephemeris of deforming body and cast to required type.
//         std::shared_ptr< SphericalHarmonicsGravityField > sphericalHarmonicsGravityField =
//             std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( deformingBody->getGravityFieldModel( ) );

//         std::shared_ptr< RotationalEphemeris> rotationalEphemeris = deformingBody->getRotationalEphemeris( );
//         if( sphericalHarmonicsGravityField == nullptr )
//         {
//             throw std::runtime_error(
//                         std::string( "Error, spherical harmonic gravity field model not set when ")
//                         + " creating Maxwell gravity deformation model of " + nameOfDeformingBody );
//         }
//         else
//         {
//             if( rotationalEphemeris == nullptr )
//             {
//                 throw std::runtime_error( "Warning when creating Maxwell gravity deformation of " + nameOfDeformingBody 
//                     + "no rotation model found for " + nameOfDeformingBody );  
//             }

//             // Create gravity deformation object.
//             deformationModel = std::make_shared< MaxwellGravityDeformationModel >(
//                     std::bind( &Body::getPositionByReference, deformingBody, std::placeholders::_1 ),
//                     maxwellDeformationSettings->maxwellRelaxationTime_,
//                     maxwellDeformationSettings->globalRelaxationTime_,
//                     sphericalHarmonicsGravityField->getGravitationalParameter( ),
//                     perturbingBody->getGravitationalParameter( ),
//                     sphericalHarmonicsGravityField->getReferenceRadius( ),
//                     maxwellDeformationSettings->rotationRate_,
//                     maxwellDeformationSettings->loveNumber_,
//                     std::bind( &SphericalHarmonicsGravityField::getCosineCoefficientsBlock,
//                                 sphericalHarmonicsGravityField,
//                                 maxwellDeformationSettings->maximumDegree_,
//                                 maxwellDeformationSettings->maximumOrder_ ),
//                     std::bind( &SphericalHarmonicsGravityField::getSineCoefficientsBlock,
//                                 sphericalHarmonicsGravityField,
//                                 maxwellDeformationSettings->maximumDegree_,
//                                 maxwellDeformationSettings->maximumOrder_ ),
//                     std::bind( &Body::getPositionByReference, perturbingBody, std::placeholders::_1 ),
//                     std::bind( &Body::getCurrentRotationToGlobalFrame, deformingBody ) );
//         }
//     }
//     return deformationModel;
// };


// //! Function to create acceleration model object.
// /*!
//  *  Function to create acceleration model object.
//  *  Type of requested model is checked and corresponding factory function is called.
//  *  \param bodyUndergoingAcceleration Pointer to object of body that is being accelerated.
//  *  \param bodyExertingAcceleration Pointer to object of body that is exerting acceleration,
//  *  \param accelerationSettings Settings for acceleration model that is to be created.
//  *  \param nameOfBodyUndergoingAcceleration Name of object of body that is being accelerated.
//  *  \param nameOfBodyExertingAcceleration Name of object of body that is exerting the acceleration.
//  *  \param centralBody Pointer to central body in frame centered at which acceleration is to be
//  *  calculated (optional, only relevant for third body accelerations).
//  *  \param nameOfCentralBody Name of central body in frame cenetered at which acceleration is to
//  *  be calculated (optional, only relevant for third body accelerations).
//  *  \param bodies List of pointers to bodies required for the creation of the acceleration model
//  *  objects.
//  *  \return Acceleration model pointer.
//  */
// std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
// createAccelerationModel(
//         const std::shared_ptr< Body > bodyUndergoingAcceleration,
//         const std::shared_ptr< Body > bodyExertingAcceleration,
//         const std::shared_ptr< AccelerationSettings > accelerationSettings,
//         const std::string& nameOfBodyUndergoingAcceleration,
//         const std::string& nameOfBodyExertingAcceleration,
//         const std::shared_ptr< Body > centralBody = std::shared_ptr< Body >( ),
//         const std::string& nameOfCentralBody = "",
//         const SystemOfBodies& bodies = SystemOfBodies( ) );

// //! Function to put SelectedAccelerationMap in correct order, to ensure correct model creation
// /*!
//  * Function to put SelectedAccelerationMap in correct order, to ensure correct model creation
//  * \param selectedAccelerationPerBody List of acceleration settings per body.
//  * \return selectedAccelerationPerBody, put in order to ensure correct model creation.
//  */
// SelectedAccelerationList orderSelectedAccelerationMap( const SelectedAccelerationMap& selectedAccelerationPerBody );

// //! Function to create acceleration models from a map of bodies and acceleration model types.
// /*!
//  *  Function to create acceleration models from a map of bodies and acceleration model types.
//  *  The return type can be used to identify both the body undergoing and exerting acceleration.
//  *  \param bodies List of pointers to bodies required for the creation of the acceleration model
//  *  objects.
//  *  \param selectedAccelerationPerBody List identifying which bodies exert which type of
//  *  acceleration(s) on which bodies.
//  *  \param centralBodies Map of central bodies for each body undergoing acceleration.
//  *  \return List of acceleration model objects, in form of AccelerationMap.
//  */
// //basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
// //        const SystemOfBodies& bodies,
// //        const SelectedAccelerationMap& selectedAccelerationPerBody,
// //        const std::map< std::string, std::string >& centralBodies );
// //! Function to create a set of acceleration models from a map of bodies and acceleration model types.
// basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
//         const SystemOfBodies& bodies,
//         const SelectedAccelerationMap& selectedAccelerationPerBody,
//         const std::map< std::string, std::string >& centralBodies );

// //! Function to create acceleration models from a map of bodies and acceleration model types.
// /*!
//  *  Function to create acceleration models from a map of bodies and acceleration model types.
//  *  The return type can be used to identify both the body undergoing and exerting acceleration.
//  *  \param bodies List of pointers to bodies required for the creation of the acceleration model
//  *  objects.
//  *  \param selectedAccelerationPerBody List identifying which bodies exert which type of
//  *  acceleration(s) on which bodies.
//  *  \param propagatedBodies List of bodies that are to be propagated
//  *  \param centralBodies List of central bodies for each body undergoing acceleration (in same order as propagatedBodies).
//  *  \return List of acceleration model objects, in form of AccelerationMap.
//  */
// basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
//         const SystemOfBodies& bodies,
//         const SelectedAccelerationMap& selectedAccelerationPerBody,
//         const std::vector< std::string >& propagatedBodies,
//         const std::vector< std::string >& centralBodies );


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEGRAVITYDEFORMATIONMODELS_H