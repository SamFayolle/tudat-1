/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEGRAVITYDEFORMATIONPARTIALS_H
#define TUDAT_CREATEGRAVITYDEFORMATIONPARTIALS_H

#include <memory.h>

#include "tudat/astro/basic_astro/gravityDeformationModel.h"
#include "tudat/astro/basic_astro/gravityDeformationModelTypes.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/orbit_determination/gravity_deformation_partials/deformationPartial.h"
#include "tudat/astro/orbit_determination/gravity_deformation_partials/maxwellDeformationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"
#include "tudat/simulation/estimation_setup/createCartesianStatePartials.h"
#include "tudat/simulation/propagation_setup/environmentUpdater.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create a single gravity deformation partial derivative object.
/*!
 *  Function to create a single gravity deformation partial derivative object.
 *  \param gravityDeformationModel Gravity deformation model for which a partial derivative is to be computed.
 *  \param deformingBody Pair of name and object of body undergoing deformation
 *  \param bodies List of all body objects
 *  \param parametersToEstimate List of parameters that are to be estimated. Empty by default, only required for selected
 *  types of partials (e.g. maxwell deformation w.r.t. rotational parameters).
 *  \return Single gravity deformation partial derivative object.
 */
template< typename InitialStateParameterType = double >
std::shared_ptr< acceleration_partials::DeformationPartial > createAnalyticalGravityDeformationPartial(
        std::shared_ptr< basic_astrodynamics::GravityDeformationModel > gravityDeformationModel,
        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > deformingBody,
        const simulation_setup::SystemOfBodies& bodies = simulation_setup::SystemOfBodies( ),
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > > parametersToEstimate =
                std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >( ) )
{
    using namespace gravitation;
    using namespace basic_astrodynamics;
    using namespace electromagnetism;
    using namespace aerodynamics;
    using namespace acceleration_partials;

    std::shared_ptr< acceleration_partials::DeformationPartial > deformationPartial;

    // Identify current gravity deformation model type
    GravityDeformationType deformationType = getGravityDeformationModelType( gravityDeformationModel );
    switch( deformationType )
    {
        case maxwell_deformation:
        {
            // Check if identifier is consistent with type.
            if( std::dynamic_pointer_cast< MaxwellGravityDeformationModel >( gravityDeformationModel ) == nullptr )
            {
                throw std::runtime_error(
                        "Deformation class type does not match deformation type (maxwell_deformation) when making gravity deformation partial" );
            }
            else
            {
                // Create partial-calculating object.
                std::map< std::pair< estimatable_parameters::EstimatebleParametersEnum, std::string >, std::shared_ptr< observation_partials::RotationMatrixPartial > >
                        rotationMatrixPartials = observation_partials::createRotationMatrixPartials( parametersToEstimate, deformingBody.first, bodies );

                std::vector< std::string > perturbingBodies = std::dynamic_pointer_cast< MaxwellGravityDeformationModel >( gravityDeformationModel )->getPerturbingBody( );

                deformationPartial = std::make_shared< acceleration_partials::MaxwellDeformationPartial >( 
                    deformingBody.first, 
                    perturbingBodies, 
                    std::dynamic_pointer_cast< MaxwellGravityDeformationModel >( gravityDeformationModel ), 
                    rotationMatrixPartials );
            }
            break;
        }
        default:
            std::string errorMessage = "Deformation model " + std::to_string( deformationType ) + " not found when making gravity deformation partial";
            throw std::runtime_error( errorMessage );
            break;
    }

    return deformationPartial;
}
 

//! This function creates gravity deformation partial objects for deformation dynamics
/*!
 *  This function creates gravity deformation partial objects for deformation dynamics from deformation models and
 *  list of bodies' states of which derivatives are needed. The return type is an StateDerivativePartialsMap,
 *  a standardized type for communicating such lists of these objects.
 *  \param deformationMap Map of maps containing list of deformation models, identifying which deformation acts on which
 *   body.
 *  \param bodies List of body objects constituting environment for calculations.
 *  \param parametersToEstimate List of parameters which are to be estimated.
 *  \return List of gravity-deformation-partial-calculating objects in StateDerivativePartialsMap type.
 */
template< typename InitialStateParameterType >
orbit_determination::StateDerivativePartialsMap createGravityDeformationPartialsMap(
        const basic_astrodynamics::GravityDeformationModelMap& deformationMap,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > > parametersToEstimate )
{
    // Declare return map.
    orbit_determination::StateDerivativePartialsMap gravityDeformationPartialsList;

    std::vector< std::shared_ptr<
            estimatable_parameters::EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
            initialDynamicalParameters = estimatable_parameters::getListOfGravityDeformationStateParametersToEstimate( parametersToEstimate );

    gravityDeformationPartialsList.resize( initialDynamicalParameters.size( ) );

    // Iterate over list of bodies of which the partials of the gravity deformation acting on them are required.
    for( basic_astrodynamics::GravityDeformationModelMap::const_iterator deformationIterator = deformationMap.begin( ); deformationIterator != deformationMap.end( );
         deformationIterator++ )
    {
        for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
        {
            if( initialDynamicalParameters.at( i )->getParameterName( ).second.first == deformationIterator->first )
            {
                if( ( initialDynamicalParameters.at( i )->getParameterName( ).first ==
                      estimatable_parameters::initial_gravity_deformation_state ) )
                {
                    // Get object for body undergoing deformation
                    const std::string deformingBody = deformationIterator->first;
                    std::shared_ptr< simulation_setup::Body > deformingBodyObject = bodies.at( deformingBody );

                    // Retrieve list of deformations acting on current body.
                    std::vector< std::shared_ptr< GravityDeformationModel > > deformationVector = deformationMap.at( deformingBody );

                    // Declare list of deformation partials of current body.
                    std::vector< std::shared_ptr< orbit_determination::StateDerivativePartial > > deformationPartialVector;

                    // Iterate over all deformation models and generate their partial-calculating objects.
                    for( unsigned int k = 0 ; k < deformationVector.size( ) ; k++ )
                    {
                        // Create single partial object
                        std::shared_ptr< acceleration_partials::DeformationPartial > currentDeformationPartial =
                                createAnalyticalGravityDeformationPartial( deformationVector[ k ],
                                                                std::make_pair( deformingBody, deformingBodyObject ),
                                                                bodies,
                                                                parametersToEstimate );

                        deformationPartialVector.push_back( currentDeformationPartial );
                        
                    }

                    // Add partials of current body's deformation to list.
                    gravityDeformationPartialsList[ i ] = deformationPartialVector;
                }
            }
        }
    }
    return gravityDeformationPartialsList;
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_CREATEGRAVITYDEFORMATIONPARTIALS_H
