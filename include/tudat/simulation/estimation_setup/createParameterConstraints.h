/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEPARAMETERCONSTRAINTS_H
#define TUDAT_CREATEPARAMETERCONSTRAINTS_H

#include "tudat/astro/basic_astro/accelerationModel.h"

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialRotationalState.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialMassState.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantDragCoefficient.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationalOrientation.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/empiricalAccelerationCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravitationalParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/observationBiasParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/groundStationPosition.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicCosineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicSineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/radiationPressureCoefficient.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/ppnParameters.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/equivalencePrincipleViolationParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/tidalLoveNumber.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/directTidalTimeLag.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/inverseTidalQualityFactor.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/meanMomentOfInertiaParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/desaturationDeltaV.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/periodicSpinVariation.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/polarMotionAmplitude.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/coreFactor.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/freeCoreNutationRate.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/desaturationDeltaV.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/longitudeLibrationAmplitude.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantThrust.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/yarkovskyParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/referencePointPosition.h"
#include "tudat/astro/relativity/metric.h"
#include "tudat/astro/basic_astro/accelerationModelTypes.h"
#include "tudat/simulation/estimation_setup/estimatableParameterSettings.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/parameterConstraints.h"

namespace tudat
{

namespace simulation_setup
{

//template< typename InitialStateParameterType = double >
std::vector< estimatable_parameters::EstimatebleParametersEnum > getParametersTypesPerConstraint(
        const estimatable_parameters::ParameterConstraintsEnum constraintType );


template< typename InitialStateParameterType = double >
std::set< estimatable_parameters::ParameterConstraintsEnum > getPotentialConstraints(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parameters )
{
    using namespace estimatable_parameters;

    std::set< EstimatebleParametersEnum > listParameterTypes = parameters->getListParameterTypes( );

    std::set< estimatable_parameters::ParameterConstraintsEnum > relevantConstraints;

    if ( listParameterTypes.count( single_degree_variable_tidal_love_number ) > 0 &&
            listParameterTypes.count( inverse_tidal_quality_factor ) > 0 )
    {
        relevantConstraints.insert( tidal_quality_factor_single_love_number_constraint );
    }
    if ( listParameterTypes.count( single_degree_variable_tidal_love_number ) > 0 &&
         listParameterTypes.count( direct_dissipation_tidal_time_lag ) > 0 )
    {
        relevantConstraints.insert( tidal_time_lag_single_love_number_constraint );
    }
    if ( listParameterTypes.count( full_degree_tidal_love_number ) > 0 &&
         listParameterTypes.count( inverse_tidal_quality_factor ) > 0 )
    {
        relevantConstraints.insert( tidal_quality_factor_full_love_number_constraint );
    }
    if ( listParameterTypes.count( full_degree_tidal_love_number ) > 0 &&
         listParameterTypes.count( direct_dissipation_tidal_time_lag ) > 0 )
    {
        relevantConstraints.insert( tidal_time_lag_full_love_number_constraint );
    }

    return relevantConstraints;
}

//template< typename InitialStateParameterType = double >
void constraintSingleOrderLoveNumber(
    const estimatable_parameters::ParameterConstraintsEnum constraint,
    const std::shared_ptr< estimatable_parameters::TidalLoveNumber< Eigen::VectorXd > > loveNumber,
    std::vector< std::pair< unsigned int, unsigned int > >& constraintIndices,
    Eigen::VectorXd& k2Value );


//! Function to ....
/*!
 *  Function to .....
 *  \param
 *  \return
 */
template< typename InitialStateParameterType = double >
void createParameterConstraints(
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parameters )

{
    using namespace tudat::estimatable_parameters;

    std::cout << "TEST, in createParameterConstraints function" << "\n\n";

//    std::cout << "list parameter types" << "\n\n";
//    std::set< EstimatebleParametersEnum > listParameterTypes = parameters->getListParameterTypes( );
//    for ( auto itr : listParameterTypes )
//    {
//        std::cout << itr << "\n\n";
//    }
//    std::cout << "end list parameter types" << "\n\n";

    // Retrieve list of potential constraints
    std::set< ParameterConstraintsEnum > potentialConstraints = getPotentialConstraints( parameters );

    // Retrieve list of potentially constrained parameters
    std::set< EstimatebleParametersEnum > potentiallyConstrainedParameters;
    std::cout << "list potential constraints: " << "\n\n";
    for ( auto constraint : potentialConstraints )
    {
        std::cout << constraint << "\n\n";
        if ( getParametersTypesPerConstraint( constraint ).size( ) < 2 )
        {
            throw std::runtime_error( "Error when creating constraint of type " + std::to_string( constraint )
                                      + ", the number of parameters involved should be > 2." );
        }
        for ( const auto param : getParametersTypesPerConstraint( constraint ) )
        {
            potentiallyConstrainedParameters.insert( param );
        }
    }
    std::cout << "\n\n";


    // Retrieve indices and bodies for potential constrained parameters
    std::map< EstimatebleParametersEnum, std::pair< std::vector< std::pair< std::string, std::string > >, std::vector< std::pair< int, int > > > >
            indicesPossiblyConstrainedParameters;
    for ( auto parameterType : potentiallyConstrainedParameters )
    {
        if ( isParameterPossiblyConstrained( parameterType ) )
        {
            std::pair< std::vector< std::pair< std::string, std::string > >, std::vector< std::pair< int, int > > > indicesParameter;

            switch( parameterType )
            {
                case single_degree_variable_tidal_love_number:
                case full_degree_tidal_love_number:
                {
                    std::pair< std::vector< std::pair< std::string, std::string > >, std::vector< std::pair< int, int > > > indicesAllDegrees =
                            parameters->getIndicesForParameterType( parameterType );
                    for ( unsigned int i = 0 ; i < indicesAllDegrees.second.size( ) ; i++ )
                    {
                        std::shared_ptr< TidalLoveNumber< Eigen::VectorXd > > loveNumber = std::dynamic_pointer_cast< TidalLoveNumber< Eigen::VectorXd > >(
                                parameters->getVectorParameters( ).at( indicesAllDegrees.second.at( i ).first ) );
                        if ( loveNumber->getDegree( ) == 2 && loveNumber->useComplexComponents( ) )
                        {
                            indicesParameter.first.push_back( indicesAllDegrees.first.at( i ) );
                            indicesParameter.second.push_back( indicesAllDegrees.second.at( i ) );
                        }
                    }
                    break;
                }
//                case full_degree_tidal_love_number:
//                {
//                    std::pair< std::vector< std::pair< std::string, std::string > >, std::vector< std::pair< int, int > > > indicesAllDegrees =
//                            parameters->getIndicesForParameterType( parameterType );
//                    for ( unsigned int i = 0 ; i < indicesAllDegrees.second.size( ) ; i++ )
//                    {
//                        std::shared_ptr< TidalLoveNumber< Eigen::VectorXd > > loveNumber = std::dynamic_pointer_cast< TidalLoveNumber< Eigen::VectorXd > >(
//                                parameters->getVectorParameters( ).at( indicesAllDegrees.second.at( i ).first ) );
//                        if ( loveNumber->getDegree( ) == 2 && loveNumber->useComplexComponents( ) )
//                        {
//                            indicesParameter.first.push_back( indicesAllDegrees.first.at( i ) );
//                            indicesParameter.second.push_back( indicesAllDegrees.second.at( i ) );
//                        }
//                    }
//                    break;
//                }
                default:
                {
                    indicesParameter = parameters->getIndicesForParameterType( parameterType );
                    break;
                }
            }

            if ( !indicesParameter.first.empty( ) )
            {
                indicesPossiblyConstrainedParameters[ parameterType ] = indicesParameter;
            }
        }
    }

    // print statements
    std::cout << "print indicesPossiblyConstrainedParameters" << "\n\n";
    for ( auto itr : indicesPossiblyConstrainedParameters )
    {
        std::cout << "param: " << itr.first << "\n\n";
        for ( unsigned int i = 0 ; i < itr.second.first.size( ) ; i++ )
        {
            std::cout << "bodies: " << itr.second.first.at( i ).first << " - " << itr.second.first.at( i ).second << "\n\n";
            std::cout << "indices: " << itr.second.second.at( i ).first << " - " << itr.second.second.at( i ).second << "\n\n";
        }
    }
    // end print statements

    // Check
    if ( !indicesPossiblyConstrainedParameters.empty( ) )
    {
        std::vector< std::shared_ptr< ParameterConstraints > > parameterConstraints;

        for ( auto constraint : potentialConstraints )
        {
            std::vector< EstimatebleParametersEnum > parametersInConstraint = getParametersTypesPerConstraint( constraint );
            unsigned int nbParametersInConstraint = parametersInConstraint.size( );

            // Initialise common parameters identifiers
            std::vector< std::pair< std::string, std::string > > detectedBodies;
            if ( indicesPossiblyConstrainedParameters.count( parametersInConstraint[ 0 ] ) > 0 )
            {
                detectedBodies = indicesPossiblyConstrainedParameters.at( parametersInConstraint[ 0 ] ).first;
            }

            // PRINTING STATEMENTS
            std::cout << "initial common parameters: " << "\n\n";
            for ( unsigned int i = 0 ; i < detectedBodies.size( ) ; i++ )
            {
                std::cout << detectedBodies[ i ].first << " - " << detectedBodies[ i ].second << "\n\n";
            }
            std::cout << "end parsing initial common parameters." << "\n\n";
            // END PRINTING STATEMENTS

            // TEST
            for ( unsigned int k = 0 ; k < detectedBodies.size( ) ; k++ )
            {
                std::pair< std::string, std::string > constrainedBody = detectedBodies.at( k );

                std::vector< std::pair< unsigned int, unsigned int > > constraintIndices =
                        { indicesPossiblyConstrainedParameters.at( parametersInConstraint[ 0 ] ).second.at( k ) };
                std::vector< std::pair< std::string, std::string > > constrainedBodies = { constrainedBody };
                std::vector< Eigen::VectorXd > constrainedParametersValues;


                bool constraintConfirmed = true;
                for ( unsigned int i = 1 ; i < nbParametersInConstraint ; i++ )
                {
                    // Retrieve list of bodies detected for other parameter
                    std::vector< std::pair< std::string, std::string > > constrainedBodiesOtherParam;
                    if ( indicesPossiblyConstrainedParameters.count( parametersInConstraint[ i ] ) > 0 )
                    {
                        constrainedBodiesOtherParam = indicesPossiblyConstrainedParameters.at( parametersInConstraint[ i ] ).first;
                    }

                    auto it = find( constrainedBodiesOtherParam.begin( ), constrainedBodiesOtherParam.end( ), constrainedBody );
                    if (it != constrainedBodiesOtherParam.end( ) )
                    {
                        int index = it - constrainedBodiesOtherParam.begin( );
                        constraintIndices.push_back( indicesPossiblyConstrainedParameters.at( parametersInConstraint[ i ] ).second.at( index ) );
                        constrainedBodies.push_back( constrainedBody );
                    }
                    else
                    {
                        constraintConfirmed = false;
                    }
                }

                if ( constraintConfirmed )
                {
                    // PRINTING STATEMENTS
                    std::cout << "constraintIndices: " << "\n\n";
                    for ( unsigned int k = 0 ; k < constraintIndices.size( ) ; k++ )
                    {
                        std::cout << constraintIndices.at( k ).first << " (" << constraintIndices.at(k).second << ")" << " - ";
                    }
                    std::cout << "\n\n";
                    // END PRINTING STATEMENTS


                    // Create parameter constraint object
                    switch( constraint )
                    {
                        case tidal_quality_factor_single_love_number_constraint:
                        case tidal_quality_factor_full_love_number_constraint:
                        {
                            std::cout << "common check for tidal_quality_factor_single_love_number_constraint and "
                                         "tidal_quality_factor_full_love_number_constraint" << "\n\n";

                            std::shared_ptr< TidalLoveNumber< Eigen::VectorXd > > loveNumber = std::dynamic_pointer_cast< TidalLoveNumber< Eigen::VectorXd > >(
                                    parameters->getVectorParameters( ).at( constraintIndices.at( 0 ).first ) );
                            std::shared_ptr< InverseTidalQualityFactor > invQ = std::dynamic_pointer_cast< InverseTidalQualityFactor >(
                                    parameters->getDoubleParameters( ).at( constraintIndices.at( 1 ).first ) );

                            Eigen::VectorXd k2Value = loveNumber->getParameterValue( );
                            std::cout << "k2Value: " << k2Value.transpose( ) << "\n\n";

                            // If constrained parameters include single love numbers, check orders to be estimated.
                            if ( constraint == tidal_quality_factor_single_love_number_constraint )
                            {
                                std::cout << "special check for tidal_quality_factor_single_love_number_constraint" << "\n\n";

                                std::cout << "before call to constraintSingleOrderLoveNumber, k2 value: " << k2Value.transpose( ) << "\n\n";
                                constraintSingleOrderLoveNumber( constraint, loveNumber, constraintIndices, k2Value );
                                std::cout << "after call to constraintSingleOrderLoveNumber - constraint indices: " << "\n\n";
                                for ( unsigned int k = 0 ; k < constraintIndices.size( ) ; k++ )
                                {
                                    std::cout << constraintIndices.at( k ).first << " (" << constraintIndices.at(k).second << ")" << " - ";
                                }
                                std::cout << "\n\n";

                                std::cout << "after call to constraintSingleOrderLoveNumber, k2 value: " << k2Value.transpose( ) << "\n\n";

//                                std::vector< int > orders = loveNumber->getOrders( );
//                                if ( orders.size( ) > 1 )
//                                {
//                                    if ( find( orders.begin( ), orders.end( ), 2 ) == orders.end( ) )
//                                    {
//                                        throw std::runtime_error( "Error for constraint of type " + std::to_string( constraint ) +
//                                                                  ", multiple orders for k2. k22 should be used for the constraint, but order 2 is not estimated." );
//                                    }
//                                    else
//                                    {
//                                        std::cerr << "Warning, for  constraint of type " << constraint << ", multiple orders for k2 "
//                                                      "are considered, but only k22 will be included in the constraint." << "\n\n";
//                                    }
//                                    constraintIndices.clear( );
//                                    k2Value.resize( 2 );
//
//                                    unsigned int indexParam = constraintIndices.at( 0 ).first;
//                                    unsigned int indexK2 = 0;
//                                    for ( unsigned int i = 0 ; i < orders.size( ) ; i++ )
//                                    {
//                                        if ( orders.at( i ) == 2 )
//                                        {
//                                            constraintIndices.push_back( std::make_pair( indexParam, 2 ) );
//                                            k2Value = k2Value.segment( indexK2, 2 );
//                                        }
//                                        indexParam += 2;
//                                        indexK2 += 2;
//                                    }
//
//                                    // Add other tidal parameter indices
//                                    constraintIndices.push_back( constraintIndices.at( 1 ) );
//                                }
                            }

                            double invQValue = invQ->getParameterValue( );
                            double theoreticalInvQValue = k2Value[ 1 ] / std::sqrt( k2Value[ 0 ] * k2Value[ 0 ] + k2Value[ 1 ] * k2Value[ 1 ] );
                            std::cout << "invQValue: " << invQValue << "\n\n";
                            std::cout << "theoreticalInvQValue: " << theoreticalInvQValue << "\n\n";

                            if ( theoreticalInvQValue != invQValue )
                            {
                                invQ->setParameterValue( theoreticalInvQValue );
                                std::string warning = "Warning when creating constraint " + std::to_string( constraint ) + ", the values for k2 and invQ for body "
                                        + loveNumber->getParameterName( ).second.first;
                                if ( !loveNumber->getParameterName( ).second.second.empty( ) )
                                {
                                    warning += " (tides caused by body " + loveNumber->getParameterName( ).second.second + ")";
                                }
                                warning += " are inconsistent. InvQ overwritten.";
                                std::cerr << warning << "\n\n";
                            }
                            std::cout << "updated invQ value: " << invQ->getParameterValue( ) << "\n\n";

                            parameterConstraints.push_back( std::make_shared< TidalQualityFactorLoveNumberConstraints >(
                                    constraint, constrainedBodies, constraintIndices,
                                    parameters->template getFullParameterValues< InitialStateParameterType >( ) ) );



//                            unsigned int indexLoveNumberParam = constraintIndices.at( 0 ).first;
//                            unsigned int indexInvQParam = constraintIndices.at( 1 ).first;
//
//                            std::shared_ptr< TidalLoveNumber< Eigen::VectorXd > > loveNumber = std::dynamic_pointer_cast< TidalLoveNumber< Eigen::VectorXd > >(
//                                    parameters->getVectorParameters( ).at( indexLoveNumberParam ) );
//                            Eigen::VectorXd k2Value = loveNumber->getParameterValue( );
//                            std::cout << "k2Value: " << k2Value.transpose( ) << "\n\n";
//
//                            std::vector< int > orders = loveNumber->getOrders( );
//
//                            std::shared_ptr< InverseTidalQualityFactor > invQ = std::dynamic_pointer_cast< InverseTidalQualityFactor >(
//                                    parameters->getDoubleParameters( ).at( indexInvQParam ) );
//                            double invQValue = invQ->getParameterValue( );
//
//                            std::vector< std::pair< unsigned int, unsigned int > > refinedConstraintIndices;
//                            Eigen::Vector2d effectiveK2Value = Eigen::Vector2d::Zero( );
//                            if ( orders.size( ) == 1 )
//                            {
//                                refinedConstraintIndices = constraintIndices;
//                                effectiveK2Value = k2Value.segment( 0, 2 );
//                            }
//                            else
//                            {
//                                std::cerr << "Warning, for  constraint of type " << constraint << ", multiple orders for k2 "
//                                                                                                  "are considered, but only k22 will be included in the constraint." << "\n\n";
//                                unsigned int indexParam = constraintIndices.at( 0 ).first;
//                                unsigned int indexK2 = 0;
//                                for ( unsigned int i = 0 ; i < orders.size( ) ; i++ )
//                                {
//                                    if ( orders.at( i ) == 2 )
//                                    {
//                                        refinedConstraintIndices.push_back( std::make_pair( indexParam, 2 ) );
//                                        effectiveK2Value = k2Value.segment( indexK2, 2 );
//                                    }
//                                    indexParam += 2;
//                                    indexK2 += 2;
//                                }
//
//                                // Add invQ indices
//                                refinedConstraintIndices.push_back( constraintIndices.at( 1 ) );
//                            }
//
//
//
//                            std::cout << "effectiveK2Value: " << effectiveK2Value.transpose( ) << "\n\n";
//                            double theoreticalInvQValue = effectiveK2Value[ 1 ]; // /
////                            std::sqrt( effectiveK2Value[ 0 ] * effectiveK2Value[ 0 ] + effectiveK2Value[ 1 ] * effectiveK2Value[ 1 ] );
//                            std::cout << "invQValue: " << invQValue << "\n\n";
//                            std::cout << "theoreticalInvQValue: " << theoreticalInvQValue << "\n\n";
//
//                            if ( theoreticalInvQValue != invQValue )
//                            {
//                                invQ->setParameterValue( theoreticalInvQValue );
//                                std::cerr << "Warning when creating tidal_quality_factor_single_love_number_constraint, the values for k2 "
//                                             "and invQ for body " << loveNumber->getParameterName( ).second.first << " (tides caused by body " <<
//                                          loveNumber->getParameterName( ).second.second << ") are inconsistent. InvQ overwritten." << "\n\n";
//                            }
//                            std::cout << "updated invQ value: " << invQ->getParameterValue( ) << "\n\n";
//
//                            parameterConstraints.push_back( std::make_shared< TidalQualityFactorLoveNumberConstraints >(
//                                    constraint, parameterId, refinedConstraintIndices,
//                                    parameters->template getFullParameterValues< InitialStateParameterType >( ) ) );

                            break;
                        }
                        case tidal_time_lag_single_love_number_constraint:
                        case tidal_time_lag_full_love_number_constraint:
                        {
                            std::cout << "common check for tidal_time_lag_single_love_number_constraint and "
                                         "tidal_time_lag_full_love_number_constraint" << "\n\n";

                            std::shared_ptr< TidalLoveNumber< Eigen::VectorXd > > loveNumber = std::dynamic_pointer_cast< TidalLoveNumber< Eigen::VectorXd > >(
                                    parameters->getVectorParameters( ).at( constraintIndices.at( 0 ).first ) );
                            std::shared_ptr< DirectTidalTimeLag > timeLag = std::dynamic_pointer_cast< DirectTidalTimeLag >(
                                    parameters->getDoubleParameters( ).at( constraintIndices.at( 1 ).first ) );

                            Eigen::VectorXd k2Value = loveNumber->getParameterValue( );
                            std::cout << "k2Value: " << k2Value.transpose( ) << "\n\n";

                            // If constrained parameters include single love numbers, check orders to be estimated.
                            if ( constraint == tidal_time_lag_single_love_number_constraint )
                            {
                                std::cout << "special check for tidal_time_lag_single_love_number_constraint" << "\n\n";

                                std::cout << "before call to constraintSingleOrderLoveNumber, k2 value: " << k2Value.transpose( ) << "\n\n";
                                constraintSingleOrderLoveNumber( constraint, loveNumber, constraintIndices, k2Value );
                                std::cout << "after call to constraintSingleOrderLoveNumber - constraint indices: " << "\n\n";
                                for ( unsigned int k = 0 ; k < constraintIndices.size( ) ; k++ )
                                {
                                    std::cout << constraintIndices.at( k ).first << " (" << constraintIndices.at(k).second << ")" << " - ";
                                }
                                std::cout << "\n\n";

                                std::cout << "after call to constraintSingleOrderLoveNumber, k2 value: " << k2Value.transpose( ) << "\n\n";
                            }

                            double timeLagValue = timeLag->getParameterValue( );
                            double theoreticalInvQValue = k2Value[ 1 ] /
                                    std::sqrt( k2Value[ 0 ] * k2Value[ 0 ] + k2Value[ 1 ] * k2Value[ 1 ] );
                            std::cout << "theoreticalInvQValue: " << theoreticalInvQValue << "\n\n";
                            std::cout << "timeLag->getTidalPeriod( ): " <<  timeLag->getTidalPeriod( ) << "\n\n";
                            double theoreticalTimeLagValue = timeLag->getTidalPeriod( ) * std::atan( theoreticalInvQValue ) / ( 2.0 * mathematical_constants::PI );
                            std::cout << "timeLagValue: " << timeLagValue << "\n\n";
                            std::cout << "theoreticalTimeLagValue: " << theoreticalTimeLagValue << "\n\n";

                            if ( theoreticalTimeLagValue != timeLagValue )
                            {
                                timeLag->setParameterValue( theoreticalTimeLagValue );
                                std::string warning = "Warning when creating constraint " + std::to_string( constraint ) + ", the values for k2 and time lag for body "
                                                      + loveNumber->getParameterName( ).second.first;
                                if ( !loveNumber->getParameterName( ).second.second.empty( ) )
                                {
                                    warning += " (tides caused by body " + loveNumber->getParameterName( ).second.second + ")";
                                }
                                warning += " are inconsistent. Time lag overwritten.";
                                std::cerr << warning << "\n\n";
                            }
                            std::cout << "updated time lag value: " << timeLag->getParameterValue( ) << "\n\n";

                            parameterConstraints.push_back( std::make_shared< TidalTimeLagLoveNumberConstraints >(
                                    constraint, constrainedBodies, constraintIndices,
                                    parameters->template getFullParameterValues< InitialStateParameterType >( ), timeLag->getTidalPeriod( ) ) );
                            break;
                        }
                        default:
                            break;
                    }
                }
            }
            // END TEST





//            for ( unsigned int i = 1 ; i < nbParametersInConstraint ; i++ )
//            {
//                std::vector< std::pair< std::string, std::string > > otherParameter;
//                if ( indicesPossiblyConstrainedParameters.count( parametersInConstraint[ i ] ) > 0 )
//                {
//                    otherParameter = indicesPossiblyConstrainedParameters.at( parametersInConstraint[ i ] ).first;
//                }
//
//                std::vector< int > parametersToRemove;
//                for ( unsigned int j = 0 ; j < detectedBodies.size( ) ; j++ )
//                {
//                    // Verify if element initially in common parameters is absent from other parameter involved in the constraint
//                    if ( count( otherParameter.begin( ), otherParameter.end( ), detectedBodies[ j ] ) == 0 )
//                    {
//                        parametersToRemove.push_back( j );
//                    }
//                }
//                for ( unsigned int j = 0 ; j < parametersToRemove.size( ) ; j++ )
//                {
//                    detectedBodies.erase( detectedBodies.begin( ) + parametersToRemove.at( j ) - j); // remove element from common parameters
//                }
//            }
//
//            // PRINTING STATEMENTS
//            std::cout << "common parameters left: " << "\n\n";
//            for ( unsigned int i = 0 ; i < detectedBodies.size( ) ; i++ )
//            {
//                std::cout << detectedBodies[ i ].first << " - " << detectedBodies[ i ].second << "\n\n";
//            }
//            std::cout << "end parsing common parameters." << "\n\n";
//            // END PRINTING STATEMENTS
//
//            // Retrieve indices of parameters in constraint
//            for ( auto parameterId : detectedBodies )
//            {
//                std::vector< std::pair< unsigned int, unsigned int > > constraintIndices;
//                for ( auto param : parametersInConstraint )
//                {
//                    std::pair< std::vector< std::pair< std::string, std::string > >, std::vector< std::pair< int, int > > > indices =
//                            indicesPossiblyConstrainedParameters.at( param );
//                    for ( unsigned int j = 0 ; j < indices.first.size( ) ; j++ )
//                    {
//                        if ( indices.first[ j ] == parameterId )
//                        {
//                            constraintIndices.push_back( indices.second.at( j ) );
//                        }
//                    }
//                }
//
//                // PRINTING STATEMENTS
//                std::cout << "constraintIndices: " << "\n\n";
//                for ( unsigned int k = 0 ; k < constraintIndices.size( ) ; k++ )
//                {
//                    std::cout << constraintIndices.at( k ).first << " (" << constraintIndices.at(k).second << ")" << " - ";
//                }
//                std::cout << "\n\n";
//                // END PRINTING STATEMENTS
//
//
//                // Create parameter constraint object
//                switch( constraint )
//                {
//                    case tidal_quality_factor_single_love_number_constraint:
//                    {
//                        unsigned int indexLoveNumberParam = constraintIndices.at( 0 ).first;
//                        unsigned int indexInvQParam = constraintIndices.at( 1 ).first;
//
//                        std::shared_ptr< TidalLoveNumber< Eigen::VectorXd > > loveNumber = std::dynamic_pointer_cast< TidalLoveNumber< Eigen::VectorXd > >(
//                                parameters->getVectorParameters( ).at( indexLoveNumberParam ) );
//                        std::vector< int > orders = loveNumber->getOrders( );
//                        Eigen::VectorXd k2Value = loveNumber->getParameterValue( );
//                        std::cout << "k2Value: " << k2Value.transpose( ) << "\n\n";
//
//                        std::shared_ptr< InverseTidalQualityFactor > invQ = std::dynamic_pointer_cast< InverseTidalQualityFactor >(
//                                parameters->getDoubleParameters( ).at( indexInvQParam ) );
//                        double invQValue = invQ->getParameterValue( );
//
//                        std::vector< std::pair< unsigned int, unsigned int > > refinedConstraintIndices;
//                        Eigen::Vector2d effectiveK2Value = Eigen::Vector2d::Zero( );
//                        if ( orders.size( ) == 1 )
//                        {
//                            refinedConstraintIndices = constraintIndices;
//                            effectiveK2Value = k2Value.segment( 0, 2 );
//                        }
//                        else
//                        {
//                            std::cerr << "Warning, for  constraint of type " << constraint << ", multiple orders for k2 "
//                                         "are considered, but only k22 will be included in the constraint." << "\n\n";
//                            unsigned int indexParam = constraintIndices.at( 0 ).first;
//                            unsigned int indexK2 = 0;
//                            for ( unsigned int i = 0 ; i < orders.size( ) ; i++ )
//                            {
//                                if ( orders.at( i ) == 2 )
//                                {
//                                    refinedConstraintIndices.push_back( std::make_pair( indexParam, 2 ) );
//                                    effectiveK2Value = k2Value.segment( indexK2, 2 );
//                                }
//                                indexParam += 2;
//                                indexK2 += 2;
//                            }
//
//                            // Add invQ indices
//                            refinedConstraintIndices.push_back( constraintIndices.at( 1 ) );
//                        }
//
//
//
//                        std::cout << "effectiveK2Value: " << effectiveK2Value.transpose( ) << "\n\n";
//                        double theoreticalInvQValue = effectiveK2Value[ 1 ]; // /
////                            std::sqrt( effectiveK2Value[ 0 ] * effectiveK2Value[ 0 ] + effectiveK2Value[ 1 ] * effectiveK2Value[ 1 ] );
//                        std::cout << "invQValue: " << invQValue << "\n\n";
//                        std::cout << "theoreticalInvQValue: " << theoreticalInvQValue << "\n\n";
//
//                        if ( theoreticalInvQValue != invQValue )
//                        {
//                            invQ->setParameterValue( theoreticalInvQValue );
//                            std::cerr << "Warning when creating tidal_quality_factor_single_love_number_constraint, the values for k2 "
//                                         "and invQ for body " << loveNumber->getParameterName( ).second.first << " (tides caused by body " <<
//                                      loveNumber->getParameterName( ).second.second << ") are inconsistent. InvQ overwritten." << "\n\n";
//                        }
//                        std::cout << "updated invQ value: " << invQ->getParameterValue( ) << "\n\n";
//
//                        parameterConstraints.push_back( std::make_shared< TidalQualityFactorLoveNumberConstraints >(
//                                tidal_quality_factor_single_love_number_constraint, parameterId, refinedConstraintIndices,
//                                parameters->template getFullParameterValues< InitialStateParameterType >( ) ) );
//                        break;
//                    }
//                    case tidal_time_lag_single_love_number_constraint:
//                    {
//                        unsigned int indexLoveNumberParam = constraintIndices.at( 0 ).first;
//                        unsigned int indexTimeLag = constraintIndices.at( 1 ).first;
//
//                        std::shared_ptr< TidalLoveNumber< Eigen::VectorXd > > loveNumber = std::dynamic_pointer_cast< TidalLoveNumber< Eigen::VectorXd > >(
//                                parameters->getVectorParameters( ).at( indexLoveNumberParam ) );
//                        std::vector< int > orders = loveNumber->getOrders( );
//                        Eigen::VectorXd k2Value = loveNumber->getParameterValue( );
//                        std::cout << "k2Value: " << k2Value.transpose( ) << "\n\n";
//
//                        std::shared_ptr< DirectTidalTimeLag > timeLag = std::dynamic_pointer_cast< DirectTidalTimeLag >(
//                                parameters->getDoubleParameters( ).at( indexTimeLag ) );
//                        double timeLagValue = timeLag->getParameterValue( );
//
//                        std::vector< std::pair< unsigned int, unsigned int > > refinedConstraintIndices;
//
//                        Eigen::Vector2d effectiveK2Value = Eigen::Vector2d::Zero( );
//                        if ( orders.size( ) == 1 )
//                        {
//                            refinedConstraintIndices = constraintIndices;
//                            effectiveK2Value = k2Value.segment( 0, 2 );
//                        }
//                        else
//                        {
//                            unsigned int indexParam = constraintIndices.at( 0 ).first;
//                            unsigned int indexK2 = 0;
//                            for ( unsigned int i = 0 ; i < orders.size( ) ; i++ )
//                            {
//                                if ( orders.at( i ) == 0 || orders.at( i ) == 2 )
//                                {
//                                    refinedConstraintIndices.push_back( std::make_pair( indexParam, 2 ) );
//                                    effectiveK2Value = k2Value.segment( indexK2, 2 );
//                                }
//                                indexParam += 2;
//                                indexK2 += 2;
//                            }
//
//                            // Add time lag indices
//                            refinedConstraintIndices.push_back( constraintIndices.at( 1 ) );
//                        }
//
//
//
//                        std::cout << "effectiveK2Value: " << effectiveK2Value.transpose( ) << "\n\n";
//                        double theoreticalInvQValue = effectiveK2Value[ 1 ] /
//                                                      std::sqrt( effectiveK2Value[ 0 ] * effectiveK2Value[ 0 ] + effectiveK2Value[ 1 ] * effectiveK2Value[ 1 ] );
//                        double theoreticalTimeLagValue = timeLag->getTidalPeriod( ) * std::atan( theoreticalInvQValue ) / ( 2.0 * mathematical_constants::PI );
//                        std::cout << "timeLagValue: " << timeLagValue << "\n\n";
//                        std::cout << "theoreticalTimeLagValue: " << theoreticalTimeLagValue << "\n\n";
//
//                        if ( theoreticalTimeLagValue != timeLagValue )
//                        {
//                            timeLag->setParameterValue( theoreticalTimeLagValue );
//                            std::cerr << "Warning when creating tidal_time_lag_single_love_number_constraint, the values for k2 "
//                                         "and time lag for body " << loveNumber->getParameterName( ).second.first << " (tides caused by body " <<
//                                      loveNumber->getParameterName( ).second.second << ") are inconsistent. Time lag overwritten." << "\n\n";
//                        }
//                        std::cout << "updated time lag value: " << timeLag->getParameterValue( ) << "\n\n";
//
//                        parameterConstraints.push_back( std::make_shared< TidalTimeLagLoveNumberConstraints >(
//                                tidal_time_lag_single_love_number_constraint, parameterId, refinedConstraintIndices,
//                                parameters->template getFullParameterValues< InitialStateParameterType >( ) ) );
//                        break;
//                    }
//                    case tidal_quality_factor_full_love_number_constraint:
//                    {
//                        unsigned int indexLoveNumberParam = constraintIndices.at( 0 ).first;
//                        unsigned int indexInvQParam = constraintIndices.at( 1 ).first;
//
//                        std::shared_ptr< TidalLoveNumber< Eigen::VectorXd > > loveNumber = std::dynamic_pointer_cast< TidalLoveNumber< Eigen::VectorXd > >(
//                                parameters->getVectorParameters( ).at( indexLoveNumberParam ) );
//                        Eigen::VectorXd k2Value = loveNumber->getParameterValue( );
//                        std::cout << "k2Value: " << k2Value.transpose( ) << "\n\n";
//
//                        std::shared_ptr< InverseTidalQualityFactor > invQ = std::dynamic_pointer_cast< InverseTidalQualityFactor >(
//                                parameters->getDoubleParameters( ).at( indexInvQParam ) );
//                        double invQValue = invQ->getParameterValue( );
//                        double theoreticalInvQValue = k2Value[ 1 ] / std::sqrt( k2Value[ 0 ] * k2Value[ 0 ] + k2Value[ 1 ] * k2Value[ 1 ] );
//                        std::cout << "invQValue: " << invQValue << "\n\n";
//                        std::cout << "theoreticalInvQValue: " << theoreticalInvQValue << "\n\n";
//
//                        if ( theoreticalInvQValue != invQValue )
//                        {
//                            invQ->setParameterValue( theoreticalInvQValue );
//                            std::cerr << "Warning when creating tidal_quality_factor_full_love_number_constraint, the values for k2"
//                                         " and invQ for body " << loveNumber->getParameterName( ).second.first << " (tides caused by body " <<
//                                      loveNumber->getParameterName( ).second.second << ") are inconsistent. InvQ overwritten." << "\n\n";
//                        }
//                        std::cout << "updated invQ value: " << invQ->getParameterValue( ) << "\n\n";
//
//                        parameterConstraints.push_back( std::make_shared< TidalQualityFactorLoveNumberConstraints >(
//                                tidal_quality_factor_full_love_number_constraint, parameterId, constraintIndices,
//                                parameters->template getFullParameterValues< InitialStateParameterType >( ) ) );
//                        break;
//                    }
//                    case tidal_time_lag_full_love_number_constraint:
//                    {
//                        unsigned int indexLoveNumberParam = constraintIndices.at( 0 ).first;
//                        unsigned int indexTimeLag = constraintIndices.at( 1 ).first;
//
//                        std::shared_ptr< TidalLoveNumber< Eigen::VectorXd > > loveNumber = std::dynamic_pointer_cast< TidalLoveNumber< Eigen::VectorXd > >(
//                                parameters->getVectorParameters( ).at( indexLoveNumberParam ) );
//                        Eigen::VectorXd k2Value = loveNumber->getParameterValue( );
//                        std::cout << "k2Value: " << k2Value.transpose( ) << "\n\n";
//
//                        std::shared_ptr< DirectTidalTimeLag > timeLag = std::dynamic_pointer_cast< DirectTidalTimeLag >(
//                                parameters->getDoubleParameters( ).at( indexTimeLag ) );
//                        double timeLagValue = timeLag->getParameterValue( );
//                        double theoreticalInvQValue = k2Value[ 1 ] / std::sqrt( k2Value[ 0 ] * k2Value[ 0 ] + k2Value[ 1 ] * k2Value[ 1 ] );
//                        double theoreticalTimeLagValue = timeLag->getTidalPeriod( ) * std::atan( theoreticalInvQValue ) / ( 2.0 * mathematical_constants::PI );
//                        std::cout << "timeLagValue: " << timeLagValue << "\n\n";
//                        std::cout << "theoreticalTimeLagValue: " << theoreticalTimeLagValue << "\n\n";
//
//                        if ( theoreticalTimeLagValue != timeLagValue )
//                        {
//                            timeLag->setParameterValue( theoreticalTimeLagValue );
//                            std::cerr << "Warning when creating tidal_time_lag_full_love_number_constraint, the values for k2 "
//                                         "and time lag for body " << loveNumber->getParameterName( ).second.first << " (tides caused by body " <<
//                                      loveNumber->getParameterName( ).second.second << ") are inconsistent. Time lag overwritten." << "\n\n";
//                        }
//                        std::cout << "updated time lag value: " << timeLag->getParameterValue( ) << "\n\n";
//
//                        parameterConstraints.push_back( std::make_shared< TidalTimeLagLoveNumberConstraints >(
//                                tidal_time_lag_full_love_number_constraint, parameterId, constraintIndices,
//                                parameters->template getFullParameterValues< InitialStateParameterType >( ) ) );
//                        break;
//                    }
//
//                    default:
//                        break;
//                }
//
//            }
        }

        // Add constraints to estimatable parameters set
        parameters->resetConstraints( parameterConstraints );
    }

}


template< typename InitialStateParameterType = double >
std::vector< std::map< estimatable_parameters::EstimatebleParametersEnum, std::pair<
        std::pair< std::string, std::string >, std::pair< unsigned int, unsigned int > > > > getActiveCustomConstraints(
                std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parameters,
                const std::vector< estimatable_parameters::EstimatebleParametersEnum >& parametersInvolved,
                std::vector< std::pair< std::string, std::string > > bodiesInvolved,
                const bool applyPerBody,
                const bool applyToAllBodies )
{
    // Check that the number of parameters involved is >= 2
    if ( parametersInvolved.size( ) < 2 )
    {
        throw std::runtime_error( "Error when creating custom constraint, at least two parameters should be involved." );
    }

    if ( applyToAllBodies && !applyPerBody )
    {
        throw std::runtime_error( "Error when creating custom constraint, inconsistent input. When applyToAllBodies is true, the "
                                  "constraint should be applied between parameters of a same body (i.e., applyPerBody should be true). " );
    }
    if ( ( bodiesInvolved.size( ) != 0 && applyToAllBodies ) || ( bodiesInvolved.size( ) == 0 && !applyToAllBodies ) )
    {
        throw std::runtime_error( "Error when creating custom constraint, inconsistent input. If no body is specified, the constraint should be applied "
                                  "to all bodies (boolean applyToAllBodies to true). Conversely, if applyToAllBodies is true, then no body needs to be "
                                  "specified." );
    }
    if ( !applyPerBody && ( bodiesInvolved.size( ) != parametersInvolved.size( ) ) )
    {
        throw std::runtime_error( "Error when creating custom constraint, inconsistent input. If the constraint involves parameters related "
                                  "to different bodies, the sizes (and orders) of bodiesInvolved and parametersInvolved should match." );
    }

    // Check if relevant parameters are estimated for the requested bodies
    std::vector< std::map< estimatable_parameters::EstimatebleParametersEnum, std::pair<
            std::pair< std::string, std::string >, std::pair< unsigned int, unsigned int > > > > activeConstraintsMetadata;

    if ( applyPerBody ) // if constraint has to be applied between parameters related to be the same body
    {
        if ( applyToAllBodies )
        {
            // Initialise bodies involved in constraint by listing all bodies for which the first constraint parameter is estimated.
            bodiesInvolved = parameters->getIndicesForParameterType( parametersInvolved.at( 0 ) ).first;
        }
        for ( auto parameterBodyId : bodiesInvolved )
        {
            bool constraintDetected = true;
            std::map< estimatable_parameters::EstimatebleParametersEnum, std::pair<
                    std::pair< std::string, std::string >, std::pair< unsigned int, unsigned int > > > singleConstraintMetadata;

            for ( auto param : parametersInvolved )
            {
                std::pair< std::vector< std::pair< std::string, std::string > >, std::vector< std::pair< int, int > > > indicesParameter =
                        parameters->getIndicesForParameterType( param );
                auto it = find( indicesParameter.first.begin( ), indicesParameter.first.end( ), parameterBodyId );
                if (it != indicesParameter.first.end( ) )
                {
                    int index = it - indicesParameter.first.begin( );
                    singleConstraintMetadata[ param ] = std::make_pair( parameterBodyId, indicesParameter.second.at( index ) );
                }
                else
                {
                    constraintDetected = false;
                    if ( !applyToAllBodies )
                    {
                        std::cerr << "Warning: custom constraint cannot be created: parameter " << std::to_string( param )
                                  << " is not estimated for body " << parameterBodyId.first <<
                                  " (secondary body " << parameterBodyId.second <<")\n\n";
                    }
                }
            }
            if ( constraintDetected )
            {
                activeConstraintsMetadata.push_back( singleConstraintMetadata );
            }
        }
    }
    else // if constraint has to be applied between parameters related to different bodies
    {
        bool constraintDetected = true;
        std::map< estimatable_parameters::EstimatebleParametersEnum, std::pair<
                std::pair< std::string, std::string >, std::pair< unsigned int, unsigned int > > > singleConstraintMetadata;

        for ( unsigned int k = 0 ; k < parametersInvolved.size( ) ; k++ )
        {
            std::pair< std::vector< std::pair< std::string, std::string > >, std::vector< std::pair< int, int > > > indicesParameter =
                    parameters->getIndicesForParameterType( parametersInvolved.at( k ) );

            auto it = find( indicesParameter.first.begin( ), indicesParameter.first.end( ), bodiesInvolved.at( k ) );
            if (it != indicesParameter.first.end( ) )
            {
                int index = it - indicesParameter.first.begin( );
                singleConstraintMetadata[ parametersInvolved.at( k ) ] = std::make_pair( bodiesInvolved.at( k ), indicesParameter.second.at( index ) );
            }
            else
            {
                constraintDetected = false;
                std::cerr << "Warning: custom constraint cannot be created: parameter " << std::to_string( parametersInvolved.at( k ) )
                          << " is not estimated for body " << bodiesInvolved.at( k ).first <<
                          " (secondary body " << bodiesInvolved.at( k ).second <<")\n\n";
            }
        }
        if ( constraintDetected )
        {
            activeConstraintsMetadata.push_back( singleConstraintMetadata );
        }
    }


    // PRINTING STATEMENTS
    std::cout << "activeConstraintsMetadata: " << "\n\n";
    for ( auto constraintMetadata : activeConstraintsMetadata )
    {
        std::cout << "detected custom constraint between \n\n";
        for ( auto itr : constraintMetadata )
        {
            std::cout << itr.first << " for " << itr.second.first.first << " - " << itr.second.first.second <<
                      " | indices: " << itr.second.second.first << " - " << itr.second.second.second << "\n\n";
        }
    }
    // END PRINTING STATEMENTS

    return activeConstraintsMetadata;
}


template< typename InitialStateParameterType = double >
void createCustomParameterConstraints(
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parameters,
        const std::vector< std::function< Eigen::MatrixXd( Eigen::VectorXd ) > >& linearisedFactorsConstraint,
        const std::vector< estimatable_parameters::EstimatebleParametersEnum >& parametersInvolved,
        std::vector< std::pair< std::string, std::string > > bodiesInvolved = { },
        const bool applyPerBody = false,
        const bool applyToAllBodies = false )

{
    using namespace tudat::estimatable_parameters;
    std::cout << "TEST, in createCustomParameterConstraint function" << "\n\n";

    // Check if the number of parameters involved in the constraint matches with the number of linearisation factors provided as input
    if ( linearisedFactorsConstraint.size ( ) != parametersInvolved.size( ) )
    {
        throw std::runtime_error( "Number of linearised factors for custom constraint should be consistent with number of "
                                  "parameters involved in said constraint." );
    }

    // Check if relevant parameters are estimated for the requested bodies
    std::vector< std::map< EstimatebleParametersEnum, std::pair<
            std::pair< std::string, std::string >, std::pair< unsigned int, unsigned int > > > > activeConstraintsMetadata =
                    getActiveCustomConstraints( parameters, parametersInvolved, bodiesInvolved, applyPerBody, applyToAllBodies );

    // Define custom constraints
    std::vector< std::shared_ptr< estimatable_parameters::ParameterConstraints > > customConstraints;

    for ( auto constraint : activeConstraintsMetadata )
    {
        std::vector< std::pair< std::string, std::string > > constraintBodies;
        std::vector< std::pair< unsigned int, unsigned int > > constraintIndices;
        std::vector< std::string > involvedParameters;
        for ( auto itr : constraint )
        {
            involvedParameters.push_back( getParameterTypeString( itr.first ) );
            constraintBodies.push_back( itr.second.first );
            constraintIndices.push_back( itr.second.second );
        }
        customConstraints.push_back( std::make_shared< CustomConstraint >(
                involvedParameters, constraintBodies, linearisedFactorsConstraint, constraintIndices, parameters->template getFullParameterValues< InitialStateParameterType >( ) ) );
    }

    parameters->addConstraints( customConstraints );

}


template< typename InitialStateParameterType = double >
void createConstantCustomParameterConstraints(
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parameters,
        const std::vector< Eigen::MatrixXd >& linearisedFactorsConstraint,
        const std::vector< estimatable_parameters::EstimatebleParametersEnum >& parametersInvolved,
        std::vector< std::pair< std::string, std::string > > bodiesInvolved = { },
        const bool applyPerBody = false,
        const bool applyToAllBodies = false )

{
    using namespace tudat::estimatable_parameters;
    std::cout << "TEST, in createConstantCustomParameterConstraints function" << "\n\n";

    // Check if the number of parameters involved in the constraint matches with the number of linearisation factors provided as input
    if ( linearisedFactorsConstraint.size ( ) != parametersInvolved.size( ) )
    {
        throw std::runtime_error( "Number of linearised factors for custom constraint should be consistent with number of "
                                  "parameters involved in said constraint." );
    }

    // Check if relevant parameters are estimated for the requested bodies
    std::vector< std::map< EstimatebleParametersEnum, std::pair<
            std::pair< std::string, std::string >, std::pair< unsigned int, unsigned int > > > > activeConstraintsMetadata =
            getActiveCustomConstraints( parameters, parametersInvolved, bodiesInvolved, applyPerBody, applyToAllBodies );

    // Define custom constraints
    std::vector< std::shared_ptr< estimatable_parameters::ParameterConstraints > > customConstraints;

    for ( auto constraint : activeConstraintsMetadata )
    {
        std::vector< std::pair< std::string, std::string > > constraintBodies;
        std::vector< std::pair< unsigned int, unsigned int > > constraintIndices;
        std::vector< std::string > involvedParameters;

        for ( auto itr : constraint )
        {
            involvedParameters.push_back( getParameterTypeString( itr.first ) );
            constraintBodies.push_back( itr.second.first );
            constraintIndices.push_back( itr.second.second );
        }
        customConstraints.push_back( std::make_shared< CustomConstantConstraint >(
                involvedParameters, constraintBodies, linearisedFactorsConstraint, constraintIndices, parameters->template getFullParameterValues< InitialStateParameterType >( ) ) );
    }

    parameters->addConstraints( customConstraints );

}


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEPARAMETERCONSTRAINTS_H
