/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/estimation_setup/createParameterConstraints.h"

namespace tudat
{

namespace simulation_setup
{

std::vector< estimatable_parameters::EstimatebleParametersEnum > getParametersTypesPerConstraint(
        const estimatable_parameters::ParameterConstraintsEnum constraintType )
{
    using namespace estimatable_parameters;

    std::vector< EstimatebleParametersEnum > parametersList;

    switch( constraintType )
    {
        case tidal_quality_factor_single_love_number_constraint:
            parametersList = { single_degree_variable_tidal_love_number, inverse_tidal_quality_factor };
            break;
        case tidal_time_lag_single_love_number_constraint:
            parametersList = { single_degree_variable_tidal_love_number, direct_dissipation_tidal_time_lag };
            break;
        case tidal_quality_factor_full_love_number_constraint:
            parametersList = { full_degree_tidal_love_number, inverse_tidal_quality_factor };
            break;
        case tidal_time_lag_full_love_number_constraint:
            parametersList = { full_degree_tidal_love_number, direct_dissipation_tidal_time_lag };
            break;
        default:
            break;
    }

    return parametersList;
}


//template< typename InitialStateParameterType >
//std::set< estimatable_parameters::ParameterConstraintsEnum > getPotentialConstraints(
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parameters )
//{
//    using namespace estimatable_parameters;
//
//    std::set< EstimatebleParametersEnum > listParameterTypes = parameters->getListParameterTypes( );
//
//    std::set< estimatable_parameters::ParameterConstraintsEnum > relevantConstraints;
//
//    if ( listParameterTypes.count( single_degree_variable_tidal_love_number ) > 0 &&
//         listParameterTypes.count( inverse_tidal_quality_factor ) > 0 )
//    {
//        relevantConstraints.insert( tidal_quality_factor_single_love_number_constraint );
//    }
//    if ( listParameterTypes.count( single_degree_variable_tidal_love_number ) > 0 &&
//         listParameterTypes.count( direct_dissipation_tidal_time_lag ) > 0 )
//    {
//        relevantConstraints.insert( tidal_time_lag_single_love_number_constraint );
//    }
//    if ( listParameterTypes.count( full_degree_tidal_love_number ) > 0 &&
//         listParameterTypes.count( inverse_tidal_quality_factor ) > 0 )
//    {
//        relevantConstraints.insert( tidal_quality_factor_full_love_number_constraint );
//    }
//    if ( listParameterTypes.count( full_degree_tidal_love_number ) > 0 &&
//         listParameterTypes.count( direct_dissipation_tidal_time_lag ) > 0 )
//    {
//        relevantConstraints.insert( tidal_time_lag_full_love_number_constraint );
//    }
//
//    return relevantConstraints;
//}


//template< typename InitialStateParameterType >
//std::vector< std::map< estimatable_parameters::EstimatebleParametersEnum, std::pair<
//        std::pair< std::string, std::string >, std::pair< unsigned int, unsigned int > > > > getActiveCustomConstraints(
//        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parameters,
//        const std::vector< estimatable_parameters::EstimatebleParametersEnum >& parametersInvolved,
//        std::vector< std::pair< std::string, std::string > > bodiesInvolved,
//        const bool applyPerBody,
//        const bool applyToAllBodies )
//{
//    // Check that the number of parameters involved is >= 2
//    if ( parametersInvolved.size( ) < 2 )
//    {
//        throw std::runtime_error( "Error when creating custom constraint, at least two parameters should be involved." );
//    }
//
//    if ( applyToAllBodies && !applyPerBody )
//    {
//        throw std::runtime_error( "Error when creating custom constraint, inconsistent input. When applyToAllBodies is true, the "
//                                  "constraint should be applied between parameters of a same body (i.e., applyPerBody should be true). " );
//    }
//    if ( ( bodiesInvolved.size( ) != 0 && applyToAllBodies ) || ( bodiesInvolved.size( ) == 0 && !applyToAllBodies ) )
//    {
//        throw std::runtime_error( "Error when creating custom constraint, inconsistent input. If no body is specified, the constraint should be applied "
//                                  "to all bodies (boolean applyToAllBodies to true). Conversely, if applyToAllBodies is true, then no body needs to be "
//                                  "specified." );
//    }
//    if ( !applyPerBody && ( bodiesInvolved.size( ) != parametersInvolved.size( ) ) )
//    {
//        throw std::runtime_error( "Error when creating custom constraint, inconsistent input. If the constraint involves parameters related "
//                                  "to different bodies, the sizes (and orders) of bodiesInvolved and parametersInvolved should match." );
//    }
//
//    // Check if relevant parameters are estimated for the requested bodies
//    std::vector< std::map< estimatable_parameters::EstimatebleParametersEnum, std::pair<
//            std::pair< std::string, std::string >, std::pair< unsigned int, unsigned int > > > > activeConstraintsMetadata;
//
//    if ( applyPerBody ) // if constraint has to be applied between parameters related to be the same body
//    {
//        if ( applyToAllBodies )
//        {
//            // Initialise bodies involved in constraint by listing all bodies for which the first constraint parameter is estimated.
//            bodiesInvolved = parameters->getIndicesForParameterType( parametersInvolved.at( 0 ) ).first;
//        }
//        for ( auto parameterBodyId : bodiesInvolved )
//        {
//            bool constraintDetected = true;
//            std::map< estimatable_parameters::EstimatebleParametersEnum, std::pair<
//                    std::pair< std::string, std::string >, std::pair< unsigned int, unsigned int > > > singleConstraintMetadata;
//
//            for ( auto param : parametersInvolved )
//            {
//                std::pair< std::vector< std::pair< std::string, std::string > >, std::vector< std::pair< int, int > > > indicesParameter =
//                        parameters->getIndicesForParameterType( param );
//                auto it = find( indicesParameter.first.begin( ), indicesParameter.first.end( ), parameterBodyId );
//                if (it != indicesParameter.first.end( ) )
//                {
//                    int index = it - indicesParameter.first.begin( );
//                    singleConstraintMetadata[ param ] = std::make_pair( parameterBodyId, indicesParameter.second.at( index ) );
//                }
//                else
//                {
//                    constraintDetected = false;
//                    if ( !applyToAllBodies )
//                    {
//                        std::cerr << "Warning: custom constraint cannot be created: parameter " << std::to_string( param )
//                                  << " is not estimated for body " << parameterBodyId.first <<
//                                  " (secondary body " << parameterBodyId.second <<")\n\n";
//                    }
//                }
//            }
//            if ( constraintDetected )
//            {
//                activeConstraintsMetadata.push_back( singleConstraintMetadata );
//            }
//        }
//    }
//    else // if constraint has to be applied between parameters related to different bodies
//    {
//        bool constraintDetected = true;
//        std::map< estimatable_parameters::EstimatebleParametersEnum, std::pair<
//                std::pair< std::string, std::string >, std::pair< unsigned int, unsigned int > > > singleConstraintMetadata;
//
//        for ( unsigned int k = 0 ; k < parametersInvolved.size( ) ; k++ )
//        {
//            std::pair< std::vector< std::pair< std::string, std::string > >, std::vector< std::pair< int, int > > > indicesParameter =
//                    parameters->getIndicesForParameterType( parametersInvolved.at( k ) );
//
//            auto it = find( indicesParameter.first.begin( ), indicesParameter.first.end( ), bodiesInvolved.at( k ) );
//            if (it != indicesParameter.first.end( ) )
//            {
//                int index = it - indicesParameter.first.begin( );
//                singleConstraintMetadata[ parametersInvolved.at( k ) ] = std::make_pair( bodiesInvolved.at( k ), indicesParameter.second.at( index ) );
//            }
//            else
//            {
//                constraintDetected = false;
//                std::cerr << "Warning: custom constraint cannot be created: parameter " << std::to_string( parametersInvolved.at( k ) )
//                          << " is not estimated for body " << bodiesInvolved.at( k ).first <<
//                          " (secondary body " << bodiesInvolved.at( k ).second <<")\n\n";
//            }
//        }
//        if ( constraintDetected )
//        {
//            activeConstraintsMetadata.push_back( singleConstraintMetadata );
//        }
//    }
//
//
//    // PRINTING STATEMENTS
//    std::cout << "activeConstraintsMetadata: " << "\n\n";
//    for ( auto constraintMetadata : activeConstraintsMetadata )
//    {
//        std::cout << "detected custom constraint between \n\n";
//        for ( auto itr : constraintMetadata )
//        {
//            std::cout << itr.first << " for " << itr.second.first.first << " - " << itr.second.first.second <<
//                      " | indices: " << itr.second.second.first << " - " << itr.second.second.second << "\n\n";
//        }
//    }
//    // END PRINTING STATEMENTS
//
//    return activeConstraintsMetadata;
//}


void constraintSingleOrderLoveNumber(
        const estimatable_parameters::ParameterConstraintsEnum constraint,
        const std::shared_ptr< estimatable_parameters::TidalLoveNumber< Eigen::VectorXd > > loveNumber,
        std::vector< std::pair< unsigned int, unsigned int > >& constraintIndices,
        Eigen::VectorXd& k2Value )
{
    std::vector< int > orders = loveNumber->getOrders( );
    if ( orders.size( ) > 1 )
    {
        if ( find( orders.begin( ), orders.end( ), 2 ) == orders.end( ) )
        {
            throw std::runtime_error( "Error for constraint of type " + std::to_string( constraint ) + "(body " + loveNumber->getParameterName( ).second.first +
                                      ", multiple orders for k2. k22 should be used for the constraint, but order 2 is not estimated." );
        }
        else
        {
            std::cerr << "Warning, for  constraint of type " << constraint << "(body " << loveNumber->getParameterName( ).second.first << "), multiple orders for k2 "
                                                                              "are considered, but only k22 will be included in the constraint." << "\n\n";
        }

        std::cout << "inside constraintSingleOrderLoveNumber, constraintIndices: " << "\n\n";
        for ( auto ind : constraintIndices )
        {
            std::cout << ind.first << " - " << ind.second << "\n\n";
        }

        std::pair< unsigned int, unsigned int > savedIndexOtherTidalParameter = constraintIndices.at( 1 );
        Eigen::VectorXd savedOldK2Value = k2Value;

        unsigned int indexParam = constraintIndices.at( 0 ).first;
        unsigned int indexK2 = 0;

        std::cout << "indexParam: " << indexParam << "\n\n";
        std::cout << "indexK2: " << indexK2 << "\n\n";

        constraintIndices.clear( );
        k2Value.resize( 2 );
        for ( unsigned int i = 0 ; i < orders.size( ) ; i++ )
        {
            if ( orders.at( i ) == 2 )
            {
                constraintIndices.push_back( std::make_pair( indexParam, 2 ) );
                k2Value = savedOldK2Value.segment( indexK2, 2 );
            }
            indexParam += 2;
            indexK2 += 2;
        }

        // Add other tidal parameter indices
        constraintIndices.push_back( savedIndexOtherTidalParameter );
    }
}




}

}