/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PARAMETERCONSTRAINTS_H
#define TUDAT_PARAMETERCONSTRAINTS_H

#include <iostream>
#include <vector>
#include <string>
#include <vector>
#include <map>



#include <memory>
#include <Eigen/Geometry>

//#include "tudat/astro/propagators/singleStateTypeDerivative.h"
//#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
//#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/basics/utilities.h"

namespace tudat
{

namespace estimatable_parameters
{

enum ParameterConstraintsEnum
{
    tidal_quality_factor_single_love_number_constraint,
    tidal_time_lag_single_love_number_constraint,
    tidal_quality_factor_full_love_number_constraint,
    tidal_time_lag_full_love_number_constraint,
    custom_constraint
};

std::string getConstraintTypeString( const ParameterConstraintsEnum constraintType );
//{
//    std::string constraintDescription;
//    switch( constraintType )
//    {
//        case tidal_quality_factor_single_love_number_constraint:
//            constraintDescription = "constraint between tidal quality factor and single Love number ";
//            break;
//        case tidal_time_lag_single_love_number_constraint:
//            constraintDescription = "constraint between tidal time lag and single Love number ";
//            break;
//        case tidal_quality_factor_full_love_number_constraint:
//            constraintDescription = "constraint between tidal quality factor and full Love number ";
//            break;
//        case tidal_time_lag_full_love_number_constraint:
//            constraintDescription = "constraint between tidal time lag and full Love number ";
//            break;
//        case custom_constraint:
//            constraintDescription = "custom constraint ";
//            break;
//        default:
//            std::string errorMessage = "Error when getting constraint string, did not recognize constraint " +
//                                       std::to_string( constraintType );
//            throw std::runtime_error( errorMessage );
//    }
//    return constraintDescription;
//}



//! Container class for all parameters that are to be estimated.
/*!
 *  Container class for all parameters that are to be estimated. Class is templated with the scalar type used for the
 *  estimation of any initial dynamical states that may be included
 */
class ParameterConstraints
{
public:

    //! Constructor of parameter set.
    /*!
     *  Constructor of parameter set.
     *  \param estimatedDoubleParameters List of double parameters that are estimated.
     *  \param estimatedVectorParameters List of vector parameters that are estimated.
     *  \param estimateInitialStateParameters List of initial dynamical states that are to be estimated.
     */
    ParameterConstraints(
            const ParameterConstraintsEnum constraintType,
            const std::vector< std::pair< std::string, std::string > >& bodiesInvolved,
            const std::vector< std::pair< unsigned int, unsigned int > >& constraintIndices,
            const Eigen::VectorXd& parameterValues ):
            constraintType_( constraintType ), bodiesInvolved_( bodiesInvolved ), constraintIndices_( constraintIndices ), parameterValues_( parameterValues )
//            parameters_( parameters ), singleConstraintSize_( 0 ), fullConstraintsSize_( 0 )
    {
//        constraints_ = Eigen::MatrixXd::Zero( fullConstraintsSize_, parameters_->getParameterSetSize( ) );

    }

    //! Function to update the parameter constraints at each estimation iteration
    /*!
     *
     * parameters updated parameter set object
     */
    virtual void updateConstraints( const Eigen::VectorXd& newParameterValues )
    {
        std::cout << "in updateConstraints (base class)" << "\n\n";
        parameterValues_ = newParameterValues;
    }


    //! Function to return the indices of the parameters involved in the constraint
    /*!
     * Function to return the indices of the parameters involved in the constraint
     * \return vector containing the indices of each parameter involved
     */
     std::vector< std::pair< unsigned int, unsigned int > > getConstraintIndices( ) const
    {
         return constraintIndices_;
    }

    //! Function to return the constraint size
    /*!
     * Function to return the constraint size
     * \return constraintSize_
     */
    unsigned int getConstraintSize( ) const
    {
        return constraintSize_;
    }

    Eigen::MatrixXd getConstraint( ) const
    {
        return constraints_;
    }

//    virtual std::string getConstraintDescription( ) const
//    {
//        std::string constraintDescription = getConstraintTypeString( constraintType_ ) + "for body " + bodiesInvolved_.at( 0 ).first;
//        if( bodiesInvolved_.at( 0 ).second == "" )
//        {
//            constraintDescription += ".";
//        }
//        else
//        {
//            constraintDescription += ", " + bodiesInvolved_.at( 0 ).second + ".";
//        }
//        return constraintDescription;
//    }

    virtual std::string getConstraintDescription( ) const
    {
        std::string constraintDescription = getConstraintTypeString( constraintType_ ) + " (body " + bodiesInvolved_.at( 0 ).first;
        if ( bodiesInvolved_.at( 0 ).second == "" )
        {
            constraintDescription += ")";
        }
        else
        {
            constraintDescription += ", " + bodiesInvolved_.at( 0 ).second + ")";
        }

        return constraintDescription;
    }

protected:

    ParameterConstraintsEnum constraintType_;

    std::vector< std::pair< std::string, std::string > > bodiesInvolved_;

    std::vector< std::pair< unsigned int, unsigned int > > constraintIndices_;

    Eigen::VectorXd parameterValues_;

    unsigned int constraintSize_;

    Eigen::MatrixXd constraints_;

};


class TidalQualityFactorLoveNumberConstraints: public ParameterConstraints
{

public:

//! Constructor
/*!
 * Constructor
 * \param
 */
TidalQualityFactorLoveNumberConstraints(
        const ParameterConstraintsEnum constraintType,
        const std::vector< std::pair< std::string, std::string > >& involvedBodies,
        const std::vector< std::pair< unsigned int, unsigned int > >& constraintIndices,
        const Eigen::VectorXd& parameterValues ):
        ParameterConstraints( constraintType, involvedBodies, constraintIndices, parameterValues )
{
    std::cout << "size constraint indices: " << constraintIndices.size( ) << "\n\n";

    // Set constraint size
    this->constraintSize_ = 1;

    // Check number of indices provided
    if ( constraintIndices.size( ) != 2 )
    {
        throw std::runtime_error( "Error when creating TidalQualityFactorLoveNumberConstraints, the number of parameter indices provided should be 2." );
    }

    if ( ( constraintType != tidal_quality_factor_single_love_number_constraint  ) && ( constraintType != tidal_quality_factor_full_love_number_constraint ) )
    {
        throw std::runtime_error( "Error when creating TidalQualityFactorLoveNumberConstraints, the constraint must be of type "
                                  "tidal_quality_factor_single_love_number_constraint or tidal_quality_factor_full_love_number_constraint." );
    }

    // Set indices
//    for ( unsigned int i = 0 ; i < this->constraintIndices_.size( ) - 1 ; i++ )
//    {
//        indicesRealLoveNumbers_.push_back( this->constraintIndices_[ i ].first );
//        indicesImagLoveNumbers_.push_back( this->constraintIndices_[ i ].first + 1 );
//    }
//    nbK2Values_ = indicesRealLoveNumbers_.size( );
    indexRealLoveNumber_ = this->constraintIndices_[ 0 ].first;
    indexImagLoveNumber_ = this->constraintIndices_[ 0 ].first + 1;
    indexInvQ_ = this->constraintIndices_[ constraintIndices.size( ) - 1 ].first;

    std::cout << "indexRealLoveNumber_: " << indexRealLoveNumber_ << "\n\n";
    std::cout << "indexImagLoveNumber_: " << indexImagLoveNumber_ << "\n\n";
//    for ( unsigned int i = 0 ; i < indicesRealLoveNumbers_.size( ) ; i++ )
//    {
//        std::cout << "indexRealLoveNumber_: " << indicesRealLoveNumbers_.at( i ) << "\n\n";
//        std::cout << "indexImagLoveNumber_: " << indicesImagLoveNumbers_.at( i ) << "\n\n";
//    }
    std::cout << "indexInvQ_: " << indexInvQ_ << "\n\n";

    // Initialise constraints
    updateConstraints( this->parameterValues_ );

    std::cout << "constraints: " << this->constraints_ << "\n\n";

}

//! Destructor
~TidalQualityFactorLoveNumberConstraints( ){ }


//! Function to update the constraint(s) based on the parameter values
void updateConstraints( const Eigen::VectorXd& newParameterValues )
{
    this->parameterValues_ = newParameterValues;

    std::cout << "inside update parameter values: " << "\n\n";
    std::cout << this->parameterValues_.transpose( ) << "\n\n";
    std::cout << "indexRealLoveNumber_: " << indexRealLoveNumber_ << "\n\n";
    std::cout << "indexImagLoveNumber_: " << indexImagLoveNumber_ << "\n\n";

    // Re-compute k2
    k2_ = std::complex< double >( this->parameterValues_[ indexRealLoveNumber_ ], this->parameterValues_[ indexImagLoveNumber_ ] );
//    updateEffectiveK2Value( this->parameterValues_ );

    std::cout << "k2_: " << k2_ << "\n\n";

    std::cout << "in updateConstraints (TidalQualityFactorLoveNumberConstraints)" << "\n\n";

    double normLoveNumber = std::sqrt( k2_.real( ) * k2_.real( ) + k2_.imag( ) * k2_.imag( ) );
    double invQ = this->parameterValues_[ indexInvQ_ ];

    std::cout << "realLoveNumber: " << k2_.real( ) << "\n\n";
    std::cout << "imagLoveNumber: " << k2_.imag( ) << "\n\n";
    std::cout << "invQ: " << invQ << "\n\n";
    std::cout << "normLoveNumber: " << normLoveNumber << "\n\n";

    // update constraints
    this->constraints_ = Eigen::MatrixXd::Zero( this->constraintSize_, this->parameterValues_.size( ) );

//    for ( unsigned int i = 0 ; i < nbK2Values_ ; i++ )
//    {
        this->constraints_( 0, indexRealLoveNumber_ ) =
//              ( - k2_.real( ) *  k2_.imag( ) ) / std::pow( normLoveNumber, 3 );
                0.0;
        std::cout << "constraint Re(k2): " << this->constraints_( 0, indexRealLoveNumber_ ) << "\n\n";
        this->constraints_( 0, indexImagLoveNumber_ ) =
//                ( k2_.real( ) * k2_.real( ) ) / std::pow( normLoveNumber, 3 );
                1.0;
        std::cout << "constraint Im(k2): " << this->constraints_( 0, indexImagLoveNumber_ ) << "\n\n";
//    }
//    this->constraints_( 0, indexRealLoveNumber_ ) = 0.0; // - ( effectiveK2_.real( ) *  effectiveK2_.imag( ) ) / std::pow( normLoveNumber, 3 );
//    std::cout << "constraint Re(k2): " << this->constraints_( 0, indexRealLoveNumber_ ) << "\n\n";
//    this->constraints_( 0, indexImagLoveNumber_ ) = 1.0; //( effectiveK2_.real( ) * effectiveK2_.real( ) ) / std::pow( normLoveNumber, 3 );
//    std::cout << "constraint Im(k2): " << this->constraints_( 0, indexImagLoveNumber_ ) << "\n\n";

    this->constraints_( 0, indexInvQ_ ) = - 1.0;

    std::cout << "new constraint: " << this->constraints_ << "\n\n";

}

//void updateEffectiveK2Value( const Eigen::VectorXd& newParameterValues )
//{
//    effectiveK2_ = std::complex< double >( 0.0, 0.0 );
//    for ( unsigned int i = 0 ; i < nbK2Values_ ; i++ )
//    {
//        effectiveK2_.real( effectiveK2_.real( ) + this->parameterValues_[ indicesRealLoveNumbers_.at( i ) ] );
//        effectiveK2_.imag( effectiveK2_.imag( ) + this->parameterValues_[ indicesImagLoveNumbers_.at( i ) ] );
//    }
//    effectiveK2_.real( effectiveK2_.real( ) / nbK2Values_ );
//    effectiveK2_.imag( effectiveK2_.imag( ) / nbK2Values_ );
//
//    std::cout << "effective k2: " << effectiveK2_ << "\n\n";
//}

private:

    unsigned int indexRealLoveNumber_;
    unsigned int indexImagLoveNumber_;
//    std::vector< unsigned int > indicesRealLoveNumbers_;
//    std::vector< unsigned int > indicesImagLoveNumbers_;
//    unsigned int nbK2Values_;

//    std::complex< double > effectiveK2_;
    std::complex< double > k2_;

    unsigned int indexInvQ_;

};


class TidalTimeLagLoveNumberConstraints: public ParameterConstraints
{

public:

//! Constructor
/*!
* Constructor
* \param
*/
TidalTimeLagLoveNumberConstraints(
        const ParameterConstraintsEnum constraintType,
        const std::vector< std::pair< std::string, std::string > >& involvedBodies,
        const std::vector< std::pair< unsigned int, unsigned int > >& constraintIndices,
        const Eigen::VectorXd& parameterValues ):
        ParameterConstraints( constraintType, involvedBodies, constraintIndices, parameterValues )
{
    std::cout << "in constructor TidalTimeLagLoveNumberConstraints" << "\n\n";

    // Set constraint size
    this->constraintSize_ = 1;

    // Check number of indices provided
    if ( constraintIndices.size( ) != 2 )
    {
        throw std::runtime_error( "Error when creating TidalTimeLagLoveNumberConstraint, the number of parameter indices provided should be 2." );
    }

    if ( ( constraintType != tidal_time_lag_single_love_number_constraint  ) && ( constraintType != tidal_time_lag_full_love_number_constraint ) )
    {
        throw std::runtime_error( "Error when creating TidalTimeLagLoveNumberConstraints, the constraint must be of type "
                                  "tidal_time_lag_single_love_number_constraint or tidal_time_lag_full_love_number_constraint." );
    }

    // Set indices
    indexRealLoveNumber_ = this->constraintIndices_[ 0 ].first;
    indexImagLoveNumber_ = this->constraintIndices_[ 0 ].first + 1;
//    for ( unsigned int i = 0 ; i < this->constraintIndices_.size( ) - 1 ; i++ )
//    {
//        indicesRealLoveNumbers_.push_back( this->constraintIndices_[ i ].first );
//        indicesImagLoveNumbers_.push_back( this->constraintIndices_[ i ].first + 1 );
//    }
//    nbK2Values_ = indicesRealLoveNumbers_.size( );

    indexTimeLag_ = this->constraintIndices_[ constraintIndices.size( ) - 1 ].first;

    std::cout << "indexRealLoveNumber_: " << indexRealLoveNumber_ << "\n\n";
    std::cout << "indexImagLoveNumber_: " << indexImagLoveNumber_ << "\n\n";
//    for ( unsigned int i = 0 ; i < indicesRealLoveNumbers_.size( ) ; i++ )
//    {
//        std::cout << "indexRealLoveNumber_: " << indicesRealLoveNumbers_.at( i ) << "\n\n";
//        std::cout << "indexImagLoveNumber_: " << indicesImagLoveNumbers_.at( i ) << "\n\n";
//    }
    std::cout << "indexTimeLag_: " << indexTimeLag_ << "\n\n";

    // Initialise constraints
    updateConstraints( this->parameterValues_ );

    std::cout << "constraints: " << this->constraints_ << "\n\n";
}

//! Destructor
~TidalTimeLagLoveNumberConstraints( ){ }


//! Function to update the constraint(s) based on the parameter values
void updateConstraints( const Eigen::VectorXd& newParameterValues )
{
    this->parameterValues_ = newParameterValues;
    std::cout << "in updateConstraints (TidalTimeLagLoveNumberConstraints)" << "\n\n";

    std::cout << "inside update parameter values: " << "\n\n";
    std::cout << this->parameterValues_.transpose( ) << "\n\n";
    std::cout << "indexRealLoveNumber_: " << indexRealLoveNumber_ << "\n\n";
    std::cout << "indexImagLoveNumber_: " << indexImagLoveNumber_ << "\n\n";

    for ( unsigned int i = 0 ; i < this->parameterValues_.size( ) ; i++ )
    {
        std::cout << "test ind " << i << ": " << this->parameterValues_[ i ] << "\n\n";
    }


//    // Re-compute effective value for k2
//    updateEffectiveK2Value( this->parameterValues_ );
    // Re-compute k2
    k2_ = std::complex< double >( this->parameterValues_[ indexRealLoveNumber_ ], this->parameterValues_[ indexImagLoveNumber_ ] );

    std::cout << "k2_: " << k2_ << "\n\n";

    double normLoveNumber = std::sqrt( k2_.real( ) * k2_.real( ) + k2_.imag( ) * k2_.imag( ) );
    double invQ = k2_.imag( ) / normLoveNumber;

    double timeLag = this->parameterValues_[ indexTimeLag_ ];

    std::cout << "realLoveNumber: " << k2_.real( ) << "\n\n";
    std::cout << "imagLoveNumber: " << k2_.imag( ) << "\n\n";
    std::cout << "invQ: " << invQ << "\n\n";
    std::cout << "timeLag: " << timeLag << "\n\n";
    std::cout << "normLoveNumber: " << normLoveNumber << "\n\n";


    // update constraints
    this->constraints_ = Eigen::MatrixXd::Zero( this->constraintSize_, this->parameterValues_.size( ) );

    double multiplyingFactor = ( timeLag / std::atan( invQ ) )
            * 1.0 / ( normLoveNumber * ( k2_.real( ) * k2_.real( ) + 2.0 * k2_.imag( ) * k2_.imag( ) ) );
//    std::cout << "multiplyingFactor: " << multiplyingFactor << "\n\n";
//    std::cout << "tidal period: " << ( timeLag / std::atan( invQ ) ) * 2.0 * mathematical_constants::PI << "\n\n";

//    for ( unsigned int i = 0 ; i < nbK2Values_ ; i++ )
//    {
        this->constraints_( 0, indexRealLoveNumber_ ) =
                ( - k2_.real( ) * k2_.imag( ) ) * multiplyingFactor;
        this->constraints_( 0, indexImagLoveNumber_ ) =
                ( k2_.real( ) * k2_.real( ) ) * multiplyingFactor;
//    }
    this->constraints_( 0, indexTimeLag_ ) = - 1.0;

    std::cout << "new constraint: " << this->constraints_ << "\n\n";

}

//void updateEffectiveK2Value( const Eigen::VectorXd& newParameterValues )
//{
//    effectiveK2_ = std::complex< double >( 0.0, 0.0 );
//    for ( unsigned int i = 0 ; i < nbK2Values_ ; i++ )
//    {
//        effectiveK2_.real( effectiveK2_.real( ) + this->parameterValues_[ indicesRealLoveNumbers_.at( i ) ] );
//        effectiveK2_.imag( effectiveK2_.imag( ) + this->parameterValues_[ indicesImagLoveNumbers_.at( i ) ] );
//    }
//    effectiveK2_.real( effectiveK2_.real( ) / indicesRealLoveNumbers_.size( ) );
//    effectiveK2_.imag( effectiveK2_.imag( ) / indicesRealLoveNumbers_.size( ) );
//
//    std::cout << "effective k2: " << effectiveK2_ << "\n\n";
//}

private:

    unsigned int indexRealLoveNumber_;
    unsigned int indexImagLoveNumber_;
//    std::vector< unsigned int > indicesRealLoveNumbers_;
//    std::vector< unsigned int > indicesImagLoveNumbers_;
//    unsigned int nbK2Values_;
    std::complex< double > k2_;

//    std::complex< double > effectiveK2_;

    unsigned int indexTimeLag_;

};


class CustomConstraint: public ParameterConstraints
{

public:

//! Constructor
/*!
* Constructor
* \param
*/
CustomConstraint(
        const std::vector< std::string >& involvedParameters,
        const std::vector< std::pair< std::string, std::string > >& involvedBodies,
//        const std::pair< std::string, std::string >& constraintBodies,
        const std::vector< std::function< Eigen::MatrixXd ( Eigen::VectorXd ) > >& constraintFunctions,
        const std::vector< std::pair< unsigned int, unsigned int > >& constraintIndices,
        const Eigen::VectorXd& parameterValues ):
        ParameterConstraints( custom_constraint, involvedBodies, constraintIndices, parameterValues ),
        involvedParameters_( involvedParameters ), constraintFunctions_( constraintFunctions )
{
    // Check consistency input parameters
    nbParametersInvolved_ = constraintIndices.size( );

    if ( constraintIndices.size( ) != constraintFunctions.size( ) )
    {
        throw std::runtime_error("Error in CustomConstraint object, the numbers of parameter indices and functions provided do not match.");
    }

    if ( nbParametersInvolved_ < 2 )
    {
        throw std::runtime_error( "Error in CustomConstraint object, a minimum of two parameters must be involved in the constraint." );
    }

    // Initialise values of parameters involved in the constraint
    sizeParametersInvolved_ = 0;
    for ( unsigned int i = 0 ; i < nbParametersInvolved_ ; i++ )
    {
        sizeParametersInvolved_ += constraintIndices.at( i ).second;
    }
    updateParametersInvolvedValues( );


    // Check consistency of the constraint functions
    this->constraintSize_ = constraintFunctions.at( 0 )( parametersInvolvedValues_ ).rows( );
    for ( unsigned int i = 1 ; i < nbParametersInvolved_ ; i++ )
    {
        if ( this->constraintSize_ != constraintFunctions.at( i )( parametersInvolvedValues_ ).rows( ) )
        {
            throw std::runtime_error( "Error in CustomConstraint object, the size of the constraint function is not consistent for all parameters involved" );
        }
    }

    // Initialise constraints
    updateConstraints( this->parameterValues_ );

//    std::cout << "constraints: " << this->constraints_ << "\n\n";

}


//! Destructor
~CustomConstraint( ){ }


//! Function to update the constraint(s) based on the parameter values
void updateConstraints( const Eigen::VectorXd& newParameterValues )
{
    this->parameterValues_ = newParameterValues;
    updateParametersInvolvedValues( );
    std::cout << "in updateConstraints (CustomConstraint)" << "\n\n";

    Eigen::MatrixXd newConstraint = Eigen::MatrixXd::Zero( this->constraintSize_, this->parameterValues_.size( ) );
    for ( unsigned int i = 0 ; i < nbParametersInvolved_ ; i++ )
    {
        Eigen::MatrixXd constraintFactor = constraintFunctions_.at( i )( parametersInvolvedValues_ );
        std::cout << "new constraint: " << constraintFactor << "\n\n";
        unsigned int indexParam = this->constraintIndices_.at( i ).first;
        newConstraint.block( 0, indexParam, this->constraintSize_, constraintFactor.cols( )  ) = constraintFactor;
    }

    this->constraints_ = newConstraint;
    std::cout << "new constraint: " << this->constraints_ << "\n\n";
}

std::string getConstraintDescription( ) const
{
    std::string constraintDescription = getConstraintTypeString( constraintType_ ) + "between ";
    for ( unsigned int k = 0 ; k < involvedParameters_.size( ) ; k++ )
    {
        constraintDescription += involvedParameters_.at( k ) + " (body " + this->bodiesInvolved_.at( k ).first;
        if ( bodiesInvolved_.at( k ).second == "" )
        {
            constraintDescription += ")";
        }
        else
        {
            constraintDescription += ", " + this->bodiesInvolved_.at( k ).second + ")";
        }
        if ( k <  involvedParameters_.size( ) - 1 )
        {
            constraintDescription += " and ";
        }
        else
        {
            constraintDescription += ".";
        }
    }

    return constraintDescription;
}

private:

    void updateParametersInvolvedValues( )
    {
        parametersInvolvedValues_ = Eigen::VectorXd::Zero( sizeParametersInvolved_ );
        unsigned int indexParam = 0;
        for ( unsigned int i = 0 ; i < this->constraintIndices_.size( ) ; i++ )
        {
            unsigned int sizeParam = this->constraintIndices_.at( i ).second;
            std::cout << "size param: " << sizeParam << "\n\n";
            parametersInvolvedValues_.segment( indexParam, sizeParam ) = this->parameterValues_.segment( this->constraintIndices_.at( i ).first, sizeParam );
            indexParam += sizeParam;
        }

        std::cout << "parametersInvolvedValues_: " << parametersInvolvedValues_ << "\n\n";

    }

    std::vector< std::string > involvedParameters_;

    std::vector< std::function< Eigen::MatrixXd ( Eigen::VectorXd ) > > constraintFunctions_;

    unsigned int nbParametersInvolved_;

    unsigned int sizeParametersInvolved_;

    Eigen::VectorXd parametersInvolvedValues_;


};


class CustomConstantConstraint: public ParameterConstraints
{

public:

//! Constructor
/*!
* Constructor
* \param
*/
CustomConstantConstraint(
        const std::vector< std::string >& involvedParameters,
        const std::vector< std::pair< std::string, std::string > >& involvedBodies,
//        const std::pair< std::string, std::string >& constraintBodies,
        const std::vector< Eigen::MatrixXd >& constantConstraint,
        const std::vector< std::pair< unsigned int, unsigned int > >& constraintIndices,
        const Eigen::VectorXd& parameterValues ):
        ParameterConstraints( custom_constraint, involvedBodies, constraintIndices, parameterValues ),
        involvedParameters_( involvedParameters ), constantConstraint_( constantConstraint )
{
    // Check consistency input parameters
    nbParametersInvolved_ = constraintIndices.size( );

    if ( constraintIndices.size( ) != constantConstraint.size( ) )
    {
        throw std::runtime_error("Error in CustomConstantConstraint object, the numbers of parameter indices and functions provided do not match.");
    }

    if ( nbParametersInvolved_ < 2 )
    {
        throw std::runtime_error( "Error in CustomConstantConstraint object, a minimum of two parameters must be involved in the constraint." );
    }

    // Initialise values of parameters involved in the constraint
    sizeParametersInvolved_ = 0;
    for ( unsigned int i = 0 ; i < nbParametersInvolved_ ; i++ )
    {
        sizeParametersInvolved_ += constraintIndices.at( i ).second;
    }
    updateParametersInvolvedValues( );

    // Check consistency of the constraint functions
    this->constraintSize_ = constantConstraint.at( 0 ).rows( );
    for ( unsigned int i = 1 ; i < nbParametersInvolved_ ; i++ )
    {
        if ( this->constraintSize_ != constantConstraint.at( i ).rows( ) )
        {
            throw std::runtime_error( "Error in CustomConstantConstraint object, the size of the constraint function is not consistent for all parameters involved" );
        }
    }

    // Initialise constraints
    updateConstraints( this->parameterValues_ );

}

//! Destructor
    ~CustomConstantConstraint( ){ }


//! Function to update the constraint(s) based on the parameter values
    void updateConstraints( const Eigen::VectorXd& newParameterValues )
    {
        this->parameterValues_ = newParameterValues;
        updateParametersInvolvedValues( );
        std::cout << "in updateConstraints (CustomConstraint)" << "\n\n";

        Eigen::MatrixXd newConstraint = Eigen::MatrixXd::Zero( this->constraintSize_, this->parameterValues_.size( ) );
        for ( unsigned int i = 0 ; i < nbParametersInvolved_ ; i++ )
        {
            Eigen::MatrixXd constraintFactor = constantConstraint_.at( i );
            std::cout << "new constraint: " << constraintFactor << "\n\n";
            unsigned int indexParam = this->constraintIndices_.at( i ).first;
            newConstraint.block( 0, indexParam, this->constraintSize_, constraintFactor.cols( )  ) = constraintFactor;
        }

        this->constraints_ = newConstraint;
        std::cout << "new constraint: " << this->constraints_ << "\n\n";
    }

    std::string getConstraintDescription( ) const
    {
        std::string constraintDescription = getConstraintTypeString( constraintType_ ) + "between ";
        for ( unsigned int k = 0 ; k < involvedParameters_.size( ) ; k++ )
        {
            constraintDescription += involvedParameters_.at( k ) + " (body " + this->bodiesInvolved_.at( k ).first;
            if ( this->bodiesInvolved_.at( k ).second == "" )
            {
                constraintDescription += ")";
            }
            else
            {
                constraintDescription += ", " + this->bodiesInvolved_.at( k ).second + ")";
            }
            if ( k <  involvedParameters_.size( ) - 1 )
            {
                constraintDescription += " and ";
            }
            else
            {
                constraintDescription += ".";
            }
        }

        return constraintDescription;
    }

private:

    void updateParametersInvolvedValues( )
    {
        parametersInvolvedValues_ = Eigen::VectorXd::Zero( sizeParametersInvolved_ );
        unsigned int indexParam = 0;
        for ( unsigned int i = 0 ; i < this->constraintIndices_.size( ) ; i++ )
        {
            unsigned int sizeParam = this->constraintIndices_.at( i ).second;
            std::cout << "size param: " << sizeParam << "\n\n";
            parametersInvolvedValues_.segment( indexParam, sizeParam ) = this->parameterValues_.segment( this->constraintIndices_.at( i ).first, sizeParam );
            indexParam += sizeParam;
        }

        std::cout << "parametersInvolvedValues_: " << parametersInvolvedValues_ << "\n\n";

    }

    std::vector< std::string > involvedParameters_;

    std::vector< Eigen::MatrixXd > constantConstraint_;

    unsigned int nbParametersInvolved_;

    unsigned int sizeParametersInvolved_;

    Eigen::VectorXd parametersInvolvedValues_;




};



} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_PARAMETERCONSTRAINTS_H
