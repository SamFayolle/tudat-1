/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DEFORMATIONPARTIALS_H
#define TUDAT_DEFORMATIONPARTIALS_H

#include <string>
#include <map>
#include <Eigen/Core>

#include "tudat/astro/basic_astro/gravityDeformationModel.h"
#include "tudat/astro/basic_astro/gravityDeformationModelTypes.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/stateDerivativePartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Base class for objects calculating partial derivatives of accelerations w.r.t. states, model parameters.
/*!
 *  Base class for objects calculating partial derivatives of accelerations  w.r.t. states, model parameters. Such
 *  calculations are used in orbit determination, for the computation of the state transition; sensitivity matrices.
 *  Derived classes implement derivative-calculating models for specific acceleration models, so that the calculation
 *  of all partials of a single type acceleration model is encompassed in a single derived class.
 */
class DeformationPartial : public orbit_determination::StateDerivativePartial
{
public:
    //! Base class constructor.
    /*!
     *  Constructor of base class, sets the base class member variables identifying the body undergoing and exerting the
     *  acceleration.
     *  \param acceleratedBody Body undergoing acceleration.
     *  \param acceleratingBody Body exerting acceleration.
     *  \param accelerationType Type of acceleration w.r.t. which partial is taken.
     */
    DeformationPartial( const std::string& deformingBody,
                        const std::vector< std::string >& perturbingBodies, // WARNING: SHOULD BE A LIST
                        const basic_astrodynamics::GravityDeformationType deformationType ):
        StateDerivativePartial( propagators::gravity_deformation_state, std::make_pair( deformingBody, "" ) ),
        deformingBody_( deformingBody ), perturbingBodies_( perturbingBodies ), deformationType_( deformationType )
    { }

    //! Virtual destructor.
    virtual ~DeformationPartial( ) { }

    //! Function to retrieve the function that returns the partial derivative w.r.t. a propagated state.
    /*!
     * Function to retrieve the function that returns the partial derivative w.r.t. a propagated state.
     * \param stateReferencePoint Reference point (id) for propagated state (i.e. body name for translational dynamics).
     * \param integratedStateType Type of propagated state.
     * \return Pair with function, returning partial derivative, and number of columns in partial vector,
     */
    std::pair< std::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int > getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        // Initialize to empty function; 0 parameter size.
        std::pair< std::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int > partialFunction =
                std::make_pair( std::function< void( Eigen::Block< Eigen::MatrixXd > ) >( ), 0 );

        // Check if state dependency exists
        switch( integratedStateType )
        {
            case propagators::gravity_deformation_state: 
            {
                // Check if reference id is consistent.
                if( stateReferencePoint.second != "" )
                {
                    throw std::runtime_error(
                            "Error when getting state derivative partial deformation model, cannot have reference point on body for "
                            "dynamics" );
                }
                // Check if propagated body corresponds to deforming, perturbing, or relevant third body.
                else if( stateReferencePoint.first == deformingBody_ )
                {
                    partialFunction =
                            std::make_pair( std::bind( &DeformationPartial::wrtStateOfDeformingBody, this, std::placeholders::_1 ), 5 );
                }
                // else if( std::find( perturbingBodies_.begin( ), perturbingBodies_.end( ), stateReferencePoint.first ) != perturbingBodies_.end( ) )
                // {
                //     // TO BE MODIFIED: Add index perturbing body maybe?
                //     partialFunction =
                //             std::make_pair( std::bind( &DeformationPartial::wrtStateOfDeformingBody, this, std::placeholders::_1 ), 5 );
                // }
                // else if( isDeformationPartialWrtAdditionalBodyNonnullptr( stateReferencePoint.first ) )
                // {
                //     partialFunction = std::make_pair( std::bind( &DeformationPartial::wrtStateOfAdditionalBody,
                //                                                  this,
                //                                                  std::placeholders::_1,
                //                                                  stateReferencePoint.first ), 5 );
                // }
                break;
            }
            case propagators::translational_state:
            {
                // Check if reference id is consistent.
                if( stateReferencePoint.second != "" )
                {
                    throw std::runtime_error(
                            "Error when getting state derivative partial deformation model, cannot have reference point on body for body "
                            "translational state" );
                }
                else if( isStateDerivativeDependentOnIntegratedAdditionalStateTypes( stateReferencePoint, integratedStateType ) )
                {
                    partialFunction = std::make_pair( std::bind( &DeformationPartial::wrtNonDeformationStateOfAdditionalBody,
                                                                 this,
                                                                 std::placeholders::_1,
                                                                 stateReferencePoint,
                                                                 integratedStateType,
                                                                 true ), 6 );
                }
                break;
            }
            case propagators::rotational_state: {
                // Check if reference id is consistent.
                if( stateReferencePoint.second != "" )
                {
                    throw std::runtime_error(
                            "Error when getting state derivative partial deformation model, cannot have reference point on body for body "
                            "rotational state" );
                }
                else if( isStateDerivativeDependentOnIntegratedAdditionalStateTypes( stateReferencePoint, integratedStateType ) )
                {
                    partialFunction = std::make_pair( std::bind( &DeformationPartial::wrtNonDeformationStateOfAdditionalBody,
                                                                 this,
                                                                 std::placeholders::_1,
                                                                 stateReferencePoint,
                                                                 integratedStateType,
                                                                 true ), 7 );
                }
                break;
            }
            case propagators::body_mass_state: {
                // Check if reference id is consistent.
                if( stateReferencePoint.second != "" )
                {
                    throw std::runtime_error(
                            "Error when getting state derivative partial deformation model, cannot have reference point on body for body "
                            "mass" );
                }
                else if( isStateDerivativeDependentOnIntegratedAdditionalStateTypes( stateReferencePoint, integratedStateType ) )
                {
                    partialFunction = std::make_pair( std::bind( &DeformationPartial::wrtNonDeformationStateOfAdditionalBody,
                                                                 this,
                                                                 std::placeholders::_1,
                                                                 stateReferencePoint,
                                                                 integratedStateType,
                                                                 true ), 1 );
                }
                break;
            }
            case propagators::custom_state: {
                break;
            }
            default:
                std::string errorMessage = "Error when getting state derivative partial deformation model, dynamics type " +
                        std::to_string( integratedStateType ) + "not recognized";
                throw std::runtime_error( errorMessage );
                break;
        }

        return partialFunction;
    }

    //! Function for determining if the acceleration is dependent on a non-deformational integrated state (default none).
    /*!
     *  Function for determining if the acceleration is dependent on a non-deformational integrated state.
     *  No dependency is implemented is returned in this base class function, but may be overriden by derived class.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which dependency is to be determined.
     *  \return True if dependency exists (non-zero partial), false otherwise.
     */
    virtual bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        return false;
    }


    //! Pure virtual function for calculating the partial of the acceleration w.r.t. the position of the accelerated body.
    /*!
     *  Pure virtual function for calculating the partial of the acceleration w.r.t. the position of the accelerated body and
     *  adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  undergoing acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    virtual void wrtDeformationOfDeformingBody( 
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const bool addContribution = 1,
        const int startRow = 0,
         const int startColumn = 0 ) = 0;


    //! Function for calculating the partial of the acceleration w.r.t. the Cartesian state of the body undergoing acceleration.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the Cartesian state of the body
     *  undergoing acceleration  and adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian state of body
     *  undergoing acceleration where current partial is to be added.
     */
    void wrtStateOfDeformingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix )
    {
        wrtDeformationOfDeformingBody( partialMatrix, true, 0, 0 );
        // wrtVelocityOfAcceleratedBody( partialMatrix, true, 0, 3 );
    }


    //! Function for calculating the partial of the acceleration w.r.t. a non-translational integrated state
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. a non-translational integrated state
     *  and adding it to the existing partial block. Function may be overridden in derived class, default dependency is
     *  none.
     *  \param partialMatrix Block of partial derivatives of where current partial is to be added.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which partial is to be computed.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     */
    virtual void wrtNonDeformationStateOfAdditionalBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                                           const std::pair< std::string, std::string >& stateReferencePoint,
                                                           const propagators::IntegratedStateType integratedStateType,
                                                           const bool addContribution = true )
    { }

    //! Function to retrieve the name of the body undergoing acceleration.
    /*!
     *  Function to retrieve the name of the body undergoing acceleration.
     *  \return Name of the body undergoing acceleration.
     */
    std::string getDeformingBody( )
    {
        return deformingBody_;
    }

    //! Function to retrieve the name of the body exerting acceleration.
    /*!
     *  Function to retrieve the name of the body exerting acceleration.
     *  \return Name of the body exerting acceleration.
     */
    std::vector< std::string > getPerturbingBodies( )
    {
        return perturbingBodies_;
    }

    //! Function to retrieve the type of acceleration w.r.t. which partial is taken..
    /*!
     *  Function to retrieve the type of acceleration w.r.t. which partial is taken..
     *  \return Type of acceleration w.r.t. which partial is taken..
     */
    basic_astrodynamics::GravityDeformationType getDeformationType( )
    {
        return deformationType_;
    }

protected:
    //! Name of the body undergoing acceleration.
    std::string deformingBody_;

    //! Name of the body exerting acceleration.
    std::vector< std::string > perturbingBodies_;

    //! Type of deformation the body is undergoing
    basic_astrodynamics::GravityDeformationType deformationType_;
};

}  // namespace acceleration_partials

}  // namespace tudat

#endif  // TUDAT_DEFORMATIONPARTIALS_H
