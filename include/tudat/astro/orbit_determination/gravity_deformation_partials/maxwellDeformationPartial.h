/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MAXWELLDEFORMATIONPARTIAL_H
#define TUDAT_MAXWELLDEFORMATIONPARTIAL_H

// #include "tudat/math/basic/coordinateConversions.h"

#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"

// #include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/astro/orbit_determination/gravity_deformation_partials/deformationPartial.h"
// #include "tudat/astro/orbit_determination/acceleration_partials/sphericalHarmonicAccelerationPartial.h"
// #include "tudat/astro/orbit_determination/acceleration_partials/tidalLoveNumberPartialInterface.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"

namespace tudat
{

namespace acceleration_partials
{

Eigen::Matrix3d computeSphericalJacobian(const Eigen::Vector3d& p);

void computeFullSphericalStatePartials(const Eigen::Vector3d& p_b, 
                                       const Eigen::Vector3d& p_b_dot,
                                       const Eigen::Vector3d& omega_b,
                                       Eigen::Matrix<double,6,6>& dX_dS);

//! Class for calculating partial derivatives of a spherical harmonic gravitational acceleration.
/*!
 *  Class for calculating partial derivatives of a spherical harmonic gravitational acceleration, as calculated by the
 *  SphericalHarmonicsGravitationalAccelerationModel class.
 */
class MaxwellDeformationPartial : public DeformationPartial
{
public:
    //! Contructor.
    /*!
     *  Constructor, requires input on the acceleration model as of which partials are to be computed.
     *  If any partials of parameters of the rotation model of the body exerting acceleration are to be calculated,
     *  RotationMatrixPartial objects must be pre-constructed and passed here as a map, with one object for each parameter
     *  wrt which a partial is to be taken.
     *  \param acceleratedBody Name of body undergoing acceleration.
     *  \param acceleratingBody Name of body exerting acceleration.
     *  \param accelerationModel Spherical harmonic gravity acceleration model from which acceleration is calculated wrt
     *  which the object being constructed is to calculate partials.
     *  \param rotationMatrixPartials Map of RotationMatrixPartial, one for each paramater representing a property of the
     *  rotation of the body exerting the acceleration wrt which an acceleration partial will be calculated.
     *  \param tidalLoveNumberPartialInterfaces List of objects to compute partials of tidal gravity field variations, one
     *  per corresponding variation object in acceleratedBody.
     */
    MaxwellDeformationPartial(
            const std::string& deformingBody,
            const std::vector< std::string >& perturbingBodies, // WARNING: SHOULD BE A LIST
            const std::shared_ptr< basic_astrodynamics::MaxwellGravityDeformationModel > deformationModel,
            const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartials =
                    observation_partials::RotationMatrixPartialNamedList( ) );

    //! Destructor
    ~MaxwellDeformationPartial( ) { }

    //! Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration
     *  and adding it to the existing partial block
     *  Update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  undergoing acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtDeformationOfDeformingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                       const bool addContribution = 1,
                                       const int startRow = 0,
                                       const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 5, 5 ) += currentPartialWrtGravity_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 5, 5 ) -= currentPartialWrtGravity_;
        }
    }


    //! Function for determining if the acceleration is dependent on a non-translational integrated state.
    /*!
     *  Function for determining if the acceleration is dependent on a non-translational integrated state.
     *  No dependency is implemented, but a warning is provided if partial w.r.t. mass of body exerting acceleration
     *  (and undergoing acceleration if mutual attraction is used) is requested.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which dependency is to be determined.
     *  \return True if dependency exists (non-zero partial), false otherwise.
     */
    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes( const std::pair< std::string, std::string >& stateReferencePoint,
                                                                     const propagators::IntegratedStateType integratedStateType )
    {
        bool doesDependencyExist = false;
        if( ( stateReferencePoint.first == deformingBody_ || ( std::find( perturbingBodies_.begin( ), perturbingBodies_.end( ), stateReferencePoint.first ) != perturbingBodies_.end( ) ) ) &&
              integratedStateType == propagators::translational_state )
        {
            doesDependencyExist = true;
        }
        if ( ( stateReferencePoint.first == deformingBody_ ) && ( integratedStateType == propagators::rotational_state ) )
        {
            doesDependencyExist = true;
        }

        if( ( stateReferencePoint.first == deformingBody_ || ( std::find( perturbingBodies_.begin( ), perturbingBodies_.end( ), stateReferencePoint.first ) != perturbingBodies_.end( ) ) ) &&
              integratedStateType == propagators::body_mass_state )
        {
            throw std::runtime_error( "Warning, dependency of gravity deformation on body masses not yet implemented" );
        }
        return doesDependencyExist;
    }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency, 1 otherwise).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter );

    //! Function for calculating the partial of the acceleration w.r.t. a non-translational integrated state
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. a non-translational integrated state
     *  and adding it to the existing partial block. Function calls constituent spherical harmonic model functions
     *  \param partialMatrix Block of partial derivatives of where current partial is to be added.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which partial is to be computed.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     */
    void wrtNonDeformationStateOfAdditionalBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                                   const std::pair< std::string, std::string >& stateReferencePoint,
                                                   const propagators::IntegratedStateType integratedStateType,
                                                   const bool addContribution = true )
    { 
        // wrt rotational state of deforming body
        if( stateReferencePoint.first == deformingBody_ && integratedStateType == propagators::rotational_state )
        {
            Eigen::MatrixXd tempMatrix = Eigen::MatrixXd::Zero( 5, 7 );
            wrtRotationModelParameter( tempMatrix, estimatable_parameters::initial_rotational_body_state, "" );
            partialMatrix.block( 0, 0, 5, 7 ) = ( addContribution ? 1.0 : -1.0 ) * tempMatrix;
        }
        // wrt translational state of deforming body
        else if ( stateReferencePoint.first == deformingBody_ && integratedStateType == propagators::translational_state )
        {
            partialMatrix.block( 0, 0, 5, 6 ) += ( addContribution ? 1.0 : -1.0 ) * currentPartialWrtDeformingState_;
        }
        // wrt translational state of perturbing body
        else if ( ( std::find( perturbingBodies_.begin( ), perturbingBodies_.end( ), stateReferencePoint.first ) != perturbingBodies_.end( ) ) 
            && integratedStateType == propagators::translational_state )
        {
            partialMatrix.block( 0, 0, 5, 6 ) += ( addContribution ? 1.0 : -1.0 ) * currentPartialWrtPerturbingState_;
        }
    }

    //! Function to create a function returning the current partial w.r.t. a gravitational parameter.
    /*!
     * Function to create a function returning the current partial w.r.t. a gravitational parameter.
     * \param parameterId Identifier of parameter for which the partial is to be created.
     * \return Pair with partial function and paramater partial size. The partial function is non-empty only
     * if the parameterId input represents the gravitational parameter of acceleratingBody_ (or acceleratedBody_ if
     * accelerationUsesMutualAttraction_ is true).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getGravitationalParameterPartialFunction(
            const estimatable_parameters::EstimatebleParameterIdentifier& parameterId );

    //! Function for updating the partial object to current state and time.
    /*!
     *  Function for updating the partial object to current state and time. Calculates the variables that are
     *  used for the calculation of multple partials, to prevent multiple calculations of same function.
     *  \param currentTime Time to which object is to be updated (note that most update functions are time-independent,
     *  since the 'current' state of the bodies is typically updated globally by the NBodyStateDerivative class).
     */
    virtual void update( const double currentTime = TUDAT_NAN );


    void wrtGravitationalParameterOfDeformingBody
    ( Eigen::MatrixXd& partialMatrix, const int addPartial = 0 )
    {
        if( deformationModel_->getGravitationalParameterDeformingBody( ) != 0.0 )
        {
            Eigen::MatrixXd tempPartial = - 1.0 / ( 
                deformationModel_->getGlobalRelaxationTime( ) * deformationModel_->getGravitationalParameterDeformingBody( ) ) 
                    * ( equilibriumCoefficients_ 
                        + deformationModel_->getMaxwellRelaxationTime( ) * derivativeEquilibriumCoefficients_ );
            if( addPartial == 0 )
            {
                partialMatrix = tempPartial;
            }
            else if( addPartial == 1 )
            {
                partialMatrix += tempPartial;;
            }
            else if( addPartial == -1 )
            {
                partialMatrix -= tempPartial;;
            }
            else
            {
                throw std::runtime_error( "Error when adding partial of maxwell gravity deformation w.r.t. deforming body's mu, input is inconsistent" );
            }
        }
        else
        {
            throw std::runtime_error( "Error cannot compute partial of maxwell gravity deformation w.r.t deforming body's mu for zero value" );
        }
    }

    // TO BE MODIFIED - MULTIPLE PERTURBING BODIES
    void wrtGravitationalParameterOfPerturbingBody
    ( Eigen::MatrixXd& partialMatrix, const int addPartial = 0 )
    {
        if( deformationModel_->getGravitationalParameterDeformingBody( ) != 0.0 )
        {
            Eigen::MatrixXd tempPartial = 1.0 / ( 
                deformationModel_->getGlobalRelaxationTime( ) * deformationModel_->getGravitationalParameterPerturbingBody( )[ 0 ] ) // TO BE MODIFIED -> INDEX FIXED TO 0 FOR NOW
                    * ( equilibriumCoefficients_ 
                        + deformationModel_->getMaxwellRelaxationTime( ) * derivativeEquilibriumCoefficients_ );
            if( addPartial == 0 )
            {
                partialMatrix = tempPartial;
            }
            else if( addPartial == 1 )
            {
                partialMatrix += tempPartial;
            }
            else if( addPartial == -1 )
            {
                partialMatrix -= tempPartial; 
            }
            else
            {
                throw std::runtime_error( "Error when adding partial of maxwell gravity deformation w.r.t. perturbing body's mu, input is inconsistent" );
            }
        }
        else
        {
            throw std::runtime_error( "Error cannot compute partial of maxwell gravity deformation w.r.t perturbing body's mu for zero value" );
        }
    }

    //! Function to retrieve partial of acceleration wrt the position of body undergoing acceleration, in inertial coordinates.
    /*!
     * Function to retrieve the current partial of the acceleration wrt the position of the body undergoing the acceleration,
     * in inertial coordinates
     * \return Current partial of the acceleration wrt the position of the body undergoing the acceleration, in inertial coordinates.
     */
    Eigen::Matrix3d getCurrentPartialWrtGravity( )
    {
        return currentPartialWrtGravity_;
    }


    //! Function to calculate an acceleration partial wrt a rotational parameter.
    /*!
     *  Function to calculate an acceleration partial wrt a rotational parameter of the rotation model of the body
     *  exerting the acceleration.
     *  \param accelerationPartial Matrix of partials of spherical harmonic acceleration wrt a rotational parameter
     *  that is set by this function (returned by reference)
     *  \param parameterType Type of parameter wrt which a partial is to be calculated.
     *  An entry of the requested type must be present in the rotationMatrixPartials_ map.
     *  \param secondaryIdentifier Identifier required to unambiguously define the parameter (in addition to information in
     *  parameterType.
     */
    void wrtRotationModelParameter( Eigen::MatrixXd& accelerationPartial,
                                    const estimatable_parameters::EstimatebleParametersEnum parameterType,
                                    const std::string& secondaryIdentifier );

    Eigen::MatrixXd wrtOtherStateDerivative();

    Eigen::MatrixXd bodyFixedWrtRotational( );
    Eigen::MatrixXd bodyFixedWrtRotationParameter(
        const estimatable_parameters::EstimatebleParametersEnum parameterType,
        const std::string& secondaryIdentifier );
    Eigen::MatrixXd deformationWrtSphericalBodyFixedState( );
    Eigen::MatrixXd equilibriumCoefficientsWrtSphericalBodyFixedState( );
    Eigen::MatrixXd equilibriumCoefficientsDerivativeWrtSphericalBodyFixedState( );
    Eigen::MatrixXd equilibriumCoefficientsWrtAngularVelocityVector( );
    Eigen::MatrixXd equilibriumCoefficientsWrtRotationRate( );
    Eigen::MatrixXd equilibriumCoefficientsDerivativeWrtAngularVelocityVector( );
    Eigen::MatrixXd equilibriumCoefficientsDerivativeWrtAngularVelocityVectorDerivative( );
    Eigen::MatrixXd deformationWrtRotationalState( );
    Eigen::MatrixXd deformationWrtRotationRate( );
    Eigen::MatrixXd deformationWrtAngularVelocityVectorDerivative( );
    Eigen::Matrix6d sphericalWrtCartesianBodyFixedState();
    Eigen::Matrix6d bodyFixedWrtGlobalState( const double currentTime, const bool addIndirectRotationContribution );

    std::shared_ptr< basic_astrodynamics::MaxwellGravityDeformationModel > getDeformationModel( )
    { 
        return deformationModel_;
    }

protected:


    //! Function to calculate the partial of the deformation wrt a set of cosine coefficients.
    /*!
     *  Function to calculate the partial of the deformation wrt a set of cosine coefficients.
     *  The set of coefficients wrt which a partial is to be taken is provided as input.
     *  \param blockIndices List of cosine coefficient indices wrt which the partials are to be taken (first and second
     *  are degree and order for each vector entry).
     *  \param partialDerivatives Matrix of deformation partials that is set by this function (returned by reference),
     *  with each column containg the partial wrt a single coefficient (in same order as blockIndices).
     */
    void wrtCosineCoefficientBlock( const std::vector< std::pair< int, int > >& blockIndices, Eigen::MatrixXd& partialDerivatives );

    //! Function to calculate the partial of the deformation wrt a set of sine coefficients.
    /*!
     *  Function to calculate the partial of the deformation wrt a set of sine coefficients.
     *  The set of coefficients wrt which a partial is to be taken is provided as input.
     *  \param blockIndices List of sine coefficient indices wrt which the partials are to be taken (first and second
     *  are degree and order for each vector entry).
     *  \param partialDerivatives Matrix of deformation partials that is set by this function (returned by reference),
     *  with each column containg the partial wrt a single coefficient (in same order as blockIndices).
     */
    void wrtSineCoefficientBlock( const std::vector< std::pair< int, int > >& blockIndices, Eigen::MatrixXd& partialDerivatives );


    std::shared_ptr< basic_astrodynamics::MaxwellGravityDeformationModel > deformationModel_;

    bool includeOrder1_;

    // //! Cache object used for storing calculated values at current time and state for spherical harmonic gravity
    // //! calculations.
    // basic_mathematics::SphericalHarmonicsCache& sphericalHarmonicCache_;

    Eigen::Matrix3d currentRotationToInertialFrame_;

    Eigen::Matrix3d currentRotationToBodyFixedFrame_;

    //! Current body-fixed (w.r.t body exerting acceleration) position of body undergoing acceleration
    /*!
     *  Current body-fixed (w.r.t body exerting acceleration) position of body undergoing acceleration,
     *  set by update( time ) function.
     */
    Eigen::Vector3d bodyFixedPosition_;

    //! Current spherical coordinate of body undergoing acceleration
    /*!
     *  Current spherical coordinate of body undergoing acceleration, in reference frame fixed to body exerting acceleration.
     *  Order of components is radial distance (from center of body), latitude, longitude. Note that the the second entry
     *  differs from the direct output of the cartesian -> spherical coordinates, which produces a colatitude.
     */
    Eigen::Vector3d bodyFixedSphericalPosition_;

    //! The current partial of the acceleration wrt the position of the body undergoing the acceleration.
    /*!
     *  The current partial of the acceleration wrt the position of the body undergoing the acceleration.
     *  The partial wrt the position of the body exerting the acceleration is minus this value.
     *  Value is set by the update( time ) function.
     */
    Eigen::MatrixXd currentPartialWrtGravity_;

    Eigen::MatrixXd currentPartialWrtPosition_;
    Eigen::MatrixXd currentPartialWrtVelocity_;
    Eigen::MatrixXd currentPartialWrtDeformingState_;
    Eigen::MatrixXd currentPartialWrtPerturbingState_;

    //! The current partial of the acceleration wrt the position of the body undergoing the acceleration,
    //! with both acceleration and position in body-fixed frame.
    /*!
     *  The current partial of the acceleration wrt the position of the body undergoing the acceleration.
     *  with both acceleration and position in body-fixed frame. Value is set by the update( time ) function.
     */
    Eigen::Matrix3d currentBodyFixedPartialWrtPosition_;

    Eigen::VectorXd equilibriumCoefficients_;
    Eigen::VectorXd derivativeEquilibriumCoefficients_;

    //! Map of RotationMatrixPartial, one for each relevant rotation parameter
    /*!
     *  Map of RotationMatrixPartial, one for each parameter representing a property of the rotation of the
     *  body exerting the acceleration wrt which an acceleration partial will be calculated.
     *  Map is pre-created and set through the constructor.
     */
    observation_partials::RotationMatrixPartialNamedList rotationMatrixPartials_;

};

}  // namespace acceleration_partials

}  // namespace tudat

#endif  // TUDAT_MAXWELLDEFORMATIONPARTIAL_H
