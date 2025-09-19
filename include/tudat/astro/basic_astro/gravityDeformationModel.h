/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_GRAVITY_DEFORMATION_MODEL_H
#define TUDAT_GRAVITY_DEFORMATION_MODEL_H

#include <vector>
#include <map>
#include <unordered_map>

#include <memory>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/math/basic/sphericalHarmonics.h"
#include "tudat/math/basic/coordinateConversions.h"

namespace tudat
{
namespace basic_astrodynamics
{

//! Base class for gravity deformation models 
/*!
 * Base class for gravity deformation models. Derived classes should contain
 * implementations to perform calculations of gravity field deformations. 
 */
class GravityDeformationModel
{
public:

    //! Constructor.
    GravityDeformationModel( ):
        currentTime_( TUDAT_NAN ),
        currentDeformation_( Eigen::VectorXd::Zero( 5 ) ){ }

    //! Virtual destructor.
    /*!
     * Virtual destructor, necessary to ensure that derived class destructors get called correctly.
     */
    virtual ~GravityDeformationModel( ) { }


    //! Update member variables used by the gravity deformation model.
    /*!
     * Updates member variables used by the gravity deformation model. In the case of acceleration models
     * containing varying parameters, function-pointers returning such a parameter (for instance
     * the Cartesian state of a body) will be set as a member variable.
     * This function evaluates such function-pointers and updates member variables to the 'current'
     * values of these parameters. Only these current values, not the function-pointers are then
     * used by the getAcceleration() function.
     *
     * N.B.: This pure virtual function must be overridden by derived classes!
     * \param currentTime Time at which acceleration model is to be updated.
     */
    virtual void updateMembers( const double currentTime = TUDAT_NAN ) = 0;

    Eigen::VectorXd& getDeformationReference( )
    {
        return currentDeformation_;
    }

    Eigen::VectorXd getDeformation( )
    {
        return currentDeformation_;
    }

    void getAccelerationByReference( Eigen::VectorXd& deformation ) const
    {
        deformation = currentDeformation_;
    }

    void addCurrentDeformation( Eigen::VectorXd& deformation ) const
    {
        if ( deformation.size( ) != currentDeformation_.size( ) )
        {
            throw std::runtime_error( "Error when adding current gravity deformation, inconsistent sizes." );
        }
        deformation += currentDeformation_;
    }

    //! Function to reset the current time
    /*!
     * Function to reset the current time of the deformation gravity model.
     * \param currentTime Current time (default NaN).
     */
    virtual void resetCurrentTime( )
    {
        currentTime_ = TUDAT_NAN;
    }

protected:

    //! Previous time to which the gravity deformation model was updated.
    double currentTime_;

    Eigen::VectorXd currentDeformation_;

protected:

private:

};



//! Typedef for the gravity deformation model map.
typedef std::map< std::string, std::vector< std::shared_ptr< GravityDeformationModel > > > GravityDeformationModelMap;


//! Class for Maxwell rheology gravity deformation model.
/*!
 * This class implements a gravity field deformation model assuming a Maxwell rheology.
 */
class MaxwellGravityDeformationModel : public GravityDeformationModel 
{
private:

    //! Typedef for coefficient-matrix-returning function.
    typedef std::function< Eigen::MatrixXd( ) > CoefficientMatrixReturningFunction;


public:

    //! Typedef for a position-returning function.
    typedef std::function< void( Eigen::Vector6d& ) > StateFunction;

    //! Constructor taking position-functions for bodies, and constant parameters of spherical
    //! harmonics expansion.
    /*!
     * Constructor taking a pointer to a function returning the position of the body subject to
     * gravitational acceleration, constant gravitational parameter and equatorial radius of the
     * body exerting the acceleration, constant coefficient matrices for the spherical harmonics
     * expansion, and a pointer to a function returning the position of the body exerting the
     * gravitational acceleration (typically the central body). This constructor uses the
     * Boost::lambda library to create a function on-the-fly that returns the constant
     * gravitational parameter, equatorial radius and coefficient matrices provided. The
     * constructor also updates all the internal members. The position of the body exerting the
     * gravitational acceleration is an optional parameter; the default position is the origin.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     * \param aGravitationalParameter A (constant) gravitational parameter [m^2 s^-3].
     * \param anEquatorialRadius A (constant) equatorial radius [m].
     * \param aCosineHarmonicCoefficientMatrix A (constant) cosine harmonic coefficient matrix.
     * \param aSineHarmonicCoefficientMatrix A (constant) sine harmonic coefficient matrix.
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration (default = (0,0,0)).
     * \param rotationFromBodyFixedToIntegrationFrameFunction Function providing the rotation from
     * body-fixes from to the frame in which the numerical integration is performed.
     * \param isMutualAttractionUsed Variable denoting whether attraction from body undergoing acceleration on
     * body exerting acceleration is included (i.e. whether aGravitationalParameter refers to the property
     * of the body exerting the acceleration, if variable is false, or the sum of the gravitational parameters,
     * if the variable is true.
     * \param sphericalHarmonicsCache Cache object for computing/retrieving repeated terms in spherical harmonics potential
     *          gradient calculation.
     */
    MaxwellGravityDeformationModel(
            const StateFunction stateOfDeformingBodyFunction,
            const std::string perturbingBody,
            const double maxwellRelaxationTime,
            const double globalRelaxationTime,
            const double gravitationalParameterDeformingBody,
            const double gravitationalParameterPerturbingBody,
            const double referenceRadius,
            const std::function< Eigen::Vector3d( ) > angularVelocityDeformingBody,
            const std::function< Eigen::Vector3d( ) > angularVelocityDerivativeDeformingBody,
            const double k2,
            CoefficientMatrixReturningFunction cosineCoefficients,
            CoefficientMatrixReturningFunction sineCoefficients,
            const StateFunction stateOfPerturbingBodyFunction =
            [ ]( Eigen::Vector6d& input ){ input = Eigen::Vector6d::Zero( ); },
            const std::function< Eigen::Quaterniond( ) >
            rotationFromBodyFixedToIntegrationFrameFunction =
            [ ]( ){ return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); },
            const std::function< Eigen::Matrix3d( ) >
                rotationToLocalFrameDerivativeFunction = [ ]( ){ return Eigen::Matrix3d::Zero( ); },
            const Eigen::VectorXd staticCoefficients = Eigen::VectorXd::Zero( 3 ),
            const bool includeOrder1 = true ) //,
            /*std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache =
            std::make_shared< basic_mathematics::SphericalHarmonicsCache >( ) )*/ : 
            GravityDeformationModel( ),
            stateOfDeformingBodyFunction_( stateOfDeformingBodyFunction ),
            perturbingBody_( perturbingBody ),
            maxwellRelaxationTime_( maxwellRelaxationTime ),
            globalRelaxationTime_( globalRelaxationTime ),
            gravitationalParameterDeformingBody_( gravitationalParameterDeformingBody ),
            gravitationalParameterPerturbingBody_( gravitationalParameterPerturbingBody ),
            referenceRadius_( referenceRadius ),
            angularVelocityDeformingBody_( angularVelocityDeformingBody ),
            angularVelocityDerivativeDeformingBody_( angularVelocityDerivativeDeformingBody ),
            k2_( k2 ),
            // rotationRateDerivative_( 0.0 ),
            getCosineHarmonicsCoefficients( cosineCoefficients ),
            getSineHarmonicsCoefficients( sineCoefficients ),
          /*getCosineHarmonicsCoefficients( [ = ]( ){ return aCosineHarmonicCoefficientMatrix; } ),
          getSineHarmonicsCoefficients( [ = ]( ){ return aSineHarmonicCoefficientMatrix; } ),*/
          stateOfPerturbingBodyFunction_( stateOfPerturbingBodyFunction ),
          rotationFromBodyFixedToIntegrationFrameFunction_(
              rotationFromBodyFixedToIntegrationFrameFunction ),
          rotationToBodyFixedDerivativeFunction_( rotationToLocalFrameDerivativeFunction ),
          staticCoefficients_( staticCoefficients ),
          includeOrder1_( includeOrder1 ),
          /*sphericalHarmonicsCache_( sphericalHarmonicsCache ),*/
          saveSphericalHarmonicTermsSeparately_( false )
    {

        // Initialise nominal coefficient values
        nominalCoefficients_ = Eigen::VectorXd::Zero( 5 );
        nominalCoefficients_[ 0 ] = getCosineHarmonicsCoefficients( )( 2, 0 );
        nominalCoefficients_[ 1 ] = getCosineHarmonicsCoefficients( )( 2, 1 );
        nominalCoefficients_[ 2 ] = getCosineHarmonicsCoefficients( )( 2, 2 );
        nominalCoefficients_[ 3 ] = getSineHarmonicsCoefficients( )( 2, 1 );
        nominalCoefficients_[ 4 ] = getSineHarmonicsCoefficients( )( 2, 2 );
        std::cout << "original nominal coefficients: " << nominalCoefficients_.transpose( ) << std::endl;

        equilibriumCoefficients_ = Eigen::VectorXd::Zero( 5 ); 
        derivativeEquilibriumCoefficients_ = Eigen::VectorXd::Zero( 5 );

        // Tranform to **unnormalised** coefficients
        staticCoefficients_[ 0 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
        staticCoefficients_[ 1 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
        staticCoefficients_[ 2 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
        staticCoefficients_[ 3 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
        staticCoefficients_[ 4 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );

        // // Update rotation and states
        // rotationToIntegrationFrame_ = rotationFromBodyFixedToIntegrationFrameFunction_( );
        // stateOfDeformingBodyFunction_( stateOfDeformingBody_ );
        // stateOfPerturbingBodyFunction_( stateOfPerturbingBody_ );
        // // std::cout << "positionOfDeformingBody_ " << positionOfDeformingBody_.transpose( ) << std::endl;
        // // std::cout << "positionOfPerturbingBody_ " << positionOfPerturbingBody_.transpose( ) << std::endl;

        // // Compute relative state
        // currentInertialRelativeState_ = stateOfPerturbingBody_ - stateOfDeformingBody_;
        // currentRelativePosition_ = rotationToIntegrationFrame_.inverse( ) * ( currentInertialRelativeState_.segment( 0, 3 ) );

        // Eigen::Vector3d currentSphericalPositionPerturbingBody = coordinate_conversions::convertCartesianToSpherical( 
        //     currentRelativePosition_ );
        // currentLongitude_ = currentSphericalPositionPerturbingBody[ 2 ];
        // double latitude = mathematical_constants::PI / 2.0 - currentSphericalPositionPerturbingBody.y( );

        // // rotationRateDerivative_ = 0.0;

        // // Initialise equilibrium coefficients
        // updateEquilibriumDeformation( );
        // std::cout << "original equilibriumCoefficients: " << equilibriumCoefficients_.transpose( ) << std::endl;

    }



    // //! Get gravitational acceleration in body-fixed frame of body undergoing acceleration.
    // /*!
    //  * Returns the gravitational acceleration in body-fixed frame of body undergoing acceleration computed
    //  * computed by the updateMembers function.
    //  * \return Computed gravitational acceleration vector.
    //  */
    // Eigen::Vector3d getAccelerationInBodyFixedFrame( )
    // {
    //     return currentAccelerationInBodyFixedFrame_;
    // }

    //! Update class members.
    /*!
     * Updates all the base class members to their current values and also updates the class
     * members of this class.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        std::cout.precision( 20 );

        // std::cout << "in updateMembers -- " << std::endl;

        // if( !( this->currentTime_ == currentTime ) )
        // {
            // std::cout << "in updateMembers -- " << std::endl;
            // std::cout << "currentTime " << currentTime << std::endl;

            // Update gravity coefficients
            cosineHarmonicCoefficients = getCosineHarmonicsCoefficients( );
            sineHarmonicCoefficients = getSineHarmonicsCoefficients( );

            nominalCoefficients_[ 0 ] = cosineHarmonicCoefficients( 2, 0 );
            nominalCoefficients_[ 1 ] = cosineHarmonicCoefficients( 2, 1 );
            nominalCoefficients_[ 2 ] = cosineHarmonicCoefficients( 2, 2 );
            nominalCoefficients_[ 3 ] = sineHarmonicCoefficients( 2, 1 );   
            nominalCoefficients_[ 4 ] = sineHarmonicCoefficients( 2, 2 );            
            // std::cout << "in update, normalised nominal coefs with static " << nominalCoefficients_[ 0 ] << " " << nominalCoefficients_[ 1 ]
            // << " " << nominalCoefficients_[ 2 ] << std::endl;

            // Tranform to **unnormalised** coefficients
            nominalCoefficients_[ 0 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
            nominalCoefficients_[ 1 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
            nominalCoefficients_[ 2 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
            nominalCoefficients_[ 3 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
            nominalCoefficients_[ 4 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );

            // std::cout << "in update, unnormalised nominal coefs with static " << nominalCoefficients_[ 0 ] << " " << nominalCoefficients_[ 1 ]
            // << " " << nominalCoefficients_[ 2 ] << std::endl;

            // Remove static field contribution
            nominalCoefficients_ = nominalCoefficients_ - staticCoefficients_;
            // std::cout << "in update, unnormalised nominal coefs without static " << nominalCoefficients_[ 0 ] << " " << nominalCoefficients_[ 1 ]
            // << " " << nominalCoefficients_[ 2 ] << std::endl;

            // Update rotation and positions
            rotationToIntegrationFrame_ = rotationFromBodyFixedToIntegrationFrameFunction_( );
            stateOfDeformingBodyFunction_( stateOfDeformingBody_ );
            stateOfPerturbingBodyFunction_( stateOfPerturbingBody_ );
            // this->updateBaseMembers( );

            // Compute relative inertial state
            currentInertialRelativeState_ = stateOfPerturbingBody_ - stateOfDeformingBody_;
            // std::cout << "currentInertialRelativeState_ " << currentInertialRelativeState_.transpose( ) << std::endl;

            Eigen::Matrix3d currentRotationToLocalFrameDerivative = rotationToBodyFixedDerivativeFunction_( );
            // std::cout << "currentRotationToLocalFrameDerivative " << std::endl;
            // std::cout << currentRotationToLocalFrameDerivative << std::endl;

            // std::cout << "rotationToIntegrationFrame_ " << std::endl;
            // std::cout << rotationToIntegrationFrame_.toRotationMatrix( ) << std::endl;

            // std::cout << "test 1" << std::endl;
            // std::cout << rotationToIntegrationFrame_.inverse( ) * currentInertialRelativeState_.segment( 3, 3 ) << std::endl;
            // std::cout << "test 2" << std::endl;
            // std::cout << currentRotationToLocalFrameDerivative * currentInertialRelativeState_.segment( 0, 3 ) << std::endl;


            // Compute current relative state in body-fixed frame
            currentRelativePosition_ = rotationToIntegrationFrame_.inverse( ) * currentInertialRelativeState_.segment( 0, 3 );
            currentRelativeVelocity_ = rotationToIntegrationFrame_.inverse( ) * currentInertialRelativeState_.segment( 3, 3 )
                + currentRotationToLocalFrameDerivative * currentInertialRelativeState_.segment( 0, 3 );

            // std::cout << "rotationToIntegrationFrame_.inverse( ) " << ( rotationToIntegrationFrame_.inverse( ) ).toRotationMatrix() << std::endl;

            // Compute spherical coordinates of perturbing body in body-fixed frame
            Eigen::Vector3d currentSphericalPositionPerturbingBody = 
                coordinate_conversions::convertCartesianToSpherical( currentRelativePosition_ );
            currentLongitude_ = currentSphericalPositionPerturbingBody[ 2 ]; 
            currentLatitude_ = mathematical_constants::PI / 2.0 - currentSphericalPositionPerturbingBody.y( );

            // MOVE DERIVATIVE CALCULATION TO UPDATE EQUILIBRIUM(?)
            // Compute current derivative of the perturbing body's body-fixed longitude
            currentLongitudeDerivative_ = (
                ( currentRelativeVelocity_[ 1 ] * currentRelativePosition_[ 0 ] 
                - currentRelativeVelocity_[ 0 ] * currentRelativePosition_[ 1 ] ) 
                / ( currentRelativePosition_[ 0 ] * currentRelativePosition_[ 0 ] + currentRelativePosition_[ 1 ] * currentRelativePosition_[ 1 ] ) );

            // Compute current derivative of the perturbing body's body-fixed latitude
            double currentDistance = currentRelativePosition_.segment( 0, 3 ).norm( );
            double currentDistanceDerivative = 
                ( currentRelativePosition_[ 0 ] * currentRelativeVelocity_[ 0 ] 
                + currentRelativePosition_[ 1 ] * currentRelativeVelocity_[ 1 ] 
                + currentRelativePosition_[ 2 ] * currentRelativeVelocity_[ 2 ] ) / currentDistance;
            currentLatitudeDerivative_ = 
                ( currentRelativeVelocity_[2] * currentDistance - currentRelativePosition_[2] * currentDistanceDerivative ) /
                ( currentDistance * std::sqrt( currentRelativePosition_[ 0 ] * currentRelativePosition_[ 0 ] + currentRelativePosition_[ 1 ] * currentRelativePosition_[ 1 ] ) );

            // std::cout << "in updateMembers time " << currentTime << " longitude " << currentLongitude_ << " distance " <<
                //  currentRelativePosition_.segment( 0, 3 ).norm( ) << " longitude derivative " << currentLongitudeDerivative_ << std::endl;

            // std::cout << "rotation " << rotationToIntegrationFrame_.toRotationMatrix( ) << std::endl;;  
            // std::cout << "position deforming " << positionOfDeformingBody_.transpose( ) << std::endl;
            // std::cout << "position perturbing " << positionOfPerturbingBody_.transpose( ) << std::endl;

            updateEquilibriumDeformation( );
            // updatePropagatedCoefficients( );

            // // THE CURRENT COEFFICIENTS SHOULD BE UPDATED AT THIS POINT
            // std::cout << "currentLongitude_ " << currentLongitude_ << std::endl;
            // std::cout << "currentLatitude_ " << currentLatitude_ << std::endl;
            //  std::cout << "equilibriumCoefficients_: " << equilibriumCoefficients_.transpose( ) << std::endl;
            //  std::cout << "in update equilibriumCoefficients: " << equilibriumCoefficients_[0] << " " <<
            //  equilibriumCoefficients_[1] << " " << equilibriumCoefficients_[2] << std::endl;
            //  std::cout << "derivativeEquilibriumCoefficients_ " << derivativeEquilibriumCoefficients_.transpose( ) << std::endl; 

            currentDeformation_ = ( 1.0 / globalRelaxationTime_ ) * ( 
                equilibriumCoefficients_ - nominalCoefficients_  + maxwellRelaxationTime_ * derivativeEquilibriumCoefficients_ ); // + ( 1.0 / 1.0e6 ) * ( equilibriumCoefficients_ - nominalCoefficients_ );
            // currentDeformation_[ 0 ] = 0.0;
            // currentDeformation_[ 1 ] = 0.0; 
            // currentDeformation_[ 2 ] = 0.0; 

            // currentDeformation_ = ( 1.0 / globalRelaxationTime_ ) * ( equilibriumCoefficients_ - nuCoefficients_ );

            // std::cout << currentLongitude_ * 180.0 / mathematical_constants::PI << " " << currentRelativePosition_.segment( 0, 3 ).norm( ) << " " <<
                //  equilibriumCoefficients_.transpose( ) << std::endl;
            // std::cout << "currentDeformation_: " << currentDeformation_.transpose( ) << std::endl;
            // std::cout << "nominalCoefficients_: " << nominalCoefficients_.transpose( ) << std::endl;

        // }
    }


    void updateEquilibriumDeformation( const double currentTime = TUDAT_NAN )
    {
        std::cout.precision(20);

        double relativeDistance = currentRelativePosition_.segment( 0, 3 ).norm( );
        // std::cout << "relativeDistance " << relativeDistance << std::endl;
        double radiusRatioPowerThree = referenceRadius_ * referenceRadius_ * referenceRadius_ / ( 
            relativeDistance * relativeDistance * relativeDistance );
        // std::cout << "radiusRatioPowerThree " << radiusRatioPowerThree << std::endl;
        // std::cout << "currentLongitude_ " << currentLongitude_ << std::endl;

        double gravitationalParametersRatio = gravitationalParameterPerturbingBody_ / gravitationalParameterDeformingBody_ ;
        // // std::cout << "gravitationalParametersRatio " << gravitationalParametersRatio << std::endl;

        // std::cout << "in update equilibrium " << currentTime << std::endl;
        // std::cout << "rotation " << rotationToIntegrationFrame_.toRotationMatrix( ) << std::endl;;  
        // std::cout << "position deforming " << positionOfDeformingBody_.transpose( ) << std::endl;
        // std::cout << "position perturbing " << positionOfPerturbingBody_.transpose( ) << std::endl;

        double rotationRate = angularVelocityDeformingBody_( ).norm( );
        double rotationRateDerivative = angularVelocityDerivativeDeformingBody_( )[ 2 ]; // SHOULD BE MODIFIED
        // std::cout << "rotationRateDerivative " << rotationRateDerivative << std::endl;

        // // Alternative computation of cos(2gamma), sin(2gamma)
        // Eigen::Vector3d e_A = ( Eigen::Vector3d( ) << 1.0, 0.0, 0.0 ).finished( );
        // Eigen::Vector3d e_B = ( Eigen::Vector3d( ) << 0.0, 1.0, 0.0 ).finished( );

        // Eigen::Vector3d bodyFixedUnitVector = currentRelativePosition_ / currentRelativePosition_.norm( );

        // double cos2gamma = e_A.dot( bodyFixedUnitVector ) * e_A.dot( bodyFixedUnitVector )
        //     - e_B.dot( bodyFixedUnitVector ) * e_B.dot( bodyFixedUnitVector );
        // double sin2gamma = - 2.0 * e_A.dot( bodyFixedUnitVector ) * e_B.dot( bodyFixedUnitVector );

        // std::cout << "cos2gamma " << cos2gamma << " " << std::cos( 2.0 * currentLongitude_ ) << 
        // " sin2gamma " << sin2gamma << " " << - std::sin( 2.0 * currentLongitude_ ) << std::endl;

        equilibriumCoefficients_[ 0 ] = k2_ * ( 
            // - rotationRate * rotationRate * referenceRadius_* referenceRadius_ * referenceRadius_ 
            //  / ( 3.0 * gravitationalParameterDeformingBody_ ) 
            + 0.5 * gravitationalParametersRatio * radiusRatioPowerThree 
            * ( 3.0 * std::sin( currentLatitude_ ) * std::sin( currentLatitude_ ) - 1.0 ) ); 
        equilibriumCoefficients_[ 2 ] = k2_ / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
            ( 1.0 - std::sin( currentLatitude_ ) * std::sin( currentLatitude_ ) ) * std::cos( 2.0 * currentLongitude_ );
        equilibriumCoefficients_[ 4 ] = k2_ / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
            ( 1.0 - std::sin( currentLatitude_ ) * std::sin( currentLatitude_ ) ) * std::sin( 2.0 * currentLongitude_ );
        if ( includeOrder1_ ) 
        {
            equilibriumCoefficients_[ 1 ] = - k2_ * gravitationalParametersRatio * radiusRatioPowerThree 
            * ( - std::cos( currentLatitude_ ) * std::sin( currentLatitude_ ) ) * std::cos( currentLongitude_ );
            equilibriumCoefficients_[ 3 ] = - k2_ * gravitationalParametersRatio * radiusRatioPowerThree 
            * ( - std::cos( currentLatitude_ ) * std::sin( currentLatitude_ ) ) * std::sin( currentLongitude_ );
        }


        // std::cout << "equilibriumCoefficients_ " << equilibriumCoefficients_.transpose( ) << std::endl;
        // std::cout << "position_norm " << relativeDistance << " radius " << referenceRadius_ << 
        // " mu ratio " << gravitationalParametersRatio << " cos 2 gamma " << std::cos( 2.0 * currentLongitude_ ) 
        // << " sin 2 gamma " << std::sin( 2.0 * currentLongitude_ ) << std::endl;
        // std::cout << "currentRelativePosition_ " << currentRelativePosition_.transpose( ) << std::endl;

        // equilibriumCoefficients_ += staticCoefficients_;

        double currentDistanceDerivative = 
            ( currentRelativePosition_[ 0 ] * currentRelativeVelocity_[ 0 ] 
            + currentRelativePosition_[ 1 ] * currentRelativeVelocity_[ 1 ] 
            + currentRelativePosition_[ 2 ] * currentRelativeVelocity_[ 2 ] ) / relativeDistance;

        derivativeEquilibriumCoefficients_[ 0 ] = - k2_ * ( 
            //2.0 * rotationRate * referenceRadius_ * referenceRadius_ * referenceRadius_ / ( 3.0 * gravitationalParameterDeformingBody_ ) * rotationRateDerivative
            + 1.0 / 2.0 * gravitationalParametersRatio * radiusRatioPowerThree 
            * 3.0 * currentDistanceDerivative / relativeDistance * ( 3.0 * std::sin( currentLatitude_ ) * std::sin( currentLatitude_ ) - 1.0 )
            + 1.0 / 2.0 * gravitationalParametersRatio * radiusRatioPowerThree 
            * ( 6.0 * currentLatitudeDerivative_ * std::sin( currentLatitude_ ) * std::cos( currentLatitude_ ) )  );

        derivativeEquilibriumCoefficients_[ 2 ] = - k2_ / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * (
            3.0 * currentDistanceDerivative / relativeDistance 
            * ( 1.0 - std::sin( currentLatitude_ ) * std::sin( currentLatitude_ ) ) * std::cos( 2.0 * currentLongitude_ )
            + 2.0 * currentLongitudeDerivative_ * std::sin( 2.0 * currentLongitude_ ) * ( 1.0 - std::sin( currentLatitude_ ) * std::sin( currentLatitude_ ) ) 
            + 2.0 * std::cos( 2.0 * currentLongitude_ ) * currentLatitudeDerivative_ * std::sin( currentLatitude_ ) * std::cos( currentLatitude_ ) );

        derivativeEquilibriumCoefficients_[ 4 ] = - k2_ / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * (
            3.0 * currentDistanceDerivative / relativeDistance
            * ( 1.0 - std::sin( currentLatitude_ ) * std::sin( currentLatitude_ ) ) * std::sin( 2.0 * currentLongitude_ )
            - 2.0 * currentLongitudeDerivative_ * std::cos( 2.0 * currentLongitude_ ) * ( 1.0 - std::sin( currentLatitude_ ) * std::sin( currentLatitude_ ) ) 
            + 2.0 * std::sin( 2.0 * currentLongitude_ ) * currentLatitudeDerivative_ * std::sin( currentLatitude_ ) * std::cos( currentLatitude_ ) );

        if ( includeOrder1_ ) 
        {
            derivativeEquilibriumCoefficients_[ 1 ] = - k2_ * gravitationalParametersRatio * radiusRatioPowerThree * (
            3.0 * currentDistanceDerivative / relativeDistance 
            * std::sin( currentLatitude_ ) * std::cos( currentLatitude_ ) * std::cos( currentLongitude_ )
            + currentLongitudeDerivative_ * std::sin( currentLatitude_ ) * std::cos( currentLatitude_ ) * std::sin( currentLongitude_ )
            - currentLatitudeDerivative_ * 
            ( std::cos( currentLatitude_ ) * std::cos( currentLatitude_ ) - std::sin( currentLatitude_ ) * std::sin( currentLatitude_ ) ) 
            * std::cos( currentLongitude_ ) );

            derivativeEquilibriumCoefficients_[ 3 ] = - k2_ * gravitationalParametersRatio * radiusRatioPowerThree * (
            3.0 * currentDistanceDerivative / relativeDistance 
            * std::sin( currentLatitude_ ) * std::cos( currentLatitude_ ) * std::sin( currentLongitude_ )
            - currentLongitudeDerivative_ * std::sin( currentLatitude_ ) * std::cos( currentLatitude_ ) * std::cos( currentLongitude_ )
            - currentLatitudeDerivative_ * 
            ( std::cos( currentLatitude_ ) * std::cos( currentLatitude_ ) - std::sin( currentLatitude_ ) * std::sin( currentLatitude_ ) ) 
            * std::sin( currentLongitude_ ) );
        }

        // std::cout << "equilibriumCoefficients_ " << equilibriumCoefficients_.transpose( ) << std::endl;
        // std::cout << "derivativeEquilibriumCoefficients_ " << derivativeEquilibriumCoefficients_.transpose( ) << std::endl;
    }



    Eigen::VectorXd getCurrentCoefficients( )
    {
        return nominalCoefficients_;
    }

    // Eigen::VectorXd computeNominalCoefficients( Eigen::VectorXd propagatedCoefficients )
    // {
    //     std::cout << "propagatedCoefficients " << propagatedCoefficients.transpose( ) << std::endl;
    //     std::cout << " in computeNominalCoefficients " << 
    //         ( ( 1 - maxwellRelaxationTime_ / globalRelaxationTime_ ) * propagatedCoefficients 
    //         + maxwellRelaxationTime_ / globalRelaxationTime_ * equilibriumCoefficients_ ).transpose( ) << std::endl;
    //     updateEquilibriumDeformation( );
    //     return ( 1 - maxwellRelaxationTime_ / globalRelaxationTime_ ) * propagatedCoefficients 
    //         + maxwellRelaxationTime_ / globalRelaxationTime_ * equilibriumCoefficients_;
    // }


    //! Function to retrieve the spherical harmonics cache for this acceleration.
    /*!
     *  Function to retrieve the spherical harmonics cache for this acceleration.
     *  \return Spherical harmonics cache for this acceleration
     */
    std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > getSphericalHarmonicsCache( )
    {
        return sphericalHarmonicsCache_;
    }

    //! Function to return current position vector from body exerting acceleration to body undergoing acceleration, in frame
    //! fixed to body undergoing acceleration
    /*!
     * Function to return current position vector from body exerting acceleration to body undergoing acceleration, in frame
     * fixed to bodyundergoing acceleration
     * \return Current position vector from body exerting acceleration to body undergoing acceleration, in frame
     * fixed to bodyundergoing acceleration
     */
    Eigen::Vector3d getCurrentRelativePosition( )
    {
        return currentRelativePosition_;
    }

    //! Function to return current position vector from body exerting acceleration to body undergoing acceleration, in inertial
    //! frame
    /*!
     * Function to return current position vector from body exerting acceleration to body undergoing acceleration, in inertial
     * frame
     * \return Current position vector from body exerting acceleration to body undergoing acceleration, in inertial frame
     */
    Eigen::Vector6d getCurrentInertialRelativeState( )
    {
        return currentInertialRelativeState_;
    }

    //! Function to retrieve the spherical harmonics reference radius.
    /*!
     *  Function to retrieve the spherical harmonics reference radius.
     *  \return Spherical harmonics reference radius.
     */
    double getReferenceRadius( )
    {
        return referenceRadius_;
    }

    //! Matrix of cosine coefficients.
    /*!
     * Matrix containing coefficients of cosine terms for spherical harmonics expansion.
     */
    CoefficientMatrixReturningFunction getCosineHarmonicCoefficientsFunction( )
    {
        return getCosineHarmonicsCoefficients;
    }

    //! Matrix of sine coefficients.
    /*!
     * Matrix containing coefficients of sine terms for spherical harmonics expansion.
     */
    CoefficientMatrixReturningFunction getSineHarmonicCoefficientsFunction( )
    {
        return getSineHarmonicsCoefficients;
    }

    //! Function to retrieve the current rotation from body-fixed frame to integration frame, in the form of a quaternion.
    /*!
     *  Function to retrieve the current rotation from body-fixed frame to integration frame, in the form of a quaternion.
     *  \return current rotation from body-fixed frame to integration frame, in the form of a quaternion.
     */
    Eigen::Quaterniond getCurrentRotationToIntegrationFrame( )
    {
        return rotationToIntegrationFrame_;
    }

    //! Function to retrieve the current rotation from body-fixed frame to integration frame, as a rotation matrix.
    /*!
     *  Function to retrieve the current rotation from body-fixed frame to integration frame, as a rotation matrix.
     *  \return current rotation from body-fixed frame to integration frame, as a rotation matrix.
     */
    Eigen::Matrix3d getCurrentRotationToIntegrationFrameMatrix( )
    {
        return rotationToIntegrationFrame_.toRotationMatrix( );
    }

    // //! Function to set whether each of the separate spherical harmonic terms should be saved
    // /*!
    //  * Function to set whether each of the separate spherical harmonic terms should be saved (in accelerationPerTerm_ member
    //  * variable of this class)
    //  * \param saveSphericalHarmonicTermsSeparately Boolean denoting whether each of the separate spherical harmonic terms should
    //  * be saved (in accelerationPerTerm_ member variable of this class
    //  */
    // void setSaveSphericalHarmonicTermsSeparately( const bool saveSphericalHarmonicTermsSeparately )
    // {
    //     saveSphericalHarmonicTermsSeparately_ = saveSphericalHarmonicTermsSeparately;
    // }

    // //! Function to retrieve the contributions of separate degrees/ordesr to the acceleration, concatenated in a single vector
    // /*!
    //  * Function to retrieve the contributions of specific separate degree/order to the acceleration, concatenated in a single
    //  * vector
    //  * \param coefficientIndices List of degree/order at which the contributions to the full acceleration are to be retrieved
    //  * \return Contributions of separate degrees/ordesr to the acceleration, concatenated in a single vector
    //  */
    // Eigen::VectorXd getConcatenatedAccelerationComponents( const std::vector< std::pair< int, int > >& coefficientIndices )
    // {
    //     if( !saveSphericalHarmonicTermsSeparately_ )
    //     {
    //         throw std::runtime_error( "Error when retrieving component accelerations from spherial harmonic acceleration, components not saved" );
    //     }

    //     Eigen::VectorXd returnVector = Eigen::VectorXd( 3 * coefficientIndices.size( ) );
    //     for( unsigned int i = 0; i < coefficientIndices.size( ); i++ )
    //     {
    //         if( accelerationPerTerm_.count( coefficientIndices.at( i ) ) != 0 )
    //         {
    //             returnVector.segment( i * 3, 3 ) = accelerationPerTerm_.at( coefficientIndices.at( i ) );
    //         }
    //         else
    //         {
    //             throw std::runtime_error( "Error when retrieving spherical harmonic acceleration at degree/order: " +
    //                                       std::to_string( coefficientIndices.at( i ).first ) + "/" +
    //                                       std::to_string( coefficientIndices.at( i ).second ) +
    //                                       ". This degree/order combination is not within the selected range of the current acceleration model." );
    //         }

    //     }
    //     return returnVector;
    // }

    // Eigen::VectorXd getConcatenatedAccelerationComponentNorms( const std::vector< std::pair< int, int > >& coefficientIndices )
    // {
    //     if( !saveSphericalHarmonicTermsSeparately_ )
    //     {
    //         throw std::runtime_error( "Error when retrieving component accelerations from spherial harmonic acceleration, components not saved" );
    //     }

    //     Eigen::VectorXd returnVector = Eigen::VectorXd( coefficientIndices.size( ) );
    //     for( unsigned int i = 0; i < coefficientIndices.size( ); i++ )
    //     {
    //         if( accelerationPerTerm_.count( coefficientIndices.at( i ) ) != 0 )
    //         {
    //             returnVector( i ) = accelerationPerTerm_.at( coefficientIndices.at( i ) ).norm( );
    //         }
    //         else
    //         {
    //             throw std::runtime_error( "Error when retrieving spherical harmonic acceleration at degree/order: " +
    //                                       std::to_string( coefficientIndices.at( i ).first ) + "/" +
    //                                       std::to_string( coefficientIndices.at( i ).second ) +
    //                                       ". This degree/order combination is not within the selected range of the current acceleration model." );
    //         }
    //     }
    //     return returnVector;
    // }

    //! Function to retrieve maximum degree of gravity field expansion
    /*!
     * Function to retrieve maximum degree of gravity field expansion
     * \return Maximum degree of gravity field expansion
     */
    int getMaximumDegree( )
    {
        return maximumDegree_;
    }

    //! Function to retrieve maximum order of gravity field expansion
    /*!
     * Function to retrieve maximum order of gravity field expansion
     * \return Maximum order of gravity field expansion
     */
    int getMaximumOrder( )
    {
        return maximumOrder_;
    }

    std::string getPerturbingBody( ) const
    {
        return perturbingBody_;
    }

    double getLoveNumber( ) const
    { 
        return k2_;
    }

    double getMaxwellRelaxationTime( ) const
    {
        return maxwellRelaxationTime_;
    } 

    double getGlobalRelaxationTime( ) const
    {
        return globalRelaxationTime_;
    } 

    double getGravitationalParameterDeformingBody( ) const
    {
        return gravitationalParameterDeformingBody_;
    } 

    double getGravitationalParameterPerturbingBody( ) const
    {
        return gravitationalParameterPerturbingBody_;
    } 

    double getReferenceRadius( ) const
    {  
        return referenceRadius_;
    }

    StateFunction getStateOfDeformingBodyFunction( )
    {
        return stateOfDeformingBodyFunction_;
    }

    StateFunction getStateOfPerturbingBodyFunction( )
    {
        return stateOfPerturbingBodyFunction_;
    }

    // void resetRotationRateDerivative( const double rotationRateDerivative )
    // {
    //     rotationRateDerivative_ = rotationRateDerivative;
    // }


protected:

private:

    Eigen::Vector6d stateOfDeformingBody_;

    const std::string perturbingBody_;

    Eigen::Vector6d stateOfPerturbingBody_;

    const double maxwellRelaxationTime_;

    const double globalRelaxationTime_;

    const double gravitationalParameterDeformingBody_;

    const double gravitationalParameterPerturbingBody_;

    const double referenceRadius_;

    //! Love number k2
    const double k2_;

    // double rotationRateDerivative_;

    Eigen::VectorXd equilibriumCoefficients_;

    Eigen::VectorXd derivativeEquilibriumCoefficients_;

    Eigen::VectorXd nominalCoefficients_;

    Eigen::VectorXd staticCoefficients_;

    //! Matrix of cosine coefficients.
    /*!
     * Matrix containing coefficients of cosine terms for spherical harmonics expansion.
     */
    Eigen::MatrixXd cosineHarmonicCoefficients;

    //! Matrix of sine coefficients.
    /*!
     * Matrix containing coefficients of sine terms for spherical harmonics expansion.
     */
    Eigen::MatrixXd sineHarmonicCoefficients;

    //! Pointer to function returning cosine harmonics coefficients matrix.
    /*!
     * Pointer to function that returns the current coefficients of the cosine terms of the
     * spherical harmonics expansion.
     */
    const CoefficientMatrixReturningFunction getCosineHarmonicsCoefficients;

    //! Pointer to function returning sine harmonics coefficients matrix.
    /*!
     * Pointer to function that returns the current coefficients of the sine terms of the
     * spherical harmonics expansion.
     */
    const CoefficientMatrixReturningFunction getSineHarmonicsCoefficients;

    //! Function returning the current rotation from body-fixed frame to integration frame.
    std::function< Eigen::Quaterniond( ) > rotationFromBodyFixedToIntegrationFrameFunction_;

    std::function< Eigen::Matrix3d( ) > rotationToBodyFixedDerivativeFunction_;

    //! Current rotation from body-fixed frame to integration frame.
    Eigen::Quaterniond rotationToIntegrationFrame_;

    //! Current position vector from body exerting acceleration to body undergoing acceleration, in frame fixed to body
    //! undergoing acceleration
    Eigen::Vector3d currentRelativePosition_;

    Eigen::Vector3d currentRelativeVelocity_;

    //! Current position vector from body exerting acceleration to body undergoing acceleration, in inertial frame
    Eigen::Vector3d currentInertialRelativePosition_;

    Eigen::Vector6d currentInertialRelativeState_;

    //!  Spherical harmonics cache for this acceleration
    std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache_;

    std::map< std::pair< int, int >, Eigen::VectorXd > deformationPerTerm_;

    //! List of contributions to accelerations at given degrees/orders, represented by first/second entry of map key pair.
    std::map< std::pair< int, int >, Eigen::Vector3d > accelerationPerTerm_;

    //! Boolean that denotes whether each of the separate spherical harmonic terms should be saved (in accelerationPerTerm_)
    bool saveSphericalHarmonicTermsSeparately_;

    //! Maximum degree of gravity field expansion
    int maximumDegree_;

    //! Maximum order of gravity field expansion
    int maximumOrder_;

    //! Function returning the state of the body undergoing deformation
    StateFunction stateOfDeformingBodyFunction_;

    //! Function returning the state of the body causing the deformation
    StateFunction stateOfPerturbingBodyFunction_;

    //! Current body-fixed longitude of the perturbing body
    double currentLongitude_;

    //! Current body-fixed latitude of the perturbing body
    double currentLatitude_;

    double currentLongitudeDerivative_;
    double currentLatitudeDerivative_;

    const bool includeOrder1_;

    std::function< Eigen::Vector3d( ) > angularVelocityDeformingBody_;
    std::function< Eigen::Vector3d( ) > angularVelocityDerivativeDeformingBody_;

};

} // namespace basic_astrodynamics
} // namespace tudat

#endif // TUDAT_GRAVITY_DEFORMATION_MODEL_H