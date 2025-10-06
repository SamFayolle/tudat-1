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

#include "tudat/astro/basic_astro/gravityDeformationModel.h"

namespace tudat
{

namespace basic_astrodynamics
{

void MaxwellGravityDeformationModel::updateMembers( const double currentTime )
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

        for ( unsigned int k = 0 ; k < perturbingBody_.size( ) ; k++ )
        {
            stateOfPerturbingBodyFunction_.at( k )( stateOfPerturbingBody_.at( k ) );

            // Compute relative inertial state
            currentInertialRelativeState_.at( k ) = stateOfPerturbingBody_.at( k ) - stateOfDeformingBody_;
        }
        // this->updateBaseMembers( );

        // // Compute relative inertial state
        // currentInertialRelativeState_ = stateOfPerturbingBody_ - stateOfDeformingBody_;
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


        for ( unsigned int k = 0 ; k < perturbingBody_.size( ) ; k++ )
        {
            // Compute current relative state in body-fixed frame
            currentRelativePosition_.at( k ) = rotationToIntegrationFrame_.inverse( ) * 
                currentInertialRelativeState_.at( k ).segment( 0, 3 );
            currentRelativeVelocity_.at( k ) = rotationToIntegrationFrame_.inverse( ) * 
                currentInertialRelativeState_.at( k ).segment( 3, 3 )
                + currentRotationToLocalFrameDerivative * currentInertialRelativeState_.at( k ).segment( 0, 3 );

            // std::cout << "rotationToIntegrationFrame_.inverse( ) " << ( rotationToIntegrationFrame_.inverse( ) ).toRotationMatrix() << std::endl;

            // Compute spherical coordinates of perturbing body in body-fixed frame
            Eigen::Vector3d currentSphericalPositionPerturbingBody = 
                coordinate_conversions::convertCartesianToSpherical( currentRelativePosition_.at( k ) );
            currentLongitude_.at( k ) = currentSphericalPositionPerturbingBody[ 2 ]; 
            currentLatitude_.at( k ) = mathematical_constants::PI / 2.0 - currentSphericalPositionPerturbingBody.y( );

            // MOVE DERIVATIVE CALCULATION TO UPDATE EQUILIBRIUM(?)
            // Compute current derivative of the perturbing body's body-fixed longitude
            currentLongitudeDerivative_.at( k ) = (
                ( currentRelativeVelocity_.at( k )[ 1 ] * currentRelativePosition_.at( k )[ 0 ] 
                - currentRelativeVelocity_.at( k )[ 0 ] * currentRelativePosition_.at( k )[ 1 ] ) 
                / ( currentRelativePosition_.at( k )[ 0 ] * currentRelativePosition_.at( k )[ 0 ] + currentRelativePosition_.at( k )[ 1 ] * currentRelativePosition_.at( k )[ 1 ] ) );

            // Compute current derivative of the perturbing body's body-fixed latitude
            double currentDistance = currentRelativePosition_.at( k ).segment( 0, 3 ).norm( );
            double currentDistanceDerivative = 
                ( currentRelativePosition_.at( k )[ 0 ] * currentRelativeVelocity_.at( k )[ 0 ] 
                + currentRelativePosition_.at( k )[ 1 ] * currentRelativeVelocity_.at( k )[ 1 ] 
                + currentRelativePosition_.at( k )[ 2 ] * currentRelativeVelocity_.at( k )[ 2 ] ) / currentDistance;
            currentLatitudeDerivative_.at( k ) =
                ( currentRelativeVelocity_.at( k )[2] * currentDistance - currentRelativePosition_.at( k )[2] * currentDistanceDerivative ) /
                ( currentDistance * std::sqrt( currentRelativePosition_.at( k )[ 0 ] * currentRelativePosition_.at( k )[ 0 ] + currentRelativePosition_.at( k )[ 1 ] * currentRelativePosition_.at( k )[ 1 ] ) );
        }
        

        //  if (negativeSignLatitude_)
        //  {
        //     currentLatitudeDerivative_ = - currentLatitudeDerivative_;
        //     currentLatitude_ = - currentLatitude_;
        //  }   
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


void MaxwellGravityDeformationModel::updateEquilibriumDeformation( const double currentTime )
{
    std::cout.precision(20);

    // Reset equilibrium coefficients and derivatives to zero
    equilibriumCoefficients_ = Eigen::VectorXd::Zero( 5 ); 
    derivativeEquilibriumCoefficients_ = Eigen::VectorXd::Zero( 5 );

    double rotationRate = angularVelocityDeformingBody_( ).norm( );
    double rotationRateDerivative = angularVelocityDerivativeDeformingBody_( )[ 2 ]; // SHOULD BE MODIFIED
    // std::cout << "rotationRateDerivative " << rotationRateDerivative << std::endl;

    for ( unsigned int k = 0 ; k < perturbingBody_.size( ) ; k++ )
    {
         double relativeDistance = currentRelativePosition_.at( k ).segment( 0, 3 ).norm( );
        // std::cout << "relativeDistance " << relativeDistance << std::endl;
        double radiusRatioPowerThree = referenceRadius_ * referenceRadius_ * referenceRadius_ / ( 
            relativeDistance * relativeDistance * relativeDistance );
        // std::cout << "radiusRatioPowerThree " << radiusRatioPowerThree << std::endl;
        // std::cout << "currentLongitude_ " << currentLongitude_ << std::endl;

        double gravitationalParametersRatio = gravitationalParameterPerturbingBody_.at( k ) / gravitationalParameterDeformingBody_ ;
        // // std::cout << "gravitationalParametersRatio " << gravitationalParametersRatio << std::endl;

        // std::cout << "in update equilibrium " << currentTime << std::endl;
        // std::cout << "rotation " << rotationToIntegrationFrame_.toRotationMatrix( ) << std::endl;;  
        // std::cout << "position deforming " << positionOfDeformingBody_.transpose( ) << std::endl;
        // std::cout << "position perturbing " << positionOfPerturbingBody_.transpose( ) << std::endl;

        

        // // Alternative computation of cos(2gamma), sin(2gamma)
        // Eigen::Vector3d e_A = ( Eigen::Vector3d( ) << 1.0, 0.0, 0.0 ).finished( );
        // Eigen::Vector3d e_B = ( Eigen::Vector3d( ) << 0.0, 1.0, 0.0 ).finished( );

        // Eigen::Vector3d bodyFixedUnitVector = currentRelativePosition_ / currentRelativePosition_.norm( );

        // double cos2gamma = e_A.dot( bodyFixedUnitVector ) * e_A.dot( bodyFixedUnitVector )
        //     - e_B.dot( bodyFixedUnitVector ) * e_B.dot( bodyFixedUnitVector );
        // double sin2gamma = - 2.0 * e_A.dot( bodyFixedUnitVector ) * e_B.dot( bodyFixedUnitVector );

        // std::cout << "cos2gamma " << cos2gamma << " " << std::cos( 2.0 * currentLongitude_ ) << 
        // " sin2gamma " << sin2gamma << " " << - std::sin( 2.0 * currentLongitude_ ) << std::endl;

        equilibriumCoefficients_[ 0 ] += k2_ * ( 
            // - rotationRate * rotationRate * referenceRadius_* referenceRadius_ * referenceRadius_ 
            //     / ( 3.0 * gravitationalParameterDeformingBody_ ) 
            + 0.5 * gravitationalParametersRatio * radiusRatioPowerThree 
            * ( 3.0 * std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) - 1.0 ) ); 
        equilibriumCoefficients_[ 2 ] += k2_ / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
            ( 1.0 - std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) * std::cos( 2.0 * currentLongitude_.at( k ) );
        equilibriumCoefficients_[ 4 ] += k2_ / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
            ( 1.0 - std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) * std::sin( 2.0 * currentLongitude_.at( k ) );
        if ( includeOrder1_ ) 
        {
            equilibriumCoefficients_[ 1 ] += - k2_ * gravitationalParametersRatio * radiusRatioPowerThree 
            * ( - std::cos( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) * std::cos( currentLongitude_.at( k ) );
            equilibriumCoefficients_[ 3 ] += - k2_ * gravitationalParametersRatio * radiusRatioPowerThree 
            * ( - std::cos( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) * std::sin( currentLongitude_.at( k ) );
        }


        // std::cout << "equilibriumCoefficients_ " << equilibriumCoefficients_.transpose( ) << std::endl;
        // std::cout << "position_norm " << relativeDistance << " radius " << referenceRadius_ << 
        // " mu ratio " << gravitationalParametersRatio << " cos 2 gamma " << std::cos( 2.0 * currentLongitude_ ) 
        // << " sin 2 gamma " << std::sin( 2.0 * currentLongitude_ ) << std::endl;
        // std::cout << "currentRelativePosition_ " << currentRelativePosition_.transpose( ) << std::endl;

        // equilibriumCoefficients_ += staticCoefficients_;

        double currentDistanceDerivative = 
            ( currentRelativePosition_.at( k )[ 0 ] * currentRelativeVelocity_.at( k )[ 0 ] 
            + currentRelativePosition_.at( k )[ 1 ] * currentRelativeVelocity_.at( k )[ 1 ] 
            + currentRelativePosition_.at( k )[ 2 ] * currentRelativeVelocity_.at( k )[ 2 ] ) / relativeDistance;

        derivativeEquilibriumCoefficients_[ 0 ] += - k2_ * ( 
            //2.0 * rotationRate * referenceRadius_ * referenceRadius_ * referenceRadius_ / ( 3.0 * gravitationalParameterDeformingBody_ ) * rotationRateDerivative
            + 1.0 / 2.0 * gravitationalParametersRatio * radiusRatioPowerThree 
            * 3.0 * currentDistanceDerivative / relativeDistance * ( 3.0 * std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) - 1.0 )
            - 1.0 / 2.0 * gravitationalParametersRatio * radiusRatioPowerThree 
            * ( 6.0 * currentLatitudeDerivative_.at( k ) * std::sin( currentLatitude_.at( k ) ) * std::cos( currentLatitude_.at( k ) ) )  );

        derivativeEquilibriumCoefficients_[ 2 ] += - k2_ / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * (
            3.0 * currentDistanceDerivative / relativeDistance 
            * ( 1.0 - std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) * std::cos( 2.0 * currentLongitude_.at( k ) )
            + 2.0 * currentLongitudeDerivative_.at( k ) * std::sin( 2.0 * currentLongitude_.at( k ) ) * ( 1.0 - std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) 
            + 2.0 * std::cos( 2.0 * currentLongitude_.at( k ) ) * currentLatitudeDerivative_.at( k ) * std::sin( currentLatitude_.at( k ) ) * std::cos( currentLatitude_.at( k ) ) );

        derivativeEquilibriumCoefficients_[ 4 ] += - k2_ / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * (
            3.0 * currentDistanceDerivative / relativeDistance
            * ( 1.0 - std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) * std::sin( 2.0 * currentLongitude_.at( k ) )
            - 2.0 * currentLongitudeDerivative_.at( k ) * std::cos( 2.0 * currentLongitude_.at( k ) ) * ( 1.0 - std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) 
            + 2.0 * std::sin( 2.0 * currentLongitude_.at( k ) ) * currentLatitudeDerivative_.at( k ) * std::sin( currentLatitude_.at( k ) ) * std::cos( currentLatitude_.at( k ) ) );

        if ( includeOrder1_ ) 
        {
            derivativeEquilibriumCoefficients_[ 1 ] += - k2_ * gravitationalParametersRatio * radiusRatioPowerThree * (
            3.0 * currentDistanceDerivative / relativeDistance 
            * std::sin( currentLatitude_.at( k ) ) * std::cos( currentLatitude_.at( k ) ) * std::cos( currentLongitude_.at( k ) )
            + currentLongitudeDerivative_.at( k ) * std::sin( currentLatitude_.at( k ) ) * std::cos( currentLatitude_.at( k ) ) * std::sin( currentLongitude_.at( k ) )
            - currentLatitudeDerivative_.at( k ) * 
            ( std::cos( currentLatitude_.at( k ) ) * std::cos( currentLatitude_.at( k ) ) - std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) 
            * std::cos( currentLongitude_.at( k ) ) );

            derivativeEquilibriumCoefficients_[ 3 ] += - k2_ * gravitationalParametersRatio * radiusRatioPowerThree * (
            3.0 * currentDistanceDerivative / relativeDistance 
            * std::sin( currentLatitude_.at( k ) ) * std::cos( currentLatitude_.at( k ) ) * std::sin( currentLongitude_.at( k ) )
            - currentLongitudeDerivative_.at( k ) * std::sin( currentLatitude_.at( k ) ) * std::cos( currentLatitude_.at( k ) ) * std::cos( currentLongitude_.at( k ) )
            - currentLatitudeDerivative_.at( k ) * 
            ( std::cos( currentLatitude_.at( k ) ) * std::cos( currentLatitude_.at( k ) ) - std::sin( currentLatitude_.at( k ) ) * std::sin( currentLatitude_.at( k ) ) ) 
            * std::sin( currentLongitude_.at( k ) ) );

            // if (negativeSignLatitude_)
            // {
            //     derivativeEquilibriumCoefficients_[ 1 ] = - derivativeEquilibriumCoefficients_[ 1 ];
            //     derivativeEquilibriumCoefficients_[ 3 ] = - derivativeEquilibriumCoefficients_[ 3 ];
            // }
        }
    }

   

    // std::cout << "equilibriumCoefficients_ " << equilibriumCoefficients_.transpose( ) << std::endl;
    // std::cout << "derivativeEquilibriumCoefficients_ " << derivativeEquilibriumCoefficients_.transpose( ) << std::endl;
}

} // namespace basic_astrodynamics

} // namespace tudat