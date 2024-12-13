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

#include <Eigen/Core>

#include "tudat/math/basic/mathematicalConstants.h"

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
        currentDeformation_( Eigen::VectorXd::Constant( TUDAT_NAN ) ){ }

    //! Virtual destructor.
    /*!
     * Virtual destructor, necessary to ensure that derived class destructors get called correctly.
     */
    virtual ~GravityDeformationModel( ) { }

//    //! Get acceleration.
//    /*!
//     * Returns the acceleration. No arguments are passed to this function for generality.
//     * Instead, all data required for computation is to be obtained from pointers to functions/
//     * classes/structs, etc which are to be set in a derived class and evaluated by the
//     * updateMembers() function below.
//     * \return Acceleration.
//     * \sa updateMembers().
//     */
//    virtual AccelerationDataType getAcceleration( ) = 0;

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


//! Update the members of a gravity deformation model and evaluate the deformation.
/*!
 * Updates the member variables of a gravity deformation model and subsequently evaluates the
 * deformation. 
 * \param deformationModel Gravity deformation model that is to be evaluated.
 * \param currentTime Time at which the gravity deformation model is to be updated.
 * \return Deformation that is obtained following the member update.
 */
Eigen::VectorXd updateAndGetDeformation(
        const std::shared_ptr< GravityDeformationModel > deformationModel,
        const double currentTime = TUDAT_NAN )
{
    // Update members.
    deformationModel->updateMembers( currentTime );

    // Evaluate and return deformation.
    return deformationModel->getDeformation( );
}


// //! Typedef defining a list of accelerations acting on a single body, key is the name of each
// //! body exerting a acceletation, value is a list of accelerations exerted by that body.
// typedef std::unordered_map< std::string, std::vector<
// std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > > >
// SingleBodyAccelerationMap;


// //! Typedef defining a list of accelerations acting on a set of bodies, key is the name of each
// //! body undergoing an acceletation, value is SingleBodyAccelerationMap, defining all accelerations
// //! acting on it.
// typedef std::unordered_map< std::string, SingleBodyAccelerationMap > AccelerationMap;


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
            const StateFunction positionOfDeformingBodyFunction,
            const double maxwellRelaxationTime,
            const double globalRelaxationTime,
            const double gravitationalParameterDeformingBody,
            const double gravitationalParamemeterPerturbingBody,
            const double referenceRadius,
            const double rotationRate,
            const double k2,
            Eigen::VectorXd& transientCosineCoefficients,
            Eigen::VectorXd& transientSineCoefficients,
            const StateFunction positionOfPerturbingBodyFunction =
            [ ]( Eigen::Vector3d& input ){ input = Eigen::Vector3d::Zero( ); },
            const std::function< Eigen::Quaterniond( ) >
            rotationFromBodyFixedToIntegrationFrameFunction =
            [ ]( ){ return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); } ) //,
            /*std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache =
            std::make_shared< basic_mathematics::SphericalHarmonicsCache >( ) )*/ : 
            GravityDeformationModel( ),
            positionOfDeformingBodyFunction_( positionOfDeformingBodyFunction ),
            maxwellRelaxationTime_( maxwellRelaxationTime ),
            globalRelaxationTime_( globalRelaxationTime ),
            gravitationalParameterDeformingBody_( gravitationalParameterDeformingBody ),
            gravitationalParamemeterPerturbingBody_( gravitationalParameterPerturbingBody ),
            referenceRadius_( referenceRadius ),
            rotationRate_( rotationRate ),
            k2_( k2 ),
          getCosineHarmonicsCoefficients( [ = ]( ){ return aCosineHarmonicCoefficientMatrix; } ),
          getSineHarmonicsCoefficients( [ = ]( ){ return aSineHarmonicCoefficientMatrix; } ),
          positionOfPerturbingBodyFunction_( positionOfPerturbingBodyFunction ),
          rotationFromBodyFixedToIntegrationFrameFunction_(
              rotationFromBodyFixedToIntegrationFrameFunction ),
          sphericalHarmonicsCache_( sphericalHarmonicsCache ),
          saveSphericalHarmonicTermsSeparately_( false )
    {
        // maximumDegree_ = static_cast< int >( getCosineHarmonicsCoefficients( ).rows( ) ) - 1 ;
        // maximumOrder_ = static_cast< int >( getCosineHarmonicsCoefficients( ).cols( ) )- 1 ;
        // sphericalHarmonicsCache_->resetMaximumDegreeAndOrder(
        //             std::max< int >( maximumDegree_,
        //                              sphericalHarmonicsCache_->getMaximumDegree( ) ) + 1,
        //             std::max< int >( maximumOrder_,
        //                              sphericalHarmonicsCache_->getMaximumOrder( ) ) + 1 );
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
        if( !( this->currentTime_ == currentTime ) )
        {

            // NEED FOR THIS
            cosineHarmonicCoefficients = getCosineHarmonicsCoefficients( );
            sineHarmonicCoefficients = getSineHarmonicsCoefficients( );

            rotationToIntegrationFrame_ = rotationFromBodyFixedToIntegrationFrameFunction_( );
            // this->updateBaseMembers( );

            currentInertialRelativePosition_ = positionOfDeformingBodyFunction_( ) - positionOfPerturbingBodyFunction_( );
                    // this->positionOfBodySubjectToAcceleration - this->positionOfBodyExertingAcceleration ;

            currentRelativePosition_ = rotationToIntegrationFrame_.inverse( ) * ( currentInertialRelativePosition_ );

            updateEquilibriumDeformation();

            // THE CURRENT COEFFICIENTS SHOULD BE UPDATED AT THIS POINT

            // Current (transient) deformation
            currentDeformation_ = ( 1.0 / globalRelaxationTime_ ) * ( equilibriumCoefficients_ - currentCoefficients_ );


            

            //         computeGeodesyNormalizedGravitationalAccelerationSum(
            //             currentRelativePosition_,
            //             gravitationalParameter,
            //             equatorialRadius,
            //             cosineHarmonicCoefficients,
            //             sineHarmonicCoefficients, sphericalHarmonicsCache_,
            //             accelerationPerTerm_,
            //             saveSphericalHarmonicTermsSeparately_,
            //             rotationToIntegrationFrame_.toRotationMatrix( ) );
            // currentAccelerationInBodyFixedFrame_ = rotationToIntegrationFrame_.inverse( ) * currentAcceleration_;

            // if ( this->updatePotential_ )
            // {
            //     this->currentPotential_ = gravitation::calculateSphericalHarmonicGravitationalPotential(
            //             currentRelativePosition_,
            //             gravitationalParameter,
            //             equatorialRadius,
            //             cosineHarmonicCoefficients,
            //             sineHarmonicCoefficients,
            //             sphericalHarmonicsCache_ );
            // }
        }
    }

    void updateEquilibriumDeformation( const double currentTime = TUDAT_NAN )
    {
        double relativeDistance = currentRelativePosition_.norm( );
        double radiusRatioPowerThree = referenceRadius_ * referenceRadius_ * referenceRadius_ / ( 
            relativeDistance * relativeDistance * relativeDistance );

        double gravitationalParametersRatio = gravitationalParameterPerturbingBody_ / gravitationalParameterDeformingBody_ ;

        equilibriumCoefficients_[ 0 ] = - k2_ * ( 
            rotationRate_ * rotationRate_ * referenceRadius_* referenceRadius_ * referenceRadius_ ) 
            / ( ( 3.0 * gravitationalParameterDeformingBody_ )
            + 0.5 * gravitationalParametersRatio * radiusRatioPowerThree );
        equilibriumCoefficients_[ 1 ] = k2_ / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
            std::cos( 2.0 * currentLongitude_ );
        equilibriumCoefficients_[ 2 ] = - k2_ / 4.0 * gravitationalParametersRatio * radiusRatioPowerThree * 
            std::sin( 2.0 * currentLongitude_ );
    }

    void computeFullCoefficients( )
    {
        currentFullCoefficients_ = ( 1.0 - maxwellRelaxationTime_ / globalRelaxationTime_ ) * currentCoefficients_
        + ( maxwellRelaxationTime_ / globalRelaxationTime_ ) * equilibriumCoefficients_;
    }

    Eigen::VectorXd getFullCoefficients( )
    {
        return currentFullCoefficients_;
    }

    
   

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
    Eigen::Vector3d getCurrentInertialRelativePosition( )
    {
        return currentInertialRelativePosition_;
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

    //! Function to set whether each of the separate spherical harmonic terms should be saved
    /*!
     * Function to set whether each of the separate spherical harmonic terms should be saved (in accelerationPerTerm_ member
     * variable of this class)
     * \param saveSphericalHarmonicTermsSeparately Boolean denoting whether each of the separate spherical harmonic terms should
     * be saved (in accelerationPerTerm_ member variable of this class
     */
    void setSaveSphericalHarmonicTermsSeparately( const bool saveSphericalHarmonicTermsSeparately )
    {
        saveSphericalHarmonicTermsSeparately_ = saveSphericalHarmonicTermsSeparately;
    }

    //! Function to retrieve the contributions of separate degrees/ordesr to the acceleration, concatenated in a single vector
    /*!
     * Function to retrieve the contributions of specific separate degree/order to the acceleration, concatenated in a single
     * vector
     * \param coefficientIndices List of degree/order at which the contributions to the full acceleration are to be retrieved
     * \return Contributions of separate degrees/ordesr to the acceleration, concatenated in a single vector
     */
    Eigen::VectorXd getConcatenatedAccelerationComponents( const std::vector< std::pair< int, int > >& coefficientIndices )
    {
        if( !saveSphericalHarmonicTermsSeparately_ )
        {
            throw std::runtime_error( "Error when retrieving component accelerations from spherial harmonic acceleration, components not saved" );
        }

        Eigen::VectorXd returnVector = Eigen::VectorXd( 3 * coefficientIndices.size( ) );
        for( unsigned int i = 0; i < coefficientIndices.size( ); i++ )
        {
            if( accelerationPerTerm_.count( coefficientIndices.at( i ) ) != 0 )
            {
                returnVector.segment( i * 3, 3 ) = accelerationPerTerm_.at( coefficientIndices.at( i ) );
            }
            else
            {
                throw std::runtime_error( "Error when retrieving spherical harmonic acceleration at degree/order: " +
                                          std::to_string( coefficientIndices.at( i ).first ) + "/" +
                                          std::to_string( coefficientIndices.at( i ).second ) +
                                          ". This degree/order combination is not within the selected range of the current acceleration model." );
            }

        }
        return returnVector;
    }

    Eigen::VectorXd getConcatenatedAccelerationComponentNorms( const std::vector< std::pair< int, int > >& coefficientIndices )
    {
        if( !saveSphericalHarmonicTermsSeparately_ )
        {
            throw std::runtime_error( "Error when retrieving component accelerations from spherial harmonic acceleration, components not saved" );
        }

        Eigen::VectorXd returnVector = Eigen::VectorXd( coefficientIndices.size( ) );
        for( unsigned int i = 0; i < coefficientIndices.size( ); i++ )
        {
            if( accelerationPerTerm_.count( coefficientIndices.at( i ) ) != 0 )
            {
                returnVector( i ) = accelerationPerTerm_.at( coefficientIndices.at( i ) ).norm( );
            }
            else
            {
                throw std::runtime_error( "Error when retrieving spherical harmonic acceleration at degree/order: " +
                                          std::to_string( coefficientIndices.at( i ).first ) + "/" +
                                          std::to_string( coefficientIndices.at( i ).second ) +
                                          ". This degree/order combination is not within the selected range of the current acceleration model." );
            }
        }
        return returnVector;
    }

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


protected:

private:

    const double maxwellRelaxationTime_;

    const double globalRelaxationTime_;
            
    const double gravitationalParameterDeformingBody_
            
    const double gravitationalParameterPerturbingBody_;
            
    const double referenceRadius_;

    const double rotationRate_;

    //! Love number k2
    const double k2_;

    Eigen::VectorXd equilibriumCoefficients_;

    Eigen::VectorXd currentCoefficients_;

    Eigen::VectorXd currentFullCoefficients_;

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

    //! Current rotation from body-fixed frame to integration frame.
    Eigen::Quaterniond rotationToIntegrationFrame_;

    //! Current position vector from body exerting acceleration to body undergoing acceleration, in frame fixed to body
    //! undergoing acceleration
    Eigen::Vector3d currentRelativePosition_;

    //! Current position vector from body exerting acceleration to body undergoing acceleration, in inertial frame
    Eigen::Vector3d currentInertialRelativePosition_;

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
    StateFunction positionOfDeformingBodyFunction_;

    //! Function returning the state of the body causing the deformation
    StateFunction positionOfPerturbingBodyFunction_;

    //! Current body-fixed longitude of the perturbing body
    double currentLongitude_;

};

} // namespace basic_astrodynamics
} // namespace tudat

#endif // TUDAT_GRAVITY_DEFORMATION_MODEL_H
