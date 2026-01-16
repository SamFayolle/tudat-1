/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/rotational_dynamics_partials/inertialTorquePartial.h"
#include "tudat/math/basic/legendrePolynomials.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function for calculating the partial of the torque w.r.t. a non-rotational integrated state
void InertialTorquePartial::wrtNonRotationalStateOfAdditionalBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType )
{
    if( stateReferencePoint.first == bodyUndergoingTorque_  && integratedStateType == propagators::gravity_deformation_state )
    {
        wrtGravityDeformation( partialMatrix );
    }
}

Eigen::MatrixXd InertialTorquePartial::wrtOtherStateDerivative()
{
    Eigen::MatrixXd partialWrtDeformationDerivative = Eigen::MatrixXd::Zero( 3, 5 );

    // std::cout << "currentInertiaTensorNormalizationFactor_ " << currentInertiaTensorNormalizationFactor_ << std::endl;

    Eigen::VectorXd partialWrtDerivativeC20 = - currentInertiaTensorNormalizationFactor_ * UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C20 * currentAngularVelocityVector_;
    Eigen::VectorXd partialWrtDerivativeC21 = - currentInertiaTensorNormalizationFactor_ * UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C21 * currentAngularVelocityVector_;
    Eigen::VectorXd partialWrtDerivativeC22 = - currentInertiaTensorNormalizationFactor_ * UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C22 * currentAngularVelocityVector_;
    Eigen::VectorXd partialWrtDerivativeS21 = - currentInertiaTensorNormalizationFactor_ * UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S21 * currentAngularVelocityVector_;
    Eigen::VectorXd partialWrtDerivativeS22 = - currentInertiaTensorNormalizationFactor_ * UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S22 * currentAngularVelocityVector_;

    partialWrtDeformationDerivative.block(0, 0, 3, 1) = partialWrtDerivativeC20;
    partialWrtDeformationDerivative.block(0, 1, 3, 1) = partialWrtDerivativeC21;
    partialWrtDeformationDerivative.block(0, 2, 3, 1) = partialWrtDerivativeC22;
    partialWrtDeformationDerivative.block(0, 3, 3, 1) = partialWrtDerivativeS21;
    partialWrtDeformationDerivative.block(0, 4, 3, 1) = partialWrtDerivativeS22;

    // std::cout << "partialWrtDeformationDerivative" << std::endl;
    // std::cout << partialWrtDeformationDerivative << std::endl;

    partialWrtDeformationDerivative = currentInertiaTensor_.inverse( ) * partialWrtDeformationDerivative;

    return partialWrtDeformationDerivative;
}

//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > InertialTorquePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    using namespace estimatable_parameters;

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunction =
            std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );

    if( parameter->getParameterName( ).second.first == bodyUndergoingTorque_ )
    {
        switch( parameter->getParameterName( ).first )
        {
            case gravitational_parameter: {
                if( bodyGravitationalParameterFunction_ == nullptr )
                {
                    throw std::runtime_error(
                            "Error when getting partial of inertial torque w.r.t. gravitational parameter, gravitational parameter "
                            "function not found." );
                }
                partialFunction =
                        std::make_pair( std::bind( &InertialTorquePartial::wrtGravitationalParameter, this, std::placeholders::_1 ), 1 );

                break;
            }
            case mean_moment_of_inertia: {
                if( getInertiaTensorNormalizationFactor_ == nullptr )
                {
                    throw std::runtime_error(
                            "Error when getting partial of inertial torque w.r.t. mean moment of inertia, inertia tensor normalization "
                            "function not found." );
                }
                partialFunction =
                        std::make_pair( std::bind( &InertialTorquePartial::wrtMeanMomentOfInertia, this, std::placeholders::_1 ), 1 );

                break;
            }
            default:
                break;
        }
    }
    return partialFunction;
}

//! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > InertialTorquePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    using namespace estimatable_parameters;

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunction =
            std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );

    if( parameter->getParameterName( ).second.first == bodyUndergoingTorque_ )
    {
        switch( parameter->getParameterName( ).first )
        {
            case spherical_harmonics_cosine_coefficient_block: {
                // Cast parameter object to required type.
                std::shared_ptr< SphericalHarmonicsCosineCoefficients > coefficientsParameter =
                        std::dynamic_pointer_cast< SphericalHarmonicsCosineCoefficients >( parameter );

                int c20Index, c21Index, c22Index;
                coefficientsParameter->getDegreeTwoEntries( c20Index, c21Index, c22Index );

                if( c20Index >= 0 || c21Index >= 0 || c22Index >= 0 )
                {
                    if( getInertiaTensorNormalizationFactor_ == nullptr )
                    {
                        throw std::runtime_error(
                                "Error when getting partial of inertial torque w.r.t. cosine sh parameters, inertia tensor normalization "
                                "function not found." );
                    }
                    partialFunction =
                            std::make_pair( std::bind( &InertialTorquePartial::wrtCosineSphericalHarmonicCoefficientsOfCentralBody,
                                                       this,
                                                       std::placeholders::_1,
                                                       c20Index,
                                                       c21Index,
                                                       c22Index ),
                                            coefficientsParameter->getParameterSize( ) );
                }

                break;
            }
            case spherical_harmonics_sine_coefficient_block: {
                // Cast parameter object to required type.

                std::shared_ptr< SphericalHarmonicsSineCoefficients > coefficientsParameter =
                        std::dynamic_pointer_cast< SphericalHarmonicsSineCoefficients >( parameter );

                int s21Index, s22Index;
                coefficientsParameter->getDegreeTwoEntries( s21Index, s22Index );

                if( s21Index >= 0 || s22Index >= 0 )
                {
                    if( ( getInertiaTensorNormalizationFactor_ == nullptr ) )
                    {
                        throw std::runtime_error(
                                "Error when getting partial of inertial torque w.r.t. sine sh parameters, inertia tensor normalization "
                                "function not found." );
                    }
                    partialFunction = std::make_pair( std::bind( &InertialTorquePartial::wrtSineSphericalHarmonicCoefficientsOfCentralBody,
                                                                 this,
                                                                 std::placeholders::_1,
                                                                 s21Index,
                                                                 s22Index ),
                                                      coefficientsParameter->getParameterSize( ) );
                }

                break;
            }
            default:
                break;
        }
    }
    return partialFunction;
}

//! Function to compute partial of torque w.r.t. mean moment of inertia
void InertialTorquePartial::wrtMeanMomentOfInertia( Eigen::MatrixXd& momentOfInertiaPartial )
{
    momentOfInertiaPartial.block( 0, 0, 3, 1 ) = -currentAngularVelocityCrossProductMatrix_ * currentInertiaTensorNormalizationFactor_ *
            UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_MEAN_MOMENT * currentAngularVelocityVector_;
}

//! Function to compute partial of torque w.r.t. gravitational parameter
void InertialTorquePartial::wrtGravitationalParameter( Eigen::MatrixXd& gravitationalParameterPartial )
{
    gravitationalParameterPartial.block( 0, 0, 3, 1 ) = -currentAngularVelocityCrossProductMatrix_ * currentInertiaTensor_ *
            currentAngularVelocityVector_ / currentGravitationalParameter_;
}

//! Function to compute partial of torque w.r.t. spherical harmonic cosine coefficients
void InertialTorquePartial::wrtCosineSphericalHarmonicCoefficientsOfCentralBody( Eigen::MatrixXd& sphericalHarmonicCoefficientPartial,
                                                                                 const int c20Index,
                                                                                 const int c21Index,
                                                                                 const int c22Index )
{
    sphericalHarmonicCoefficientPartial.setZero( );

    if( c20Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial.block( 0, c20Index, 3, 1 ) =
                -basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 ) * currentAngularVelocityCrossProductMatrix_ *
                currentInertiaTensorNormalizationFactor_ * UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C20 * currentAngularVelocityVector_;
    }

    if( c21Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial.block( 0, c21Index, 3, 1 ) =
                -basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 ) * currentAngularVelocityCrossProductMatrix_ *
                currentInertiaTensorNormalizationFactor_ * UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C21 * currentAngularVelocityVector_;
    }

    if( c22Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial.block( 0, c22Index, 3, 1 ) =
                -basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 ) * currentAngularVelocityCrossProductMatrix_ *
                currentInertiaTensorNormalizationFactor_ * UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C22 * currentAngularVelocityVector_;
    }
}

//! Function to compute partial of torque w.r.t. spherical harmonic sine coefficients
void InertialTorquePartial::wrtSineSphericalHarmonicCoefficientsOfCentralBody( Eigen::MatrixXd& sphericalHarmonicCoefficientPartial,
                                                                               const int s21Index,
                                                                               const int s22Index )
{
    sphericalHarmonicCoefficientPartial.setZero( );

    if( s21Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial.block( 0, s21Index, 3, 1 ) =
                -basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 ) * currentAngularVelocityCrossProductMatrix_ *
                currentInertiaTensorNormalizationFactor_ * UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S21 * currentAngularVelocityVector_;
    }

    if( s22Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial.block( 0, s22Index, 3, 1 ) =
                -basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 ) * currentAngularVelocityCrossProductMatrix_ *
                currentInertiaTensorNormalizationFactor_ * UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S22 * currentAngularVelocityVector_;
    }
}

//! Function to compute partial of torque w.r.t. gravity field deformation
void InertialTorquePartial::wrtGravityDeformation(
        Eigen::Block< Eigen::MatrixXd >& deformationPartial )
{
    // deformationPartial.setZero( );

    Eigen::MatrixXd cosinePartials = Eigen::MatrixXd::Zero( 3, 3 );
    wrtCosineSphericalHarmonicCoefficientsOfCentralBody( cosinePartials, 0, 1, 2 );
    Eigen::MatrixXd sinePartials = Eigen::MatrixXd::Zero( 3, 2 );
    wrtSineSphericalHarmonicCoefficientsOfCentralBody( sinePartials, 0, 1 );

    // std::cout << "cosine partials" << std::endl;
    // std::cout << cosinePartials << std::endl;
    // std::cout << "sine partials" << std::endl;
    // std::cout << sinePartials << std::endl;
    // std::cout << "(2,0) " << basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 ) << std::endl;
    // std::cout << "(2,1) " << basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 ) << std::endl;
    // std::cout << "(2,2) " << basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 ) << std::endl;

    cosinePartials.block(0, 0, 3, 1) /= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
    cosinePartials.block(0, 1, 3, 1) /= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
    cosinePartials.block(0, 2, 3, 1) /= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );

    sinePartials.block(0, 0, 3, 1) /= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
    sinePartials.block(0, 1, 3, 1) /= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );

    // // in case no order 1
    // cosinePartials.block(0, 1, 3, 1) = Eigen::MatrixXd::Zero( 3, 1 );
    // sinePartials.block(0, 0, 3, 1) = Eigen::MatrixXd::Zero( 3, 1 );

    deformationPartial.block(0, 0, 3, 3) += cosinePartials;
    deformationPartial.block(0, 3, 3, 2) += sinePartials;
    
}

}  // namespace acceleration_partials

}  // namespace tudat
