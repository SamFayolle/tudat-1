/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/gravity_deformation_partials/maxwellDeformationPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicCosineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicSineCoefficients.h"

namespace tudat
{

namespace acceleration_partials
{

// Compute spherical Jacobian
Eigen::Matrix3d computeSphericalJacobian(const Eigen::Vector3d& p)
{
    double x = p(0);
    double y = p(1);
    double z = p(2);
    double r = p.norm();
    double s = std::sqrt(x*x + y*y);

    Eigen::Matrix3d jacobian;
    
    // Radial distance
    jacobian(0,0) = x/r; 
    jacobian(0,1) = y/r; 
    jacobian(0,2) = z/r;

    // Colatitude
    if(s < 1e-12) // pole limit case
    { 
        jacobian(1,0) = 0.0; 
        jacobian(1,1) = 0.0; 
        jacobian(1,2) = 0.0;  
    }
    else 
    { 
        jacobian(1,0) = x*z/(r*r*s); 
        jacobian(1,1) = y*z/(r*r*s); 
        jacobian(1,2) = - s/(r*r); 
    }
    
    // Longitude
    if(s < 1e-12) // pole limit case
    { 
        jacobian(2,0) = 0.0; 
        jacobian(2,1) = 0.0; 
        jacobian(2,2) = 0.0;  
    }
    else 
    { 
        jacobian(2,0) = -y/(s*s); 
        jacobian(2,1) = x/(s*s); 
        jacobian(2,2)=0; 
    }
    return jacobian;
}


void computeFullSphericalStatePartials(const Eigen::Vector3d& bodyFixedPosition, 
                                       const Eigen::Vector3d& bodyFixedVelocity,
                                       const Eigen::Vector3d& omega_b,
                                       Eigen::Matrix< double, 6, 6 >& dX_dS )
{
    Eigen::Matrix3d jacobian = computeSphericalJacobian(bodyFixedPosition); 
    
    // Partials spherical position wrt body-fixed position  
    dX_dS.block< 3, 3 >( 0, 0 ) = jacobian; 
    
    // Partials spherical position wrt body-fixed velocity
    dX_dS.block< 3, 3 >( 0, 3 ) = Eigen::Matrix3d::Zero(); 
    
    // Partials spherical velocity wrt body-fixed position
    Eigen::Matrix3d dxdot_dp; 
    
    // Partial distance derivative wrt position 
    double distance = bodyFixedPosition.segment( 0, 3 ).norm( ); 
    double distanceDerivative = ( bodyFixedPosition[ 0 ] * bodyFixedVelocity[ 0 ] + bodyFixedPosition[ 1 ] * bodyFixedVelocity[ 1 ] + bodyFixedPosition[ 2 ] * bodyFixedVelocity[ 2 ] ) 
        / distance; 
    dxdot_dp(0,0) = ( bodyFixedVelocity[0] - ( bodyFixedPosition[0] / distance ) * distanceDerivative ) / distance; 
    dxdot_dp(0,1) = ( bodyFixedVelocity[1] - ( bodyFixedPosition[1] / distance ) * distanceDerivative ) / distance; 
    dxdot_dp(0,2) = ( bodyFixedVelocity[2] - ( bodyFixedPosition[2] / distance ) * distanceDerivative ) / distance; 
    
    // Partial ( -latitude ) derivative wrt position 
    double x = bodyFixedPosition.x(), y = bodyFixedPosition.y(), z = bodyFixedPosition.z(); 
    double vx = bodyFixedVelocity.x(), vy = bodyFixedVelocity.y(), vz = bodyFixedVelocity.z(); 
    double r2 = x*x + y*y + z*z, r = std::sqrt(r2); 
    double s2 = x*x + y*y, s = std::sqrt(s2); 
    double r_dot = (x*vx + y*vy + z*vz)/r; 
    
    // Define numerator and denominator 
    double N = vz * r - z * r_dot; 
    double D = r * s; 
    
    // partials of numerator 
    double dNdx = (vz * x - z * vx)/r + z * x * r_dot / (r*r); 
    double dNdy = (vz * y - z * vy)/r + z * y * r_dot / (r*r); 
    double dNdz = vz * z / r - (r_dot + z * (vz - z * r_dot / r) / r); 
    
    // partials of denominator 
    double dDdx = s * (x/r) + r * (x/s); 
    double dDdy = s * (y/r) + r * (y/s); 
    double dDdz = s * (z/r); 
    dxdot_dp(1,0) = ( N * dDdx - dNdx * D ) / (D*D); 
    dxdot_dp(1,1) = ( N * dDdy - dNdy * D ) / (D*D); 
    dxdot_dp(1,2) = ( N * dDdz - dNdz * D ) / (D*D); 
    
    // Partial longitude derivative wrt position 
    N = x*vy - y*vx; 
    
    Eigen::Vector3d dLambdaDot; 
    if ( s2 >= 1e-12 ) 
    { 
        dxdot_dp(2,0) = ( vy * s2 - 2.0 * x * N ) / ( s2 * s2 ); 
        dxdot_dp(2,1) = ( -vx * s2 - 2.0 * y * N ) / ( s2 * s2 ); 
        dxdot_dp(2,2) = 0.0; 
    } 
    
    dX_dS.block< 3, 3 >( 3, 0 ) = dxdot_dp; 
     
    // partials spherical velocity wrt body-fixed velocity
    dX_dS.block< 3, 3 >( 3, 3 ) = jacobian;

}


//! Contructor.
MaxwellDeformationPartial::MaxwellDeformationPartial(
        const std::string& deformingBody,
        const std::vector< std::string >& perturbingBodies,
        const std::shared_ptr< basic_astrodynamics::MaxwellGravityDeformationModel > deformationModel,
        const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartials ):
    DeformationPartial( deformingBody, perturbingBodies, basic_astrodynamics::maxwell_deformation ),
    deformationModel_( deformationModel ), 
    rotationMatrixPartials_( rotationMatrixPartials ) 
{
    includeOrder1_ = deformationModel->isOrder1Included( );
}

//! Function to create a function returning a partial w.r.t. a double parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > MaxwellDeformationPartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    // Declare return variables, default number of rows = 0 (i.e. no dependency)
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;

    // Check gravitational parameters.
    if( parameter->getParameterName( ).first == estimatable_parameters::gravitational_parameter )
    {
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionPair =
                getGravitationalParameterPartialFunction( parameter->getParameterName( ) );
        partialFunction = partialFunctionPair.first;
        numberOfRows = partialFunctionPair.second;
    }
    // Check rotational properties of the deforming body
    else if( parameter->getParameterName( ).second.first == deformingBody_ )
    {
        // Check if partial is a rotational property of body exerting acceleration.
        if( estimatable_parameters::isParameterRotationMatrixProperty( parameter->getParameterName( ).first ) )
        {
            // Check if required rotation matrix partial exists.
            if( rotationMatrixPartials_.count(
                        std::make_pair( parameter->getParameterName( ).first, parameter->getSecondaryIdentifier( ) ) ) != 0 )
            {
                // Get partial function.
                partialFunction = std::bind( &MaxwellDeformationPartial::wrtRotationModelParameter,
                                             this,
                                             std::placeholders::_1,
                                             parameter->getParameterName( ).first,
                                             parameter->getSecondaryIdentifier( ) );
                numberOfRows = 1;
            }
            else
            {
                std::string errorMessage = "Error, not taking partial of maxwell deformation wrt rotational parameter" +
                        std::to_string( parameter->getParameterName( ).first ) + " of " + parameter->getParameterName( ).second.first;
                throw std::runtime_error( errorMessage );
            }
        }

        // // Check if partial is a tidal property of body exerting acceleration.
        // else if( estimatable_parameters::isParameterTidalProperty( parameter->getParameterName( ).first ) )
        // {
        //     // Check input consistency
        //     std::shared_ptr< estimatable_parameters::TidalLoveNumber< double > > tidalLoveNumber =
        //             std::dynamic_pointer_cast< estimatable_parameters::TidalLoveNumber< double > >( parameter );
        //     if( tidalLoveNumber == nullptr )
        //     {
        //         throw std::runtime_error( "Error when getting tidal Love number vector parameter, object is nullptr" );
        //     }

        //     // Get degree and order(s) of tidal variations
        //     int degree = tidalLoveNumber->getDegree( );
        //     std::vector< int > orders = tidalLoveNumber->getOrders( );
        //     int sumOrders = tidalLoveNumber->getSumOrders( );

        //     std::pair< int, std::pair< int, int > > currentTidalPartialOutput;
        //     for( unsigned int i = 0; i < tidalLoveNumberPartialInterfaces_.size( ); i++ )
        //     {
        //         // Check dependency on current partial object
        //         currentTidalPartialOutput =
        //                 tidalLoveNumberPartialInterfaces_.at( i )->setParameterPartialFunction( parameter, maximumDegree_, maximumOrder_ );

        //         // Check consistency
        //         if( numberOfRows != 0 && currentTidalPartialOutput.first > 0 )
        //         {
        //             throw std::runtime_error( "Error when getting double tidal parameter partial, multiple dependencies found " +
        //                                       std::to_string( numberOfRows ) + ", " + std::to_string( currentTidalPartialOutput.first ) );
        //         }
        //         else
        //         {
        //             // If tidal dependency esists, set partial function
        //             if( currentTidalPartialOutput.first > 0 )
        //             {
        //                 std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > coefficientPartialFunction =
        //                         std::bind( &orbit_determination::TidalLoveNumberPartialInterface::getCurrentDoubleParameterPartial,
        //                                    tidalLoveNumberPartialInterfaces_.at( i ),
        //                                    parameter,
        //                                    currentTidalPartialOutput.second );
        //                 partialFunction = std::bind( &MaxwellDeformationPartial::wrtTidalModelParameter,
        //                                              this,
        //                                              coefficientPartialFunction,
        //                                              degree,
        //                                              orders,
        //                                              sumOrders,
        //                                              parameter->getParameterSize( ),
        //                                              std::placeholders::_1 );
        //                 numberOfRows = currentTidalPartialOutput.first;
        //             }
        //         }
        //     }
        // }
    }

    // Return partial function and partial size.
    return std::make_pair( partialFunction, numberOfRows );
}

//! Function to create a function returning a partial w.r.t. a vector parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > MaxwellDeformationPartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    using namespace tudat::estimatable_parameters;

    // Declare return variables, default number of rows = 0 (i.e. no dependency)
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;

    // Check properties of body exerting acceleration.
    if( parameter->getParameterName( ).second.first == deformingBody_ )
    {
        // Check if partial is a rotational property of body exerting acceleration.
        if( estimatable_parameters::isParameterRotationMatrixProperty( parameter->getParameterName( ).first ) )
        {
            // Check if required rotation matrix partial exists.
            if( rotationMatrixPartials_.count(
                        std::make_pair( parameter->getParameterName( ).first, parameter->getSecondaryIdentifier( ) ) ) != 0 )
            {
                // Get partial function.
                partialFunction = std::bind( &MaxwellDeformationPartial::wrtRotationModelParameter,
                                             this,
                                             std::placeholders::_1,
                                             parameter->getParameterName( ).first,
                                             parameter->getSecondaryIdentifier( ) );
                numberOfRows = parameter->getParameterSize( );
            }
            else
            {
                std::string errorMessage = "Error, not taking partial of sh acceleration wrt rotational parameter" +
                        std::to_string( parameter->getParameterName( ).first ) + " of " + parameter->getParameterName( ).second.first;
                throw std::runtime_error( errorMessage );
            }
        }
        // // Check if partial is a tidal property of body exerting acceleration.
        // else if( estimatable_parameters::isParameterTidalProperty( parameter->getParameterName( ).first ) )
        // {
        //     if( parameter->getParameterName( ).first != mode_coupled_tidal_love_numbers )
        //     {
        //         // Check input consistency
        //         std::shared_ptr< estimatable_parameters::TidalLoveNumber< Eigen::VectorXd > > tidalLoveNumber =
        //                 std::dynamic_pointer_cast< estimatable_parameters::TidalLoveNumber< Eigen::VectorXd > >( parameter );
        //         if( tidalLoveNumber == nullptr )
        //         {
        //             throw std::runtime_error( "Error when getting tidal Love number vector parameter, object is nullptr" );
        //         }

        //         // Get degree and order(s) of tidal variations
        //         int degree = tidalLoveNumber->getDegree( );
        //         std::vector< int > orders = tidalLoveNumber->getOrders( );
        //         int sumOrders = tidalLoveNumber->getSumOrders( );

        //         std::pair< int, std::pair< int, int > > currentTidalPartialOutput;
        //         for( unsigned int i = 0; i < tidalLoveNumberPartialInterfaces_.size( ); i++ )
        //         {
        //             // Check dependency on current partial object
        //             currentTidalPartialOutput = tidalLoveNumberPartialInterfaces_.at( i )->setParameterPartialFunction(
        //                     parameter, maximumDegree_, maximumOrder_ );
        //             if( numberOfRows != 0 && currentTidalPartialOutput.first > 0 )
        //             {
        //                 std::cout << i << std::endl;
        //                 throw std::runtime_error( "Error when getting vector tidal parameter partial B, inconsistent output" +
        //                                           std::to_string( numberOfRows ) + ", " +
        //                                           std::to_string( currentTidalPartialOutput.first ) );
        //             }
        //             else
        //             {
        //                 // If tidal dependency esists, set partial function
        //                 if( currentTidalPartialOutput.first > 0 )
        //                 {
        //                     std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > coefficientPartialFunction =
        //                             std::bind( &orbit_determination::TidalLoveNumberPartialInterface::getCurrentVectorParameterPartial,
        //                                        tidalLoveNumberPartialInterfaces_.at( i ),
        //                                        parameter,
        //                                        currentTidalPartialOutput.second );
        //                     partialFunction = std::bind( &MaxwellDeformationPartial::wrtTidalModelParameter,
        //                                                  this,
        //                                                  coefficientPartialFunction,
        //                                                  degree,
        //                                                  orders,
        //                                                  sumOrders,
        //                                                  parameter->getParameterSize( ),
        //                                                  std::placeholders::_1 );

        //                     numberOfRows = currentTidalPartialOutput.first;
        //                 }
        //             }
        //         }
        //     }
        //     else
        //     {
        //         // Check input consistency
        //         std::shared_ptr< estimatable_parameters::ModeCoupledTidalLoveNumber > tidalLoveNumber =
        //                 std::dynamic_pointer_cast< estimatable_parameters::ModeCoupledTidalLoveNumber >( parameter );
        //         if( tidalLoveNumber == nullptr )
        //         {
        //             throw std::runtime_error( "Error when getting mode coupled tidal Love number vector parameter, object is nullptr" );
        //         }

        //         std::pair< int, std::pair< int, int > > currentTidalPartialOutput;
        //         for( unsigned int i = 0; i < tidalLoveNumberPartialInterfaces_.size( ); i++ )
        //         {
        //             // Check dependency on current partial object
        //             currentTidalPartialOutput = tidalLoveNumberPartialInterfaces_.at( i )->setParameterPartialFunction(
        //                     parameter, maximumDegree_, maximumOrder_ );
        //             if( numberOfRows != 0 && currentTidalPartialOutput.first > 0 )
        //             {
        //                 throw std::runtime_error( "Error when getting vector tidal parameter partial A, inconsistent output" +
        //                                           std::to_string( numberOfRows ) + ", " +
        //                                           std::to_string( currentTidalPartialOutput.first ) );
        //             }
        //             else
        //             {
        //                 // If tidal dependency esists, set partial function
        //                 if( currentTidalPartialOutput.first > 0 )
        //                 {
        //                     std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > coefficientPartialFunction =
        //                             std::bind( &orbit_determination::TidalLoveNumberPartialInterface::getCurrentVectorParameterPartial,
        //                                        tidalLoveNumberPartialInterfaces_.at( i ),
        //                                        parameter,
        //                                        currentTidalPartialOutput.second );
        //                     partialFunction = std::bind( &MaxwellDeformationPartial::wrtModeCoupledLoveNumbers,
        //                                                  this,
        //                                                  coefficientPartialFunction,
        //                                                  tidalLoveNumber->getResponseIndices( ),
        //                                                  tidalLoveNumber->getResponseDegreeOrders( ),
        //                                                  tidalLoveNumber->getParameterSize( ),
        //                                                  std::placeholders::_1 );

        //                     numberOfRows = currentTidalPartialOutput.first;
        //                 }
        //             }
        //         }
        //     }
        // }
        // else if( estimatable_parameters::isParameterNonTidalGravityFieldVariationProperty( parameter->getParameterName( ).first ) )
        // {
        //     switch( parameter->getParameterName( ).first )
        //     {
        //         case polynomial_gravity_field_variation_amplitudes: {
        //             std::shared_ptr< PolynomialGravityFieldVariationsParameters > polynomialVariationParameter =
        //                     std::dynamic_pointer_cast< PolynomialGravityFieldVariationsParameters >( parameter );
        //             std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > indexAndPowerPerCosineBlockIndex =
        //                     polynomialVariationParameter->getIndexAndPowerPerCosineBlockIndex( );
        //             std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > indexAndPowerPerSineBlockIndex =
        //                     polynomialVariationParameter->getIndexAndPowerPerSineBlockIndex( );

        //             partialFunction = std::bind( &MaxwellDeformationPartial::wrtPolynomialGravityFieldVariations,
        //                                          this,
        //                                          utilities::createVectorFromMapKeys( indexAndPowerPerCosineBlockIndex ),
        //                                          utilities::createVectorFromMapKeys( indexAndPowerPerSineBlockIndex ),
        //                                          utilities::createVectorFromMapValues( indexAndPowerPerCosineBlockIndex ),
        //                                          utilities::createVectorFromMapValues( indexAndPowerPerSineBlockIndex ),
        //                                          polynomialVariationParameter->getPolynomialVariationModel( )->getReferenceEpoch( ),
        //                                          std::placeholders::_1 );

        //             numberOfRows = parameter->getParameterSize( );
        //             break;
        //         }
        //         case periodic_gravity_field_variation_amplitudes: {
        //             std::shared_ptr< PeriodicGravityFieldVariationsParameters > periodicVariationParameter =
        //                     std::dynamic_pointer_cast< PeriodicGravityFieldVariationsParameters >( parameter );
        //             std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > indexAndPowerPerCosineBlockIndex =
        //                     periodicVariationParameter->getIndexAndPowerPerCosineBlockIndex( );
        //             std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > indexAndPowerPerSineBlockIndex =
        //                     periodicVariationParameter->getIndexAndPowerPerSineBlockIndex( );

        //             partialFunction = std::bind( &MaxwellDeformationPartial::wrtPeriodicGravityFieldVariations,
        //                                          this,
        //                                          utilities::createVectorFromMapKeys( indexAndPowerPerCosineBlockIndex ),
        //                                          utilities::createVectorFromMapKeys( indexAndPowerPerSineBlockIndex ),
        //                                          utilities::createVectorFromMapValues( indexAndPowerPerCosineBlockIndex ),
        //                                          utilities::createVectorFromMapValues( indexAndPowerPerSineBlockIndex ),
        //                                          periodicVariationParameter->getPeriodicVariationModel( )->getFrequencies( ),
        //                                          periodicVariationParameter->getPeriodicVariationModel( )->getReferenceEpoch( ),
        //                                          std::placeholders::_1 );

        //             numberOfRows = parameter->getParameterSize( );
        //             break;
        //         }
        //         default:
        //             break;
        //     }
        // }
        // Check non-rotational parameters.
        else
        {
            switch( parameter->getParameterName( ).first )
            {
                case spherical_harmonics_cosine_coefficient_block: 
                {
                    std::shared_ptr< SphericalHarmonicsCosineCoefficients > coefficientsParameter =
                            std::dynamic_pointer_cast< SphericalHarmonicsCosineCoefficients >( parameter );

                    partialFunction = std::bind( &MaxwellDeformationPartial::wrtCosineCoefficientBlock,
                                                 this,
                                                 coefficientsParameter->getBlockIndices( ),
                                                 std::placeholders::_1 );
                    numberOfRows = coefficientsParameter->getParameterSize( );

                    break;
                }
                case spherical_harmonics_sine_coefficient_block: 
                {
                    std::shared_ptr< SphericalHarmonicsSineCoefficients > coefficientsParameter =
                            std::dynamic_pointer_cast< SphericalHarmonicsSineCoefficients >( parameter );

                    partialFunction = std::bind( &MaxwellDeformationPartial::wrtSineCoefficientBlock,
                                                 this,
                                                 coefficientsParameter->getBlockIndices( ),
                                                 std::placeholders::_1 );
                    numberOfRows = coefficientsParameter->getParameterSize( );
                    break;
                }
                default:
                    break;
            }
        }
    }

    // Return partial function and partial size.
    return std::make_pair( partialFunction, numberOfRows );
}

//! Function to create a function returning the current partial w.r.t. a gravitational parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > MaxwellDeformationPartial::getGravitationalParameterPartialFunction(
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterId )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfColumns = 0;

    if( parameterId.first == estimatable_parameters::gravitational_parameter )
    {
        // Check for dependency
        if( parameterId.second.first == deformingBody_ )
        {
            partialFunction = std::bind( 
                &MaxwellDeformationPartial::wrtGravitationalParameterOfDeformingBody, this, std::placeholders::_1, 0 );
            numberOfColumns = 1;
        }

        if( std::find( perturbingBodies_.begin( ), perturbingBodies_.end( ), parameterId.second.first ) != perturbingBodies_.end( ) )
        {
            partialFunction = std::bind(
                &MaxwellDeformationPartial::wrtGravitationalParameterOfPerturbingBody, this, std::placeholders::_1, 0 );
            numberOfColumns = 1;
        }
    }
    return std::make_pair( partialFunction, numberOfColumns );
}

Eigen::MatrixXd MaxwellDeformationPartial::equilibriumCoefficientsWrtSphericalBodyFixedState( )
{
    double kf = deformationModel_->getLoveNumber( );
    double muRatio = deformationModel_->getGravitationalParametersRatio();
    double r = deformationModel_->getCurrentRelativePosition( ).at( 0 ).segment( 0, 3 ).norm( ); // TO BE MODIFIED -> SET TO INDEX 0 FOR NOW
    double currentLatitude = deformationModel_->getCurrentLatitude( ).at( 0 ); // TO BE MODIFIED -> SET TO INDEX 0 FOR NOW
    double currentLongitude = deformationModel_->getCurrentLongitude( ).at( 0 ); // TO BE MODIFIED -> SET TO INDEX 0 FOR NOW
    double radiusOverDistance = deformationModel_->getReferenceRadius( ) / r;
    double radiusOverDistancePowerThree = radiusOverDistance * radiusOverDistance * radiusOverDistance;
    
    Eigen::MatrixXd partial = Eigen::MatrixXd::Zero( 5, 6 );

    double c20eq_wrt_r = - 3.0 * kf / 2.0 * muRatio * radiusOverDistancePowerThree / r 
        * ( 3.0 * std::sin( currentLatitude ) * std::sin( currentLatitude ) - 1.0 );
    double c21eq_wrt_r = - 3.0 * kf * muRatio * radiusOverDistancePowerThree / r * std::cos( currentLongitude ) 
        * std::cos( currentLatitude ) * std::sin( currentLatitude );
    double c22eq_wrt_r = - 3.0 * kf / 4.0 * muRatio * radiusOverDistancePowerThree / r * std::cos( 2.0 * currentLongitude )
        * ( 1.0 - std::sin( currentLatitude ) * std::sin( currentLatitude ) );
    double s21eq_wrt_r = - 3.0 * kf * muRatio * radiusOverDistancePowerThree / r * std::sin( currentLongitude ) 
        * std::cos( currentLatitude ) * std::sin( currentLatitude );
    double s22eq_wrt_r = - 3.0 * kf / 4.0 * muRatio * radiusOverDistancePowerThree / r * std::sin( 2.0 * currentLongitude )
        * ( 1.0 - std::sin( currentLatitude ) * std::sin( currentLatitude ) );

    double c20eq_wrt_longitude = 0.0;
    double c21eq_wrt_longitude = kf * muRatio * radiusOverDistancePowerThree * ( - std::sin( currentLongitude ) )
        * std::cos( currentLatitude ) * std::sin( currentLatitude );
    double c22eq_wrt_longitude = kf / 4.0 * muRatio * radiusOverDistancePowerThree * ( - 2.0 * std::sin( 2.0 * currentLongitude ) )
        * ( 1.0 - std::sin( currentLatitude ) * std::sin( currentLatitude ) );
    double s21eq_wrt_longitude = kf * muRatio * radiusOverDistancePowerThree * ( std::cos( currentLongitude ) ) 
        * std::cos( currentLatitude ) * std::sin( currentLatitude );
    double s22eq_wrt_longitude = kf / 4.0 * muRatio * radiusOverDistancePowerThree * ( 2.0 * std::cos( 2.0 * currentLongitude ) )
        * ( 1.0 - std::sin( currentLatitude ) * std::sin( currentLatitude ) );

    double c20eq_wrt_latitude = kf / 2.0 * muRatio * radiusOverDistancePowerThree 
        * ( 3.0 * 2.0 * std::sin( currentLatitude ) * std::cos( currentLatitude ) );
    double c21eq_wrt_latitude = kf * muRatio * radiusOverDistancePowerThree 
        * ( std::cos( currentLatitude ) * std::cos( currentLatitude ) - std::sin( currentLatitude ) * std::sin( currentLatitude ) ) * std::cos( currentLongitude );
    double c22eq_wrt_latitude = kf / 4.0 * muRatio * radiusOverDistancePowerThree 
        * ( - 2.0 * std::sin( currentLatitude ) * std::cos( currentLatitude ) ) * std::cos( 2.0 * currentLongitude );
    double s21eq_wrt_latitude = kf * muRatio * radiusOverDistancePowerThree 
        * ( std::cos( currentLatitude ) * std::cos( currentLatitude ) - std::sin( currentLatitude ) * std::sin( currentLatitude ) ) * std::sin( currentLongitude );
    double s22eq_wrt_latitude = kf / 4.0 * muRatio * radiusOverDistancePowerThree 
        * ( - 2.0 * std::sin( currentLatitude ) * std::cos( currentLatitude ) ) * std::sin( 2.0 * currentLongitude );

    // Set partials wrt radial distance
    partial( 0, 0 ) = c20eq_wrt_r;
    partial( 2, 0 ) = c22eq_wrt_r;
    partial( 4, 0 ) = s22eq_wrt_r;
    if ( includeOrder1_ )
    {
        partial( 1, 0 ) = c21eq_wrt_r;
        partial( 3, 0 ) = s21eq_wrt_r;
    }

   // Set partials wrt longitude
   partial( 0, 2 ) = c20eq_wrt_longitude;
   partial( 2, 2 ) = c22eq_wrt_longitude;
   partial( 4, 2 ) = s22eq_wrt_longitude;
   if ( includeOrder1_ )
   {
        partial( 1, 2 ) = c21eq_wrt_longitude;
        partial( 3, 2 ) = s21eq_wrt_longitude;
   }

   // Set partials wrt - latitude
   partial( 0, 1 ) = - c20eq_wrt_latitude;
   partial( 2, 1 ) = - c22eq_wrt_latitude;
   partial( 4, 1 ) = - s22eq_wrt_latitude;
   if ( includeOrder1_ )
    {
        partial( 1, 1 ) = - c21eq_wrt_latitude;
        partial( 3, 1 ) = - s21eq_wrt_latitude;
    }

//    partial = Eigen::MatrixXd::Zero( 5, 6 );

   return partial;

}

Eigen::MatrixXd MaxwellDeformationPartial::equilibriumCoefficientsDerivativeWrtSphericalBodyFixedState( )
{
    double kf = deformationModel_->getLoveNumber( );
    double muRatio = deformationModel_->getGravitationalParametersRatio();
    double r = deformationModel_->getCurrentRelativePosition( ).at( 0 ).segment( 0, 3 ).norm( ); // TO BE MODIFIED -> SET TO INDEX 0 FOR NOW
    double currentLatitude = deformationModel_->getCurrentLatitude( ).at( 0 ); // TO BE MODIFIED -> SET TO INDEX 0 FOR NOW
    double currentLongitude = deformationModel_->getCurrentLongitude( ).at( 0 ); // TO BE MODIFIED -> SET TO INDEX 0 FOR NOW
    double radiusOverDistance = deformationModel_->getReferenceRadius( ) / r;
    double radiusOverDistancePowerThree = radiusOverDistance * radiusOverDistance * radiusOverDistance;

    double distanceDerivative = deformationModel_->getCurrentDistanceDerivative( ).at( 0 ); // TO BE MODIFIED -> SET TO INDEX 0 FOR NOW
    double latitudeDerivative = deformationModel_->getCurrentLatitudeDerivative( ).at( 0 ); // TO BE MODIFIED -> SET TO INDEX 0 FOR NOW
    double longitudeDerivative = deformationModel_->getCurrentLongitudeDerivative( ).at( 0 ); // TO BE MODIFIED -> SET TO INDEX 0 FOR NOW

    double rSquared = r * r;
    double sinPhiCosPhi = std::sin( currentLatitude ) * std::cos( currentLatitude );
    double sinPhiSquared = std::sin( currentLatitude ) * std::sin( currentLatitude );
    double cosPhi2MinusSinPhi2 = 
        std::cos( currentLatitude ) * std::cos( currentLatitude ) - std::sin( currentLatitude ) * std::sin( currentLatitude );
    
    Eigen::MatrixXd partial = Eigen::MatrixXd::Zero( 5, 6 );

    double factor = - kf * muRatio * radiusOverDistancePowerThree;

    // partials equilibrium coefficients derivatives wrt radial distance
    double dc20eq_r = factor / 2.0 * ( - 12.0 * distanceDerivative / rSquared * ( 3.0 * sinPhiSquared - 1.0 )
                                       - 3.0 / r * ( - 6.0 * latitudeDerivative * sinPhiCosPhi ) );
    double dc21eq_r = factor * ( 
        - 12.0 * distanceDerivative / rSquared * sinPhiCosPhi * std::cos( currentLongitude )
        - 3.0 / r * ( 
            longitudeDerivative * sinPhiCosPhi * std::sin( currentLongitude )
            - latitudeDerivative * cosPhi2MinusSinPhi2 * std::cos( currentLongitude ) ) );
    double ds21eq_r = factor * ( 
        - 12.0 * distanceDerivative / rSquared * sinPhiCosPhi * std::sin( currentLongitude )
        - 3.0 / r * ( 
            - longitudeDerivative * sinPhiCosPhi * std::cos( currentLongitude )
            - latitudeDerivative * cosPhi2MinusSinPhi2 * std::sin( currentLongitude ) ) );
    double dc22eq_r = factor / 4.0 * (
        - 12.0 * distanceDerivative / rSquared * ( 1.0 - sinPhiSquared ) * std::cos( 2.0 * currentLongitude )
        - 3.0 / r * (
            2.0 * longitudeDerivative * ( 1.0 - sinPhiSquared ) * std::sin( 2.0 * currentLongitude )
            + 2.0 * latitudeDerivative * sinPhiCosPhi * std::cos( 2.0 * currentLongitude ) ) );
    double ds22eq_r = factor / 4.0 * (
        - 12.0 * distanceDerivative / rSquared * ( 1.0 - sinPhiSquared ) * std::sin( 2.0 * currentLongitude )
        - 3.0 / r * (
            - 2.0 * longitudeDerivative * ( 1.0 - sinPhiSquared ) * std::cos( 2.0 * currentLongitude )
            + 2.0 * latitudeDerivative * sinPhiCosPhi * std::sin( 2.0 * currentLongitude ) ) );

    // partials equilibrium coefficients derivatives wrt longitude
    double dc20eq_longitude = 0.0;
    double dc21eq_longitude = factor * ( 
        - 3.0 * distanceDerivative / r * sinPhiCosPhi * std::sin( currentLongitude )
        + longitudeDerivative * sinPhiCosPhi * std::cos( currentLongitude )
        + latitudeDerivative * cosPhi2MinusSinPhi2 * std::sin( currentLongitude ) );
    double ds21eq_longitude = factor * ( 
        3.0 * distanceDerivative / r * sinPhiCosPhi * std::cos( currentLongitude )
        + longitudeDerivative * sinPhiCosPhi * std::sin( currentLongitude )
        - latitudeDerivative * cosPhi2MinusSinPhi2 * std::cos( currentLongitude ) );
    double dc22eq_longitude = factor / 4.0 * ( 
        - 6.0 * distanceDerivative / r * ( 1.0 - sinPhiSquared ) * std::sin( 2.0 * currentLongitude )
        + 4.0 * longitudeDerivative * ( 1.0 - sinPhiSquared ) * std::cos( 2.0 * currentLongitude )
        - 4.0 * latitudeDerivative * sinPhiCosPhi * std::sin( 2.0 * currentLongitude ) );
    double ds22eq_longitude = factor / 4.0 * ( 
        6.0 * distanceDerivative / r * ( 1.0 - sinPhiSquared ) * std::cos( 2.0 * currentLongitude )
        + 4.0 * longitudeDerivative * ( 1.0 - sinPhiSquared ) * std::sin( 2.0 * currentLongitude )
        + 4.0 * latitudeDerivative * sinPhiCosPhi * std::cos( 2.0 * currentLongitude ) );

    // partials equilibrium coefficients derivatives wrt latitude    
    double dc20eq_latitude = factor / 2.0 * ( 
        + 18.0 * distanceDerivative / r * sinPhiCosPhi 
        - 6.0 * latitudeDerivative * cosPhi2MinusSinPhi2 );
    double dc21eq_latitude = factor * ( 
        3.0 * distanceDerivative / r * cosPhi2MinusSinPhi2 * std::cos( currentLongitude )
        + longitudeDerivative * cosPhi2MinusSinPhi2 * std::sin( currentLongitude )
        + 4.0 * latitudeDerivative * sinPhiCosPhi * std::cos( currentLongitude ) );
    double ds21eq_latitude = factor * ( 
        3.0 * distanceDerivative / r * cosPhi2MinusSinPhi2 * std::sin( currentLongitude )
        - longitudeDerivative * cosPhi2MinusSinPhi2 * std::cos( currentLongitude )
        + 4.0 * latitudeDerivative * sinPhiCosPhi * std::sin( currentLongitude ) );
    double dc22eq_latitude = factor / 4.0 * ( 
        - 6.0 * distanceDerivative / r * sinPhiCosPhi * std::cos( 2.0 * currentLongitude )
        - 4.0 * longitudeDerivative * sinPhiCosPhi * std::sin( 2.0 * currentLongitude )
        + 2.0 * latitudeDerivative * cosPhi2MinusSinPhi2 * std::cos( 2.0 * currentLongitude ) );
    double ds22eq_latitude = factor / 4.0 * ( 
        - 6.0 * distanceDerivative / r * sinPhiCosPhi * std::sin( 2.0 * currentLongitude )
        + 4.0 * longitudeDerivative * sinPhiCosPhi * std::cos( 2.0 * currentLongitude )
        + 2.0 * latitudeDerivative * cosPhi2MinusSinPhi2 * std::sin( 2.0 * currentLongitude ) );

    // partials equilibrium coefficients derivatives wrt radial distance derivative    
    double dc20eq_dr = - kf / 2.0 * muRatio * radiusOverDistancePowerThree * 3.0 / r * ( 3.0 * sinPhiSquared - 1.0 );
    double dc21eq_dr = - kf * muRatio * radiusOverDistancePowerThree * 3.0 / r * sinPhiCosPhi * std::cos( currentLongitude );
    double ds21eq_dr = - kf * muRatio * radiusOverDistancePowerThree * 3.0 / r * sinPhiCosPhi * std::sin( currentLongitude );
    double dc22eq_dr = - kf / 4.0 * muRatio * radiusOverDistancePowerThree * 3.0 / r * ( 1.0 - sinPhiSquared ) * std::cos( 2.0 * currentLongitude );
    double ds22eq_dr = - kf / 4.0 * muRatio * radiusOverDistancePowerThree * 3.0 / r * ( 1.0 - sinPhiSquared ) * std::sin( 2.0 * currentLongitude );

    // partials equilibrium coefficients derivatives wrt longitude derivative
    double dc20eq_dlongitude = 0.0;
    double dc21eq_dlongitude = - kf * muRatio * radiusOverDistancePowerThree * sinPhiCosPhi * ( std::sin( currentLongitude ) );
    double ds21eq_dlongitude = - kf * muRatio * radiusOverDistancePowerThree * sinPhiCosPhi * ( - std::cos( currentLongitude ) );
    double dc22eq_dlongitude = - kf / 4.0 * muRatio * radiusOverDistancePowerThree * ( 1.0 - sinPhiSquared ) * ( 2.0 * std::sin( 2.0 * currentLongitude ) );
    double ds22eq_dlongitude = - kf / 4.0 * muRatio * radiusOverDistancePowerThree * ( 1.0 - sinPhiSquared ) * ( - 2.0 * std::cos( 2.0 * currentLongitude ) );

    // partials equilibrium coefficients derivatives wrt latitude derivative
    double dc20eq_dlatitude = - kf / 2.0 * muRatio * radiusOverDistancePowerThree * ( - 6.0 * sinPhiCosPhi );
    double dc21eq_dlatitude = - kf * muRatio * radiusOverDistancePowerThree * ( - cosPhi2MinusSinPhi2 ) * std::cos( currentLongitude );
    double ds21eq_dlatitude = - kf * muRatio * radiusOverDistancePowerThree * ( - cosPhi2MinusSinPhi2 ) * std::sin( currentLongitude );
    double dc22eq_dlatitude = - kf / 4.0 * muRatio * radiusOverDistancePowerThree * 2.0 * sinPhiCosPhi * std::cos( 2.0 * currentLongitude );
    double ds22eq_dlatitude = - kf / 4.0 * muRatio * radiusOverDistancePowerThree * 2.0 * sinPhiCosPhi * std::sin( 2.0 * currentLongitude );


    // Set partials wrt radial distance
    partial( 0, 0 ) = dc20eq_r;
    partial( 2, 0 ) = dc22eq_r;
    partial( 4, 0 ) = ds22eq_r;
    if ( includeOrder1_ )
    {
        partial( 1, 0 ) = dc21eq_r;
        partial( 3, 0 ) = ds21eq_r;
    }

   // Set partials wrt longitude
   partial( 0, 2 ) = dc20eq_longitude;
   partial( 2, 2 ) = dc22eq_longitude;
   partial( 4, 2 ) = ds22eq_longitude;
   if ( includeOrder1_ )
    {
        partial( 1, 2 ) = dc21eq_longitude;
        partial( 3, 2 ) = ds21eq_longitude;
    }

   // Set partials wrt - latitude 
   partial( 0, 1 ) = - dc20eq_latitude;
   partial( 2, 1 ) = - dc22eq_latitude;
   partial( 4, 1 ) = - ds22eq_latitude;
   if ( includeOrder1_ )
    {
        partial( 1, 1 ) = - dc21eq_latitude;
        partial( 3, 1 ) = - ds21eq_latitude;
    }

   // Set partials wrt time derivative radial distance
    partial( 0, 3 ) = dc20eq_dr;
    partial( 2, 3 ) = dc22eq_dr;
    partial( 4, 3 ) = ds22eq_dr;
    if ( includeOrder1_ )
    {
        partial( 1, 3 ) = dc21eq_dr;
        partial( 3, 3 ) = ds21eq_dr;
    }

   // Set partials wrt time derivative longitude
   partial( 0, 5 ) = dc20eq_dlongitude;
   partial( 2, 5 ) = dc22eq_dlongitude;
   partial( 4, 5 ) = ds22eq_dlongitude;
   if ( includeOrder1_ )
    {
        partial( 1, 5 ) = dc21eq_dlongitude;
        partial( 3, 5 ) = ds21eq_dlongitude;
    }

   // Set partials wrt time derivative longitude
   partial( 0, 4 ) = - dc20eq_dlatitude;
   partial( 2, 4 ) = - dc22eq_dlatitude;
   partial( 4, 4 ) = - ds22eq_dlatitude;
   if ( includeOrder1_ )
    {
        partial( 1, 4 ) = - dc21eq_dlatitude;
        partial( 3, 4 ) = - ds21eq_dlatitude;
    }

//    partial = Eigen::MatrixXd::Zero( 5, 6 );

   return partial;

}

Eigen::MatrixXd MaxwellDeformationPartial::equilibriumCoefficientsWrtAngularVelocityVector( )
{
    double kf = deformationModel_->getLoveNumber( );
    double muPerturbingBody = deformationModel_->getGravitationalParameterDeformingBody( );
    double radius = deformationModel_->getReferenceRadius( );
    Eigen::Vector3d angularVelocityVector = deformationModel_->getCurrentAngularVelocityDeformingBody( ); 
    
    Eigen::MatrixXd partial = Eigen::MatrixXd::Zero( 5, 3 );

    // only C20eq depends on the angular velocity vector (and only if the centrifugal perturbing potential is included)
    if ( deformationModel_->isCentrifugalPotentialIncluded( ) )
    {
        partial.block( 0, 0, 1, 3 ) = - 2.0 * kf / ( 3.0 * muPerturbingBody ) * radius * radius * radius * angularVelocityVector.transpose();
    }

   return partial;
}

Eigen::MatrixXd MaxwellDeformationPartial::equilibriumCoefficientsWrtRotationRate( )
{
    double kf = deformationModel_->getLoveNumber( );
    double muPerturbingBody = deformationModel_->getGravitationalParameterDeformingBody( );
    double radius = deformationModel_->getReferenceRadius( );
    double rotationRate = deformationModel_->getCurrentAngularVelocityDeformingBody( ).norm( ); 
    
    Eigen::MatrixXd partial = Eigen::MatrixXd::Zero( 5, 1 );

    // only C20eq depends on the rotation rate (and only if the centrifugal perturbing potential is included)
    if ( deformationModel_->isCentrifugalPotentialIncluded( ) )
    {
        partial( 0, 0 ) = - 2.0 * kf / ( 3.0 * muPerturbingBody ) * radius * radius * radius * rotationRate;
    }

   return partial;
}

Eigen::MatrixXd MaxwellDeformationPartial::equilibriumCoefficientsDerivativeWrtAngularVelocityVector( )
{
    double kf = deformationModel_->getLoveNumber( );
    double muPerturbingBody = deformationModel_->getGravitationalParameterDeformingBody( );
    double radius = deformationModel_->getReferenceRadius( );
    Eigen::Vector3d angularVelocityVectorDerivative = deformationModel_->getCurrentAngularVelocityDerivativeDeformingBody( ); 
    
    Eigen::MatrixXd partial = Eigen::MatrixXd::Zero( 5, 3 );

    // only dC20eq depends on the angular velocity vector (and only if the centrifugal perturbing potential is included)
    if ( deformationModel_->isCentrifugalPotentialIncluded( ) )
    {
        partial.block( 0, 0, 1, 3 ) = - 2.0 * kf / ( 3.0 * muPerturbingBody ) * radius * radius * radius * angularVelocityVectorDerivative.transpose();
    }

   return partial;
}

Eigen::MatrixXd MaxwellDeformationPartial::equilibriumCoefficientsDerivativeWrtAngularVelocityVectorDerivative( )
{
    double kf = deformationModel_->getLoveNumber( );
    double muPerturbingBody = deformationModel_->getGravitationalParameterDeformingBody( );
    double radius = deformationModel_->getReferenceRadius( );
    Eigen::Vector3d angularVelocityVector = deformationModel_->getCurrentAngularVelocityDeformingBody( ); 
    
    Eigen::MatrixXd partial = Eigen::MatrixXd::Zero( 5, 3 );

    // only dC20eq depends on the angular velocity vector derivative (and only if the centrifugal perturbing potential is included)
    if ( deformationModel_->isCentrifugalPotentialIncluded( ) )
    {
        partial.block( 0, 0, 1, 3 ) = - 2.0 * kf / ( 3.0 * muPerturbingBody ) * radius * radius * radius * angularVelocityVector.transpose();
    }

   return partial;
}

Eigen::MatrixXd MaxwellDeformationPartial::deformationWrtSphericalBodyFixedState( )
{
    double invGlobalRelaxationTime = 1.0 / deformationModel_->getGlobalRelaxationTime( );
    double maxwellRelaxationTime = deformationModel_->getMaxwellRelaxationTime( );

    // Retrieve partials equilibrium coefficients and derivatives wrt spherical body-fixed state
    Eigen::MatrixXd equilibriumCoefficientsWrtSpherical = equilibriumCoefficientsWrtSphericalBodyFixedState();
    Eigen::MatrixXd equilibriumCoefficientsDerivativeWrtSpherical = equilibriumCoefficientsDerivativeWrtSphericalBodyFixedState();

    // Compute deformation partials
    Eigen::MatrixXd partials = Eigen::MatrixXd::Zero( 5, 6 );
    partials = invGlobalRelaxationTime * ( equilibriumCoefficientsWrtSpherical + maxwellRelaxationTime * equilibriumCoefficientsDerivativeWrtSpherical  );

    return partials;
}

Eigen::Matrix6d MaxwellDeformationPartial::sphericalWrtCartesianBodyFixedState( )
{
    Eigen::Matrix6d partial = Eigen::Matrix6d::Zero( );

    Eigen::Vector3d cartesianBodyFixedPosition = deformationModel_->getCurrentRelativePosition( )[0]; // TO BE MODIFIED -> INDEX 0
    Eigen::Vector3d cartesianBodyFixedVelocity = deformationModel_->getCurrentRelativeVelocity( )[0]; // TO BE MODIFIED -> INDEX 0
    Eigen::Vector3d angularVelocityBodyFixedFrame = deformationModel_->getCurrentAngularVelocityDeformingBody(); 

    // Eigen::Matrix< double, 6, 6 > partial;
    computeFullSphericalStatePartials( cartesianBodyFixedPosition, cartesianBodyFixedVelocity, angularVelocityBodyFixedFrame, partial );
    
    // std::cout << "spherical wrt cartesian =\n" << partial << "\n";

    return partial;
}

Eigen::Matrix6d MaxwellDeformationPartial::bodyFixedWrtGlobalState( 
    const double currentTime, const bool addIndirectRotationContribution )
{
    Eigen::Matrix3d rotationMatrixInertialToBodyFixed = deformationModel_->getCurrentRotationToIntegrationFrameMatrix( ).transpose( );

    Eigen::Matrix6d partial = Eigen::Matrix6d::Zero( );

    partial.block( 0, 0, 3, 3 ) = rotationMatrixInertialToBodyFixed;
    partial.block( 0, 3, 3, 3 ) = Eigen::Matrix3d::Zero();

    partial.block( 3, 0, 3, 3 ) = deformationModel_->getCurrentRotationToLocalFrameDerivative( );
    partial.block( 3, 3, 3, 3 ) = rotationMatrixInertialToBodyFixed;

    // If rotation matrix depends on translational state, add correction partials
    if( addIndirectRotationContribution && rotationMatrixPartials_.count( std::make_pair( estimatable_parameters::initial_body_state, "" ) ) > 0 )
    {
        // Compute rotation matrix partials
        std::vector< Eigen::Matrix3d > rotationMatrixPartials =
                rotationMatrixPartials_.at( std::make_pair( estimatable_parameters::initial_body_state, "" ) )
                        ->calculatePartialOfRotationMatrixToBaseFrameWrParameter( currentTime );

        Eigen::Vector6d currentInertialState = deformationModel_->getStateOfDeformingBody( ); 
        Eigen::Vector3d currentInertialPosition = currentInertialState.segment( 0, 3 ); //Eigen::Vector3d::Zero();;
        Eigen::Vector3d currentInertialVelocity = currentInertialState.segment( 3, 3 ); //Eigen::Vector3d::Zero();
        // std::cout << "currentInertialState " << currentInertialState.transpose( ) << std::endl;

        std::vector< Eigen::Matrix3d > rotationMatrixDerivativePartials =
                rotationMatrixPartials_.at( std::make_pair( estimatable_parameters::initial_body_state, "" ) )
                        ->calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter( currentTime );

        // Add correction terms to position and velocity partials
        for( unsigned int i = 0; i < 6; i++ )
        {
            // d R / d x,y,z,vx,vy,vz
            Eigen::Matrix3d currentRotationMatrixPartial = rotationMatrixPartials[ i ];

            Eigen::Vector3d currentPositionPartial = rotationMatrixPartials[ i ].transpose() * currentInertialPosition;
            partial.block( 0, i, 3, 1 ) += currentPositionPartial;

            Eigen::Vector3d currentVelocityPartial = rotationMatrixPartials[ i ].transpose() * currentInertialVelocity;
            partial.block( 3, i, 3, 1 ) += currentVelocityPartial;

            currentVelocityPartial = rotationMatrixDerivativePartials[ i ].transpose() * currentInertialPosition;
            partial.block( 3, i, 3, 1 ) += currentVelocityPartial;

            // std::cout << i << " : " << std::endl;
            // // std::cout << rotationMatrixPartials[ i ] << std::endl;
            // std::cout << "pos " << currentPositionPartial.transpose() << std::endl; 
            // std::cout << "vel " << currentVelocityPartial.transpose() << std::endl; 
        }
        // std::cout << "partialBodyFixedWrtGlobalState - after" << std::endl;
        // std::cout << partialBodyFixedWrtGlobalState << std::endl;
    }

    return partial;
}


Eigen::MatrixXd MaxwellDeformationPartial::deformationWrtRotationalState( )
{
    double invGlobalRelaxationTime = 1.0 / deformationModel_->getGlobalRelaxationTime( );
    double maxwellRelaxationTime = deformationModel_->getMaxwellRelaxationTime( );

    // Retrieve partials of equilibrium coefficients and derivatives wrt angular velocity vector
    Eigen::MatrixXd equilibriumCoefficientsPartials = equilibriumCoefficientsWrtAngularVelocityVector();
    Eigen::MatrixXd equilibriumCoefficientsDerivativePartials = equilibriumCoefficientsDerivativeWrtAngularVelocityVector();

    // Compute deformation partials
    Eigen::MatrixXd partials = Eigen::MatrixXd::Zero( 5, 7 );
    partials.block( 0, 4, 5, 3 ) = invGlobalRelaxationTime * ( equilibriumCoefficientsPartials + maxwellRelaxationTime * equilibriumCoefficientsDerivativePartials );

    return partials;
}

Eigen::MatrixXd MaxwellDeformationPartial::deformationWrtRotationRate( )
{
    double invGlobalRelaxationTime = 1.0 / deformationModel_->getGlobalRelaxationTime( );

    // Retrieve partials of equilibrium coefficients wrt (constant) rotation rate
    Eigen::MatrixXd equilibriumCoefficientsPartials = equilibriumCoefficientsWrtRotationRate( );

    // Compute deformation partials
    Eigen::MatrixXd partials = invGlobalRelaxationTime * equilibriumCoefficientsPartials;

    return partials;
}

Eigen::MatrixXd MaxwellDeformationPartial::deformationWrtAngularVelocityVectorDerivative( )
{
    double invGlobalRelaxationTime = 1.0 / deformationModel_->getGlobalRelaxationTime( );
    double maxwellRelaxationTime = deformationModel_->getMaxwellRelaxationTime( );

    // Retrieve partials of equilibrium coefficients' derivatives wrt angular velocity vector derivative
    Eigen::MatrixXd equilibriumCoefficientsDerivativePartials = equilibriumCoefficientsDerivativeWrtAngularVelocityVectorDerivative();

    // Compute deformation partials
    Eigen::MatrixXd partials = invGlobalRelaxationTime * maxwellRelaxationTime * equilibriumCoefficientsDerivativePartials;

    return partials;
}

Eigen::MatrixXd MaxwellDeformationPartial::wrtOtherStateDerivative()
{
    Eigen::MatrixXd partialWrtAngularVelocityDerivative = deformationWrtAngularVelocityVectorDerivative( );

    // std::cout << "partialWrtAngularVelocityDerivative" << std::endl;
    // std::cout << partialWrtAngularVelocityDerivative << std::endl;

    return partialWrtAngularVelocityDerivative;
}

//! Function for updating the partial object to current state and time.
void MaxwellDeformationPartial::update( const double currentTime )
{
    using namespace tudat::coordinate_conversions;

    if( !( currentTime_ == currentTime ) )
    {
        // Update deformation model
        deformationModel_->updateMembers( currentTime );

        equilibriumCoefficients_ = deformationModel_->getEquilibriumCoefficients( );
        derivativeEquilibriumCoefficients_ = deformationModel_->getDerivativeEquilibriumCoefficients( );

        // Compute partial w.r.t. gravity field // TO BE MODIFIED -> DOES NOT NEED TO BE UPDATED
        currentPartialWrtGravity_ = Eigen::MatrixXd::Zero( 5, 5 );
        for ( unsigned int i = 0 ; i < 5 ; i++ )
        {
            currentPartialWrtGravity_( i, i ) = - 1.0 / deformationModel_->getGlobalRelaxationTime( );
        }
        // std::cout << "currentPartialWrtGravity_" << std::endl;
        // std::cout << currentPartialWrtGravity_ << std::endl;
        

        // Compute partial w.r.t. translational state 
        currentPartialWrtPosition_ = Eigen::MatrixXd::Zero( 5, 3 );
        currentPartialWrtVelocity_ = Eigen::MatrixXd::Zero( 5, 3 );

        
        // Eigen::Matrix6d partialBodyFixedWrtGlobalState = bodyFixedWrtGlobalState( currentTime );
        // std::cout << "partialBodyFixedWrtGlobalState" << std::endl;
        // std::cout << partialBodyFixedWrtGlobalState << std::endl;

        // std::cout << "wrt state " << wrtSphericalBodyFixedState( ) << std::endl;
        // std::cout << "partial " << sphericalWrtCartesianBodyFixedState( ) * bodyFixedWrtGlobalState( ) << std::endl;
        // std::cout << "currentPartialWrtState_ " << currentPartialWrtState_ << std::endl;

        

        // // Compute partial w.r.t. position in inertial frame
        // currentPartialWrtPosition_ +=
        //         currentRotationToBodyFixedFrame_.inverse( ) * currentBodyFixedPartialWrtPosition_ * currentRotationToBodyFixedFrame_;

        // // If rotation matrix depends on translational state, add correction partials
        // if( rotationMatrixPartials_.count( std::make_pair( estimatable_parameters::initial_body_state, "" ) ) > 0 )
        // {
        //     // Compute rotation matrix partials
        //     std::vector< Eigen::Matrix3d > rotationMatrixPartials =
        //             rotationMatrixPartials_.at( std::make_pair( estimatable_parameters::initial_body_state, "" ) )
        //                     ->calculatePartialOfRotationMatrixToBaseFrameWrParameter( currentTime );

        //     Eigen::Vector6d currentInertialState = deformationModel_->getStateOfDeformingBody( ); 
        //     Eigen::Vector3d currentInertialPosition = currentInertialState.segment( 0, 3 ); //Eigen::Vector3d::Zero();;
        //     Eigen::Vector3d currentInertialVelocity = currentInertialState.segment( 3, 3 ); //Eigen::Vector3d::Zero();
        //     // std::cout << "currentInertialState " << currentInertialState.transpose( ) << std::endl;

        //     std::vector< Eigen::Matrix3d > rotationMatrixDerivativePartials =
        //             rotationMatrixPartials_.at( std::make_pair( estimatable_parameters::initial_body_state, "" ) )
        //                     ->calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter( currentTime );

        //     // Add correction terms to position and velocity partials
        //     for( unsigned int i = 0; i < 6; i++ )
        //     {
        //         // d R / d x,y,z,vx,vy,vz
        //         Eigen::Matrix3d currentRotationMatrixPartial = rotationMatrixPartials[ i ];

        //         Eigen::Vector3d currentPositionPartial = rotationMatrixPartials[ i ].transpose() * currentInertialPosition;
        //         partialBodyFixedWrtGlobalState.block( 0, i, 3, 1 ) += currentPositionPartial;

        //         Eigen::Vector3d currentVelocityPartial = rotationMatrixPartials[ i ].transpose() * currentInertialVelocity;
        //         partialBodyFixedWrtGlobalState.block( 3, i, 3, 1 ) += currentVelocityPartial;

        //         currentVelocityPartial = rotationMatrixDerivativePartials[ i ].transpose() * currentInertialPosition;
        //         partialBodyFixedWrtGlobalState.block( 3, i, 3, 1 ) += currentVelocityPartial;

        //         // std::cout << i << " : " << std::endl;
        //         // // std::cout << rotationMatrixPartials[ i ] << std::endl;
        //         // std::cout << "pos " << currentPositionPartial.transpose() << std::endl; 
        //         // std::cout << "vel " << currentVelocityPartial.transpose() << std::endl; 
        //     }
        //     // std::cout << "partialBodyFixedWrtGlobalState - after" << std::endl;
        //     // std::cout << partialBodyFixedWrtGlobalState << std::endl;
        // }

        currentPartialWrtDeformingState_ = - deformationWrtSphericalBodyFixedState( ) * sphericalWrtCartesianBodyFixedState( ) 
            * bodyFixedWrtGlobalState( currentTime, true );

        currentPartialWrtPerturbingState_ = deformationWrtSphericalBodyFixedState( ) * sphericalWrtCartesianBodyFixedState( ) 
            * bodyFixedWrtGlobalState( currentTime, false );

        currentTime_ = currentTime;

        // // // Update tidal interfaces
        // // for( unsigned int i = 0; i < tidalLoveNumberPartialInterfaces_.size( ); i++ )
        // // {
        // //     tidalLoveNumberPartialInterfaces_.at( i )->update( currentTime );
        // // }
    }
}


//! Function to calculate the partial of the acceleration wrt a set of cosine coefficients.
void MaxwellDeformationPartial::wrtCosineCoefficientBlock( const std::vector< std::pair< int, int > >& blockIndices,
                                                           Eigen::MatrixXd& partialDerivatives )
{
    partialDerivatives.setZero( );

    int c20Index = -1;
    int c21Index = -1;
    int c22Index = -1;

    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        if( blockIndices.at( i ).first == 2 && blockIndices.at( i ).second == 0 )
        {
            c20Index = i;
        }

        if( blockIndices.at( i ).first == 2 && blockIndices.at( i ).second == 1 )
        {
            c21Index = i;
        }

        if( blockIndices.at( i ).first == 2 && blockIndices.at( i ).second == 2 )
        {
            c22Index = i;
        }
    }

    partialDerivatives( 0, c20Index ) =  - 1.0 / deformationModel_->getGlobalRelaxationTime( ) 
        * basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
    partialDerivatives( 1, c21Index )  = - 1.0 / deformationModel_->getGlobalRelaxationTime( ) 
        * basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
    partialDerivatives( 2, c22Index )  = - 1.0 / deformationModel_->getGlobalRelaxationTime( ) 
        * basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
}

//! Function to calculate the partial of the acceleration wrt a set of sine coefficients.
void MaxwellDeformationPartial::wrtSineCoefficientBlock( const std::vector< std::pair< int, int > >& blockIndices,
                                                         Eigen::MatrixXd& partialDerivatives )
{
    partialDerivatives.setZero( );

    int s21Index = -1;
    int s22Index = -1;

    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        if( blockIndices.at( i ).first == 2 && blockIndices.at( i ).second == 1 )
        {
            s21Index = i;
        }

        if( blockIndices.at( i ).first == 2 && blockIndices.at( i ).second == 2 )
        {
            s22Index = i;
        }
    }

    partialDerivatives( 3, s21Index )  = - 1.0 / deformationModel_->getGlobalRelaxationTime( ) 
        * basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
    partialDerivatives( 4, s22Index )  = - 1.0 / deformationModel_->getGlobalRelaxationTime( ) 
        * basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
}

Eigen::MatrixXd MaxwellDeformationPartial::bodyFixedWrtRotational( )
{
    // std::cout << "in body fixed wrt rotation" << std::endl;

    // Get rotation matrix partial(s) wrt requested parameter
    std::vector< Eigen::Matrix3d > rotationMatrixPartials =
            rotationMatrixPartials_.at( std::make_pair( estimatable_parameters::initial_rotational_body_state, "" ) )
                    ->calculatePartialOfRotationMatrixToBaseFrameWrParameter( currentTime_ );

    // for ( unsigned int i = 0 ; i < 4 ; i++ )
    // {
    //     std::cout << "partial R" << std::endl;
    //     std::cout << rotationMatrixPartials[i] << std::endl;
    // }

    Eigen::Vector6d currentInertialState = - deformationModel_->getStateOfDeformingBody( ); 
    Eigen::Vector3d currentInertialPosition = currentInertialState.segment( 0, 3 ); //Eigen::Vector3d::Zero();;
    Eigen::Vector3d currentInertialVelocity = currentInertialState.segment( 3, 3 );
    Eigen::Vector3d angularVelocityBodyFixedFrame = deformationModel_->getCurrentAngularVelocityDeformingBody();
    Eigen::Vector3d currentBodyFixedPosition = deformationModel_->getCurrentRelativePosition( )[0]; // TO BE MODIFIED -> INDEX 0

    Eigen::MatrixXd tempPartial = Eigen::MatrixXd::Zero( 6, 7 );

    // Add correction terms to position and velocity partials
    for( unsigned int i = 0; i < 4; i++ )
    {
        // d R / d q0, q1, q2, q3
        Eigen::Matrix3d currentRotationMatrixPartial = rotationMatrixPartials[ i ].transpose();

        Eigen::Vector3d temp = currentRotationMatrixPartial * currentInertialPosition;
        tempPartial.block( 0, i, 3, 1) = temp;
    
        Eigen::Vector3d term = currentRotationMatrixPartial * currentInertialVelocity;
        Eigen::Vector3d crossTerm = angularVelocityBodyFixedFrame.cross( currentRotationMatrixPartial * currentInertialPosition );
        tempPartial.block(3, i, 3, 1) = term - crossTerm; 
    }

    // std::cout << "body fixed wrt quaternion" << std::endl;
    // std::cout << tempPartial.block( 0, 0, 3, 4) << std::endl;

    // std::cout << "currentBodyFixedPosition " << currentBodyFixedPosition.transpose( ) << std::endl;
    // std::cout << "currentInertialPosition " << currentInertialPosition.transpose( ) << std::endl;
    // std::cout << "rotation matrix from body-fixed to inertial " << std::endl;
    // std::cout << deformationModel_->getCurrentRotationToIntegrationFrameMatrix( ) << std::endl;

    Eigen::Matrix3d drb_domega = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d drdotb_domega = linear_algebra::getCrossProductMatrix(currentBodyFixedPosition); // [r_b]_x
    tempPartial.block(3, 4, 3, 3) = drdotb_domega;

    return tempPartial;
}

Eigen::MatrixXd MaxwellDeformationPartial::bodyFixedWrtRotationParameter(
    const estimatable_parameters::EstimatebleParametersEnum parameterType,
    const std::string& secondaryIdentifier )
{

    // Compute rotation matrix partials
    std::vector< Eigen::Matrix3d > rotationMatrixPartials =
        rotationMatrixPartials_.at( std::make_pair( parameterType, secondaryIdentifier ) )
                ->calculatePartialOfRotationMatrixToBaseFrameWrParameter( currentTime_ );

    // Compute rotation matrix derivative partials
    std::vector< Eigen::Matrix3d > rotationMatrixDerivativePartials =
            rotationMatrixPartials_.at( std::make_pair( parameterType, secondaryIdentifier ) )
                    ->calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter( currentTime_ );

    Eigen::Vector6d currentInertialState = - deformationModel_->getStateOfDeformingBody( ); 
    Eigen::Vector3d currentInertialPosition = currentInertialState.segment( 0, 3 ); 
    Eigen::Vector3d currentInertialVelocity = currentInertialState.segment( 3, 3 );

    
    int parameterSize = rotationMatrixPartials.size( );           
    Eigen::Matrix< double, 6, Eigen::Dynamic > tempPartial = Eigen::Matrix< double, 6, Eigen::Dynamic >::Zero( 6, parameterSize );
                        
        for( int i = 0 ; i < parameterSize ; i++ )
        {
            Eigen::Matrix3d currentRotationMatrixPartial = rotationMatrixPartials[ i ].transpose();
            Eigen::Matrix3d currentRotationMatrixDerivativePartial = rotationMatrixDerivativePartials[ i ].transpose();

            tempPartial.block( 0, i, 3, 1) = currentRotationMatrixPartial * currentInertialPosition;
            tempPartial.block( 3, i, 3, 1) = currentRotationMatrixDerivativePartial * currentInertialPosition 
                + currentRotationMatrixPartial * currentInertialVelocity;

            // std::cout << "in deformation partial: currentRotationMatrixPartial" << std::endl;
            // std::cout << rotationMatrixPartials[ i ] << std::endl;
            // std::cout << "in deformation partial: currentRotationMatrixDerivativePartial" << std::endl;
            // std::cout << currentRotationMatrixDerivativePartial << std::endl;
        }

        // std::cout << "tempPartial" << std::endl;
        // std::cout << tempPartial.transpose() << std::endl;

        return tempPartial;
}

//! Function to calculate an acceleration partial wrt a rotational parameter.
void MaxwellDeformationPartial::wrtRotationModelParameter( 
    Eigen::MatrixXd& accelerationPartial,
    const estimatable_parameters::EstimatebleParametersEnum parameterType,
    const std::string& secondaryIdentifier )
{
    // std::cout << "wrtSphericalBodyFixedState " << std::endl;
    // std::cout << wrtSphericalBodyFixedState( ) << std::endl;
    // std::cout << "sphericalWrtCartesianBodyFixedState " << std::endl;
    // std::cout << sphericalWrtCartesianBodyFixedState( ) << std::endl;
    // std::cout << "bodyFixedWrtRotational " << std::endl;
    // std::cout << bodyFixedWrtRotational() << std::endl;

    if ( parameterType == estimatable_parameters::initial_rotational_body_state )
    {
        accelerationPartial = deformationWrtSphericalBodyFixedState( ) * sphericalWrtCartesianBodyFixedState( ) * bodyFixedWrtRotational( )
            + deformationWrtRotationalState( );
    }
    else
    {
        accelerationPartial = deformationWrtSphericalBodyFixedState( ) 
            * sphericalWrtCartesianBodyFixedState( ) * bodyFixedWrtRotationParameter( parameterType, secondaryIdentifier );
        if ( parameterType == estimatable_parameters::constant_rotation_rate )
        {
            accelerationPartial += deformationWrtRotationRate( );
        }
    }

    

    // std::cout << " deformation wrt (q,w) " << std::endl;
    // std::cout << accelerationPartial << std::endl;

    // // Iterate for each single parameter entry partial.
    // for( unsigned int i = 0; i < rotationMatrixPartials.size( ); i++ )
    // {
    //     // Calculate acceleration partial for current parameter entry.
    //     accelerationPartial.block( 0, i, 3, 1 ) =
    //             rotationMatrixPartials[ i ] * (currentRotationToBodyFixedFrame_)*accelerationModel_->getAcceleration( ) +
    //             currentRotationToInertialFrame_ * currentBodyFixedPartialWrtPosition_ * rotationMatrixPartials[ i ].transpose( ) *
    //                     distanceVector;
    // }
    // std::cout << "end wrtRotationModelParameter" << std::endl;
}


}  // namespace acceleration_partials

}  // namespace tudat
