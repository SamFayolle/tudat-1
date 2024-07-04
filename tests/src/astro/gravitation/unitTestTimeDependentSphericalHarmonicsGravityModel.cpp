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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/lambda/lambda.hpp>
#include <boost/make_shared.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"
#include "tudat/math/basic/sphericalHarmonics.h"
#include "tudat/astro/gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/environment_setup/createGravityFieldVariations.h"
#include "tudat/simulation/environment_setup.h"
#include "tudat/simulation/estimation_setup.h"

namespace tudat
{
namespace unit_tests
{
using namespace tudat::gravitation;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::estimatable_parameters;

BOOST_AUTO_TEST_SUITE( test_TimeDependentSphericalHarmonicsGravity )

void addGravityFieldVariationsMoons( BodyListSettings& bodySettings )
{
    // (Complex) Love numbers
    std::map< int, std::vector< std::complex< double > > > loveNumbersIo, loveNumbersEuropa, loveNumbersGanymede, loveNumbersCallisto;

    double k2Moons = 0.2;
    double k2ImIo = 0.01;
    double k2ImMoons = 0.0;

    std::vector< std::complex< double > > loveNumbersIoList;
    loveNumbersIoList.push_back( std::complex< double >( 0.0, 0.0 ) );
    loveNumbersIoList.push_back( std::complex< double >( 0.0, 0.0 ) );
    loveNumbersIoList.push_back( std::complex< double >( k2Moons, k2ImIo ) );
    loveNumbersIo[ 2 ] = loveNumbersIoList;

    std::vector< std::complex< double > > loveNumbersEuropaList;
    loveNumbersEuropaList.push_back( std::complex< double >( 0.0, 0.0 ) );
    loveNumbersEuropaList.push_back( std::complex< double >( 0.0, 0.0 ) );
    loveNumbersEuropaList.push_back( std::complex< double >( k2Moons, k2ImMoons ) );
    loveNumbersEuropa[ 2 ] = loveNumbersEuropaList;

    std::vector< std::complex< double > > loveNumbersGanymedeList;
    loveNumbersGanymedeList.push_back( std::complex< double >( 0.0, 0.0 ) );
    loveNumbersGanymedeList.push_back( std::complex< double >( 0.0, 0.0 ) );
    loveNumbersGanymedeList.push_back( std::complex< double >( k2Moons, k2ImMoons ) );
    loveNumbersGanymede[ 2 ] = loveNumbersGanymedeList;

    std::vector< std::complex< double > > loveNumbersCallistoList;
    loveNumbersCallistoList.push_back( std::complex< double >( 0.0, 0.0 ) );
    loveNumbersCallistoList.push_back( std::complex< double >( 0.0, 0.0 ) );
    loveNumbersCallistoList.push_back( std::complex< double >( k2Moons, k2ImMoons ) );
    loveNumbersCallisto[ 2 ] = loveNumbersCallistoList;

    // Gravity field variations
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationsIo, gravityFieldVariationsEuropa, gravityFieldVariationsGanymede,
            gravityFieldVariationsCallisto;

    std::vector< std::string > deformingBodiesMoons = { "Jupiter" };
    // For Io
    gravityFieldVariationsIo.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >( deformingBodiesMoons, loveNumbersIo ) );
    // For Europa
    gravityFieldVariationsEuropa.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >( deformingBodiesMoons, loveNumbersEuropa ) );
    // For Ganymede
    gravityFieldVariationsGanymede.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >( deformingBodiesMoons, loveNumbersGanymede ) );
    // For Callisto
    gravityFieldVariationsCallisto.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >( deformingBodiesMoons, loveNumbersCallisto ) );

    // Define gravity field variations settings
    bodySettings.at( "Io" )->gravityFieldVariationSettings = gravityFieldVariationsIo;
    bodySettings.at( "Europa" )->gravityFieldVariationSettings = gravityFieldVariationsEuropa;
    bodySettings.at( "Ganymede" )->gravityFieldVariationSettings = gravityFieldVariationsGanymede;
    bodySettings.at( "Callisto" )->gravityFieldVariationSettings = gravityFieldVariationsCallisto;

}

void addGravityFieldVariationsJupiter( BodyListSettings& bodySettings )
{
    // (Complex) Love numbers
    std::map< int, std::vector< std::complex< double > > > loveNumbersJupiter;

    double k2Jupiter = 0.5;
    double k2ImJupiter = -1.102e-5;

    std::vector< std::complex< double > > loveNumbersJupiterList;
    loveNumbersJupiterList.push_back( std::complex< double >( 0.0, 0.0 ) );
    loveNumbersJupiterList.push_back( std::complex< double >( 0.0, 0.0 ) );
    loveNumbersJupiterList.push_back( std::complex< double >( k2Jupiter, k2ImJupiter ) );
    loveNumbersJupiter[ 2 ] = loveNumbersJupiterList;

    // Gravity field variations
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationsJupiter;
    gravityFieldVariationsJupiter.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
            std::vector< std::string >( { "Io" } ), loveNumbersJupiter ) );
    gravityFieldVariationsJupiter.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
            std::vector< std::string >( { "Europa" } ), loveNumbersJupiter ) );
    gravityFieldVariationsJupiter.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
            std::vector< std::string >( { "Ganymede" } ), loveNumbersJupiter ) );
    gravityFieldVariationsJupiter.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
            std::vector< std::string >( { "Callisto" } ), loveNumbersJupiter ) );

    // Define gravity field variations settings
    bodySettings.at( "Jupiter" )->gravityFieldVariationSettings = gravityFieldVariationsJupiter;
}


BOOST_AUTO_TEST_CASE( testTimeVaryingGravityField )
{
    loadStandardSpiceKernels( );
    std::cout.precision( 20 );

    double initialEpoch = 32.0 * physical_constants::JULIAN_YEAR;
    double finalEpoch = initialEpoch + 10.0 * 86400.0;

    std::vector< std::map< double, Eigen::VectorXd > > accelerationsHistoryAllSetups;
    std::vector< std::map< double, Eigen::MatrixXd > > stateTransitionMatrixHistoryAllSetups, sensitivityMatrixHistoryAllSetups;

    for ( unsigned int setup = 0 ; setup < 5 ; setup++ )
    {
        // Define bodies settings for simulation
        std::vector< std::string > bodiesToCreate = { "Io", "Europa", "Ganymede", "Callisto", "Jupiter" };

        // Create body objects.
        BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, initialEpoch, finalEpoch, "Jupiter", "J2000" );
        bodySettings.at( "Io" )->rotationModelSettings = std::make_shared< SynchronousRotationModelSettings >(
                "Jupiter", "J2000", "IAU_Io" );
        bodySettings.at( "Io" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
                initialEpoch, finalEpoch, 3600.0, "Jupiter", "J2000" );

        bodySettings.at( "Europa" )->rotationModelSettings = std::make_shared< SynchronousRotationModelSettings >(
                "Jupiter", "J2000", "IAU_Europa" );
        bodySettings.at( "Europa" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
                initialEpoch, finalEpoch, 3600.0, "Jupiter", "J2000" );

        bodySettings.at( "Ganymede" )->rotationModelSettings = std::make_shared< SynchronousRotationModelSettings >(
                "Jupiter", "J2000", "IAU_Ganymede" );
        bodySettings.at( "Ganymede" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
                initialEpoch, finalEpoch, 3600.0, "Jupiter", "J2000" );

        bodySettings.at( "Callisto" )->rotationModelSettings = std::make_shared< SynchronousRotationModelSettings >(
                "Jupiter", "J2000", "IAU_Callisto" );
        bodySettings.at( "Callisto" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
                initialEpoch, finalEpoch, 3600.0, "Jupiter", "J2000" );

        if ( setup == 1 || setup == 3 || setup == 4 )
        {
            addGravityFieldVariationsMoons( bodySettings );
            addGravityFieldVariationsJupiter( bodySettings );
        }
        else if( setup == 2 )
        {
            addGravityFieldVariationsJupiter( bodySettings );
        }

        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        // Define bodies to propagate
        std::vector< std::string > bodiesToPropagate = { "Io", "Europa", "Ganymede", "Callisto" };
        std::vector< std::string > centralBodies = { "Jupiter", "Jupiter", "Jupiter", "Jupiter" };

        // Define accelerations
        int maximumSatelliteDegree = 2;
        int maximumSatelliteOrder = 2;
        int maximumJupiterDegree = 10;
        int maximumJupiterOrder = 10;

        SelectedAccelerationMap accelerationSettings;
        for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
        {
            if ( setup == 0 || setup == 1 || setup == 4 )
            {
                accelerationSettings[ bodiesToPropagate.at( i ) ][ "Jupiter" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >(
                        maximumJupiterDegree, maximumJupiterOrder, maximumSatelliteDegree, maximumSatelliteOrder, false ) );
            }
            else
            {
                accelerationSettings[ bodiesToPropagate.at( i ) ][ "Jupiter" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >(
                        maximumJupiterDegree, maximumJupiterOrder, maximumSatelliteDegree, maximumSatelliteOrder, true ) );
            }


            for( unsigned int j = 0; j < bodiesToPropagate.size( ); j++ )
            {
                if( bodiesToPropagate.at( j ) != bodiesToPropagate.at( i ) )
                {
                    if ( setup == 0 || setup == 1 || setup == 3 || setup == 4 )
                    {
                        accelerationSettings[ bodiesToPropagate.at( i ) ][ bodiesToPropagate.at( j ) ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >(
                                maximumSatelliteDegree, maximumSatelliteOrder, false, false ) );
                    }
                    else
                    {
                        accelerationSettings[ bodiesToPropagate.at( i ) ][ bodiesToPropagate.at( j ) ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >(
                                maximumSatelliteDegree, maximumSatelliteOrder, false, true ) );
                    }

                }
            }
        }
        basic_astrodynamics::AccelerationMap accelerations = createAccelerationModelsMap(
                bodies, accelerationSettings, bodiesToPropagate, centralBodies );

        // Define propagator settings
        Eigen::VectorXd initialStates = getInitialStatesOfBodies( bodiesToPropagate, centralBodies, bodies, initialEpoch );

        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesSettings;
        dependentVariablesSettings.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                basic_astrodynamics::spherical_harmonic_gravity, "Io", "Europa" ) );
        dependentVariablesSettings.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                basic_astrodynamics::spherical_harmonic_gravity, "Io", "Ganymede" ) );
        dependentVariablesSettings.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                basic_astrodynamics::spherical_harmonic_gravity, "Io", "Callisto" ) );
        std::shared_ptr< DependentVariableSaveSettings > dependentVariables = std::make_shared< DependentVariableSaveSettings >( dependentVariablesSettings );

        std::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< > >(
                        centralBodies, accelerations, bodiesToPropagate, initialStates, finalEpoch, cowell, dependentVariables );

        // Define integrator settings
        double timeStep = 2700.0;
        std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
                std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                        initialEpoch, timeStep, numerical_integrators::RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78, timeStep, timeStep, 1.0e3, 1.0e3 );

        // Propagate dynamics
        SingleArcDynamicsSimulator< > dynamicsSimulator = SingleArcDynamicsSimulator< >( bodies, integratorSettings, propagatorSettings, true, false, false );

        accelerationsHistoryAllSetups.push_back( dynamicsSimulator.getDependentVariableHistory( ) );
//        std::map< double, Eigen::VectorXd > accelerationsHistory = dynamicsSimulator.getDependentVariableHistory( );


        // Create parameters to estimate
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
        for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
        {
            parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                    bodiesToPropagate.at( i ), initialStates.segment( i * 6, 6 ), "Jupiter" ) );
        }

        if ( setup < 4 )
        {
            for ( unsigned int i = 0 ; i < bodiesToPropagate.size( ) ; i++ )
            {
                parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                        2, 0, maximumSatelliteDegree, maximumSatelliteOrder, bodiesToPropagate.at( i ), spherical_harmonics_cosine_coefficient_block ) );
                parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                        2, 1,  maximumSatelliteDegree, maximumSatelliteOrder, bodiesToPropagate.at( i ), spherical_harmonics_sine_coefficient_block ) );
            }

            for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
            {
                parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( bodiesToPropagate.at( i ), gravitational_parameter ) );
            }
        }
        else
        {
            for ( unsigned int i = 0 ; i < bodiesToPropagate.size( ) ; i++ )
            {
                parameterNames.push_back( std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
                        bodiesToPropagate.at( i ), 2, std::vector< int >{ 2 }, "Jupiter", true ) );
                parameterNames.push_back( std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
                        "Jupiter", 2, std::vector< int >{ 2 }, bodiesToPropagate.at( i ), true ) );
            }
        }

        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double >( parameterNames, bodies, propagatorSettings );

        // Propagate variational equations
        SingleArcVariationalEquationsSolver< > variationalEquationsSolver = SingleArcVariationalEquationsSolver< >(
                bodies, integratorSettings, propagatorSettings, parametersToEstimate, true,
                std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), false, true, false );

        std::vector< std::map< double, Eigen::MatrixXd > > variationalEquationsSolution = variationalEquationsSolver.getNumericalVariationalEquationsSolution( );
        stateTransitionMatrixHistoryAllSetups.push_back( variationalEquationsSolution[ 0 ] );
        sensitivityMatrixHistoryAllSetups.push_back( variationalEquationsSolution[ 1 ] );

    }

    for ( double time = initialEpoch ; time < finalEpoch ; time += 20.0 * 2700.0 )
    {
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( accelerationsHistoryAllSetups[ 0 ].at( time ), accelerationsHistoryAllSetups[ 1 ].at( time ), 1.0e-15 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( accelerationsHistoryAllSetups[ 2 ].at( time ), accelerationsHistoryAllSetups[ 3 ].at( time ), 1.0e-15 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( stateTransitionMatrixHistoryAllSetups[ 0 ].at( time ), stateTransitionMatrixHistoryAllSetups[ 1 ].at( time ), 1.0e-15 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sensitivityMatrixHistoryAllSetups[ 0 ].at( time ), sensitivityMatrixHistoryAllSetups[ 1 ].at( time ), 1.0e-15 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( stateTransitionMatrixHistoryAllSetups[ 2 ].at( time ), stateTransitionMatrixHistoryAllSetups[ 3 ].at( time ), 1.0e-15 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sensitivityMatrixHistoryAllSetups[ 2 ].at( time ), sensitivityMatrixHistoryAllSetups[ 3 ].at( time ), 1.0e-15 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sensitivityMatrixHistoryAllSetups[ 4 ].at( time ), Eigen::MatrixXd::Zero( 24, 16 ), 1.0e-15 );
    }








//
//
//
//
//        SystemOfBodies perturbedBodies = createSystemOfBodies( bodySettings );
//
//        SelectedAccelerationMap accelerationSettingsPerturbed;
//        for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
//        {
//            accelerationSettingsPerturbed[ bodiesToPropagate.at( i ) ][ "Jupiter" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >(
//                    maximumJupiterDegree, maximumJupiterOrder, maximumSatelliteDegree, maximumSatelliteOrder, 0, 0, false, false, false ) );
//
//            for( unsigned int j = 0; j < bodiesToPropagate.size( ); j++ )
//            {
//                if( bodiesToPropagate.at( j ) != bodiesToPropagate.at( i ) )
//                {
//                    accelerationSettingsPerturbed[ bodiesToPropagate.at( i ) ][ bodiesToPropagate.at( j ) ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >(
//                            maximumSatelliteDegree, maximumSatelliteOrder, false ) );
//                }
//            }
//        }
//        basic_astrodynamics::AccelerationMap accelerationsPerturbed = createAccelerationModelsMap(
//                perturbedBodies, accelerationSettingsPerturbed, bodiesToPropagate, centralBodies );
//
//        // Define propagator settings
//        Eigen::VectorXd initialStatesPerturbed = getInitialStatesOfBodies( bodiesToPropagate, centralBodies, perturbedBodies, initialEpoch );
//
//        std::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettingsPerturbed =
//                std::make_shared< TranslationalStatePropagatorSettings< > >(
//                        centralBodies, accelerationsPerturbed, bodiesToPropagate, initialStatesPerturbed, finalEpoch, cowell, dependentVariables );
//
//        // Propagate dynamics
//        SingleArcDynamicsSimulator< > dynamicsSimulatorPerturbed = SingleArcDynamicsSimulator< >(
//                perturbedBodies, integratorSettings, propagatorSettingsPerturbed, true, false, false );
//
//        std::map< double, Eigen::VectorXd > accelerationsHistoryPerturbed = dynamicsSimulatorPerturbed.getDependentVariableHistory( );
//
////        for ( double time = initialEpoch ; time < finalEpoch ; time += 20.0 * timeStep )
////        {
////            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( accelerationsHistory.at( time ), accelerationsHistoryPerturbed.at( time ), 1.0e-15 );
////        }
//
////    for( auto itr : accelerationsHistory )
////    {
////        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( itr.second, accelerationsHistoryPerturbed.at( itr.first ), 1.0e-15 );
////    }
//
//        // Create parameters to estimate
//        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
//        for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
//        {
//            parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
//                    bodiesToPropagate.at( i ), initialStates.segment( i * 6, 6 ), "Jupiter" ) );
//        }
//
//        for ( unsigned int i = 0 ; i < bodiesToPropagate.size( ) ; i++ )
//        {
//            parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
//                    2, 0, maximumSatelliteDegree, maximumSatelliteOrder, bodiesToPropagate.at( i ), spherical_harmonics_cosine_coefficient_block ) );
//            parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
//                    2, 1,  maximumSatelliteDegree, maximumSatelliteOrder, bodiesToPropagate.at( i ), spherical_harmonics_sine_coefficient_block ) );
//        }
//
//        for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
//        {
//            parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( bodiesToPropagate.at( i ), gravitational_parameter ) );
//        }
//
//        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
//                createParametersToEstimate< double >( parameterNames, bodies, propagatorSettings );
//
//
//
//
//
//
//        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimatePerturbed =
//                createParametersToEstimate< double >( parameterNames, perturbedBodies, propagatorSettingsPerturbed );
//
//        SingleArcVariationalEquationsSolver< > variationalEquationsSolver = SingleArcVariationalEquationsSolver< >(
//                bodies, integratorSettings, propagatorSettings, parametersToEstimate, true,
//                std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), false, true, false );
//
//        std::vector< std::map< double, Eigen::MatrixXd > > variationalEquationsSolution = variationalEquationsSolver.getNumericalVariationalEquationsSolution( );
//        std::map< double, Eigen::MatrixXd > stateTransitionMatrices = variationalEquationsSolution[ 0 ];
//        std::map< double, Eigen::MatrixXd > sensitivityMatrices = variationalEquationsSolution[ 1 ];
//
//        SingleArcVariationalEquationsSolver< > variationalEquationsSolverPerturbed = SingleArcVariationalEquationsSolver< >(
//                perturbedBodies, integratorSettings, propagatorSettingsPerturbed, parametersToEstimatePerturbed, true,
//                std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), false, true, false );
//
//        std::vector< std::map< double, Eigen::MatrixXd > > variationalEquationsSolutionPerturbed =
//                variationalEquationsSolverPerturbed.getNumericalVariationalEquationsSolution( );
//        std::map< double, Eigen::MatrixXd > stateTransitionMatrices2 = variationalEquationsSolutionPerturbed[ 0 ];
//        std::map< double, Eigen::MatrixXd > sensitivityMatrices2 = variationalEquationsSolutionPerturbed[ 1 ];
//
//        std::cout << sensitivityMatrices.rbegin( )->second.block( 0, 0, 1, sensitivityMatrices.rbegin( )->second.cols( ) ) << "\n\n";
//        std::cout << sensitivityMatrices2.rbegin( )->second.block( 0, 0, 1, sensitivityMatrices2.rbegin( )->second.cols( ) ) << "\n\n";
//
//
//    //



}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat