/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *
 *    Notes:
 *      Test tolerance was set at 5.0e-15 (or 5.0e-7 for floats) instead of epsilon due to
 *      rounding errors in Eigen types with entries over a number of orders of magnitude,
 *      presumably causing the observed larger than epsilon relative differences between
 *      expected and computed values.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <vector>





#include <memory>
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/basic_astro/gravityDeformationModel.h"
#include "tudat/simulation/basic_astro/testAccelerationModels.h"
#include "tudat/astro/basic_astro/testBody.h"
#include "tudat/basics/basicTypedefs.h"

namespace tudat
{
namespace unit_tests
{

using basic_astrodynamics::GravityDeformationModel;


BOOST_AUTO_TEST_CASE( test_gravityDeformationModel )
{
   
}


} // namespace unit_tests
} // namespace tudat
