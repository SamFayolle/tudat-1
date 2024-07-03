/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/estimatable_parameters/parameterConstraints.h"

namespace tudat
{

namespace estimatable_parameters
{

std::string getConstraintTypeString( const ParameterConstraintsEnum constraintType )
{
    std::string constraintDescription;
    switch( constraintType )
    {
        case tidal_quality_factor_single_love_number_constraint:
            constraintDescription = "constraint between tidal quality factor and single Love number ";
            break;
        case tidal_time_lag_single_love_number_constraint:
            constraintDescription = "constraint between tidal time lag and single Love number ";
            break;
        case tidal_quality_factor_full_love_number_constraint:
            constraintDescription = "constraint between tidal quality factor and full Love number ";
            break;
        case tidal_time_lag_full_love_number_constraint:
            constraintDescription = "constraint between tidal time lag and full Love number ";
            break;
        case custom_constraint:
            constraintDescription = "custom constraint ";
            break;
        default:
            std::string errorMessage = "Error when getting constraint string, did not recognize constraint " +
                                       std::to_string( constraintType );
            throw std::runtime_error( errorMessage );
    }
    return constraintDescription;
}

}

}

