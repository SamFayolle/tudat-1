 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

# Set the source files.
set(gravity_deformation_partials_SOURCES
  "maxwellDeformationPartial.cpp"
)

# Set the header files.
set(gravity_deformation_partials_HEADERS
  "deformationPartial.h"
  "maxwellDeformationPartial.h"
)

TUDAT_ADD_LIBRARY("gravity_deformation_partials"
        "${gravity_deformation_partials_SOURCES}"
        "${gravity_deformation_partials_HEADERS}"
        )

