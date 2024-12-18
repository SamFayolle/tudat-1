/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GRAVITYDEFORMATIONSETTINGS_H
#define TUDAT_GRAVITYDEFORMATIONSETTINGS_H

#include <functional>
#include <memory>
#include <map>
#include <string>

namespace tudat
{

namespace simulation_setup
{

enum GravityDeformationType
{
    maxwell_deformation = 0
};

// Class for providing settings for gravity deformation models.
/*
*  Class for providing settings for gravity deformation models.
*/
class GravityDeformationSettings
{
public:

// Constructor, sets type of deformation.
/*
    *  Constructor, sets type of deformation.
    *  \param deformationType Type of acceleration from GravityDeformationType enum.
    */
GravityDeformationSettings( const simulation_setup::GravityDeformationType deformationType ):
    deformationType_( deformationType ){ }

// Destructor.
virtual ~GravityDeformationSettings( ){ }

// Type of acceleration from AvailableAcceleration enum.
GravityDeformationType deformationType_;

};

class MaxwellDeformationSettings: public GravityDeformationSettings
{
public:

// Constructor, sets type of acceleration.
/*
    *  Constructor, sets type of acceleration.
    *  \param accelerationType Type of acceleration from AvailableAcceleration enum.
    */
MaxwellDeformationSettings( 
    const double maxwellRelaxationTime, 
    const double globalRelaxationTime,
    const double loveNumber,
    const double rotationRate,
    const int maximumDegree,
    const int maximumOrder,
    const std::string perturbingBody ):
    GravityDeformationSettings( maxwell_deformation ), maxwellRelaxationTime_( maxwellRelaxationTime ),
    globalRelaxationTime_( globalRelaxationTime ), loveNumber_( loveNumber ), rotationRate_( rotationRate ),
    maximumDegree_( maximumDegree ), maximumOrder_( maximumOrder ), perturbingBody_( perturbingBody )
    {
        
    }

// Destructor.
virtual ~MaxwellDeformationSettings( ){ }

const double maxwellRelaxationTime_;
const double globalRelaxationTime_;
const double loveNumber_;
const double rotationRate_;
const int maximumDegree_;
const int maximumOrder_;
const std::string perturbingBody_;

};

typedef std::map< std::string, std::vector< std::shared_ptr< GravityDeformationSettings > > > SelectedGravityDeformationModelMap;


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_GRAVITYDEFORMATIONSETTINGS_H
