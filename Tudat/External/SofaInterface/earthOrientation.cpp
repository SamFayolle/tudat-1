#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "External/SofaInterface/earthOrientation.h"


namespace tudat
{

namespace sofa_interface
{

//! Function to calculate CIP and CIO locator according to requested IAU conventions
std::pair< Eigen::Vector2d, double > getPositionOfCipInGcrs(
        const double terrestrialTime, const double julianDaysEpochShift, const IAUConventions precessionNutationTheory )
{
    // Declare Sofa function return arguments (by reference)
    double xAngle, yAngle;
    double originLocator;

    // Check for IAU convention and retrieve requested values.
    switch( precessionNutationTheory )
    {
    case iau_2000_a:
        iauXys00a( julianDaysEpochShift, terrestrialTime / physical_constants::JULIAN_DAY,
                   &xAngle, &yAngle, &originLocator );
        break;

    case iau_2000_b:
        iauXys00b( julianDaysEpochShift, terrestrialTime / physical_constants::JULIAN_DAY,
                   &xAngle, &yAngle, &originLocator );
        break;

    case iau_2006:
        iauXys06a( julianDaysEpochShift, terrestrialTime / physical_constants::JULIAN_DAY,
                   &xAngle, &yAngle, &originLocator );
        break;
    default:
        std::cerr<<"Warning, precession nutation theory selection not recongnized"<<std::endl;

    }

    // Set and return requested values.
    Eigen::Vector2d cioPosition;
    cioPosition << xAngle, yAngle;
    return std::pair< Eigen::Vector2d, double >( cioPosition, originLocator );
}

double calculateGreenwichMeanSiderealTime(
        const double terrestrialTime, const double universalTime1JulianDaysSinceJ2000, const IAUConventions iauConvention )
{
    // Declare GMST variable
    double gmst;

    // Check for IAU convention and retrieve requested GMST
    switch( iauConvention )
    {
    case iau_2000_a:
        gmst = iauGmst00( basic_astrodynamics::JULIAN_DAY_ON_J2000, terrestrialTime / physical_constants::JULIAN_DAY,
                          basic_astrodynamics::JULIAN_DAY_ON_J2000, universalTime1JulianDaysSinceJ2000 / physical_constants::JULIAN_DAY);
        break;

    case iau_2000_b:
        gmst = iauGmst00( basic_astrodynamics::JULIAN_DAY_ON_J2000, terrestrialTime / physical_constants::JULIAN_DAY,
                          basic_astrodynamics::JULIAN_DAY_ON_J2000, universalTime1JulianDaysSinceJ2000 / physical_constants::JULIAN_DAY);
        break;

    case iau_2006:
        gmst = iauGmst06( basic_astrodynamics::JULIAN_DAY_ON_J2000, terrestrialTime / physical_constants::JULIAN_DAY,
                          basic_astrodynamics::JULIAN_DAY_ON_J2000, universalTime1JulianDaysSinceJ2000 / physical_constants::JULIAN_DAY );
        break;
    default:
        std::cerr<<"Warning, iau convention for GMST calculation not recongnized"<<std::endl;

    }

    return gmst;

}

//! Function to calculate ERA (earth rotation angle)
double calculateEarthRotationAngle( const double ut1, const double julianDaysEpochShift )
{
    return iauEra00Unnormalized( julianDaysEpochShift, ut1 / physical_constants::JULIAN_DAY );
}

}

}
