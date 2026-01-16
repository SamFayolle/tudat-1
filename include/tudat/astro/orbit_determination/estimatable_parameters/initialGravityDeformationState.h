/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_INITIALGRAVITYDEFORMATIONSTATE_H
#define TUDAT_INITIALGRAVITYDEFORMATIONSTATE_H

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of an initial gravity deformation state.
template< typename InitialStateParameterType = double >
class InitialGravityDeformationStateParameter : public EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > >
{
public:
    //! Constructor, sets initial value of gravity deformation state.
    /*!
     * Constructor, sets initial value of gravity deformation state.
     * \param associatedBody Body for which initial state is to be estimated.
     * \param initialGravityDeformationState Current value of initial state 
     */
    InitialGravityDeformationStateParameter( const std::string& associatedBody,
                                        const Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >& initialGravityDeformationState ):
        EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > >( initial_gravity_deformation_state, associatedBody ),
        initialGravityDeformationState_( initialGravityDeformationState )
    { }

    //! Function to get the current value of initial state.
    /*!
     * Function to get the current value of initial state.
     * \return The current value of initial state.
     */
    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > getParameterValue( )
    {
        return initialGravityDeformationState_;
    }

    //! Function to reset the current value of initial state.
    /*!
     * Function to reset the current value of initial state.
     * \param parameterValue The new value of initial state.
     */
    void setParameterValue( Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > parameterValue )
    {
        initialGravityDeformationState_ = parameterValue;
    }

    //! Function to retrieve the size of the parameter (always set to 5 for now - i.e., degree two).
    /*!
     *  Function to retrieve the size of the parameter (always set to 5 for now).
     *  \return Size of parameter value (always set to 5 for now).
     */
    int getParameterSize( )
    {
        return 5;
    }


private:
    //! Current value of initial state 
    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialGravityDeformationState_;

};


}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_INITIALGRAVITYDEFORMATIONSTATE_H
