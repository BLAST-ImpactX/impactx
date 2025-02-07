/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "Aperture.H"

#include <stdexcept>
#include <string>


std::string
impactx::elements::Aperture::shape_name (Shape const & shape)
{
    switch (shape)
    {
        case Aperture::Shape::rectangular :  // default
            return "rectangular";
        case Aperture::Shape::elliptical :
            return "elliptical";
        default:
            throw std::runtime_error("Unknown shape");
    }
}

std::string
impactx::elements::Aperture::action_name (Action const & action)
{
    switch (action)
    {
        case Aperture::Action::transmit :  // default
            return "transmit";
        case Aperture::Action::absorb :
            return "absorb";
        default:
            throw std::runtime_error("Unknown action");
    }
}
