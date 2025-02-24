/* Copyright 2022-2025 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "initialization/Algorithms.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX_ParmParse.H>

#include <algorithm>  // for std::transform
#include <stdexcept>
#include <string>


namespace impactx
{
    SpaceChargeAlgo
    get_space_charge_algo ()
    {
        amrex::ParmParse const pp_algo("algo");
        std::string space_charge;
        bool has_space_charge = pp_algo.query("space_charge", space_charge);

        if (!has_space_charge) { return SpaceChargeAlgo::False; }

        // TODO: at some point, remove backwards compatibility to pre 25.03
        std::string space_charge_lower = space_charge;
        std::transform(space_charge.begin(), space_charge.end(), space_charge_lower.begin(),
               [](unsigned char c){ return std::tolower(c); });
        if (space_charge == "1" || space_charge_lower == "true" || space_charge_lower == "on")
        {
            ablastr::warn_manager::WMRecordWarning(
                "algo.space_charge",
                "The option algo.space_charge = true is deprecated and will be removed in a future version of ImpactX. "
                "Please use algo.space_charge = 3D instead.",
                ablastr::warn_manager::WarnPriority::high
            );
            return SpaceChargeAlgo::True_3D;
        }
        if (space_charge == "0")
        {
            ablastr::warn_manager::WMRecordWarning(
                "algo.space_charge",
                "The option algo.space_charge = 0 is deprecated and will be removed in a future version of ImpactX. "
                "Please use algo.space_charge = false instead.",
                ablastr::warn_manager::WarnPriority::high
            );
            return SpaceChargeAlgo::False;
        }

        if (space_charge_lower == "false" || space_charge_lower == "off")
        {
            return SpaceChargeAlgo::False;
        }
        else if (space_charge == "3D")
        {
            return SpaceChargeAlgo::True_3D;
        }
        else if (space_charge == "2D")
        {
            return SpaceChargeAlgo::True_2D;
        }
        else
        {
            throw std::runtime_error("algo.space_charge = " + space_charge + " is not a valid option");
        }
    }

    std::string
    to_string (SpaceChargeAlgo sca)
    {
        std::string const str = amrex::getEnumNameString(sca);

        // strip True_ prefix, which we only added because var names / enum member names
        // cannot start with a number
        std::string const true_prefix = "True_";
        if (str.find(true_prefix) == 0)
        {
            return str.substr(true_prefix.length());
        }

        return str;
    }

} // namespace impactx
