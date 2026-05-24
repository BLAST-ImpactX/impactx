/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include <AMReX_MultiFab.H>

#include <unordered_map>
#include <utility>  // std::pair


namespace impactx
{
    // note:: 3D charge deposition is defined in
    // ImpactXParticleContainer::DepositCharge

    // rho 3D
    // pair: local & unique boxes
    // whole simulation index space (level 0)
    std::unordered_map<int, std::pair<amrex::MultiFab, amrex::MultiFab>>
    flatten_charge_to_2D (
        std::unordered_map<int, amrex::MultiFab> const & rho,
        amrex::Box domain_3d
    );
} // namespace impactx
