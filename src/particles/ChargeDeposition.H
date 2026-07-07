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

    /** Project a deposited 3D charge density to 2D (x-y).
     *
     * Convenience wrapper around flatten_charge_to_2D that keeps only the unique
     * (charge-conserving) projection per level.
     *
     * @param rho the deposited (raw, un-summed) 3D charge per level
     * @param domain_3d the level-0 simulation domain (cell-centered)
     * @return the 2D-projected charge per level
     */
    std::unordered_map<int, amrex::MultiFab>
    project_charge_to_2D (
        std::unordered_map<int, amrex::MultiFab> const & rho,
        amrex::Box domain_3d
    );
} // namespace impactx
