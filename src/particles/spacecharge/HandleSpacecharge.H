/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_HANDLE_SPACECHARGE_H
#define IMPACTX_HANDLE_SPACECHARGE_H

#include "initialization/Algorithms.H"
#include "initialization/AmrCoreData.H"

#include <AMReX_REAL.H>

#include <functional>
#include <memory>


namespace impactx::particles::spacecharge
{
    /** Function to handle space charge including charge deposition, Poisson solve,
     * field calculation, interpolation, and particle push.
     *
     * Returns immediately if slice_ds is zero.
     *
     * @param[in,out] amr_data AMR data like particle container and fields
     * @param[in] ResizeMesh function to call to resize the mesh
     * @param[in] slice_ds slice spacing along s
     */
    void HandleSpacecharge (
        std::unique_ptr<initialization::AmrCoreData> & amr_data,
        std::function<void()> ResizeMesh,
        amrex::Real slice_ds
    );

    /** Throw if the beam has no extent in a direction that the space charge model needs.
     *
     * The particle-in-cell models size their mesh from the beam extent, and the 2.5D and
     * Gaussian models divide by it.  A direction counts as flat when its extent is zero up
     * to the rounding of the particle coordinates.  Only the 2D model of an unbunched beam,
     * which does not resolve the longitudinal direction, accepts a beam without
     * longitudinal extent, e.g. a DC beam with all particles at the same time of flight.
     *
     * @param[in] x_min,y_min,z_min,x_max,y_max,z_max extent of the particle positions
     * @param[in] space_charge space charge model
     * @return true if the beam has no longitudinal extent (only with the 2D model)
     */
    bool CheckBeamExtent (
        amrex::ParticleReal x_min, amrex::ParticleReal y_min, amrex::ParticleReal z_min,
        amrex::ParticleReal x_max, amrex::ParticleReal y_max, amrex::ParticleReal z_max,
        SpaceChargeAlgo space_charge
    );

} // namespace impactx::particles::spacecharge

#endif // IMPACTX_HANDLE_SPACECHARGE_H
