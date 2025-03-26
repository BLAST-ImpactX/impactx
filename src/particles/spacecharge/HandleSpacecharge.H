/* Copyright 2022-2025 The Regents of the University of California, through Lawrence
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

#include "initialization/AmrCoreData.H"

#include <AMReX_REAL.H>

#include <functional>
#include <memory>


namespace impactx::particles::spacecharge
{
    /** Function to handle space charge including charge deposition, Poisson solve,
     * field calculation, interpolation, and particle push.
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

} // namespace impactx::particles::spacecharge

#endif // IMPACTX_HANDLE_SPACECHARGE_H
