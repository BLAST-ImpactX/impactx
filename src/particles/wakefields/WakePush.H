/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef WAKEPUSH_H
#define WAKEPUSH_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>


namespace impactx::particles::wakefields
{
    /** Function to update particle momentum using wake force kicks
     *
     * Gathers the per-bin forces along t (nearest grid point) and applies
     * them as momentum kicks over the slice length.
     *
     * @param[inout] pc the particle species container
     * @param[in] wake_pt per-bin longitudinal force [N];
     *            applied as pt -= slice_ds * F_t / (c * p_ref)
     * @param[in] wake_px per-bin horizontal force [N] or empty (= no kick);
     *            applied as px += slice_ds * F_x / (beta * c * p_ref)
     * @param[in] wake_py per-bin vertical force [N] or empty (= no kick);
     *            applied as py += slice_ds * F_y / (beta * c * p_ref)
     * @param[in] slice_ds slice spacing along s
     * @param[in] bin_size size of the beam in s divided by num_bins
     * @param[in] bin_min lower end of the beam in s
     */
    void WakePush (
        ImpactXParticleContainer & pc,
        amrex::Gpu::DeviceVector<amrex::Real> const & wake_pt,
        amrex::Gpu::DeviceVector<amrex::Real> const & wake_px,
        amrex::Gpu::DeviceVector<amrex::Real> const & wake_py,
        amrex::ParticleReal slice_ds,
        amrex::Real bin_size,
        amrex::Real bin_min
    );

} // namespace impactx::particles::wakefields

#endif // WAKEPUSH_H
