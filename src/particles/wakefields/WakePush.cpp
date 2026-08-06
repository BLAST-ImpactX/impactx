/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "WakePush.H"

#include <ablastr/particles/NodalFieldGather.H>

#include <AMReX_Algorithm.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BLassert.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>


namespace impactx::particles::wakefields
{
    void WakePush (
        ImpactXParticleContainer & pc,
        amrex::Gpu::DeviceVector<amrex::Real> const & wake_pt,
        amrex::Gpu::DeviceVector<amrex::Real> const & wake_px,
        amrex::Gpu::DeviceVector<amrex::Real> const & wake_py,
        amrex::ParticleReal slice_ds,
        amrex::Real bin_size,
        amrex::Real bin_min
    )
    {
        BL_PROFILE("impactx::particles::wakefields::WakePush")

        using namespace amrex::literals;

        int const num_bins = int(wake_pt.size());
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            (wake_px.empty() || int(wake_px.size()) == num_bins) &&
            (wake_py.empty() || int(wake_py.size()) == num_bins),
            "WakePush: transverse wake arrays must be empty or match wake_pt in size."
        );

        // Loop over refinement levels
        int const nLevel = pc.finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev)
        {
            // Loop over all particle boxes
            using ParIt = ImpactXParticleContainer::iterator;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (ParIt pti(pc, lev); pti.isValid(); ++pti)
            {
                const int np = pti.numParticles();

                // Physical constants and reference quantities
                amrex::ParticleReal const mc_SI = pc.GetRefParticle().mass * (ablastr::constant::SI::c);
                amrex::ParticleReal const pz_ref_SI = pc.GetRefParticle().beta_gamma() * mc_SI;
                amrex::ParticleReal const beta_ref = pc.GetRefParticle().beta();

                // Access data from StructOfArrays (soa)
                auto& soa_real = pti.GetStructOfArrays().GetRealData();

                amrex::ParticleReal* const AMREX_RESTRICT part_t = soa_real[RealSoA::t].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_pt = soa_real[RealSoA::pt].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_px = soa_real[RealSoA::px].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_py = soa_real[RealSoA::py].dataPtr();

                // Obtain constants for force normalization
                amrex::ParticleReal const push_consts_t = 1_prt / ((ablastr::constant::SI::c) * pz_ref_SI);
                amrex::ParticleReal const push_consts_xy = push_consts_t / beta_ref;

                // Per-bin forces; the transverse arrays are optional (nullptr = no kick)
                const amrex::Real* const wake_pt_ptr = wake_pt.data();
                const amrex::Real* const wake_px_ptr = wake_px.empty() ? nullptr : wake_px.data();
                const amrex::Real* const wake_py_ptr = wake_py.empty() ? nullptr : wake_py.data();

                // Gather particles and push momentum
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
                {
                    // Access SoA Real data
                    amrex::ParticleReal const & AMREX_RESTRICT t = part_t[i];

                    // Find index position along t, clamped to the wake support
                    // (guards against floating-point rounding at the beam edges)
                    int idx = static_cast<int>((t - bin_min) / bin_size);
                    idx = amrex::min(amrex::max(idx, 0), num_bins - 1);

                    // Update longitudinal momentum with the longitudinal wake force
                    part_pt[i] -= push_consts_t * slice_ds * wake_pt_ptr[idx];

                    // Update transverse momenta, if transverse forces are provided
                    if (wake_px_ptr != nullptr) {
                        part_px[i] += push_consts_xy * slice_ds * wake_px_ptr[idx];
                    }
                    if (wake_py_ptr != nullptr) {
                        part_py[i] += push_consts_xy * slice_ds * wake_py_ptr[idx];
                    }
                });
            } // End loop over all particle boxes
        } // End mesh-refinement level loop
    }
} // namespace impactx::particles::wakefields
