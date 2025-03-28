/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "ISRPush.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>
#include <AMReX_Random.H>

namespace impactx::particles::wakefields
{

    void ISRPush (
        ImpactXParticleContainer & pc,
        amrex::ParticleReal slice_ds,
        amrex::ParticleReal rc
    )
    {
        BL_PROFILE("impactx::particles::wakefields::ISRPush")

        using namespace amrex::literals;

        // Physical constants and reference quantities
        amrex::ParticleReal const mc_SI = pc.GetRefParticle().mass * (ablastr::constant::SI::c);
        amrex::ParticleReal const gamma_ref = pc.GetRefParticle().gamma();
        amrex::ParticleReal const bg_ref = pc.GetRefParticle().beta_gamma();
        amrex::ParticleReal const r_e = (ablastr::constant::SI::r_e);
        amrex::ParticleReal const hbar = (ablastr::constant::SI::hbar);

        // Obtain constants for force normalization
        amrex::ParticleReal const c1 = 2.0_prt/3.0_prt * (ablastr::constant::SI::r_e) * std::pow(bg_ref,2);
        amrex::ParticleReal const c2 = 55_prt/(24_prt*std::sqrt(3_prt)) * r_e * hbar * std::pow(bg_ref,3)/mc_SI;
        amrex::ParticleReal const rc_sqrt = std::sqrt(rc);

        amrex::ParticleReal const deterministic_coef = c1 * slice_ds / std::pow(rc,2);
        amrex::ParticleReal const stochastic_coef = std::sqrt(c2 * slice_ds) / std::pow(rc_sqrt,3);

        // Loop over refinement levels
        int const nLevel = pc.finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev)
        {
            // Loop over all particle boxes
            using ParIt = ImpactXParticleContainer::iterator;

            for (ParIt pti(pc, lev); pti.isValid(); ++pti)
            {
                const int np = pti.numParticles();

                // Access data from StructOfArrays (soa)
                auto& soa_real = pti.GetStructOfArrays().GetRealData();

                amrex::ParticleReal* const AMREX_RESTRICT part_px = soa_real[RealSoA::px].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_py = soa_real[RealSoA::py].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_pt = soa_real[RealSoA::pt].dataPtr();

                // Gather particles and push momentum
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
                {
                    // Access SoA Real data
                    amrex::ParticleReal & AMREX_RESTRICT px = part_px[i];
                    amrex::ParticleReal & AMREX_RESTRICT py = part_py[i];
                    amrex::ParticleReal & AMREX_RESTRICT pt = part_pt[i];

                    // Standard normal random variable used to kick this particle
                    amrex::ParticleReal const xi = amrex::RandomNormal(0.0,1.0);

                    // Relativistic beta*gamma for this particle
                    amrex::ParticleReal const gamma = gamma_ref - bg_ref*pt;
                    amrex::ParticleReal const bg = std::sqrt(std::pow(gamma,2)-1_prt);

                    // Value of ISR kick in the total momentum (normalized by mc)
                    amrex::ParticleReal const dp = std::pow(bg,2) * (deterministic_coef + stochastic_coef*xi);
                    amrex::ParticleReal const bg_f = bg - dp;

                    // Value of ISR kick in the total energy (normalized by mc^2)
                    amrex::ParticleReal const gamma_f = std::sqrt(1_prt + std::pow(bg_f,2));

                    // Update momentum
                    px = px * (1_prt - dp/bg);
                    py = py * (1_prt - dp/bg);
                    pt = (gamma_ref - gamma_f)/(bg_ref);

                });
            } // End loop over all particle boxes
        } // End mesh-refinement level loop
    }
} // namespace impactx::particles::wakefields
