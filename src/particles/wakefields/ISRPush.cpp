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
#include <AMReX_Gpu.H>
#include <AMReX_Random.H>
#include <AMReX_Print.H>

namespace impactx::particles::wakefields
{

    void ISRPush (
        ImpactXParticleContainer & pc,
        amrex::ParticleReal slice_ds,
        amrex::ParticleReal rc,
        [[maybe_unused]] int isr_order
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
        amrex::ParticleReal const lambda_e = hbar/mc_SI;

        // Obtain constants for force normalization
        amrex::ParticleReal const B_normal = bg_ref/std::abs(rc);
        amrex::ParticleReal const c1 = 2.0_prt/3.0_prt * r_e * slice_ds * std::pow(B_normal,2);
        amrex::ParticleReal const c2 = B_normal * lambda_e;

        // Coefficients of the Taylor expansion of polynomials g and h
        amrex::ParticleReal const g0 = 1_prt;
        amrex::ParticleReal const g1 = -55_prt*std::sqrt(3_prt)/16_prt;
        amrex::ParticleReal const g2 = 48_prt;
        // amrex::ParticleReal const g3 = -8855_prt*std::sqrt(3_prt)/32_prt;
        amrex::ParticleReal const h1 = 55_prt/(16_prt*std::sqrt(3_prt));
        amrex::ParticleReal const h2 = -28_prt;
        amrex::ParticleReal const h3 = 14245_prt*std::sqrt(3_prt)/64_prt;

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
                amrex::ParallelForRNG(np, [=] AMREX_GPU_DEVICE (int i, amrex::RandomEngine const & engine)
                {
                    // Access SoA Real data
                    amrex::ParticleReal & AMREX_RESTRICT px = part_px[i];
                    amrex::ParticleReal & AMREX_RESTRICT py = part_py[i];
                    amrex::ParticleReal & AMREX_RESTRICT pt = part_pt[i];

                    // Standard normal random variable used to kick this particle
                    amrex::ParticleReal const xi = amrex::RandomNormal(0.0,1.0,engine);

                    // Relativistic beta*gamma for this particle
                    amrex::ParticleReal const gamma = gamma_ref - bg_ref*pt;
                    amrex::ParticleReal const bg = std::sqrt(std::pow(gamma,2)-1_prt);

                    // Value of ISR kick in the total momentum (normalized by mc)
                    amrex::ParticleReal const tau = c1 * bg;
                    amrex::ParticleReal const chi = c2 * bg;

                    // Check the order of quantum corrections
                    amrex::ParticleReal g = 0_prt;
                    amrex::ParticleReal h = 0_prt;
                    if (isr_order == 1) {
                       g = g0;
                       h = h1*chi;
                    } else if (isr_order == 2) {
                       g = g0 + g1*chi;
                       h = h1*chi + h2*std::pow(chi,2);
                    } else if (isr_order == 3) {
                       g = g0 + g1*chi + g2*std::pow(chi,2);
                       h = h1*chi + h2*std::pow(chi,2) + h3*std::pow(chi,3);
                    }

                    // Value of the ISR kick in total momentum (relative to total momentum):
                    amrex::ParticleReal const dp = (-tau*g + std::sqrt(tau*h)*xi);

                    // Final value of updated particle gamma:
                    amrex::ParticleReal const bg_f = bg*(1_prt + dp);
                    amrex::ParticleReal const gamma_f = std::sqrt(1_prt + std::pow(bg_f,2));

                    // Update momentum
                    px = px * (1_prt + dp);
                    py = py * (1_prt + dp);
                    pt = (gamma_ref - gamma_f)/(bg_ref);

                });

            } // End loop over all particle boxes
        } // End mesh-refinement level loop

        amrex::Gpu::streamSynchronize();

   }
} // namespace impactx::particles::wakefields
