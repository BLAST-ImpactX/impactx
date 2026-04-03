/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "particles/ImpactXParticleContainer.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>


namespace impactx::particles
{
    void ParticleBoundary (
        ImpactXParticleContainer & pc
    )
    {
        BL_PROFILE("impactx::particles::ParticleBoundary")

        using namespace amrex::literals;

        amrex::ParmParse pp_algo("algo");
        std::string particle_bc = "periodic";
        pp_algo.queryAdd("particle_bc", particle_bc);

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

                // Access bucket length and reference particle quantities.
                amrex::ParticleReal const bucket_length = pc.GetBucketLength();
                amrex::ParticleReal const ref_beta = pc.GetRefParticle().beta();
                amrex::ParticleReal const bucket_duration = (ref_beta != 0.0)? bucket_length / ref_beta : 0.0;
                amrex::ParticleReal const bucket_half_duration = bucket_duration / 2.0;

                // Access data from StructOfArrays (soa)
                auto& soa_real = pti.GetStructOfArrays().GetRealData();

                amrex::ParticleReal* const AMREX_RESTRICT part_t = soa_real[RealSoA::t].dataPtr();

                // Gather particles and apply boundary condition
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
                {
                    // Access SoA Real data
                    amrex::ParticleReal & AMREX_RESTRICT t = part_t[i];

                    // Implement phase wrapping in t
                    amrex::ParticleReal ttest = std::fmod(t+bucket_half_duration, bucket_duration);
                    t = (bucket_duration != 0.0)? std::fmod(ttest+bucket_duration, bucket_duration)-bucket_half_duration : t;

                });
            } // End loop over all particle boxes
        } // End mesh-refinement level loop
    }
} // namespace impactx::particles
