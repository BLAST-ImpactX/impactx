/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "particles/ParticleBoundary.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <stdexcept>
#include <string>


namespace impactx::particles {
    ParticleBC
    get_particle_boundary_condition ()
    {
        auto particle_bc = ParticleBC::open;
        amrex::ParmParse("algo").queryAdd("particle_bc", particle_bc);

        return particle_bc;
    }

    void ParticleBoundary (
        ImpactXParticleContainer & pc
    )
    {
        BL_PROFILE("impactx::particles::ParticleBoundary")

        // check option and set default (open) if unset
        auto particle_bc = get_particle_boundary_condition();

        // nothing to do by default
        if (particle_bc == ParticleBC::open) {
            return;
        }

        // Access bucket length and reference particle quantities
        using namespace amrex::literals;
        amrex::ParticleReal const bucket_length = pc.GetBucketLength();
        if (bucket_length <= 0.0_prt) {
            throw std::runtime_error("ParticleBoundary: Bucket length must be set >0 for particle boundary conditions.");
        }
        amrex::ParticleReal const ref_beta = pc.GetRefParticle().beta();
        amrex::ParticleReal const bucket_duration = (ref_beta != 0.0_prt) ? bucket_length / ref_beta : 0.0_prt;
        amrex::ParticleReal const bucket_half_duration = bucket_duration * 0.5_prt;

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

                // Access data from StructOfArrays (soa)
                auto & soa_real = pti.GetStructOfArrays().GetRealData();

                amrex::ParticleReal * const AMREX_RESTRICT part_t = soa_real[RealSoA::t].dataPtr();
                amrex::ParticleReal * const AMREX_RESTRICT part_pt = soa_real[RealSoA::pt].dataPtr();
                uint64_t * const AMREX_RESTRICT part_idcpu = pti.GetStructOfArrays().GetIdCPUData().dataPtr();

                switch (particle_bc) {
                    case ParticleBC::periodic:

                        // Gather particles and apply boundary condition
                        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
                        {
                            // Access SoA Real data
                            amrex::ParticleReal & AMREX_RESTRICT t = part_t[i];

                            // Periodic particle boundary condition:
                            //   apply phase wrapping in t (modulo bucket_duration)
                            amrex::ParticleReal const ttest = std::fmod(t + bucket_half_duration, bucket_duration);
                            t = (bucket_duration != 0.0_prt) ?
                                std::fmod(ttest + bucket_duration, bucket_duration) - bucket_half_duration
                                : t;
                        });
                        break;

                    case ParticleBC::absorbing:

                        // Gather particles and apply boundary condition
                        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
                        {
                            // Access SoA Real data
                            amrex::ParticleReal const & AMREX_RESTRICT t = part_t[i];

                            // Absorbing particle boundary condition:
                            //   check particle against the boundary
                            bool const inside_boundary = (std::abs(t) < bucket_half_duration);
                            // Mark particles as lost if appropriate
                            amrex::ParticleIDWrapper{part_idcpu[i]}.make_invalid(!inside_boundary);
                        });
                        break;

                    case ParticleBC::reflecting:

                        // Gather particles and apply boundary condition
                        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
                        {
                            // Access SoA Real data
                            amrex::ParticleReal & AMREX_RESTRICT t = part_t[i];
                            amrex::ParticleReal & AMREX_RESTRICT pt = part_pt[i];

                            // Reflecting particle boundary condition.
                            // TODO:  Transform (t,pt) to (z,pz) using z-to-t transformation.
                            // The implementation below works through linear order in the phase space variables.
                            // If particle falls outside the boundary, reflect:
                            if (t > bucket_half_duration) {
                                t = bucket_duration - t;
                                pt = -pt;
                            } else if (t < -bucket_half_duration) {
                                t = -bucket_duration - t;
                                pt = -pt;
                            }
                        });
                        break;

                    default: {
                        std::string const bc_str = amrex::getEnumNameString(particle_bc);
                        std::string msg = "ParticleBoundary: Unknown particle_bc: ";
                        msg += bc_str +"\nMust be one of: ";
                        for (auto const& name : amrex::getEnumNameStrings<ParticleBC>()) {
                            msg += name + ", ";
                        }
                        msg.erase(msg.size() - 2);
                        throw std::runtime_error(msg);
                    }
                 } // End switch (particle_bc)
            } // End loop over all particle boxes
        } // End mesh-refinement level loop
    }
} // namespace impactx::particles
