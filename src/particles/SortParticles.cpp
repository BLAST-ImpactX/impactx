/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Andrew Myers
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactXParticleContainer.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_IntVect.H>
#include <AMReX_ParmParse.H>

#include <array>
#include <stdexcept>
#include <string>


namespace impactx
{
    void
    ImpactXParticleContainer::SortParticles (int step)
    {
        BL_PROFILE("ImpactXParticleContainer::SortParticles");

        amrex::ParmParse pp_algo("algo");

        // sort every sort_intervals global (slice) steps, never if <= 0
        // sorting pays off for deposition on NVIDIA and AMD GPUs, but not on CPUs
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
        int sort_intervals = 4;
#else
        int sort_intervals = -1;
#endif
        pp_algo.queryAddWithParser("sort_intervals", sort_intervals);
        if (sort_intervals <= 0 || step % sort_intervals != 0) { return; }

        // the cell-ordered sort is the faster choice for deposition on NVIDIA and AMD GPUs
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
        bool sort_particles_for_deposition = true;
#else
        bool sort_particles_for_deposition = false;
#endif
        pp_algo.queryAdd("sort_particles_for_deposition", sort_particles_for_deposition);

        if (sort_particles_for_deposition)
        {
#if !defined(AMREX_USE_CUDA) && !defined(AMREX_USE_HIP)
            throw std::runtime_error(
                "algo.sort_particles_for_deposition: sorting particles for deposition is only "
                "implemented for the CUDA and HIP backends. Set it to false to sort by bin.");
#endif


            // 0: cell-centered, 1: nodal, 2: in between, per direction
            std::array<int, AMREX_SPACEDIM> sort_idx_type{AMREX_D_DECL(0, 0, 0)};
            pp_algo.queryAdd("sort_idx_type", sort_idx_type);
            for (int const t : sort_idx_type) {
                if (t < 0 || t > 2) {
                    throw std::runtime_error(
                        "algo.sort_idx_type: values must be 0, 1 or 2, but got "
                        + std::to_string(t));
                }
            }

            SortParticlesForDeposition(amrex::IntVect(sort_idx_type));
        }
        else
        {
            std::array<int, AMREX_SPACEDIM> sort_bin_size{AMREX_D_DECL(1, 1, 1)};
            pp_algo.queryAdd("sort_bin_size", sort_bin_size);
            for (int const b : sort_bin_size) {
                if (b < 1) {
                    throw std::runtime_error(
                        "algo.sort_bin_size: values must be at least 1, but got "
                        + std::to_string(b));
                }
            }

            SortParticlesByBin(amrex::IntVect(sort_bin_size));
        }
    }

} // namespace impactx
