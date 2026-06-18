/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactXParticleContainer.H"
#include "ChargeDeposition.H"
#include "initialization/Algorithms.H"

#include <ablastr/coarsen/average.H>
#include <ablastr/utils/Communication.H>
#include <ablastr/particles/DepositCharge.H>

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParticleTile.H>

#include <array>


namespace impactx
{
    std::unordered_map<int, std::pair<amrex::MultiFab, amrex::MultiFab>>
    flatten_charge_to_2D (
        std::unordered_map<int, amrex::MultiFab> const & rho,
        amrex::Box domain_3d
    )
    {
        std::unordered_map<int, std::pair<amrex::MultiFab, amrex::MultiFab>> rho_2d;

        int const finest_level = rho.size() - 1;
        for (int lev = 0; lev <= finest_level; ++lev) {
            auto const & rho_at_level = rho.at(lev);
            amrex::Box const domain_lev = amrex::convert(
                lev == 0 ? domain_3d : rho_at_level.boxArray().minimalBox(),
                rho_at_level.ixType()
            );

            int const z_dir = 2;

            // Project the deposited 3D charge density to x-y for the 2D
            // space-charge solve. We must include *all* guard cells to conserve
            // charge:
            //  - z guards: with a single longitudinal cell and higher-order
            //    particle shapes, part of each particle's charge is deposited
            //    outside the valid z nodes: summed into the plane.
            //  - transverse (x/y) guards at internal box seams: that charge
            //    belongs to a neighbouring box's valid node and is folded back
            //    in across the seam.
            //
            // The open transverse boundary is non-periodic, so charge deposited
            // outside the domain is dropped (open BC).
            //
            // NOTE: this requires the *raw* deposited rho, i.e. DepositCharge
            // must NOT SumBoundary it for the 2D path, otherwise the valid nodes
            // are folded twice and the charge is double counted.
            auto const & ma = rho_at_level.const_arrays();

            rho_2d.erase(lev);
            rho_2d.emplace(
                lev,
                amrex::ReduceToPlaneMF2<amrex::ReduceOpSum>(
                    z_dir, domain_lev, rho_at_level,
                    [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
                    {
                        return ma[b](i,j,k);
                    },
                    rho_at_level.nGrowVect(),            // include z + transverse guards
                    amrex::Periodicity::NonPeriodic()    // open transverse boundary
                )
            );
        }

        return rho_2d;
    }

    std::unordered_map<int, amrex::MultiFab>
    project_charge_to_2D (
        std::unordered_map<int, amrex::MultiFab> const & rho,
        amrex::Box domain_3d
    )
    {
        std::unordered_map<int, amrex::MultiFab> rho_2d;
        auto rho_2d_pairs = flatten_charge_to_2D(rho, domain_3d);
        for (auto & [lev, rho_2d_pair] : rho_2d_pairs) {
            rho_2d.emplace(lev, std::move(rho_2d_pair.second));
        }
        return rho_2d;
    }

    void
    ImpactXParticleContainer::DepositCharge (
        std::unordered_map<int, amrex::MultiFab> & rho,
        amrex::Vector<amrex::IntVect> const & ref_ratio)
    {
        BL_PROFILE("ImpactXParticleContainer::DepositCharge");

        using namespace amrex::literals; // for _rt and _prt

        if (m_particle_shape.value() < 1)
            throw std::runtime_error("DepositCharge: Particle shape must be >=1");

        // The 2D space-charge model projects rho to x-y in flatten_charge_to_2D,
        // which folds the (transverse and longitudinal) guard-cell charge into
        // the valid nodes itself. We must therefore leave rho *un-summed* here:
        // a SumBoundary would fold the transverse guards a first time and the
        // 2D projection of ReduceToPlaneMF2 a second time, double counting the
        // charge at box seams.
        auto const space_charge = get_space_charge_algo();
        bool const skip_sum_boundary =
            (space_charge == SpaceChargeAlgo::True_2D ||
             space_charge == SpaceChargeAlgo::True_2p5D);

        // reset the values in rho to zero
        int const nLevel = this->finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev) {
            rho.at(lev).setVal(0.);
        }

        // loop fine-to-coarse over refinement levels
        for (int lev = nLevel; lev >= 0; --lev) {
            amrex::MultiFab & rho_at_level = rho.at(lev);

            // get simulation geometry information
            amrex::Geometry const & gm = this->Geom(lev);

            // Loop over particle tiles and deposit charge on each level
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            {
                amrex::FArrayBox local_rho_fab;

                using ParIt = ImpactXParticleContainer::iterator;
                for (ParIt pti(*this, lev); pti.isValid(); ++pti) {
                    // preparing access to particle data: SoA of Reals
                    auto & AMREX_RESTRICT soa_real = pti.GetStructOfArrays().GetRealData();
                    // after https://github.com/BLAST-WarpX/warpx/pull/2838 add const:
                    auto const wp = soa_real[RealSoA::w];
                    int const * const AMREX_RESTRICT ion_lev = nullptr;

                    // physical lower corner of the current box
                    //   Note that this includes guard cells since it is after tilebox.grow
                    amrex::Box tilebox = pti.tilebox();
                    tilebox.grow(rho_at_level.nGrowVect());
                    amrex::RealBox const grid_box{tilebox, gm.CellSize(), gm.ProbLo()};
                    amrex::Real const * const AMREX_RESTRICT xyzmin_ptr = grid_box.lo();

                    // mesh-refinement: for when we do not deposit on the same level
                    // note: would need to communicate the deposited-to boxes afterwards
                    //int const depos_lev = lev;
                    // mesh refinement ratio between lev and depos_lev
                    //auto const rel_ref_ratio = ref_ratio.at(depos_lev) / ref_ratio.at(lev);
                    amrex::ignore_unused(ref_ratio);

                    // in SI [C]
                    amrex::ParticleReal const charge = m_refpart.charge;

                    // lower end and inverse cell size of the mesh to deposit to
                    amrex::XDim3 const xyzmin = {xyzmin_ptr[0], xyzmin_ptr[1], xyzmin_ptr[2]};
                    amrex::XDim3 const dinv = {1.0_rt/gm.CellSize(0), 1.0_rt/gm.CellSize(1), 1.0_rt/gm.CellSize(2)};

                    // RZ modes (unused)
                    int const n_rz_azimuthal_modes = 0;

                    ablastr::particles::deposit_charge<ImpactXParticleContainer>
                            (pti, wp, charge, ion_lev, &rho_at_level,
                             local_rho_fab,
                             m_particle_shape.value(),
                             dinv, xyzmin, n_rz_azimuthal_modes);
                }
            }

            // TODO: Call portion's of WarpX' SyncRho from fine to coarser levels
            //   TODO: do coarsening to a temp, local lev-1 patch
            //   TODO: start communicating this patch into rho_at_level_minus_1
            //   note: this can either move parts from WarpXComm.cpp (SyncRho & AddRhoFromFineLevelandSumBoundary)
            //         or use code from SyncPhi in ABLASTR's PoissonSolve (in opposite order: FP->CP here)
            // needed for solving the levels by levels:
            // - coarser level is initial guess for finer level
            // - coarser level provides boundary values for finer level patch
            // Interpolation from phi[lev] to phi[lev+1]
            // (This provides both the boundary conditions and initial guess for phi[lev+1])

            if (lev > 0) {
                // Allocate rho_cp with the same distribution map as lev
                amrex::BoxArray ba = rho[lev].boxArray();
                const amrex::IntVect &refratio = ref_ratio[lev - 1];
                ba.coarsen(refratio);  // index space is now coarsened

                // Number of guard cells to fill on coarse patch and number of components
                const amrex::IntVect ngrow = (rho[lev].nGrowVect() + refratio - 1) / refratio;  // round up int division
                const int ncomp = 1;  // rho is a scalar
                amrex::MultiFab rho_cp(ba, rho[lev].DistributionMap(), ncomp, ngrow);

                // coarsen the data
                ablastr::coarsen::average::Coarsen(
                    rho_cp,
                    rho[lev],
                    refratio
                );

                // add to the lower level
                const amrex::Periodicity& crse_period = this->Geom(lev - 1).periodicity();

                // On a coarse level, the data in mf_comm comes from the
                // coarse patch of the fine level. They are unfiltered and uncommunicated.
                // We need to add it to the fine patch of the current level.
                amrex::MultiFab fine_lev_cp(
                    rho[lev-1].boxArray(),
                    rho[lev-1].DistributionMap(),
                    1,
                    0);
                fine_lev_cp.setVal(0.0);
                fine_lev_cp.ParallelAdd(
                    rho_cp,
                    0,
                    0,
                    1,
                    rho_cp.nGrowVect(),
                    amrex::IntVect(0),
                    crse_period
                );
                // We now need to create a mask to fix the double counting.
                auto owner_mask = amrex::OwnerMask(fine_lev_cp, crse_period);
                auto const& mma = owner_mask->const_arrays();
                auto const& sma = fine_lev_cp.const_arrays();
                auto const& dma = rho[lev-1].arrays();
                amrex::ParallelFor(
                    fine_lev_cp,
                    amrex::IntVect(0),
                    ncomp,
                    [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n)
                    {
                        if (mma[bno](i,j,k) && sma[bno](i,j,k,n) != 0.0_rt) {
                            dma[bno](i,j,k,n) += sma[bno](i,j,k,n);
                        }
                    }
                );
            } // if (lev > 0)

            // charge filters
            // note: we do this after SyncRho, because the physical (dx) size of
            //       the stencil is different for each level
            // note: we might not be able to do this after SumBoundary because we would then
            //       not have valid filtered values in the guard. We access the guard
            //       when we sum contributions from particles close to the
            //       MR border between levels.
            //ApplyFilterandSumBoundaryRho(lev, lev, rho[lev-1], 0, 1);

            // start async charge communication for this level
            //   (skipped for the 2D model: the 2D projection sums the guards)
            if (!skip_sum_boundary) {
                rho_at_level.SumBoundary_nowait(gm.periodicity());
            }
            //int const comp = 0;
            //rho_at_level.SumBoundary_nowait(comp, comp, rho_at_level.nGrowVect(), gm.periodicity());

        } // lev

        // finalize communication
        if (!skip_sum_boundary) {
            for (int lev = 0; lev <= nLevel; ++lev)
            {
                amrex::MultiFab & rho_at_level = rho.at(lev);
                rho_at_level.SumBoundary_finish();
            }
        }
    }
} // namespace impactx
