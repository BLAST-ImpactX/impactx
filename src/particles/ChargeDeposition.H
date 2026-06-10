/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_CHARGE_DEPOSITION_H
#define IMPACTX_CHARGE_DEPOSITION_H

#include <AMReX_BaseFab.H>
#include <AMReX_Box.H>
#include <AMReX_FabArray.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_iMultiFab.H>       // for amrex::OwnerMask
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>    // for amrex::ReduceToPlaneMF2
#include <AMReX_Periodicity.H>
#include <AMReX_Reduce.H>          // for amrex::ReduceOpSum
#include <AMReX_REAL.H>

#include <unordered_map>
#include <utility>  // std::pair


namespace impactx
{
    // note:: 3D charge deposition is defined in
    // ImpactXParticleContainer::DepositCharge

    /** Flatten a per-level 3D charge density to 2D (sum along z).
     *
     * @tparam T_FieldMF the mesh field array type (e.g. amrex::MultiFab for
     *                   double precision, amrex::fMultiFab for single precision)
     * @param rho the 3D charge density per level
     * @param domain_3d whole simulation index space (level 0)
     * @return per level a pair of (local boxes, unique boxes) 2D charge density
     */
    template <class T_FieldMF = amrex::MultiFab>
    std::unordered_map<int, std::pair<T_FieldMF, T_FieldMF>>
    flatten_charge_to_2D (
        std::unordered_map<int, T_FieldMF> const & rho,
        amrex::Box domain_3d
    )
    {
        std::unordered_map<int, std::pair<T_FieldMF, T_FieldMF>> rho_2d;

        int const finest_level = rho.size() - 1;
        for (int lev = 0; lev <= finest_level; ++lev) {
            auto omask = amrex::OwnerMask(rho.at(lev), amrex::Periodicity::NonPeriodic());
            auto const& oma = omask->const_arrays();
            // flattened rho
            auto const& ma = rho.at(lev).const_arrays();
            amrex::Box domain_lev = lev == 0 ? domain_3d : rho.at(lev).boxArray().minimalBox();
            rho_2d.erase(lev);
            rho_2d.emplace(
                lev,
                amrex::ReduceToPlaneMF2<amrex::ReduceOpSum>
                (2, domain_lev, rho.at(lev), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
                    -> typename T_FieldMF::value_type
                {
                    return oma[b](i,j,k) ? ma[b](i,j,k) : typename T_FieldMF::value_type(0);
                })
            );
        }

        return rho_2d;
    }
} // namespace impactx

#endif // IMPACTX_CHARGE_DEPOSITION_H
