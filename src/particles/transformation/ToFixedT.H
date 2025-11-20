/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_TO_FIXED_T_H
#define IMPACTX_TO_FIXED_T_H

#include "particles/transformation/SafeInverse.H"
#include "particles/ImpactXParticleContainer.H"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx::transformation
{
    struct ToFixedT
    {
        using PType = ImpactXParticleContainer::ParticleType;

        /** Transformation of particles from fixed location s to fixed time t.
         *
         * At fixed t, each particle is represented by phase space
         *  coordinates (x,y,z,px,py,pz).
         *
         * @param ptd Design value of pt/mc2 = -gamma.
         */
        ToFixedT (amrex::ParticleReal const ptd)
        : m_ptd(ptd)
        {
        }

        /** This is a s-to-t map, so that a variable of this type can be used like a
         *  s-to-t function.
         *
         * @param[inout] x particle position in x
         * @param[inout] y particle position in y
         * @param[inout] t particle position in t (in), in z (out)
         * @param[inout] px particle momentum in x
         * @param[inout] py particle momentum in y
         * @param[inout] pt particle momentum in t (in), in z (out)
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            amrex::ParticleReal & x,
            amrex::ParticleReal & y,
            amrex::ParticleReal & t,
            amrex::ParticleReal & px,
            amrex::ParticleReal & py,
            amrex::ParticleReal & pt) const
        {
            using namespace amrex::literals;
            using amrex::Math::powi;

            // small tolerance to avoid NaN for pz<0 and improve single precision stability
            constexpr amrex::ParticleReal tol =
#ifdef AMREX_SINGLE_PRECISION_PARTICLES
                1.0e-6_prt;  // single precision
#else
                1.0e-12_prt;  // double precision
#endif

            // compute value of reference pzd = beta*gamma
            amrex::ParticleReal const argd = -1.0_prt + powi<2>(m_ptd);
            // AMREX_ASSERT_WITH_MESSAGE(argd > 0.0_prt, "invalid pzd arg (<=0)");
            amrex::ParticleReal const pzdf = argd > 0.0_prt ? std::sqrt(argd) : tol;

            // transform momenta to dynamic units (eg, so that momenta are
            // normalized by mc):
            px = px * pzdf;
            py = py * pzdf;
            pt = pt * pzdf;

            // compute value of particle pz = beta*gamma
            amrex::ParticleReal const pt_total = m_ptd + pt;
            amrex::ParticleReal const arg = -1.0_prt + powi<2>(pt_total) - powi<2>(px) - powi<2>(py);
            // AMREX_ASSERT_WITH_MESSAGE(arg > 0.0_prt, "invalid pz arg (<=0)");
            amrex::ParticleReal const pzf = arg > 0.0_prt ? std::sqrt(arg) : tol;

            // transform position and momentum (from fixed s to fixed t)
            //   numerical stability: avoid division by very small numbers
            amrex::ParticleReal const inv_pt_total = safe_inverse(pt_total, tol);
            // Compute products once to reduce rounding errors
            amrex::ParticleReal const px_t = px * t;
            amrex::ParticleReal const py_t = py * t;
            amrex::ParticleReal const pzf_t = pzf * t;
            x = x + px_t * inv_pt_total;
            y = y + py_t * inv_pt_total;
            auto & z = t;    // We store z in the same memory slot as t.
            z = pzf_t * inv_pt_total;
            // px = px;
            // py = py;
            auto & pz = pt;  // We store pz in the same memory slot as pt.
            pz = pzf - pzdf;

            // transform momenta to static units (eg, so that momenta are
            // normalized by pzdf):
            amrex::ParticleReal const inv_pzdf = 1.0_prt / pzdf;
            px = px * inv_pzdf;
            py = py * inv_pzdf;
            pz = pz * inv_pzdf;
        }

    private:
        amrex::ParticleReal m_ptd;  ///< Design value of pt/mc2 = -gamma.
    };

} // namespace impactx::transformation

#endif // IMPACTX_TO_FIXED_T_H
