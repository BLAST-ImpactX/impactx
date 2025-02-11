/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Marco Garten, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "EnvelopeSpaceChargePush.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>       // for Real
#include <AMReX_SmallMatrix.H>

#include <cmath>

namespace impactx::spacecharge
{
        /** This function returns the linear transport map associated with a
         *  reduced 2D space charge model, based on the beam covariance matrix.
         *
         * @param[in] refpart reference particle
         * @param[in,out] cm covariance matrix
         * @param[in] current beam current [A]
         * @param[in] ds step size [m]
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        void
        envelope_space_charge2D_push (
        RefPart const & refpart,
        Map6x6 & cm,
        amrex::ParticleReal & current,
        amrex::ParticleReal & ds
    )
    {
        using namespace amrex::literals;

        // initialize the linear transport map
        Map6x6 R = Map6x6::Identity();

        // access reference particle values to find beta*gamma^2
        amrex::ParticleReal const pt_ref = refpart.pt;
        amrex::ParticleReal const betgam2 = std::pow(pt_ref, 2) - 1.0_prt;

        // evaluate the beam space charge perveance from current
        // TODO
        amrex::ParticleReal const Kpv=0.0_prt;

        // evaluate the linear transfer map
        amrex::ParticleReal const sigma2 = cm(1,1)*cm(3,3)-cm(1,3)*cm(1,3);
        amrex::ParticleReal const sigma = std::sqrt(sigma2);
        amrex::ParticleReal const D = (sigma + cm(1,1)) * (sigma + cm(3,3)) - cm(1,3)*cm(1,3);
        amrex::ParticleReal const coeff = Kpv*ds/(2.0_prt*D);

        R(2,1) = coeff * (sigma + cm(3,3));
        R(2,3) = coeff * (-cm(1,3));
        R(4,1) = R(2,3);
        R(4,3) = coeff * (sigma + cm(1,1));

        // update the beam covariance matrix
        cm = R * cm * R.transpose();

    }


} // namespace impactx::spacecharge
