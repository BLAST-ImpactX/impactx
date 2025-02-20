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

#include <AMReX_REAL.H>       // for Real
#include <AMReX_SmallMatrix.H>
#include <AMReX_Print.H>

#include <cmath>


namespace impactx::envelope::spacecharge
{
    void
    space_charge2D_push (
        RefPart const & AMREX_RESTRICT refpart,
        Map6x6 & AMREX_RESTRICT cm,
        amrex::ParticleReal current,
        amrex::ParticleReal ds
    )
    {
        using namespace amrex::literals;

        // initialize the linear transport map
        Map6x6 R = Map6x6::Identity();

        // physical constants and reference quantities
        using ablastr::constant::SI::c;
        using ablastr::constant::SI::ep0;
        using ablastr::constant::math::pi;
        amrex::ParticleReal const mass = refpart.mass;
        amrex::ParticleReal const charge = refpart.charge;
        amrex::ParticleReal const pt_ref = refpart.pt;
        amrex::ParticleReal const betgam2 = std::pow(pt_ref, 2) - 1_prt;
        amrex::ParticleReal const betgam = std::sqrt(betgam2);
        amrex::ParticleReal const betgam3 = std::pow(betgam,3);

        // evaluate the beam space charge perveance from current
        amrex::ParticleReal const IA = 4_prt * pi * ep0 * mass * std::pow(c,3) / charge;
        amrex::ParticleReal const Kpv = std::abs(current / IA) * 2_prt / betgam3;

        // evaluate the linear transfer map
        amrex::ParticleReal const sigma2 = cm(1,1) * cm(3,3) - cm(1,3) * cm(1,3);
        amrex::ParticleReal const sigma = std::sqrt(sigma2);
        amrex::ParticleReal const D = (sigma + cm(1,1)) * (sigma + cm(3,3)) - cm(1,3) * cm(1,3);
        amrex::ParticleReal const coeff = Kpv * ds / (2_prt * D);

        R(2,1) = coeff * (sigma + cm(3,3));
        R(2,3) = coeff * (-cm(1,3));
        R(4,1) = R(2,3);
        R(4,3) = coeff * (sigma + cm(1,1));

        // update the beam covariance matrix
        cm = R * cm * R.transpose();

    }


} // namespace impactx::spacecharge
