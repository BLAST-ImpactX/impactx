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
#include "Elliptic_Integral.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX_Math.H>
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
        using amrex::Math::powi;

        // skip calculations for trivial case
        if (current == 0_prt) { return; }

        // initialize the linear transport map
        Map6x6 R = Map6x6::Identity();

        // physical constants and reference quantities
        using ablastr::constant::SI::c;
        using ablastr::constant::SI::ep0;
        using ablastr::constant::math::pi;
        amrex::ParticleReal const mass = refpart.mass;
        amrex::ParticleReal const charge = refpart.charge;
        amrex::ParticleReal const pt_ref = refpart.pt;
        amrex::ParticleReal const betgam2 = powi<2>(pt_ref) - 1_prt;
        amrex::ParticleReal const betgam = std::sqrt(betgam2);
        amrex::ParticleReal const betgam3 = powi<3>(betgam);

        // evaluate the beam space charge perveance from current
        amrex::ParticleReal const IA = 4_prt * pi * ep0 * mass * powi<3>(c) / charge;
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

    void
    space_charge3D_push (
        RefPart const & AMREX_RESTRICT refpart,
        Map6x6 & AMREX_RESTRICT cm,
        amrex::ParticleReal bunch_charge,
        amrex::ParticleReal ds
    )
    {
        using namespace amrex::literals;
        using amrex::Math::powi;

        // skip calculations for trivial case
        if (bunch_charge == 0_prt) { return; }

        // initialize the linear transport map
        Map6x6 R = Map6x6::Identity();

        // physical constants and reference quantities
        using ablastr::constant::SI::c;
        using ablastr::constant::SI::ep0;
        using ablastr::constant::math::pi;
        amrex::ParticleReal const mass = refpart.mass;
        amrex::ParticleReal const charge = refpart.charge;
        amrex::ParticleReal const pt_ref = refpart.pt;
        amrex::ParticleReal const betgam2 = powi<2>(pt_ref) - 1_prt;

        // evaluate the 3D space charge intensity parameter from bunch charge
        amrex::ParticleReal const rcN = std::abs(charge * bunch_charge) / (4_prt * pi * ep0 * mass * powi<2>(c));
        amrex::ParticleReal const coeff = ds * rcN / betgam2 * (1_prt/(5_prt * std::sqrt(5_prt)));

        // set parameters for elliptic integrals
        amrex::ParticleReal const errtol = 1.0e-3;
        amrex::ParticleReal const x = cm(1,1);
        amrex::ParticleReal const y = cm(3,3);
        amrex::ParticleReal const z = betgam2 * cm(5,5);

        // the existing model assumes <xy> = <yt> = <tx> = 0
        amrex::ParticleReal const xy = x*y;
        amrex::ParticleReal const yz = y*z;
        amrex::ParticleReal const zx = z*x;
        amrex::ParticleReal const corr_xy = (xy==0.0)? 0.0 : std::abs(cm(1,3)/std::sqrt(xy));
        amrex::ParticleReal const corr_yz = (yz==0.0)? 0.0 : std::abs(cm(3,5)/std::sqrt(yz));
        amrex::ParticleReal const corr_zx = (zx==0.0)? 0.0 : std::abs(cm(5,1)/std::sqrt(zx));
        if (corr_xy > errtol || corr_yz > errtol || corr_zx > errtol) {
            ablastr::warn_manager::WMRecordWarning(
                "algo.space_charge",
                "Space charge 3D model in envelope tracking assumes <xy> = <yt> = <tx> = 0 "
                 "but nonzero correlations are present.",
                ablastr::warn_manager::WarnPriority::high
            );
        }

        // evaluate the off-identity elements of the linear transfer map
        R(2,1) = coeff * Elliptic_RD(y,z,x,errtol);
        R(4,3) = coeff * Elliptic_RD(z,x,y,errtol);
        R(6,5) = coeff * betgam2 * Elliptic_RD(x,y,z,errtol);

        // update the beam covariance matrix
        cm = R * cm * R.transpose();

    }


} // namespace impactx::spacecharge
