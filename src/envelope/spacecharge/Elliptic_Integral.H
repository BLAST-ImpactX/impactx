/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef ELLIPTIC_INTEGRAL_H
#define ELLIPTIC_INTEGRAL_H

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_Print.H>

#include <cmath>

namespace impactx::envelope::spacecharge
{

    /** Function to return the Carlson elliptic integral denoted
     *  by RD(x,y,z).  This is a slight modification of
     *  the C++ function written by John Burkardt (GNU LGPL license):
     *
     *  https://people.math.sc.edu/Burkardt/cpp_src
     *        /elliptic_integral/elliptic_integral.html
     *
     *  which in turn is based on:
     *
     * Bille Carlson,
     * Computing Elliptic Integrals by Duplication,
     * Numerische Mathematik,
     * Volume 33, 1979, pages 1-16.
     *
     * Bille Carlson, Elaine Notis,
     * Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
     * ACM Transactions on Mathematical Software,
     * Volume 7, Number 3, pages 398-403, September 1981.
     *
     * @param[in] x first argument
     * @param[in] y second argument
     * @param[in] z third argument
     * @param[in] errtol error tolerance
     * @returns the elliptic integral RD(x,y,z)
     */
    amrex::ParticleReal
    Elliptic_RD (
        amrex::ParticleReal x,
        amrex::ParticleReal y,
        amrex::ParticleReal z,
        amrex::ParticleReal errtol
    )
    {
        using namespace amrex::literals;

        amrex::ParticleReal c1;
        amrex::ParticleReal c2;
        amrex::ParticleReal c3;
        amrex::ParticleReal c4;
        amrex::ParticleReal ea;
        amrex::ParticleReal eb;
        amrex::ParticleReal ec;
        amrex::ParticleReal ed;
        amrex::ParticleReal ef;
        amrex::ParticleReal epslon;
        amrex::ParticleReal lamda;
        amrex::ParticleReal const lolim = 6.0E-51;
        amrex::ParticleReal mu;
        amrex::ParticleReal power4;
        amrex::ParticleReal sigma;
        amrex::ParticleReal s1;
        amrex::ParticleReal s2;
        amrex::ParticleReal const uplim = 1.0E+48;
        amrex::ParticleReal value;
        amrex::ParticleReal xn;
        amrex::ParticleReal xndev;
        amrex::ParticleReal xnroot;
        amrex::ParticleReal yn;
        amrex::ParticleReal yndev;
        amrex::ParticleReal ynroot;
        amrex::ParticleReal zn;
        amrex::ParticleReal zndev;
        amrex::ParticleReal znroot;
      //
      //  `lolim` and `uplim` determine the range of valid arguments.
      //  `lolim` is not less than 2 / (machine machimum) ^ (2/3).
      //  `uplim` is not greater than (0.1 * errtol / machine minimum)^(2/3), 
      //
        if (
          x < 0.0_prt ||
          y < 0.0_prt ||
          x + y < lolim ||
          z < lolim ||
          uplim < x ||
          uplim < y ||
          uplim < z )
        {
          ablastr::warn_manager::WMRecordWarning(
              "algo.space_charge",
              "Invalid input arguments in function Elliptic_RD.",
              ablastr::warn_manager::WarnPriority::high); 
          value = 0.0_prt;
        }

        xn = x;
        yn = y;
        zn = z;
        sigma = 0.0_prt;
        power4 = 1.0_prt;

        int max_iterations = 100;
        for ( int i = 0; i < max_iterations; ++i )
        {
          mu = ( xn + yn + 3.0_prt * zn ) * 0.2_prt;
          xndev = ( mu - xn ) / mu;
          yndev = ( mu - yn ) / mu;
          zndev = ( mu - zn ) / mu;
          epslon = std::max ( std::abs ( xndev ),
            std::max ( std::abs ( yndev ), std::abs ( zndev ) ) );

          if ( epslon < errtol )
          {
            c1 = 3.0_prt / 14.0_prt;
            c2 = 1.0_prt / 6.0_prt;
            c3 = 9.0_prt / 22.0_prt;
            c4 = 3.0_prt / 26.0_prt;
            ea = xndev * yndev;
            eb = zndev * zndev;
            ec = ea - eb;
            ed = ea - 6.0_prt * eb;
            ef = ed + ec + ec;
            s1 = ed * ( - c1 + 0.25_prt * c3 * ed - 1.5_prt * c4 * zndev * ef );
            s2 = zndev  * ( c2 * ef + zndev * ( - c3 * ec + zndev * c4 * ea ) );
            value = 3.0_prt * sigma  + power4 * ( 1.0_prt + s1 + s2 ) / ( mu * std::sqrt ( mu ) );
            return value;
          }

          xnroot = std::sqrt ( xn );
          ynroot = std::sqrt ( yn );
          znroot = std::sqrt ( zn );
          lamda = xnroot * ( ynroot + znroot ) + ynroot * znroot;
          sigma = sigma + power4 / ( znroot * ( zn + lamda ) );
          power4 = power4 * 0.25_prt;
          xn = ( xn + lamda ) * 0.25_prt;
          yn = ( yn + lamda ) * 0.25_prt;
          zn = ( zn + lamda ) * 0.25_prt;
        }

        ablastr::warn_manager::WMRecordWarning(
            "algo.space_charge",
            "Exceeded maximum number of iterations in function Elliptic_RD.",
            ablastr::warn_manager::WarnPriority::medium);
        value = 0.0_prt;
        return value;

      }

} // namespace impactx::envelope::spacecharge

#endif // ELLIPTIC_INTEGRAL_H
