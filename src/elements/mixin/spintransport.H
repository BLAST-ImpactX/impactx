/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_SPIN_TRANSPORT_H
#define IMPACTX_ELEMENTS_MIXIN_SPIN_TRANSPORT_H

#include "particles/ImpactXParticleContainer.H"

#include <ablastr/constant.H>

#include <AMReX_Math.H>
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>


namespace impactx::elements::mixin
{
    /** This is a helper class for Spin Transport
     *
     * If an element derives from this class, then it must implement the `spin_and_phasespace_push` method.
     *
     * It also provides a helper function `rotate_spin` for applying a spin map generated
     * by a three-component vector lambda, whose magnitude defines the angle of rotation
     * and whose direction defines the axis of rotation.
     */
    struct SpinTransport
    {
        /** Rotate the spin vector (sx,sy,sz).
         *
         * @param[in] lambdax generator x-component
         * @param[in] lambday generator y-component
         * @param[in] lambdaz generator z-component
         * @param[inout] sx spin vector x-component
         * @param[inout] sy spin vector y-component
         * @param[inout] sz spin vector z-component
         */
        template<typename T_Real=amrex::ParticleReal>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void rotate_spin (
            T_Real const & AMREX_RESTRICT lambdax,
            T_Real const & AMREX_RESTRICT lambday,
            T_Real const & AMREX_RESTRICT lambdaz,
            T_Real & AMREX_RESTRICT sx,
            T_Real & AMREX_RESTRICT sy,
            T_Real & AMREX_RESTRICT sz
        ) const
        {
            using namespace amrex::literals;  // for _prt
            using namespace std; // for cmath(float)

            // Compute quaternion coefficients
            T_Real const angle = sqrt(lambdax*lambdax + lambday*lambday + lambdaz*lambdaz);

            auto const [sin_half_angle, cos_half_angle] = amrex::Math::sincos(angle*0.5_prt);
            T_Real const w0 = cos_half_angle;
            T_Real const c_sin_half = -sin_half_angle / angle;

            auto const non_zero_angle = angle!=0_prt;
#ifdef AMREX_USE_SIMD
            T_Real w1 = 0_prt;
            amrex::simd::stdx::where(non_zero_angle, w1) = lambdax * c_sin_half;
            T_Real w2 = 0_prt;
            amrex::simd::stdx::where(non_zero_angle, w1) = lambday * c_sin_half;
            T_Real w3 = 0_prt;
            amrex::simd::stdx::where(non_zero_angle, w1) = lambdaz * c_sin_half;
#else
            T_Real w1 = non_zero_angle ? lambdax * c_sin_half : 0_prt;
            T_Real w2 = non_zero_angle ? lambday * c_sin_half : 0_prt;
            T_Real w3 = non_zero_angle ? lambdaz * c_sin_half : 0_prt;
#endif

            // Apply rotation
            T_Real const sxout = (1_prt-2_prt*(w2*w2+w3*w3))*sx + 2_prt*(w1*w2+w0*w3)*sy + 2_prt*(w1*w3-w0*w2)*sz;
            T_Real const syout = 2_prt*(w1*w2-w0*w3)*sx + (1_prt-2_prt*(w1*w1+w3*w3))*sy + 2_prt*(w0*w1+w2*w3)*sz;
            T_Real const szout = 2_prt*(w0*w2+w1*w3)*sx + 2_prt*(w2*w3-w0*w1)*sy + (1_prt-2_prt*(w1*w1+w2*w2))*sz;

            // Update spin vector
            sx = sxout;
            sy = syout;
            sz = szout;
        }
    };

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_SPIN_TRANSPORT_H
