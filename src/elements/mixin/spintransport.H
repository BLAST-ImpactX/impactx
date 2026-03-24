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
            amrex::simd::stdx::where(non_zero_angle, w2) = lambday * c_sin_half;
            T_Real w3 = 0_prt;
            amrex::simd::stdx::where(non_zero_angle, w3) = lambdaz * c_sin_half;
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

        /** Evaluate the Thomas-BMT precession vector
         *
          * Evaluate in a given field,
         *  based on the instantaneous particle orbit variables.  The model
         *  below is evaluated in a Frenet-Serret coordinate system,
         *  assuming a reference trajectory of curvature h, with (possible)
         *  bending in the x-z plane.  The field components Bx,By,Bz and
         *  are normalized by q/mc (units of 1/m), and the field components
         *  Ex,Ey,Ez are normalized by q/mc^2 (units of 1/m).
         *
         * @param[in] x horizontal coordinate (in local coordinate system)
         * @param[in] ux velocity unit direction x-component
         * @param[in] uy velocity unit direction y-component
         * @param[in] uz velocity unit direction z-component
         * @param[in] gamma relativistic gamma
         * @param[in] h design orbit curvature [1/m]
         * @param[in] gyro_anomaly gyromagnetic anomaly
         * @param[in] Bx magnetic field x-component [1/m, dynamic units]
         * @param[in] By magnetic field y-component [1/m, dynamic units]
         * @param[in] Bz magnetic field z-component [1/m, dynamic units]
         * @param[in] Ex electric field x-component [1/m, dynamic units]
         * @param[in] Ey electric field y-component [1/m, dynamic units]
         * @param[in] Ez electric field z-component [1/m, dynamic units]
         * @param[out] Omegax precession vector x-component
         * @param[out] Omegay precession vector y-component
         * @param[out] Omegaz precession vector z-component
         */
        template<typename T_Real=amrex::ParticleReal>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void tbmt_precession_vector (
            T_Real const & AMREX_RESTRICT x,
            T_Real const & AMREX_RESTRICT ux,
            T_Real const & AMREX_RESTRICT uy,
            T_Real const & AMREX_RESTRICT uz,
            T_Real const & AMREX_RESTRICT gamma,
            amrex::ParticleReal const & AMREX_RESTRICT h,
            amrex::ParticleReal const & AMREX_RESTRICT gyro_anomaly,
            T_Real const & AMREX_RESTRICT Bx,
            T_Real const & AMREX_RESTRICT By,
            T_Real const & AMREX_RESTRICT Bz,
            T_Real const & AMREX_RESTRICT Ex,
            T_Real const & AMREX_RESTRICT Ey,
            T_Real const & AMREX_RESTRICT Ez,
            T_Real & AMREX_RESTRICT Omegax,
            T_Real & AMREX_RESTRICT Omegay,
            T_Real & AMREX_RESTRICT Omegaz
        ) const
        {
            using namespace amrex::literals;  // for _prt
            using namespace std; // for cmath(float)

            // Evaluation of the full Thomas-BMT precession vector in the Frenet-Serret frame.
            T_Real const beta_gamma = sqrt(gamma*gamma - 1_prt);
            T_Real const beta_gamma_z = uz * beta_gamma;
            T_Real const leading_factor = -(1_prt + h * x) / beta_gamma_z;
            T_Real const udotB = ux*Bx + uy*By + uz*Bz;
            T_Real const ucrossE_x = uy*Ez - uz*Ey;
            T_Real const ucrossE_y = uz*Ex - ux*Ez;
            T_Real const ucrossE_z = ux*Ey - uy*Ex;
            T_Real const c1 = (1_prt + gyro_anomaly*gamma);
            T_Real const c2 = gyro_anomaly * (gamma - 1_prt) * udotB;
            T_Real const c3 = (gyro_anomaly + 1_prt / (1_prt + gamma)) * beta_gamma;
            Omegax = leading_factor * (c1 * Bx - c2 * ux - c3 * ucrossE_x);
            Omegay = leading_factor * (c1 * By - c2 * uy - c3 * ucrossE_y) + h;
            Omegaz = leading_factor * (c1 * Bz - c2 * uz - c3 * ucrossE_z);
        }

    };

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_SPIN_TRANSPORT_H
