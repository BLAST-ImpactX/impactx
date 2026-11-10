/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_HELPER_STUMPFF_H
#define IMPACTX_ELEMENTS_HELPER_STUMPFF_H

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <type_traits>


namespace impactx::elements
{
    /** Stumpff functions c_0, ..., c_3 of one argument
     *
     * @see stumpff
     */
    template<typename T_Real>
    struct Stumpff
    {
        T_Real c0;  //! cos(w)
        T_Real c1;  //! sin(w) / w
        T_Real c2;  //! (1 - cos(w)) / w^2
        T_Real c3;  //! (w - sin(w)) / w^3
    };

    /** Stumpff functions c_0, ..., c_3 from their series, up to the term (-z)^n
     *
     * @tparam n highest power of -z in the series of c_2 and c_3
     * @param z argument
     * @return c_0(z), c_1(z), c_2(z), c_3(z)
     */
    template<int n>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Stumpff<amrex::ParticleReal>
    stumpff_series (amrex::ParticleReal z)
    {
        using namespace amrex::literals; // for _rt and _prt
        using T = amrex::ParticleReal;

        // Horner scheme in -z, with the ratios of successive terms
        T c2 = 1_prt;
        T c3 = 1_prt;
        for (int j = n; j >= 1; --j)
        {
            c2 = 1_prt - z * c2 * T(1.0 / ((2 * j + 1) * (2 * j + 2)));
            c3 = 1_prt - z * c3 * T(1.0 / ((2 * j + 2) * (2 * j + 3)));
        }
        c2 *= 0.5_prt;
        c3 *= T(1.0 / 6.0);

        return {1_prt - z * c2, 1_prt - z * c3, c2, c3};
    }

    /** Stumpff functions c_0, ..., c_3
     *
     *   c_k(z) = sum_{j >= 0} (-z)^j / (2j + k)!
     *
     * For z = w^2 > 0, these are cos(w), sin(w)/w, (1 - cos(w))/w^2 and (w - sin(w))/w^3,
     * and for z = -w^2 < 0 their hyperbolic counterparts. They are the entries of the
     * transfer map of a harmonic oscillator x'' = -g x over a length s, with z = g s^2,
     * at any sign of g, including g = 0.
     *
     * For |z| <= 4, the series is summed directly, so that the small differences
     * 1 - cos(w) and w - sin(w) have no cancellation, also with approximate (fast-math)
     * trigonometric functions.
     *
     * @param z argument
     * @return c_0(z), c_1(z), c_2(z), c_3(z)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Stumpff<amrex::ParticleReal>
    stumpff (amrex::ParticleReal z)
    {
        using namespace amrex::literals; // for _rt and _prt
        using namespace std; // for cmath(float)
        using T = amrex::ParticleReal;

        // the series, up to the first term below a quarter of the rounding error
        constexpr bool is_float = std::is_same_v<T, float>;
        if (abs(z) <= 0.0625_prt)
        {
            return stumpff_series<is_float ? 2 : 5>(z);
        }
        if (abs(z) <= 4_prt)
        {
            return stumpff_series<is_float ? 6 : 10>(z);
        }

        T c0, c1, c2;
        if (z > 0_prt)
        {
            T const w = sqrt(z);
            auto const [sin_half, cos_half] = amrex::Math::sincos(0.5_prt * w);
            c0 = 1_prt - 2_prt * sin_half * sin_half;
            c1 = 2_prt * sin_half * cos_half / w;
            c2 = 2_prt * sin_half * sin_half / z;
        }
        else
        {
            T const w = sqrt(-z);
            T const sinh_half = sinh(0.5_prt * w);
            c0 = 1_prt + 2_prt * sinh_half * sinh_half;
            c1 = sinh(w) / w;
            c2 = -2_prt * sinh_half * sinh_half / z;
        }
        // |c1| <= 1/2 for z > 4 and c1 >= 1.8 for z < -4: no cancellation in 1 - c1
        return {c0, c1, c2, (1_prt - c1) / z};
    }

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_HELPER_STUMPFF_H
