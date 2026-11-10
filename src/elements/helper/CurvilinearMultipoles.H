/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_HELPER_CURVILINEAR_MULTIPOLES_H
#define IMPACTX_ELEMENTS_HELPER_CURVILINEAR_MULTIPOLES_H

#include <AMReX_ConstexprFor.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>
#include <type_traits>


/** Curvilinear multipoles of a sector bend
 *
 * In a bend of curvature h, with rho = 1 + h x, the multipoles of T. Zolkin,
 * Phys. Rev. Accel. Beams 20, 043501 (2017) are the midplane expansion of one family of
 * functions F_k(x), k = 1, 2, ...:
 *
 *   F_1 = 1,  F_2 = ln(rho)/h,  (d^2/dx^2 + h/rho d/dx) F_{k+2} = F_k,  F_k(0) = F_k'(0) = 0.
 *
 * The normal multipole of order n (n = 1: dipole) has the vertical field
 * By = sum_q (-1)^q y^2q/(2q)! F_{n-2q} on and off the midplane, and the skew multipole of
 * order n the scalar potential sum_q (-1)^q y^2q/(2q)! F_{n+1-2q}.
 *
 * With u = ln(rho) and xi = u/h, F_k = xi^(k-1) P_k(u). P_k(u) is an entire function
 * with positive Taylor coefficients (h^(k-1) F_k is the divided difference of rho^z over
 * the nodes z = 0, 0, 2, 2, 4, 4, ...):
 *
 *   p_{1,i} = delta_{i,0},  p_{k,i} = (p_{k-1,i} + z_k p_{k,i-1}) / (k - 1 + i),
 *   z_k = 2 floor((k - 1)/2),
 *
 * and the derivative follows from the same functions: rho F_k' = xi^(k-2) Q_k(u) with
 * Q_k = P_{k-1} + z_k u P_k. Evaluated in these terms, the field has no cancellation
 * between terms of order one, at any curvature, including h = 0.
 */
namespace impactx::elements::curvilinear
{
    /** Highest coefficient index (dipole: 0) */
    inline constexpr int max_index = 11;

    /** Highest function index: the skew multipole of the highest order needs F_{n+1} */
    inline constexpr int num_functions = max_index + 2;

    /** Number of tabulated Taylor coefficients per function */
    inline constexpr int max_terms = 40;

    /** The series of P_k are truncated to the precision of the particle data for
     *  |u| = |ln(1 + h x)| <= u_max, i.e., -0.22 < h x < 0.28
     */
    inline constexpr double u_max = 0.25;

    /** Node z_k of the divided difference */
    constexpr int node (int k)
    {
        return 2 * ((k - 1) / 2);
    }

    /** Taylor coefficients p[k][i] of P_k(u) */
    struct TaylorTable
    {
        double p[num_functions + 1][max_terms] = {};
    };

    constexpr TaylorTable make_taylor_table ()
    {
        TaylorTable table{};
        table.p[1][0] = 1.0;
        for (int k = 2; k <= num_functions; ++k)
        {
            for (int i = 0; i < max_terms; ++i)
            {
                double const previous = i > 0 ? table.p[k][i - 1] : 0.0;
                table.p[k][i] = (table.p[k - 1][i] + node(k) * previous) / (k - 1 + i);
            }
        }
        return table;
    }

    inline constexpr TaylorTable taylor = make_taylor_table();

    /** Number of Taylor terms of P_k for a truncation error below a quarter ulp of T_Real,
     *  for |u| <= u_max
     *
     * The coefficients are positive, so the sum of the omitted terms at |u| = u_max bounds
     * the truncation error for both signs of u, and P_k is smallest at u = -u_max.
     */
    template<typename T_Real>
    constexpr int num_terms (int k)
    {
        double const tolerance = 0.25 * static_cast<double>(std::numeric_limits<T_Real>::epsilon());

        double p_min = 0.0;
        double power = 1.0;
        for (int i = 0; i < max_terms; ++i)
        {
            p_min += taylor.p[k][i] * power;
            power *= -u_max;
        }

        for (int m = 1; m < max_terms; ++m)
        {
            double tail = 0.0;
            power = 1.0;
            for (int i = 0; i < m; ++i) { power *= u_max; }
            for (int i = m; i < max_terms; ++i)
            {
                tail += taylor.p[k][i] * power;
                power *= u_max;
            }
            if (tail <= tolerance * p_min) { return m; }
        }
        return max_terms;
    }

    static_assert(num_terms<double>(num_functions) < max_terms,
                  "Tabulate more Taylor coefficients for the series of the highest function.");

    /** 1/n! */
    constexpr double inverse_factorial (int n)
    {
        double factorial = 1.0;
        for (int i = 2; i <= n; ++i) { factorial *= i; }
        return 1.0 / factorial;
    }

    /** P_k(u) by Horner's scheme
     *
     * @tparam K function index
     * @param u logarithm ln(1 + h x)
     */
    template<int K, typename T_Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T_Real series (T_Real const & u)
    {
        constexpr int M = num_terms<amrex::ParticleReal>(K);
        T_Real sum = 0;
        amrex::constexpr_for<0, M>([&](auto j)
        {
            constexpr double coefficient = taylor.p[K][M - 1 - decltype(j)::value];
            sum = sum * u + static_cast<amrex::ParticleReal>(coefficient);
        });
        return sum;
    }

    /** Momentum kick per unit length of all multipoles except the normal dipole
     *
     * The same as @see kick, with the number of coefficients known at compile time.
     *
     * @tparam N number of coefficients, i.e., the highest multipole order
     */
    template<int N, typename T_Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void kick_order (
        T_Real const & x,
        T_Real const & y,
        amrex::ParticleReal h,
        amrex::ParticleReal const * AMREX_RESTRICT k_normal,
        amrex::ParticleReal const * AMREX_RESTRICT k_skew,
        amrex::ParticleReal scale,
        T_Real & AMREX_RESTRICT dpx,
        T_Real & AMREX_RESTRICT dpy
    )
    {
        using namespace amrex::literals; // for _prt
        using std::log1p;

        // the skew multipole of order N needs the functions up to F_{N+1}
        constexpr int K = N + 1;
        using Functions = T_Real[K + 1];

        // logarithmic radial coordinate xi = ln(1 + h x)/h, which is x for h = 0
        T_Real const X = h * x;
        T_Real const u = log1p(X);
        T_Real const xi = X != 0_prt ? x * (u / X) : x;
        T_Real const rho = 1_prt + X;
        T_Real const w = y * y;

        // P_k(u) and Q_k(u)
        Functions P = {};
        Functions Q = {};
        amrex::constexpr_for<1, K + 1>([&](auto k_)
        {
            constexpr int k = decltype(k_)::value;
            P[k] = series<k>(u);
            if constexpr (k >= 2) {
                Q[k] = P[k - 1] + static_cast<amrex::ParticleReal>(node(k)) * u * P[k];
            }
        });

        // sum_a k_n f_{a + Shift} xi^a with n = a + Offset + 2q, by Horner's scheme in xi;
        // f is P (Shift = 1) or Q (Shift = 2), and the orders below N_min are left out
        auto const horner_xi = [&](amrex::ParticleReal const * k_array, Functions const & functions,
                                   auto offset_, auto shift_, auto n_min_, auto q_) -> T_Real
        {
            constexpr int Offset = decltype(offset_)::value;
            constexpr int Shift = decltype(shift_)::value;
            constexpr int N_min = decltype(n_min_)::value;
            constexpr int q = decltype(q_)::value;
            T_Real sum = 0;
            amrex::constexpr_for<0, N>([&](auto j)
            {
                constexpr int a = N - 1 - decltype(j)::value;
                constexpr int n = a + Offset + 2 * q;
                if constexpr (n <= N && n >= N_min) {
                    sum = sum * xi + (k_array[n - 1] * scale) * functions[a + Shift];
                } else if constexpr (n <= N) {
                    sum = sum * xi;
                }
            });
            return sum;
        };

        // sum_q (-1)^q w^q/(2q + Parity)! S_q(xi), by Horner's scheme in w = y^2
        auto const horner_w = [&](auto parity_, auto const & S) -> T_Real
        {
            constexpr int Parity = decltype(parity_)::value;
            T_Real sum = 0;
            amrex::constexpr_for<0, (N + 1) / 2>([&](auto j)
            {
                constexpr int q = (N + 1) / 2 - 1 - decltype(j)::value;
                constexpr double c = (q % 2 == 0 ? 1.0 : -1.0) * inverse_factorial(2 * q + Parity);
                auto const S_q = S(std::integral_constant<int, q>{});
                sum = sum * w + static_cast<amrex::ParticleReal>(c) * S_q;
            });
            return sum;
        };

        using I0 = std::integral_constant<int, 0>;
        using I1 = std::integral_constant<int, 1>;
        using I2 = std::integral_constant<int, 2>;

        // normal, By: xi^a P_{a+1} with n = a + 1 + 2q, without the dipole
        T_Real const by_normal = horner_w(I0{}, [&](auto q) {
            return horner_xi(k_normal, P, I1{}, I1{}, I2{}, q); });
        // normal, rho Bx / y: xi^a Q_{a+2} with n = a + 2 + 2q
        T_Real const bx_normal = horner_w(I1{}, [&](auto q) {
            return horner_xi(k_normal, Q, I2{}, I2{}, I1{}, q); });
        // skew, -By / y: xi^a P_{a+1} with n = a + 2 + 2q
        T_Real const by_skew = horner_w(I1{}, [&](auto q) {
            return horner_xi(k_skew, P, I2{}, I1{}, I1{}, q); });
        // skew, rho Bx: xi^a Q_{a+2} with n = a + 1 + 2q
        T_Real const bx_skew = horner_w(I0{}, [&](auto q) {
            return horner_xi(k_skew, Q, I1{}, I2{}, I1{}, q); });

        dpx = -rho * (by_normal - y * by_skew);
        dpy = y * bx_normal + bx_skew;
    }

    /** Momentum kick per unit length of all multipoles except the normal dipole
     *
     * The kick is dpx = -(1 + h x) By and dpy = (1 + h x) Bx, with the field normalized by
     * the magnetic rigidity. The normal dipole is the curvature of the design orbit, which
     * the sector bend map of the element applies.
     *
     * @param x horizontal position in m
     * @param y vertical position in m
     * @param h curvature of the design orbit in 1/m; any sign, or zero
     * @param k_normal normal coefficients for the indices m = 0 .. ncoef-1, in 1/m^(m+1)
     *                 after multiplication with scale
     * @param k_skew skew coefficients for the indices m = 0 .. ncoef-1, in 1/m^(m+1)
     *               after multiplication with scale
     * @param ncoef number of coefficients, 1 <= ncoef <= max_index + 1
     * @param scale factor of the coefficients, e.g., the inverse magnetic rigidity
     * @param[out] dpx horizontal kick per unit length in 1/m
     * @param[out] dpy vertical kick per unit length in 1/m
     */
    template<typename T_Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void kick (
        T_Real const & x,
        T_Real const & y,
        amrex::ParticleReal h,
        amrex::ParticleReal const * AMREX_RESTRICT k_normal,
        amrex::ParticleReal const * AMREX_RESTRICT k_skew,
        int ncoef,
        amrex::ParticleReal scale,
        T_Real & AMREX_RESTRICT dpx,
        T_Real & AMREX_RESTRICT dpy
    )
    {
        dpx = 0;
        dpy = 0;

        // one evaluation for each number of coefficients, which sets the work at compile time
        amrex::constexpr_for<1, max_index + 2>([&](auto n_)
        {
            constexpr int n = decltype(n_)::value;
            if (n == ncoef) {
                kick_order<n>(x, y, h, k_normal, k_skew, scale, dpx, dpy);
            }
        });
    }

} // namespace impactx::elements::curvilinear

#endif // IMPACTX_ELEMENTS_HELPER_CURVILINEAR_MULTIPOLES_H
