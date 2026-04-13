#ifndef FF_ALGORITHM_H
#define FF_ALGORITHM_H

#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

// #include "basic_toolkit/PhysicsConstants.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/kokkos_types.h"

#include <Kokkos_Core.hpp>

#define LIBFF_USE_GSV 1

inline bool
close_to_zero(double v)
{
    return fabs(v) < 1e-13;
}

namespace FF_algorithm {
    using namespace kt;

    // default yoshida order and steps
    constexpr int default_steps = 6;
    constexpr int default_order = 4;

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    invsqrt(T const& x)
    {
        return 1.0 / sqrt(x);
    }

    KOKKOS_INLINE_FUNCTION
    constexpr double
    quiet_nan()
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // exact solution for drift spaces
    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    drift_unit(T& x,
               T const& xp,
               T& y,
               T const& yp,
               T& cdt,
               T const& dpop,
               double length,
               double reference_momentum,
               double m,
               double reference_cdt)
    {
        T uni(1.0);
        T sig((0.0 < length) - (length < 0.0));

        T vl(length);
        T vm(m);
        T vrm(reference_momentum);
        T vrc(reference_cdt);

        T dp = dpop + uni;
        T inv_npz = uni / sqrt(dp * dp - xp * xp - yp * yp);
        T lxpr = xp * vl * inv_npz;
        T lypr = yp * vl * inv_npz;
        T D2 = lxpr * lxpr + vl * vl + lypr * lypr;
        T p = dp * vrm;
        T E2 = p * p + vm * vm;
        // T beta2 = p*p / E2;
        T ibeta2 = E2 / (p * p);
        x = x + lxpr;
        y = y + lypr;
        // cdt += sqrt(D2 / beta2) - reference_cdt;
        cdt = cdt + sig * sqrt(D2 * ibeta2) - vrc;
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    slot_unit(T& x,
              T& xp,
              T& y,
              T& yp,
              T& cdt,
              T& dpop,
              double ct,
              double st,
              double pref,
              double m)
    {
        const T uni(1.0);

        const T vct(ct);
        const T vst(st);
        const T vpref(pref);
        const T vm(m);

        T r0 = x;
        T r1 = y;
        T r2 = 0.0;

        T zp = sqrt((dpop + uni) * (dpop + uni) - xp * xp - yp * yp);

        T p = (dpop + uni) * vpref;
        T e = sqrt(p * p + vm * vm);

        T b0 = xp * vpref / e;
        T b1 = yp * vpref / e;
        T b2 = zp * vpref / e;

        T bp = vst * b0 + vct * b2;

        T tau = -x * vst / bp;

        r0 = x + tau * b0;
        r1 = y + tau * b1;
        r2 = tau * b2;

        x = r0 * vct - r2 * vst;
        y = r1;

        cdt = cdt + tau;

        xp = xp * vct - zp * vst;
        yp = yp;
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    edge_unit(T const& y, T& yp, T k)
    {
        yp = yp - k * y;
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    edge_unit(T const& y, T& xp, T& yp, double kx, double ky, char)
    {
        const T vkx(kx);
        const T vky(ky);

        // this is the one that gets called from ff_sbend
        yp = yp - vkx * y;
        xp = xp + vky * y;
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    edge_unit(T const& y, T& xp, T& yp, T const& dpop, T k)
    {
        const T uni(1.0);

        T zp = sqrt((dpop + uni) * (dpop + uni) - xp * xp - yp * yp);

        T xxp = xp;
        T yyp = yp;

        yp = yp - (xxp / zp) * k * y;
        xp = xp + (yyp / zp) * k * y;
    }

    // exact solution for dipole without high order combined functions
    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    dipole_unit(T& x,
                T& xp,
                T& y,
                T& yp,
                T& cdt,
                T const& dpop,
                double l,
                double k0)
    {
        T xp1 = xp - k0 * l / (dpop + 1.0);
        T yp1 = yp;

        T poeb = (dpop + 1.0) / k0; // p/eB

        x += poeb *
             (sqrt(1 - xp1 * xp1 - yp1 * yp1) - sqrt(1 - xp * xp - yp * yp));
        y += poeb * yp *
             (atan(xp / sqrt(1 - xp * xp - yp * yp)) -
              atan(xp1 / sqrt(1 - xp1 * xp1 - yp1 * yp1)));

        xp = xp1;
    }

    // exact solution for dipoles, comes from CHEF
    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    bend_complete(T& x,
                  T& xp,
                  T& y,
                  T& yp,
                  T& cdt,
                  T const& dpop,
                  double dphi,
                  double strength,
                  double p_ref,
                  double m,
                  double cdt_ref,
                  Kokkos::complex<double> phase,
                  Kokkos::complex<double> term)
    {
        typedef Kokkos::complex<T> CT;

        const T uni(1.0);

        T p0 = p_ref;
        T p = p_ref * (dpop + 1.0);
        T E0 = sqrt(p0 * p0 + m * m);
        T E = sqrt(p * p + m * m);

        T igamma = m / E0;
        T ibeta = uni / sqrt(uni - igamma * igamma);

        T csq = pconstants::c * pconstants::c * 1e-9;
        T psq = (dpop + 1.0) * (dpop + 1.0);

        T Ef = uni / sqrt(psq + igamma * igamma * ibeta * ibeta);

        T beta1 = Ef * xp;
        T beta2 = Ef * yp;
        T beta3 = Ef * sqrt((dpop + 1.0) * (dpop + 1.0) - xp * xp - yp * yp);

        CT ui = CT(0.0, x);
        CT vui = CT(pconstants::c * beta3, pconstants::c * beta1);

        T iomega = E / (csq * strength);

        CT bi = CT(0.0, 1.0) * vui * iomega - ui;
        CT bf = bi * phase + term;

        T rho = pconstants::c * sqrt(beta1 * beta1 + beta3 * beta3) * iomega;

        T dthmphi = asin(bi.real() / rho) - asin(bf.real() / rho);

        CT expf = Kokkos::exp(CT(0.0, dthmphi));
        CT vuf = vui * expf;
        CT uf = (ui + bi) * expf - bf;

        T dtheta = dthmphi + dphi;
        T ncdt = -pconstants::c * dtheta * iomega;

        x = uf.imag();
        y += beta2 * ncdt;
        cdt += ncdt - cdt_ref;
        xp = vuf.imag() / (Ef * pconstants::c);
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    bend_unit(T& x,
              T& xp,
              T& y,
              T& yp,
              T& cdt,
              T const& dpop,
              double theta,
              double strength,
              double p_ref,
              double m,
              double cdt_ref,
              Kokkos::complex<double> phase,
              Kokkos::complex<double> term)
    {
        typedef Kokkos::complex<T> CT;

        const T uni(1.0);
        const T tc(pconstants::c);

        const T vtheta(theta);
        const T vstrength(strength);
        const T vp_ref(p_ref);
        const T vm(m);
        const T vcdt_ref(cdt_ref);

        T p0 = vp_ref;
        T p = vp_ref * (dpop + uni);
        T E0 = sqrt(p0 * p0 + vm * vm);
        T E = sqrt(p * p + vm * vm);

        T igamma = vm / E0;
        T ibeta = uni / sqrt(uni - igamma * igamma);

        T csq = T(pconstants::c * pconstants::c * 1e-9);
        T psq = (dpop + uni) * (dpop + uni);

        T Ef = uni / sqrt(psq + igamma * igamma * ibeta * ibeta);

        T beta1 = Ef * xp;
        T beta2 = Ef * yp;
        T beta3 = Ef * sqrt((dpop + uni) * (dpop + uni) - xp * xp - yp * yp);

        // CT ui  = CT(0.0, x);
        T ui_r = 0.0;
        T ui_i = x;

        // CT vui = CT(tc * beta3, tc * beta1);
        T vui_r = tc * beta3;
        T vui_i = tc * beta1;

        T iomega = E / (csq * vstrength);

        // CT bi = CT(0.0, 1.0) * vui * iomega - ui;
        T bi_r = -vui_i * iomega - ui_r;
        T bi_i = vui_r * iomega - ui_i;

        // CT bf = bi * phase + term;
        T phase_r(phase.real());
        T phase_i(phase.imag());
        T bf_r = bi_r * phase_r - bi_i * phase_i + T(term.real());
        T bf_i = bi_r * phase_i + bi_i * phase_r + T(term.imag());

        T rho = tc * sqrt(beta1 * beta1 + beta3 * beta3) * iomega;

        // T dthmphi = asin(bi.real() / rho) - asin(bf.real() / rho);
        T dthmphi = asin(bi_r / rho) - asin(bf_r / rho);

        // CT expf = Kokkos::exp( CT(0.0, dthmphi) );
        T expf_r = cos(dthmphi);
        T expf_i = sin(dthmphi);

        // CT vuf  = vui * expf;
        T vuf_i = vui_r * expf_i + vui_i * expf_r;

        // CT uf   = (ui + bi) * expf - bf;
        T uf_i = (ui_r + bi_r) * expf_i + (ui_i + bi_i) * expf_r - bf_i;

        T dtheta = dthmphi + vtheta;
        T ncdt = -tc * dtheta * iomega;

        // x = uf.imag();
        x = uf_i;

        // xp = vuf.imag() / (Ef * tc);
        xp = vuf_i / (Ef * tc);

        y = y + beta2 * ncdt;
        cdt = cdt + ncdt - vcdt_ref;
    }

    KOKKOS_INLINE_FUNCTION
    Kokkos::complex<double>
    bend_unit_phase(double theta)
    {
        return Kokkos::exp(Kokkos::complex<double>(0.0, theta));
    }

    KOKKOS_INLINE_FUNCTION
    Kokkos::complex<double>
    bend_edge_phase(double theta)
    {
        return Kokkos::exp(Kokkos::complex<double>(0.0, -theta));
    }

    KOKKOS_INLINE_FUNCTION
    Kokkos::complex<double>
    bend_unit_term(double r0, double theta)
    {
        return Kokkos::complex<double>(0.0, r0) *
               Kokkos::complex<double>(1.0 - cos(theta), -sin(theta));
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    bend_edge(T& x,
              T& xp,
              T& y,
              T& yp,
              T& cdt,
              T const& dpop,
              T e,
              Kokkos::complex<double> phase,
              T strength,
              T p_ref,
              T m)
    {
        typedef Kokkos::complex<T> CT;

        const T uni(1.0);
        const T tc(pconstants::c);

        T p0 = p_ref;
        T p = p_ref * (dpop + uni);
        T E0 = sqrt(p0 * p0 + m * m);
        T E = sqrt(p * p + m * m);

        T igamma = m / E0;
        T ibeta = uni / sqrt(uni - igamma * igamma);

        T csq = T(pconstants::c * pconstants::c * 1e-9);
        T psq = (dpop + uni) * (dpop + uni);

        T Ef = uni / sqrt(psq + igamma * igamma * ibeta * ibeta);

        T beta1 = Ef * xp;
        T beta2 = Ef * yp;
        T beta3 = Ef * sqrt((dpop + uni) * (dpop + uni) - xp * xp - yp * yp);

        // CT ui  = CT(0.0, x);
        T ui_r = 0.0;
        T ui_i = x;

        // CT vui = CT(pconstants::c * beta3, pconstants::c * beta1);
        T vui_r = tc * beta3;
        T vui_i = tc * beta1;

        T iomega = E / (csq * strength);

        // CT bi = CT(0.0, 1.0) * vui * iomega - ui;
        T bi_r = -vui_i * iomega - ui_r;
        T bi_i = vui_r * iomega - ui_i;

        // CT bf = bi * phase;
        T phase_r(phase.real());
        T phase_i(phase.imag());
        T bf_r = bi_r * phase_r - bi_i * phase_i;
        T bf_i = bi_r * phase_i + bi_i * phase_r;

        T rho = tc * sqrt(beta1 * beta1 + beta3 * beta3) * iomega;

        // T dthmphi = asin(bi.real() / rho) - asin(bf.real() / rho);
        T dthmphi = asin(bi_r / rho) - asin(bf_r / rho);

        // CT expf = Kokkos::exp( CT(0.0, dthmphi) );
        T expf_r = cos(dthmphi);
        T expf_i = sin(dthmphi);

        // CT vuf  = vui * expf;
        T vuf_i = vui_r * expf_i + vui_i * expf_r;

        // CT uf   = (ui + bi) * expf - bf;
        T uf_i = (ui_r + bi_r) * expf_i + (ui_i + bi_i) * expf_r - bf_i;

        T dtheta = dthmphi + e;
        T ncdt = -tc * dtheta * iomega;

        // x = uf.imag();
        x = uf_i;

        // xp = vuf.imag() / (Ef * pconstants::c);
        xp = vuf_i / (Ef * tc);

        y = y + beta2 * ncdt;
        cdt = cdt + ncdt;
    }

    // thin_cf_quadrupole_unit is not used.  It follows expressions in
    // CHEF's InducedKickPropagators which are now contained in
    // thin_cf_quadrupole_b{x|y} and retained here for historical reference
    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    thin_cf_quadrupole_unit(T const& x,
                            T& xp,
                            T const& y,
                            T& yp,
                            double r0,
                            double const* kL)
    {
        T vk0(kL[0]);
        T vk1(kL[1]);

        T alf = x / r0;

#define MODEL 10

#if MODEL == 10
        T Bx = r0 * (y / (r0 + x));
        T By = r0 * log(1.0 + alf);
        xp = xp - vk0 * By * (1.0 + alf);
        yp = yp + vk0 * Bx * (1.0 + alf);
#elif MODEL == 11
        T Bx = y * (1.0 - alf + alf * alf);
        T By = x * (1.0 - 0.5 * alf + alf * alf / 3.0);
        xp = xp - vk0 * By * (1.0 + alf);
        yp = yp + vk0 * Bx * (1.0 + alf);
#elif MODEL == 12
        T Bx = y * (1.0 - alf);
        T By = x * (1.0 - 0.5 * alf);
        xp = xp - vk0 * By * (1.0 + alf);
        yp = yp + vk0 * Bx * (1.0 + alf);
#elif MODEL == 13
        T Bx = y;
        T By = x;
        xp = xp - vk0 * By * (1.0 + alf);
        yp = yp + vk0 * Bx * (1.0 + alf);
#elif MODEL == 14
        T Bx = y;
        T By = x;
        xp = xp - vk0 * By;
        yp = yp + vk0 * Bx;
#endif
    }

    // the expressions in thin_cf_quadrupole_b{x|y} come from
    // Zolkin, Sector magnets or transverse electromagnetic fields in
    // cylindrical coordinates, Phys.Rev.Accel.Beams 20 (2017) no.4, 043501
    // Table XIII
    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    thin_cf_quadrupole_bx(T const& x, T const& y, T const& r0)
    {
        return r0 * y / (r0 + x);
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    thin_cf_quadrupole_by(T const& x, T const& y, T const& r0, T const& alf)
    {
        return r0 * log(T(1.0) + alf);
    }

    // the expressions in thin_cf_sectupole_b{x|y} come from
    // Zolkin, Sector magnets or transverse electromagnetic fields in
    // cylindrical coordinates, Phys.Rev.Accel.Beams 20 (2017) no.4, 043501
    // Table XIII
    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    thin_cf_sextupole_bx(T const& x, T const& y, T const& r0, T const& alf)
    {
        return x * (T(2.0) + alf) * y / (T(1.0) + alf);
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    thin_cf_sextupole_by(T const& x, T const& y, T const& r0, T const& alf)
    {
        return T(0.5) * (x * x + T(2.0) * x * r0) - y * y -
               r0 * r0 * log(T(1.0) + alf);
    }

    // combined function sbends kick up to quadrupole
    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    thin_cf_kick_1(T const& x,
                   T& xp,
                   T const& y,
                   T& yp,
                   double r0,
                   double const* kL)
    {
        const T uni(1.0);

        T vk1n(kL[0]);
        T vk1s(kL[1]);

        T vr0(r0);
        T alf = x / vr0;
        T rho = uni + alf;

        // xp = xp - vk1n * rho * thin_cf_quadrupole_by(x, y, vr0, alf);
        // yp = yp + vk1n * rho * thin_cf_quadrupole_bx(x, y, vr0);

        xp = xp - vk1n * rho * vr0 * log(rho) + vk1s * rho * y;

        yp = yp + vk1n * y + vk1s * vr0 * (rho * rho - uni) * T(0.5);
    }

    // combined function sbends kick up to sextupole
    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    thin_cf_kick_2(T const& x,
                   T& xp,
                   T const& y,
                   T& yp,
                   double r0,
                   double const* kL)
    {
        const T uni(1.0);

        T vk1n(kL[0]);
        T vk1s(kL[1]);

        T vk2n(kL[2]);
        T vk2s(kL[3]);

        T vr0(r0);
        T alf = x / vr0;
        T rho = uni + alf;

        xp = xp - vk1n * (uni + alf) * thin_cf_quadrupole_by(x, y, vr0, alf) +
             vk1s * rho * y -
             vk2n * (uni + alf) * thin_cf_sextupole_by(x, y, vr0, alf) * T(0.5);

        yp = yp + vk1n * (uni + alf) * thin_cf_quadrupole_bx(x, y, vr0) +
             vk1s * vr0 * (rho * rho - uni) * T(0.5) +
             vk2n * (uni + alf) * thin_cf_sextupole_bx(x, y, vr0, alf) * T(0.5);
    }

    // F_rn_i is the F_rho of a normal multipole at order i
    // F_yn_i is the F_ybar (ybar = y/r0) of a normal multipole at order i
    // F_rs_i is the F_rho of a skew multipole at order i
    // F_ys_i is the F_ybar (ybar = y/r0) of a skew multipole at order i
    //
    // See Zolkin, Sector magnets or transverse electromagnetic fields in
    // cylindrical coordinates, Phys.Rev.Accel.Beams 20 (2017) no.4, 043501
    // Table XII
    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    F_rn_2(T const& rho, T const& y)
    {
        return T(-1.0) * (y / rho);
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    F_rn_3(T const& rho, T const& y)
    {
        return T(-0.5) * (y / rho) * (rho * rho - T(1.0));
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    F_rn_4(T const& rho, T const& y)
    {
        return T(-1.0 / 6.0) * (y / rho) *
               (T(-3.0 / 2.0) * (rho * rho - T(1.0)) - y * y +
                T(3.0) * rho * rho * log(rho));
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    F_rn_5(T const& rho, T const& y)
    {
        return T(-1.0 / 12.0) * (y / rho) *
               ((rho * rho - T(1.0)) *
                    (T(3.0 / 4.0) * (rho * rho + T(1.0)) - y * y) -
                T(3.0) * rho * rho * log(rho));
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    F_yn_2(T const& rho, T const& y)
    {
        return log(rho);
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    F_yn_3(T const& rho, T const& y)
    {
        return T(0.5) * (T(0.5) * (rho * rho - T(1.0)) - y * y - log(rho));
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    F_yn_4(T const& rho, T const& y)
    {
        return T(-0.25) * (rho * rho - T(1.0)) +
               T(0.5) * (T(0.5) * (rho * rho + T(1.0)) - y * y) * log(rho);
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    F_yn_5(T const& rho, T const& y)
    {
        return T(0.125) * (T(0.125) * (rho * rho * rho * rho +
                                       T(4.0) * rho * rho - T(5.0)) -
                           (rho * rho - T(1.0)) * y * y + y * y * y * y -
                           (T(0.5) + rho * rho - T(2.0) * y * y) * log(rho));
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    F_rs_2(T const& rho, T const& y)
    {
        return T(0.5) / rho * (rho * rho - T(1.0));
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    F_rs_3(T const& rho, T const& y)
    {
        return T(0.5) / rho *
               (T(-0.5) * (rho * rho - T(1.0)) - y * y + rho * rho * log(rho));
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    F_rs_4(T const& rho, T const& y)
    {
        return T(0.25) / rho *
               (T(0.25) * (rho * rho * rho * rho - T(1.0)) -
                (rho * rho - T(1.0)) * y * y - rho * rho * log(rho));
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    F_rs_5(T const& rho, T const& y)
    {
        return T(0.125) / rho *
               (T(-0.625) * rho * rho * rho * rho + T(0.5) * rho * rho +
                T(0.125) + (rho * rho - T(1.0)) * y * y +
                y * y * y * y / T(3.0) +
                (T(1.0) + T(0.5) * rho * rho - T(2.0) * y * y) * rho * rho *
                    log(rho));
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    F_ys_2(T const& rho, T const& y)
    {
        return y;
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    F_ys_3(T const& rho, T const& y)
    {
        return y * log(rho);
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    F_ys_4(T const& rho, T const& y)
    {
        return T(0.5) * y *
               (T(0.5) * (rho * rho - T(1.0)) - y * y / T(3.0) - log(rho));
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION T
    F_ys_5(T const& rho, T const& y)
    {
        return T(-0.25) * y * (rho * rho - T(1.0)) +
               T(1.0 / 6.0) * y * (T(1.5) * (rho * rho + T(1.0)) - y * y) *
                   log(rho);
    }

    // combined function sbends kick up to dectapole
    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    thin_cf_kick_5(T const& x,
                   T& xp,
                   T const& y,
                   T& yp,
                   double r0,
                   double const* kL)
    {
        const T uni(1.0);

        T vr0(r0);

        T rho = uni + x / vr0;
        T ybar = y / vr0;

        // quadrupole
        if (kL[0] || kL[1]) {
            T kn(kL[0]);
            T ks(kL[1]);

            xp = xp - kn * rho * vr0 * F_yn_2(rho, ybar) +
                 ks * rho * vr0 * F_ys_2(rho, ybar);

            yp = yp - kn * rho * vr0 * F_rn_2(rho, ybar) +
                 ks * rho * vr0 * F_rs_2(rho, ybar);
        }

        // sextupole
        if (kL[2] || kL[3]) {
            T kn(kL[2]);
            T ks(kL[3]);

            xp = xp - kn * rho * vr0 * vr0 * F_yn_3(rho, ybar) +
                 ks * rho * vr0 * vr0 * F_ys_3(rho, ybar);

            yp = yp - kn * rho * vr0 * vr0 * F_rn_3(rho, ybar) +
                 ks * rho * vr0 * vr0 * F_rs_3(rho, ybar);
        }

        // octupole
        if (kL[4] || kL[5]) {
            T kn(kL[4]);
            T ks(kL[5]);

            xp = xp - kn * rho * vr0 * vr0 * vr0 * F_yn_4(rho, ybar) +
                 ks * rho * vr0 * vr0 * vr0 * F_ys_4(rho, ybar);

            yp = yp - kn * rho * vr0 * vr0 * vr0 * F_rn_4(rho, ybar) +
                 ks * rho * vr0 * vr0 * vr0 * F_rs_4(rho, ybar);
        }

        // decapole
        if (kL[6] || kL[7]) {
            T kn(kL[6]);
            T ks(kL[7]);

            xp = xp - kn * rho * vr0 * vr0 * vr0 * vr0 * F_yn_5(rho, ybar) +
                 ks * rho * vr0 * vr0 * vr0 * vr0 * F_ys_5(rho, ybar);

            yp = yp - kn * rho * vr0 * vr0 * vr0 * vr0 * F_rn_5(rho, ybar) +
                 ks * rho * vr0 * vr0 * vr0 * vr0 * F_rs_5(rho, ybar);
        }
    }

    // quadrupole with simple n kicks algorithm. non-Yoshida
    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    quadrupole_chef(T& x,
                    T& xp,
                    T& y,
                    T& yp,
                    T& cdt,
                    T const& dpop,
                    double pref,
                    double m,
                    double refcdt,
                    double len,
                    double* str,
                    int kicks)
    {
        double substep_ref_cdt = refcdt / kicks;
        double frontLength = 6.0 * (len / 4.0) / 15.0;
        double sepLength = (len - 2.0 * frontLength) / 3.0;
        double kl[2] = {str[0] * len / 4.0, str[1] * len / 4.0};

        if (kicks == 1) {
            kl[0] = str[0] * len;
            kl[1] = str[1] * len;

            drift_unit(
                x, xp, y, yp, cdt, dpop, 0.5 * len, pref, m, substep_ref_cdt);
            thin_quadrupole_unit(x, xp, y, yp, kl); // str*len
            drift_unit(
                x, xp, y, yp, cdt, dpop, 0.5 * len, pref, m, substep_ref_cdt);
        } else if ((kicks % 4) == 0) {
            int u = kicks / 4;
            double xu = u;
            frontLength /= xu;
            sepLength /= xu;
            kl[0] /= xu;
            kl[1] /= xu;

            for (int i = 0; i < u; ++i) {
                drift_unit(x,
                           xp,
                           y,
                           yp,
                           cdt,
                           dpop,
                           frontLength,
                           pref,
                           m,
                           substep_ref_cdt);
                thin_quadrupole_unit(x, xp, y, yp, kl); // quaterStrength

                for (int i = 0; i < 3; ++i) {
                    drift_unit(x,
                               xp,
                               y,
                               yp,
                               cdt,
                               dpop,
                               sepLength,
                               pref,
                               m,
                               substep_ref_cdt);
                    thin_quadrupole_unit(x, xp, y, yp, kl); // quaterStrength
                }

                drift_unit(x,
                           xp,
                           y,
                           yp,
                           cdt,
                           dpop,
                           frontLength,
                           pref,
                           m,
                           substep_ref_cdt);
            }
        } else {
            kl[0] = str[0] * len / kicks;
            kl[1] = str[1] * len / kicks;

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       len / (2.0 * kicks),
                       pref,
                       m,
                       substep_ref_cdt);
            thin_quadrupole_unit(x, xp, y, yp, kl); // str*len/kicks

            for (int i = 0; i < kicks - 1; ++i) {
                drift_unit(x,
                           xp,
                           y,
                           yp,
                           cdt,
                           dpop,
                           len / kicks,
                           pref,
                           m,
                           substep_ref_cdt);
                thin_quadrupole_unit(x, xp, y, yp, kl); // str*len/kicks
            }

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       len / (2.0 * kicks),
                       pref,
                       m,
                       substep_ref_cdt);
        }
    }

    // thin kick dipole
    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    thin_dipole_unit(T const& x, T& xp, T const& y, T& yp, double const* kL)
    {
        xp = xp - T(kL[0]);
        yp = yp + T(kL[1]);
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    thin_quadrupole_unit(T const& x, T& xp, T const& y, T& yp, double const* kL)
    {
        T vk0(kL[0]);
        T vk1(kL[1]);

        xp = xp - vk0 * x + vk1 * y;
        yp = yp + vk0 * y + vk1 * x;
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    thin_sextupole_unit(T const& x, T& xp, T const& y, T& yp, double const* kL)
    {
        T vk0(kL[0]);
        T vk1(kL[1]);

        xp = xp - T(0.5) * vk0 * (x * x - y * y) + vk1 * x * y;
        yp = yp + vk0 * x * y + T(0.5) * vk1 * (x * x - y * y);
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    thin_octupole_unit(T const& x, T& xp, T const& y, T& yp, double const* kL)
    {
        T vk0(kL[0]);
        T vk1(kL[1]);

        T n1(0.5);
        T n2(1.0 / 3.0);

        xp = xp - n1 * vk0 * (x * x * x * n2 - x * y * y) +
             n1 * vk1 * (x * x * y - y * y * y * n2);
        yp = yp - n1 * vk0 * (y * y * y * n2 - x * x * y) +
             n1 * vk1 * (x * x * x * n2 - x * y * y);
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    thin_rbend_unit(T const& x, T& xp, T const& y, T& yp, double const* kL)
    {
        thin_dipole_unit(x, xp, y, yp, kL);
        thin_quadrupole_unit(x, xp, y, yp, kL + 2);
        thin_sextupole_unit(x, xp, y, yp, 0.0, kL + 4);
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    rbend_thin_cf_kick(T const& x,
                       T& xp,
                       T const& y,
                       T& yp,
                       double r0,
                       double const* kL)
    {
        thin_quadrupole_unit(x, xp, y, yp, kL + 0);
        thin_sextupole_unit(x, xp, y, yp, 0.0, kL + 2);
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    thin_kicker_unit(T& p, double kL)
    {
        p = p + T(kL);
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    thin_kicker_unit(T& xp, T& yp, double const* kL)
    {
        xp = xp + T(kL[0]);
        yp = yp + T(kL[1]);
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    constfoc_unit(T& x, T& xp, double cs, double sn, double beta, double ibeta)
    {
        T vcs(cs), vsn(sn), vbeta(beta), vibeta(ibeta);

        T x1 = x * vcs + xp * vbeta * vsn;
        T xp1 = -x * vsn * vibeta + xp * vcs;

        x = x1;
        xp = xp1;
    }

    KOKKOS_INLINE_FUNCTION
    double
    thin_rfcavity_pnew(double pref, double m, double volt, double phi_s)
    {
        double E0 = sqrt(pref * pref + m * m);
        double E1 = E0 + volt * sin(phi_s);
        return sqrt(E1 * E1 - m * m);
    }

    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    thin_rfcavity_unit(T& px,
                       T& py,
                       T const& cdt,
                       T& dpop,
                       double w_rf,
                       double volt,
                       double phi_s,
                       double m,
                       double old_ref_p,
                       double new_ref_p,
                       double const* mhp,
                       int nh)
    {
        // 0 strength RF cavity does nothing (fast)
        if (volt == 0.0) { return; }
        double const anh_phase = 0.0;
        double mhp_harmonic_multiple = 1.0;
        double mhp_relative_strength = 1.0;
        double mhp_phase_shift = 0.0;

        T phase_slip_argument = (cdt * T(w_rf / pconstants::c)) + T(anh_phase);

        T strength_factor(0.0);

        for (int i = 0; i < nh; ++i) {
            mhp_harmonic_multiple = mhp[i * 3];
            mhp_relative_strength = mhp[i * 3 + 1];
            mhp_phase_shift = mhp[i * 3 + 2];

            strength_factor =
                strength_factor + T(mhp_relative_strength) *
                                      sin(T(mhp_harmonic_multiple) *
                                              (T(phi_s) + phase_slip_argument) +
                                          T(mhp_phase_shift));
        }

        T p = T(old_ref_p) * (dpop + T(1.0));
        T E = sqrt(p * p + T(m * m));

        E = E + T(volt) * strength_factor;

        px = px * T(old_ref_p / new_ref_p);
        py = py * T(old_ref_p / new_ref_p);

        dpop = sqrt((E - T(m)) * (E + T(m))) / T(new_ref_p) - T(1.0);
    }

    // adjust particle coordinates to use new reference energy
    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    adjust_ref_unit(T& px,
                       T& py,
                       T const& cdt,
                       T& dpop,
                       double m,
                       double old_pref,
                       double new_pref,
                       double old_E,
                       double new_E)
    {
    
        T p = T(old_pref) * (dpop + T(1.0));
        T E = sqrt(p * p + T(m * m));


        px = px * T(old_pref / new_pref);
        py = py * T(old_pref / new_pref);

        dpop = sqrt((E - T(m)) * (E + T(m))) / T(new_pref) - T(1.0);

    }
        
    // adjust particle coordinates to use new reference energy
    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    adjust_cdt_unit(T& p4, double offset)
    {
            p4 = p4 + T(offset);
    }
        
    KOKKOS_INLINE_FUNCTION
    double
    factorial(int n)
    {
        if (n == 0) return 1.0;

        double r = 1;
        for (int i = 1; i <= n; ++i)
            r *= i;
        return r;
    }

    // general thin magnets of n-th order
    // n = 1, dipole; n = 2, quadrupole; n = 3, sextupole, etc.
    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    thin_magnet_unit(T const& x,
                     T& xp,
                     T const& y,
                     T& yp,
                     double const* kL,
                     int n)
    {
        for (int k = 0; k < n; k += 4) {
            xp = xp - T(kL[0] * (n - k)) * qpow(x, n - k - 1) * qpow(y, k) /
                          T(factorial(n - k) * factorial(k));
        }

        for (int k = 2; k < n; k += 4) {
            xp = xp + T(kL[0] * (n - k)) * qpow(x, n - k - 1) * qpow(y, k) /
                          T(factorial(n - k) * factorial(k));
        }

        for (int k = 1; k < n; k += 4) {
            xp = xp + T(kL[1] * (n - k)) * qpow(x, n - k - 1) * qpow(y, k) /
                          T(factorial(n - k) * factorial(k));
        }

        for (int k = 3; k < n; k += 4) {
            xp = xp - T(kL[1] * (n - k)) * qpow(x, n - k - 1) * qpow(y, k) /
                          T(factorial(n - k) * factorial(k));
        }

        for (int k = 4; k <= n; k += 4) {
            yp = yp - T(kL[0] * k) * qpow(x, n - k) * qpow(y, k - 1) /
                          T(factorial(n - k) * factorial(k));
        }

        for (int k = 2; k <= n; k += 4) {
            yp = yp + T(kL[0] * k) * qpow(x, n - k) * qpow(y, k - 1) /
                          T(factorial(n - k) * factorial(k));
        }

        for (int k = 1; k <= n; k += 4) {
            yp = yp + T(kL[1] * k) * qpow(x, n - k) * qpow(y, k - 1) /
                          T(factorial(n - k) * factorial(k));
        }

        for (int k = 3; k <= n; k += 4) {
            yp = yp - T(kL[1] * k) * qpow(x, n - k) * qpow(y, k - 1) /
                          T(factorial(n - k) * factorial(k));
        }
    }

#if 0
    // combined thin multipole kicks of n-th order
    // n = 1, dipole; n = 2, quadrupole; n = 3, sextupole, etc.
    template <typename T>
    KOKKOS_INLINE_FUNCTION
    void thin_multipole_kick
      (T const& x, T& xp, T const& y, T& yp, double const * kL, int n)
#endif

    template <class T>
    KOKKOS_INLINE_FUNCTION void
    nllens_unit(T const& x, T& xp, T const& y, T& yp, T const&, double const* k)
    {
        const T icnll = T(k[0]);
        const T kick = T(k[1]);

        T xbar = x * icnll;
        T ybar = y * icnll;

        if (ybar == 0 && (xbar >= 1.0 || xbar <= -1.0)) {
            xp = quiet_nan();
            yp = quiet_nan();
            return;
        }

#if 0
        Kokkos::complex<double> c_i(0.0, 1.0);
        Kokkos::complex<double> c_1(1.0, 0.0);

        Kokkos::complex<double> zeta(xbar, ybar);
        Kokkos::complex<double> croot = sqrt(c_1 - zeta*zeta);

        // Kokkos::complex doesnt provide the log() method, so we
        // need to calculate the following manually:
        // Kokkos::complex<double> carcsin = -c_i * log(c_i * zeta + croot);

        // convert co to polar form
        Kokkos::complex<double> co = c_i * zeta + croot;
        double a = co.real();
        double b = co.imag();
        double r = sqrt(a*a + b*b);
        double theta = (a>0) ? atan(b/a) : atan(b/a) + mconstants::pi;
        Kokkos::complex<double> logco(log(r), theta);
        Kokkos::complex<double> carcsin = -c_i * logco;

        // keep working on the rest
        Kokkos::complex<double> dF = zeta/(croot*croot) + carcsin/(croot*croot*croot);

        T dpx = kick * dF.real();
        T dpy = -kick * dF.imag();

        xp = xp + dpx;
        yp = yp + dpy;
#endif

        T zeta_r = xbar;
        T zeta_i = ybar;

        // c = c_1 - zeta*zeta
        T c_r = T(1.0) - zeta_r * zeta_r + zeta_i * zeta_i;
        T c_i = -T(2.0) * zeta_r * zeta_i;

        // croot = sqrt(c);
        T croot_r = sqrt(T(0.5) * (c_r + sqrt(c_r * c_r + c_i * c_i)));
        T croot_i = sqrt(T(0.5) * (-c_r + sqrt(c_r * c_r + c_i * c_i)));

        // co = c_i * zeta + croot
        T co_r = croot_r - zeta_i;
        T co_i = croot_i + zeta_r;

        T r = sqrt(co_r * co_r + co_i * co_i);

        T theta = (co_r > 0) ?
                      atan(co_i / co_r) + T(0.0) :
                      atan(co_i / co_r) + T(Kokkos::numbers::pi_v<double>);

        // log(co)
        T logco_r = log(r);
        T logco_i = theta;

        // carcsin = -c_i * log(co)
        T carcsin_r = logco_i;
        T carcsin_i = -logco_r;

        // cr2 = croot * croot
        T cr2_r = croot_r * croot_r - croot_i * croot_i;
        T cr2_i = T(2.0) * croot_r * croot_i;

        // df1 = zeta/cr2
        T df1_r =
            (cr2_r * zeta_r + cr2_i * zeta_i) / (cr2_r * cr2_r + cr2_i * cr2_i);
        T df1_i =
            (cr2_r * zeta_i - cr2_i * zeta_r) / (cr2_r * cr2_r + cr2_i * cr2_i);

        // df2 = carcsin/cr2
        T df2_r = (cr2_r * carcsin_r + cr2_i * carcsin_i) /
                  (cr2_r * cr2_r + cr2_i * cr2_i);
        T df2_i = (cr2_r * carcsin_i - cr2_i * carcsin_r) /
                  (cr2_r * cr2_r + cr2_i * cr2_i);

        // df = df1 + df2
        T df_r = df1_r + df2_r;
        T df_i = df1_i + df2_i;

        T dpx = kick * df_r;
        T dpy = -kick * df_i;

        xp = xp + dpx;
        yp = yp + dpy;
    }

    // Calculations of dipedge taken from MAD-X SUBROUTINE tmfrng from
    // MAD-X source file twiss.f90.

    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    dipedge_unit(T& x,
                 T& xp,
                 T& y,
                 T& yp,
                 double re_2_1,
                 double re_4_3,
                 double const* te)
    {
	    T newx(x);
	    T newxp(xp);
	    T newy(y);
	    T newyp(yp);
	    
        // linear terms
        newxp = xp + T(re_2_1) * x;
        newyp = yp + T(re_4_3) * y;

        // quadratic terms
        newx = x + T(te[0]) * x * x + T(te[1]) * y * y;
        newxp = newxp + T(te[2]) * x * x + T(te[4]) * y * y +
             T(te[3]) * x * xp * T(2) + T(te[5]) * y * yp * T(2);
        newy = y + T(te[6]) * x * y * T(2);
        newyp = newyp + T(te[7]) * x * y * T(2) + T(te[8]) * x * yp * T(2) +
             T(te[9]) * y * xp * T(2);

	x = newx;
	xp = newxp;
	y = newy;
	yp = newyp;
    }

    template <class T>
    KOKKOS_INLINE_FUNCTION void
    solenoid_unit(T& x,
                  T& xp,
                  T& y,
                  T& yp,
                  T& cdt,
                  T const& dpop,
                  double ks,
                  double ksl,
                  double length,
                  double ref_p,
                  double mass,
                  double ref_cdt)
    {
        T uni(1.0);

        T p2 = (uni + dpop) * (uni + dpop);
        T zp = sqrt(p2 - xp * xp - yp * yp);
        T dtheta = T(ksl) / zp;

        T sn = sin(dtheta);
        T cs = cos(dtheta);

        T xpi = xp;
        T ypi = yp;

        xp = cs * xpi + sn * ypi;
        yp = (-sn) * xpi + cs * ypi;

        cs = cs - uni;

        x = x + (cs * (-ypi) + sn * xpi) / T(ks);
        y = y + ((-sn) * (-ypi) + cs * xpi) / T(ks);

        T en = sqrt(p2 * T(ref_p) * T(ref_p) + T(mass) * T(mass));
        T duration = T(length) / (zp * T(ref_p) / en);

        cdt = cdt + duration - T(ref_cdt);
    }

    template <class T>
    KOKKOS_INLINE_FUNCTION void
    solenoid_in_edge_kick(T const& x, T& xp, T const& y, T& yp, double kse)
    {
        xp = xp + T(kse) * y;
        yp = yp - T(kse) * x;
    }

    template <class T>
    KOKKOS_INLINE_FUNCTION void
    solenoid_out_edge_kick(T const& x, T& xp, T const& y, T& yp, double kse)
    {
        xp = xp - T(kse) * y;
        yp = yp + T(kse) * x;
    }

    // k[6]:
    // k[0] = beta_b, k[1] = gamma_b, k[2] = beta_e,
    // k[3] = ioe,    k[4] = l,       k[5] = radius
    template <class T>
    KOKKOS_INLINE_FUNCTION void
    elens_kick_gaussian(T const& x,
                        T& xp,
                        T const& y,
                        T& yp,
                        T const& dpop,
                        double const* k)
    {
        const double small_radius(1.0e-10);
        T r = sqrt(x * x + y * y);

        // no kick at r = 0.0
        if (r == 0.0) return;

        const T t1(1.0);
        const T t2(2.0);
        const T t4(4.0);
        const T t6(6.0);
        const T t8(8.0);

        const T beta_b = k[0];
        const T gamma_b = k[1];
        const T beta_e = k[2];
        const T ioe = k[3];
        const T l = k[4];
        const T radius = k[5];

        T betagamma_p = (t1 + dpop) * beta_b * gamma_b;
        T beta_p = betagamma_p / sqrt(t1 + betagamma_p * betagamma_p);

        T factors = -t2 * ioe * l * T(pconstants::rp) * (t1 + beta_e * beta_p) /
                    (beta_e * beta_p * beta_b * gamma_b * T(pconstants::c));

        T kick(0.0);

        if (r < small_radius) {
            kick = factors * (r / (t2 * radius * radius) -
                              (t1 / t2) * ((r * r * r) / (t4 * radius * radius *
                                                          radius * radius)) +
                              (t1 / t6) * ((r * r * r * r * r) /
                                           (t8 * radius * radius * radius *
                                            radius * radius * radius)));
        } else {
            kick = factors * (t1 - exp(-r * r / (t2 * radius * radius))) / r;
        }

        xp = xp + kick * x / r;
        yp = yp + kick * y / r;
    }

    template <class T>
    KOKKOS_INLINE_FUNCTION void
    elens_kick_uniform(T const& x,
                       T& xp,
                       T const& y,
                       T& yp,
                       T const& dpop,
                       double const* k)
    {
        // const double small_radius = 1.0e-10;
        T r = sqrt(x * x + y * y);

        // no kick at r = 0.0
        if (r == 0.0) return;

        const T beta_b = k[0];
        const T gamma_b = k[1];
        const T beta_e = k[2];
        const T ioe = k[3];
        const T l = k[4];
        const T radius = k[5];

        T betagamma_p = (T(1.0) + dpop) * beta_b * gamma_b;
        T beta_p = betagamma_p / sqrt(T(1.0) + betagamma_p * betagamma_p);

        T factors = -T(2.0) * ioe * l * T(pconstants::rp) *
                    (T(1.0) + beta_e * beta_p) /
                    (beta_e * beta_p * beta_b * gamma_b * T(pconstants::c));

        T kick = factors * r / (radius * radius);

        xp = xp + kick * x / r;
        yp = yp + kick * y / r;
    }

    // utility
    KOKKOS_INLINE_FUNCTION
    int
    full_drifts_per_step(int order)
    {
        return std::pow(3.0, (order - 2.0) / 2.0) * 2;
    }

    KOKKOS_INLINE_FUNCTION
    int
    compact_drifts_per_step(int order)
    {
        return (full_drifts_per_step(order) - 2) / 2 + 2;
    }

    // general n-th order yoshida
    template <class T>
    using yoshida_kf_t =
        void (*)(T const&, T&, T const&, T&, T const&, double const*);

    template <class T>
    using bend_yoshida_kf_t =
        void (*)(T const&, T&, T const&, T&, double, double const*);

    template <typename T, yoshida_kf_t<T> kf, int n, int components>
    struct yoshida_element {
        KOKKOS_INLINE_FUNCTION
        static void
        integral(T& x,
                 T& xp,
                 T& y,
                 T& yp,
                 T& cdt,
                 T const& dpop,
                 double pref,
                 double m,
                 double step_ref_cdt,
                 double step_length,
                 double const* step_strength,
                 int steps,
                 double c)
        {
            double substep_ref_cdt = step_ref_cdt / 3.0;

            yoshida_element<T, kf, n - 1, components>::integral(x,
                                                                xp,
                                                                y,
                                                                yp,
                                                                cdt,
                                                                dpop,
                                                                pref,
                                                                m,
                                                                substep_ref_cdt,
                                                                step_length,
                                                                step_strength,
                                                                steps,
                                                                c * x1(n));

            yoshida_element<T, kf, n - 1, components>::integral(x,
                                                                xp,
                                                                y,
                                                                yp,
                                                                cdt,
                                                                dpop,
                                                                pref,
                                                                m,
                                                                substep_ref_cdt,
                                                                step_length,
                                                                step_strength,
                                                                steps,
                                                                c * x0(n));

            yoshida_element<T, kf, n - 1, components>::integral(x,
                                                                xp,
                                                                y,
                                                                yp,
                                                                cdt,
                                                                dpop,
                                                                pref,
                                                                m,
                                                                substep_ref_cdt,
                                                                step_length,
                                                                step_strength,
                                                                steps,
                                                                c * x1(n));
        }

        KOKKOS_INLINE_FUNCTION
        static double
        x1(int nn)
        {
            return 1.0 / (2.0 - std::pow(2.0, 1.0 / (2 * nn + 1)));
        }

        KOKKOS_INLINE_FUNCTION
        static double
        x0(int nn)
        {
            return 1.0 - 2.0 * x1(nn);
        }
    };

#if 1
    template <typename T, yoshida_kf_t<T> kf, int components>
    struct yoshida_element<T, kf, 0, components> {
        KOKKOS_INLINE_FUNCTION
        static void
        integral(T& x,
                 T& xp,
                 T& y,
                 T& yp,
                 T& cdt,
                 T const& dpop,
                 double pref,
                 double m,
                 double step_ref_cdt,
                 double step_length,
                 double const* step_strength,
                 int steps,
                 double c)
        {
            double substep_ref_cdt = step_ref_cdt / 2.0;
            double kl[components * 2];

            for (int i = 0; i < components * 2; ++i)
                kl[i] = step_strength[i] * c;

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       0.5 * c * step_length,
                       pref,
                       m,
                       substep_ref_cdt);

            kf(x, xp, y, yp, cdt, kl);

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       0.5 * c * step_length,
                       pref,
                       m,
                       substep_ref_cdt);
        }
    };

#else

    template <typename T, yoshida_kf_t<T> kf, int components>
    KOKKOS_INLINE_FUNCTION struct yoshida_element<T, kf, 0, components> {
        KOKKOS_INLINE_FUNCTION
        static void
        integral(T& x,
                 T& xp,
                 T& y,
                 T& yp,
                 T& cdt,
                 T const& dpop,
                 double pref,
                 double m,
                 double step_ref_cdt,
                 double step_length,
                 double* step_strength,
                 int steps,
                 double c)
        {
            double substep_ref_cdt = step_ref_cdt;
            double kl[components * 2];

            for (int i = 0; i < components * 2; ++i)
                kl[i] = step_strength[i] * c * 0.5;

            kf(x, xp, y, yp, kl);

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       c * step_length,
                       pref,
                       m,
                       substep_ref_cdt);

            kf(x, xp, y, yp, kl);
        }
    };
#endif

    template <typename T, yoshida_kf_t<T> kf, int order, int components>
    KOKKOS_INLINE_FUNCTION void
    yoshida(T& x,
            T& xp,
            T& y,
            T& yp,
            T& cdt,
            T const& dpop,
            double pref,
            double m,
            double step_ref_cdt,
            double step_length,
            double const* step_strength,
            int steps)
    {
        const int n = (order - 2) / 2;

        for (int i = 0; i < steps; ++i) {
            yoshida_element<T, kf, n, components>::integral(x,
                                                            xp,
                                                            y,
                                                            yp,
                                                            cdt,
                                                            dpop,
                                                            pref,
                                                            m,
                                                            step_ref_cdt,
                                                            step_length,
                                                            step_strength,
                                                            steps,
                                                            1.0);
        }
    }

    // hardwired 2nd order yoshida
    template <typename T, yoshida_kf_t<T> kf, int components>
    KOKKOS_INLINE_FUNCTION void
    yoshida2(T& x,
             T& xp,
             T& y,
             T& yp,
             T& cdt,
             T const& dpop,
             double reference_momentum,
             double m,
             double substep_reference_cdt,
             double step_length,
             double const* step_strength,
             int steps)
    {
        for (int i = 0; i < steps; ++i) {
            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       0.5 * step_length,
                       reference_momentum,
                       m,
                       substep_reference_cdt);

            kf(x, xp, y, yp, cdt, step_strength);

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       0.5 * step_length,
                       reference_momentum,
                       m,
                       substep_reference_cdt);
        }
    }

    // hardwired 4th order yoshida
    template <typename T, yoshida_kf_t<T> kf, int components>
    KOKKOS_INLINE_FUNCTION void
    yoshida4(T& x,
             T& xp,
             T& y,
             T& yp,
             T& cdt,
             T const& dpop,
             double reference_momentum,
             double m,
             double step_reference_cdt,
             double step_length,
             double const* step_strength,
             int steps)
    {
        // see yoshida4.py for formulas
        const double c1 = 0.675603595979828817023843904487;
        const double c4 = c1;
        const double c2 = -0.175603595979828817023843904487;
        const double c3 = c2;
        const double d1 = 1.35120719195965763404768780897;
        const double d3 = d1;
        const double d2 = -1.70241438391931526809537561795;

        double k1[components * 2], k2[components * 2], k3[components * 2];

        double substep_reference_cdt = step_reference_cdt * 0.25;

        for (int i = 0; i < components; ++i) {
            k1[i * 2 + 0] = d1 * step_strength[i * 2 + 0];
            k1[i * 2 + 1] = d1 * step_strength[i * 2 + 1];

            k2[i * 2 + 0] = d2 * step_strength[i * 2 + 0];
            k2[i * 2 + 1] = d2 * step_strength[i * 2 + 1];

            k3[i * 2 + 0] = d3 * step_strength[i * 2 + 0];
            k3[i * 2 + 1] = d3 * step_strength[i * 2 + 1];
        }

        for (int i = 0; i < steps; ++i) {
            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       c1 * step_length,
                       reference_momentum,
                       m,
                       substep_reference_cdt);

            // kf( x, xp, y, yp, d1 * step_strength );
            kf(x, xp, y, yp, cdt, k1);

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       c2 * step_length,
                       reference_momentum,
                       m,
                       substep_reference_cdt);

            // kf( x, xp, y, yp, d2 * step_strength );
            kf(x, xp, y, yp, cdt, k2);

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       c3 * step_length,
                       reference_momentum,
                       m,
                       substep_reference_cdt);

            // kf( x, xp, y, yp, d3 * step_strength );
            kf(x, xp, y, yp, cdt, k3);

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       c4 * step_length,
                       reference_momentum,
                       m,
                       substep_reference_cdt);
        }
    }

    // hardwired 6th order yoshida
    template <typename T, yoshida_kf_t<T> kf, int components>
    KOKKOS_INLINE_FUNCTION void
    yoshida6(T& x,
             T& xp,
             T& y,
             T& yp,
             T& cdt,
             T const& dpop,
             double reference_momentum,
             double m,
             double step_reference_cdt,
             double step_length,
             double const* step_strength,
             int steps)
    {
        // see yoshida4.py for formulas
        const double c1 = 0.79361246386112147294625603763;
        const double c2 = -0.206276584816439780698721319354;
        const double c3 = -0.118008867881292655922412133144;
        const double c4 = 0.236949573653050744373598734219;

        const double d1 = 1.58722492772224294589251207526;
        const double d2 = -1.99977809735512250728995471397;
        const double d3 = -1.82324266348482825773733634154;
        const double d4 = 2.29714181079092974648453380998;

        // 10 drifts per step
        double substep_reference_cdt = step_reference_cdt * 0.1;

        double k1[components * 2], k2[components * 2], k3[components * 2],
            k4[components * 2];

        for (int i = 0; i < components; ++i) {
            k1[i * 2 + 0] = d1 * step_strength[i * 2 + 0];
            k1[i * 2 + 1] = d1 * step_strength[i * 2 + 1];

            k2[i * 2 + 0] = d2 * step_strength[i * 2 + 0];
            k2[i * 2 + 1] = d2 * step_strength[i * 2 + 1];

            k3[i * 2 + 0] = d3 * step_strength[i * 2 + 0];
            k3[i * 2 + 1] = d3 * step_strength[i * 2 + 1];

            k4[i * 2 + 0] = d4 * step_strength[i * 2 + 0];
            k4[i * 2 + 1] = d4 * step_strength[i * 2 + 1];
        }

        for (int i = 0; i < steps; ++i) {
            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       c1 * step_length,
                       reference_momentum,
                       m,
                       substep_reference_cdt);

            kf(x, xp, y, yp, cdt, k1);

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       c2 * step_length,
                       reference_momentum,
                       m,
                       substep_reference_cdt);

            kf(x, xp, y, yp, cdt, k2);

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       c2 * step_length,
                       reference_momentum,
                       m,
                       substep_reference_cdt);

            kf(x, xp, y, yp, cdt, k1);

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       c3 * step_length,
                       reference_momentum,
                       m,
                       substep_reference_cdt);

            kf(x, xp, y, yp, cdt, k3);

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       c4 * step_length,
                       reference_momentum,
                       m,
                       substep_reference_cdt);

            kf(x, xp, y, yp, cdt, k4);

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       c4 * step_length,
                       reference_momentum,
                       m,
                       substep_reference_cdt);

            kf(x, xp, y, yp, cdt, k3);

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       c3 * step_length,
                       reference_momentum,
                       m,
                       substep_reference_cdt);

            kf(x, xp, y, yp, cdt, k1);

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       c2 * step_length,
                       reference_momentum,
                       m,
                       substep_reference_cdt);

            kf(x, xp, y, yp, cdt, k2);

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       c2 * step_length,
                       reference_momentum,
                       m,
                       substep_reference_cdt);

            kf(x, xp, y, yp, cdt, k1);

            drift_unit(x,
                       xp,
                       y,
                       yp,
                       cdt,
                       dpop,
                       c1 * step_length,
                       reference_momentum,
                       m,
                       substep_reference_cdt);
        }
    }

    namespace yoshida_constants {
        // see yoshida4.py for formulas
        namespace bend_yoshida6 {
            const double c1 = 0.79361246386112147294625603763;
            const double c2 = -0.206276584816439780698721319354;
            const double c3 = -0.118008867881292655922412133144;
            const double c4 = 0.236949573653050744373598734219;

            const double d1 = 1.58722492772224294589251207526;
            const double d2 = -1.99977809735512250728995471397;
            const double d3 = -1.82324266348482825773733634154;
            const double d4 = 2.29714181079092974648453380998;
        }
    }

    template <typename T,
              void(kf)(T const& x,
                       T& xp,
                       T const& y,
                       T& yp,
                       double r0,
                       double const* kL),
              int components>
    KOKKOS_INLINE_FUNCTION void
    bend_yoshida4(T& x,
                  T& xp,
                  T& y,
                  T& yp,
                  T& cdt,
                  T const& dpop,
                  double reference_momentum,
                  double m,
                  double step_reference_cdt,
                  double step_angle,
                  double const* step_strength,
                  double r0,
                  double bend_strength,
                  int steps)
    {
#if 0
        // see yoshida4.py for formulas
        const double c1 = 0.675603595979828817023843904487;
        const double c4 = c1;
        const double c2 = -0.175603595979828817023843904487;
        const double c3 = c2;
        const double d1 = 1.35120719195965763404768780897;
        const double d3 = d1;
        const double d2 = -1.70241438391931526809537561795;

        double k1[components * 2], k2[components * 2], k3[components * 2];

        double substep_reference_cdt = step_reference_cdt * 0.25;

        for (int i=0; i<components; ++i)
        {
            k1[i*2+0] = d1 * step_strength[i*2+0];
            k1[i*2+1] = d1 * step_strength[i*2+1];

            k2[i*2+0] = d2 * step_strength[i*2+0];
            k2[i*2+1] = d2 * step_strength[i*2+1];

            k3[i*2+0] = d3 * step_strength[i*2+0];
            k3[i*2+1] = d3 * step_strength[i*2+1];
        }

        static double ssd = 0.0;
        static double sr0 = 0.0;

        static std::complex<double> step_phase[0] = bend_unit_phase(c1 * ssd);
        static std::complex<double> step_phase[1] = bend_unit_phase(c2 * ssd);
        static std::complex<double> step_phase[2] = bend_unit_phase(c3 * ssd);
        static std::complex<double> step_phase[3] = bend_unit_phase(c4 * ssd);

        static std::complex<double> step_term[0] = bend_unit_term(sr0, c1 * ssd);
        static std::complex<double> step_term[1] = bend_unit_term(sr0, c2 * ssd);
        static std::complex<double> step_term[2] = bend_unit_term(sr0, c3 * ssd);
        static std::complex<double> step_term[3] = bend_unit_term(sr0, c4 * ssd);

        double theta = step_angle;

        if (theta != ssd)
        {
            // updates both phase and term
            step_phase[0] = bend_unit_phase(c1 * theta);
            step_phase[1] = bend_unit_phase(c2 * theta);
            step_phase[2] = bend_unit_phase(c3 * theta);
            step_phase[3] = bend_unit_phase(c4 * theta);

            step_term[0] = bend_unit_term(r0, c1 * theta);
            step_term[1] = bend_unit_term(r0, c2 * theta);
            step_term[2] = bend_unit_term(r0, c3 * theta);
            step_term[3] = bend_unit_term(r0, c4 * theta);

            ssd = theta;
            sr0 = r0;
        }
        else if (r0 != sr0)
        {
            // update term only
            step_term[0] = bend_unit_term(r0, c1 * theta);
            step_term[1] = bend_unit_term(r0, c2 * theta);
            step_term[2] = bend_unit_term(r0, c3 * theta);
            step_term[3] = bend_unit_term(r0, c4 * theta);

            sr0 = r0;
        }

        for(int i = 0; i < steps; ++i)
        {
            bend_unit(x, xp, y, yp, cdt, dpop, - c1 * step_angle, bend_strength, reference_momentum,
                       m, substep_reference_cdt, step_phase[0], step_term[0]);

            kf( x, xp, y, yp, r0, k1 );

            bend_unit(x, xp, y, yp, cdt, dpop, - c2 * step_angle, bend_strength, reference_momentum,
                       m, substep_reference_cdt, step_phase[1], step_term[1]);

            kf( x, xp, y, yp, r0, k2 );

            bend_unit(x, xp, y, yp, cdt, dpop, - c3 * step_angle, bend_strength, reference_momentum,
                       m, substep_reference_cdt, step_phase[2], step_term[2]);

            kf( x, xp, y, yp, r0, k3 );

            bend_unit(x, xp, y, yp, cdt, dpop, - c4 * step_angle, bend_strength, reference_momentum,
                       m, substep_reference_cdt, step_phase[3], step_term[3]);
        }
#endif
    }

    KOKKOS_INLINE_FUNCTION
    Kokkos::complex<double>
    sbend_unit_phase(double c, double delta)
    {
        return bend_unit_phase(c * delta);
    }

    KOKKOS_INLINE_FUNCTION
    Kokkos::complex<double>
    rbend_unit_phase(double c, double delta)
    {
        return bend_unit_phase(0.0);
    }

    KOKKOS_INLINE_FUNCTION
    Kokkos::complex<double>
    sbend_unit_term(double c, double delta, double r0)
    {
        return bend_unit_term(r0, c * delta);
    }

    KOKKOS_INLINE_FUNCTION
    Kokkos::complex<double>
    rbend_unit_term(double c, double delta, double r0)
    {
        return c * delta;
    }

    KOKKOS_INLINE_FUNCTION
    double
    sbend_dphi(double c, double delta)
    {
        return -c * delta;
    }

    KOKKOS_INLINE_FUNCTION
    double
    rbend_dphi(double c, double delta)
    {
        return 0.0;
    }

#if 1
    // hardwired 6th order yoshida
    template <typename T,
              void(f_kf)(T const& x,
                         T& xp,
                         T const& y,
                         T& yp,
                         double r0,
                         double const* kL), // kick
              int components>
    KOKKOS_INLINE_FUNCTION void
    bend_yoshida6(T& x,
                  T& xp,
                  T& y,
                  T& yp,
                  T& cdt,
                  T const& dpop,
                  double reference_momentum,
                  double m,
                  double step_reference_cdt,
                  double const* step_strength,
                  double const* dphi,
                  Kokkos::complex<double> const* step_phase,
                  Kokkos::complex<double> const* step_term,
                  double r0,
                  double bend_strength,
                  int steps)
    {
        // see yoshida4.py for formulas
        namespace by6 = yoshida_constants::bend_yoshida6;

        // 10 drifts per step
        double substep_reference_cdt = step_reference_cdt * 0.1;

        double k1[components * 2], k2[components * 2], k3[components * 2],
            k4[components * 2];

        for (int i = 0; i < components; ++i) {
            k1[i * 2 + 0] = by6::d1 * step_strength[i * 2 + 0];
            k1[i * 2 + 1] = by6::d1 * step_strength[i * 2 + 1];

            k2[i * 2 + 0] = by6::d2 * step_strength[i * 2 + 0];
            k2[i * 2 + 1] = by6::d2 * step_strength[i * 2 + 1];

            k3[i * 2 + 0] = by6::d3 * step_strength[i * 2 + 0];
            k3[i * 2 + 1] = by6::d3 * step_strength[i * 2 + 1];

            k4[i * 2 + 0] = by6::d4 * step_strength[i * 2 + 0];
            k4[i * 2 + 1] = by6::d4 * step_strength[i * 2 + 1];
        }

        for (int i = 0; i < steps; ++i) {
            // drift_unit(x, xp, y, yp, cdt, dpop, c1 * step_length,
            // reference_momentum,
            //           m, substep_reference_cdt);

            bend_unit<T>(x,
                         xp,
                         y,
                         yp,
                         cdt,
                         dpop,
                         dphi[0],
                         bend_strength,
                         reference_momentum,
                         m,
                         substep_reference_cdt,
                         step_phase[0],
                         step_term[0]);

            f_kf(x, xp, y, yp, r0, k1);

            // drift_unit(x, xp, y, yp, cdt, dpop, c2 * step_length,
            // reference_momentum,
            //            m, substep_reference_cdt);

            bend_unit<T>(x,
                         xp,
                         y,
                         yp,
                         cdt,
                         dpop,
                         dphi[1],
                         bend_strength,
                         reference_momentum,
                         m,
                         substep_reference_cdt,
                         step_phase[1],
                         step_term[1]);

            f_kf(x, xp, y, yp, r0, k2);

            // drift_unit(x, xp, y, yp, cdt, dpop, c2 * step_length,
            // reference_momentum,
            //            m, substep_reference_cdt);

            bend_unit<T>(x,
                         xp,
                         y,
                         yp,
                         cdt,
                         dpop,
                         dphi[1],
                         bend_strength,
                         reference_momentum,
                         m,
                         substep_reference_cdt,
                         step_phase[1],
                         step_term[1]);

            f_kf(x, xp, y, yp, r0, k1);

            // drift_unit(x, xp, y, yp, cdt, dpop, c3 * step_length,
            // reference_momentum,
            //            m, substep_reference_cdt);

            bend_unit<T>(x,
                         xp,
                         y,
                         yp,
                         cdt,
                         dpop,
                         dphi[2],
                         bend_strength,
                         reference_momentum,
                         m,
                         substep_reference_cdt,
                         step_phase[2],
                         step_term[2]);

            f_kf(x, xp, y, yp, r0, k3);

            // drift_unit(x, xp, y, yp, cdt, dpop, c4 * step_length,
            // reference_momentum,
            //            m, substep_reference_cdt);

            bend_unit<T>(x,
                         xp,
                         y,
                         yp,
                         cdt,
                         dpop,
                         dphi[3],
                         bend_strength,
                         reference_momentum,
                         m,
                         substep_reference_cdt,
                         step_phase[3],
                         step_term[3]);

            f_kf(x, xp, y, yp, r0, k4);

            // drift_unit(x, xp, y, yp, cdt, dpop, c4 * step_length,
            // reference_momentum,
            //            m, substep_reference_cdt);

            bend_unit<T>(x,
                         xp,
                         y,
                         yp,
                         cdt,
                         dpop,
                         dphi[3],
                         bend_strength,
                         reference_momentum,
                         m,
                         substep_reference_cdt,
                         step_phase[3],
                         step_term[3]);

            f_kf(x, xp, y, yp, r0, k3);

            // drift_unit(x, xp, y, yp, cdt, dpop, c3 * step_length,
            // reference_momentum,
            //            m, substep_reference_cdt);

            bend_unit<T>(x,
                         xp,
                         y,
                         yp,
                         cdt,
                         dpop,
                         dphi[2],
                         bend_strength,
                         reference_momentum,
                         m,
                         substep_reference_cdt,
                         step_phase[2],
                         step_term[2]);

            f_kf(x, xp, y, yp, r0, k1);

            // drift_unit(x, xp, y, yp, cdt, dpop, c2 * step_length,
            // reference_momentum,
            //            m, substep_reference_cdt);

            bend_unit<T>(x,
                         xp,
                         y,
                         yp,
                         cdt,
                         dpop,
                         dphi[1],
                         bend_strength,
                         reference_momentum,
                         m,
                         substep_reference_cdt,
                         step_phase[1],
                         step_term[1]);

            f_kf(x, xp, y, yp, r0, k2);

            // drift_unit(x, xp, y, yp, cdt, dpop, c2 * step_length,
            // reference_momentum,
            //            m, substep_reference_cdt);

            bend_unit<T>(x,
                         xp,
                         y,
                         yp,
                         cdt,
                         dpop,
                         dphi[1],
                         bend_strength,
                         reference_momentum,
                         m,
                         substep_reference_cdt,
                         step_phase[1],
                         step_term[1]);

            f_kf(x, xp, y, yp, r0, k1);

            // drift_unit(x, xp, y, yp, cdt, dpop, c1 * step_length,
            // reference_momentum,
            //            m, substep_reference_cdt);

            bend_unit<T>(x,
                         xp,
                         y,
                         yp,
                         cdt,
                         dpop,
                         dphi[0],
                         bend_strength,
                         reference_momentum,
                         m,
                         substep_reference_cdt,
                         step_phase[0],
                         step_term[0]);
        }
    }

#endif

    template <typename T>
    KOKKOS_INLINE_FUNCTION void
    matrix_unit(T& x,
                T& xp,
                T& y,
                T& yp,
                T& cdt,
                T& dpop,
                double const k[6],
                double const rm[6][6],
                double const tm[6][6][6])
    {
        // Initialize coordinates with constant term
        T newx(k[0]);
        T newxp(k[1]);
        T newy(k[2]);
        T newyp(k[3]);
        T newcdt(k[4]);
        T newdpop(k[5]);

        // linear
        newx = newx + T(rm[0][0]) * x + T(rm[0][1]) * xp + T(rm[0][2]) * y +
                T(rm[0][3]) * yp + T(rm[0][4]) * cdt + T(rm[0][5]) * dpop;

        newxp = newxp + T(rm[1][0]) * x + T(rm[1][1]) * xp + T(rm[1][2]) * y +
                 T(rm[1][3]) * yp + T(rm[1][4]) * cdt + T(rm[1][5]) * dpop;

        newy = newy + T(rm[2][0]) * x + T(rm[2][1]) * xp + T(rm[2][2]) * y +
                T(rm[2][3]) * yp + T(rm[2][4]) * cdt + T(rm[2][5]) * dpop;

        newyp = newyp + T(rm[3][0]) * x + T(rm[3][1]) * xp + T(rm[3][2]) * y +
                 T(rm[3][3]) * yp + T(rm[3][4]) * cdt + T(rm[3][5]) * dpop;

        newcdt = newcdt + T(rm[4][0]) * x + T(rm[4][1]) * xp + T(rm[4][2]) * y +
                  T(rm[4][3]) * yp + T(rm[4][4]) * cdt + T(rm[4][5]) * dpop;

        newdpop = newdpop + T(rm[5][0]) * x + T(rm[5][1]) * xp + T(rm[5][2]) * y +
                   T(rm[5][3]) * yp + T(rm[5][4]) * cdt + T(rm[5][5]) * dpop;

        // tensor
        newx = newx + T(tm[0][0][0]) * x * x + T(tm[0][0][1]) * x * xp +
                T(tm[0][0][2]) * x * y + T(tm[0][0][3]) * x * yp +
                T(tm[0][0][4]) * x * cdt + T(tm[0][0][5]) * x * dpop +
                T(tm[0][1][0]) * xp * x + T(tm[0][1][1]) * xp * xp +
                T(tm[0][1][2]) * xp * y + T(tm[0][1][3]) * xp * yp +
                T(tm[0][1][4]) * xp * cdt + T(tm[0][1][5]) * xp * dpop +
                T(tm[0][2][0]) * y * x + T(tm[0][2][1]) * y * xp +
                T(tm[0][2][2]) * y * y + T(tm[0][2][3]) * y * yp +
                T(tm[0][2][4]) * y * cdt + T(tm[0][2][5]) * y * dpop +
                T(tm[0][3][0]) * yp * x + T(tm[0][3][1]) * yp * xp +
                T(tm[0][3][2]) * yp * y + T(tm[0][3][3]) * yp * yp +
                T(tm[0][3][4]) * yp * cdt + T(tm[0][3][5]) * yp * dpop +
                T(tm[0][4][0]) * cdt * x + T(tm[0][4][1]) * cdt * xp +
                T(tm[0][4][2]) * cdt * y + T(tm[0][4][3]) * cdt * yp +
                T(tm[0][4][4]) * cdt * cdt + T(tm[0][4][5]) * cdt * dpop +
                T(tm[0][5][0]) * dpop * x + T(tm[0][5][1]) * dpop * xp +
                T(tm[0][5][2]) * dpop * y + T(tm[0][5][3]) * dpop * yp +
                T(tm[0][5][4]) * dpop * cdt + T(tm[0][5][5]) * dpop * dpop;

        newxp = newxp + T(tm[1][0][0]) * x * x + T(tm[1][0][1]) * x * xp +
                 T(tm[1][0][2]) * x * y + T(tm[1][0][3]) * x * yp +
                 T(tm[1][0][4]) * x * cdt + T(tm[1][0][5]) * x * dpop +
                 T(tm[1][1][0]) * xp * x + T(tm[1][1][1]) * xp * xp +
                 T(tm[1][1][2]) * xp * y + T(tm[1][1][3]) * xp * yp +
                 T(tm[1][1][4]) * xp * cdt + T(tm[1][1][5]) * xp * dpop +
                 T(tm[1][2][0]) * y * x + T(tm[1][2][1]) * y * xp +
                 T(tm[1][2][2]) * y * y + T(tm[1][2][3]) * y * yp +
                 T(tm[1][2][4]) * y * cdt + T(tm[1][2][5]) * y * dpop +
                 T(tm[1][3][0]) * yp * x + T(tm[1][3][1]) * yp * xp +
                 T(tm[1][3][2]) * yp * y + T(tm[1][3][3]) * yp * yp +
                 T(tm[1][3][4]) * yp * cdt + T(tm[1][3][5]) * yp * dpop +
                 T(tm[1][4][0]) * cdt * x + T(tm[1][4][1]) * cdt * xp +
                 T(tm[1][4][2]) * cdt * y + T(tm[1][4][3]) * cdt * yp +
                 T(tm[1][4][4]) * cdt * cdt + T(tm[1][4][5]) * cdt * dpop +
                 T(tm[1][5][0]) * dpop * x + T(tm[1][5][1]) * dpop * xp +
                 T(tm[1][5][2]) * dpop * y + T(tm[1][5][3]) * dpop * yp +
                 T(tm[1][5][4]) * dpop * cdt + T(tm[1][5][5]) * dpop * dpop;

        newy = newy + T(tm[2][0][0]) * x * x + T(tm[2][0][1]) * x * xp +
                T(tm[2][0][2]) * x * y + T(tm[2][0][3]) * x * yp +
                T(tm[2][0][4]) * x * cdt + T(tm[2][0][5]) * x * dpop +
                T(tm[2][1][0]) * xp * x + T(tm[2][1][1]) * xp * xp +
                T(tm[2][1][2]) * xp * y + T(tm[2][1][3]) * xp * yp +
                T(tm[2][1][4]) * xp * cdt + T(tm[2][1][5]) * xp * dpop +
                T(tm[2][2][0]) * y * x + T(tm[2][2][1]) * y * xp +
                T(tm[2][2][2]) * y * y + T(tm[2][2][3]) * y * yp +
                T(tm[2][2][4]) * y * cdt + T(tm[2][2][5]) * y * dpop +
                T(tm[2][3][0]) * yp * x + T(tm[2][3][1]) * yp * xp +
                T(tm[2][3][2]) * yp * y + T(tm[2][3][3]) * yp * yp +
                T(tm[2][3][4]) * yp * cdt + T(tm[2][3][5]) * yp * dpop +
                T(tm[2][4][0]) * cdt * x + T(tm[2][4][1]) * cdt * xp +
                T(tm[2][4][2]) * cdt * y + T(tm[2][4][3]) * cdt * yp +
                T(tm[2][4][4]) * cdt * cdt + T(tm[2][4][5]) * cdt * dpop +
                T(tm[2][5][0]) * dpop * x + T(tm[2][5][1]) * dpop * xp +
                T(tm[2][5][2]) * dpop * y + T(tm[2][5][3]) * dpop * yp +
                T(tm[2][5][4]) * dpop * cdt + T(tm[2][5][5]) * dpop * dpop;

        newyp = newyp + T(tm[3][0][0]) * x * x + T(tm[3][0][1]) * x * xp +
                 T(tm[3][0][2]) * x * y + T(tm[3][0][3]) * x * yp +
                 T(tm[3][0][4]) * x * cdt + T(tm[3][0][5]) * x * dpop +
                 T(tm[3][1][0]) * xp * x + T(tm[3][1][1]) * xp * xp +
                 T(tm[3][1][2]) * xp * y + T(tm[3][1][3]) * xp * yp +
                 T(tm[3][1][4]) * xp * cdt + T(tm[3][1][5]) * xp * dpop +
                 T(tm[3][2][0]) * y * x + T(tm[3][2][1]) * y * xp +
                 T(tm[3][2][2]) * y * y + T(tm[3][2][3]) * y * yp +
                 T(tm[3][2][4]) * y * cdt + T(tm[3][2][5]) * y * dpop +
                 T(tm[3][3][0]) * yp * x + T(tm[3][3][1]) * yp * xp +
                 T(tm[3][3][2]) * yp * y + T(tm[3][3][3]) * yp * yp +
                 T(tm[3][3][4]) * yp * cdt + T(tm[3][3][5]) * yp * dpop +
                 T(tm[3][4][0]) * cdt * x + T(tm[3][4][1]) * cdt * xp +
                 T(tm[3][4][2]) * cdt * y + T(tm[3][4][3]) * cdt * yp +
                 T(tm[3][4][4]) * cdt * cdt + T(tm[3][4][5]) * cdt * dpop +
                 T(tm[3][5][0]) * dpop * x + T(tm[3][5][1]) * dpop * xp +
                 T(tm[3][5][2]) * dpop * y + T(tm[3][5][3]) * dpop * yp +
                 T(tm[3][5][4]) * dpop * cdt + T(tm[3][5][5]) * dpop * dpop;

        newcdt = newcdt + T(tm[4][0][0]) * x * x + T(tm[4][0][1]) * x * xp +
                  T(tm[4][0][2]) * x * y + T(tm[4][0][3]) * x * yp +
                  T(tm[4][0][4]) * x * cdt + T(tm[4][0][5]) * x * dpop +
                  T(tm[4][1][0]) * xp * x + T(tm[4][1][1]) * xp * xp +
                  T(tm[4][1][2]) * xp * y + T(tm[4][1][3]) * xp * yp +
                  T(tm[4][1][4]) * xp * cdt + T(tm[4][1][5]) * xp * dpop +
                  T(tm[4][2][0]) * y * x + T(tm[4][2][1]) * y * xp +
                  T(tm[4][2][2]) * y * y + T(tm[4][2][3]) * y * yp +
                  T(tm[4][2][4]) * y * cdt + T(tm[4][2][5]) * y * dpop +
                  T(tm[4][3][0]) * yp * x + T(tm[4][3][1]) * yp * xp +
                  T(tm[4][3][2]) * yp * y + T(tm[4][3][3]) * yp * yp +
                  T(tm[4][3][4]) * yp * cdt + T(tm[4][3][5]) * yp * dpop +
                  T(tm[4][4][0]) * cdt * x + T(tm[4][4][1]) * cdt * xp +
                  T(tm[4][4][2]) * cdt * y + T(tm[4][4][3]) * cdt * yp +
                  T(tm[4][4][4]) * cdt * cdt + T(tm[4][4][5]) * cdt * dpop +
                  T(tm[4][5][0]) * dpop * x + T(tm[4][5][1]) * dpop * xp +
                  T(tm[4][5][2]) * dpop * y + T(tm[4][5][3]) * dpop * yp +
                  T(tm[4][5][4]) * dpop * cdt + T(tm[4][5][5]) * dpop * dpop;

        newdpop = newdpop + T(tm[5][0][0]) * x * x + T(tm[5][0][1]) * x * xp +
                   T(tm[5][0][2]) * x * y + T(tm[5][0][3]) * x * yp +
                   T(tm[5][0][4]) * x * cdt + T(tm[5][0][5]) * x * dpop +
                   T(tm[5][1][0]) * xp * x + T(tm[5][1][1]) * xp * xp +
                   T(tm[5][1][2]) * xp * y + T(tm[5][1][3]) * xp * yp +
                   T(tm[5][1][4]) * xp * cdt + T(tm[5][1][5]) * xp * dpop +
                   T(tm[5][2][0]) * y * x + T(tm[5][2][1]) * y * xp +
                   T(tm[5][2][2]) * y * y + T(tm[5][2][3]) * y * yp +
                   T(tm[5][2][4]) * y * cdt + T(tm[5][2][5]) * y * dpop +
                   T(tm[5][3][0]) * yp * x + T(tm[5][3][1]) * yp * xp +
                   T(tm[5][3][2]) * yp * y + T(tm[5][3][3]) * yp * yp +
                   T(tm[5][3][4]) * yp * cdt + T(tm[5][3][5]) * yp * dpop +
                   T(tm[5][4][0]) * cdt * x + T(tm[5][4][1]) * cdt * xp +
                   T(tm[5][4][2]) * cdt * y + T(tm[5][4][3]) * cdt * yp +
                   T(tm[5][4][4]) * cdt * cdt + T(tm[5][4][5]) * cdt * dpop +
                   T(tm[5][5][0]) * dpop * x + T(tm[5][5][1]) * dpop * xp +
                   T(tm[5][5][2]) * dpop * y + T(tm[5][5][3]) * dpop * yp +
                   T(tm[5][5][4]) * dpop * cdt + T(tm[5][5][5]) * dpop * dpop;

        x = newx;
        xp = newxp;
        y = newy;
        yp = newyp;
        cdt = newcdt;
        dpop = newdpop;
    }
};


#endif // FF_ALGORITHM_H
