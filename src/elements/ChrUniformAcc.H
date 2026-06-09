/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_CHRACC_H
#define IMPACTX_CHRACC_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/pipeaperture.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"
#include "mixin/lineartransport.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>
#include <stdexcept>


namespace impactx::elements
{
    struct ChrAcc
    : public mixin::Named,
      public mixin::BeamOptic<ChrAcc>,
      public mixin::LinearTransport<ChrAcc>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "ChrAcc";
        using PType = ImpactXParticleContainer::ParticleType;

        /** Acceleration in a uniform field Ez, with a uniform solenoidal field Bz.
         *
         * The Hamiltonian is expanded through second order in the
         * transverse variables (x,px,y,py), with the exact pt dependence retained.
         *
         * @param ds Segment length in m.
         * @param ez electric field strength in m^(-1)
         *           = (charge * electric field Ez in V/m) / (m*c^2)
         * @param bz magnetic field strength in m^(-1)
         *           = (charge * magnetic field Bz in T) / (m*c)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        ChrAcc (
            amrex::ParticleReal ds,
            amrex::ParticleReal ez,
            amrex::ParticleReal bz,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            amrex::ParticleReal aperture_x = 0,
            amrex::ParticleReal aperture_y = 0,
            int nslice = 1,
            std::optional<std::string> name = std::nullopt
        )
        : Named(std::move(name)),
          Thick(ds, nslice),
          Alignment(dx, dy, rotation_degree),
          PipeAperture(aperture_x, aperture_y),
          m_ez(ez), m_bz(bz)
        {
        }

        /** Reverse the element (negate ds) */
        void reverse () { Thick::reverse(); }

        /** Push all particles */
        using BeamOptic::operator();

        /** Compute and cache the constants for the push.
         *
         * In particular, used to pre-compute and cache variables that are
         * independent of the individually tracked particle.
         *
         * @param refpart reference particle
         */
        void compute_constants (RefPart const & refpart)
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            Alignment::compute_constants(refpart);

            // length of the current slice
            m_slice_ds = m_ds / nslice();

            // access reference particle values (final, initial):
            m_ptf_ref = refpart.pt;
            m_pti_ref = m_ptf_ref + m_ez * m_slice_ds;
            m_bgf = std::sqrt(powi<2>(m_ptf_ref) - 1.0_prt);
            m_bgi = std::sqrt(powi<2>(m_pti_ref) - 1.0_prt);

            // compute focusing constant (1/m) and rotation angle (in rad)
            m_alpha = m_bz * 0.5_prt;
            m_alpha_iez = m_alpha / m_ez;
        }

        /** This is a chracc functor, so that a variable of this type can be used like a chracc function.
         *
         * The @see compute_constants method must be called before pushing particles through this operator.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py,
            T_Real & AMREX_RESTRICT pt,
            T_IdCpu & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;
            using namespace std; // for cmath(float)
            namespace stdx = amrex::simd::stdx;

            // initial conversion from static to dynamic units:
            px = px * m_bgi;
            py = py * m_bgi;
            pt = pt * m_bgi;

            // compute intermediate quantities related to acceleration
            T_Real const pti_tot = m_pti_ref + pt;
            T_Real const ptf_tot = m_ptf_ref + pt;
            T_Real const pti_tot2 = powi<2>(pti_tot);
            T_Real const ptf_tot2 = powi<2>(ptf_tot);

            // check whether particle lies within the domain of map definition
            auto const mask = pti_tot2 <= 1_prt || ptf_tot2 <= 1_prt;
            amrex::ParticleIDWrapper<T_IdCpu>{idcpu}.make_invalid(mask);
            {
                T_Real const pzi_tot = sqrt(powi<2>(pti_tot) - 1_prt);
                T_Real const pzf_tot = sqrt(powi<2>(ptf_tot) - 1_prt);
                T_Real const pzi_ref = sqrt(powi<2>(m_pti_ref) - 1_prt);
                T_Real const pzf_ref = sqrt(powi<2>(m_ptf_ref) - 1_prt);

                T_Real const numer = -ptf_tot + pzf_tot;
                T_Real const denom = -pti_tot + pzi_tot;

                // compute focusing constant (1/m) and rotation angle (in rad)
                T_Real const theta = m_alpha_iez * log(numer / denom);
                auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);

                // initialize output values
                T_Real xout = x;
                T_Real yout = y;
                T_Real tout = t;
                T_Real pxout = px;
                T_Real pyout = py;
                T_Real ptout = pt;

                // advance positions and momenta using map for focusing
                xout = cos_theta * x + sin_theta / m_alpha * px;
                pxout = -m_alpha * sin_theta * x + cos_theta * px;

                yout = cos_theta * y + sin_theta / m_alpha * py;
                pyout = -m_alpha * sin_theta * y + cos_theta * py;

                // the correct symplectic update for t
                tout = t + (pzf_tot - pzf_ref - pzi_tot + pzi_ref) / m_ez;
                tout += (1_prt/pzi_tot - 1_prt/pzf_tot)
                      * (powi<2>(py - m_alpha * x) +
                         powi<2>(px + m_alpha * y))
                      / (2_prt * m_ez);
                ptout = pt;

                // assign intermediate momenta
                px = pxout;
                py = pyout;
                pt = ptout;

                // advance positions and momenta using map for rotation
                x = cos_theta * xout + sin_theta * yout;
                pxout = cos_theta * px + sin_theta * py;

                y = -sin_theta * xout + cos_theta * yout;
                pyout = -sin_theta * px + cos_theta * py;

                t = tout;
                ptout = pt;

                // assign updated momenta
                px = pxout;
                py = pyout;
                pt = ptout;

            }

            // final conversion from dynamic to static units:
            px = px / m_bgf;
            py = py / m_bgf;
            pt = pt / m_bgf;
        }

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            // assign input reference particle values
            amrex::ParticleReal const x = refpart.x;
            amrex::ParticleReal const px = refpart.px;
            amrex::ParticleReal const y = refpart.y;
            amrex::ParticleReal const py = refpart.py;
            amrex::ParticleReal const z = refpart.z;
            amrex::ParticleReal const pz = refpart.pz;
            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;
            amrex::ParticleReal const s = refpart.s;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // compute initial value of beta*gamma
            amrex::ParticleReal const bgi = std::sqrt(powi<2>(pt) - 1.0_prt);

            // advance pt (uniform acceleration)
            refpart.pt = pt - m_ez*slice_ds;

            // compute final value of beta*gamma
            amrex::ParticleReal const ptf = refpart.pt;
            amrex::ParticleReal const bgf = std::sqrt(powi<2>(ptf) - 1.0_prt);

            // update t
            refpart.t = t + (bgf - bgi)/m_ez;

            // advance position (x,y,z)
            refpart.x = x + slice_ds*px/bgi;
            refpart.y = y + slice_ds*py/bgi;
            refpart.z = z + slice_ds*pz/bgi;

            // advance momentum (px,py,pz)
            refpart.px = px*bgf/bgi;
            refpart.py = py*bgf/bgi;
            refpart.pz = pz*bgf/bgi;

            // advance integrated path length
            refpart.s = s + slice_ds;
        }

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
         *
         * Jacobian at the reference trajectory of the ChrAcc push (see
         * @see operator()). The nonlinear push composes, in order:
         *   (1) static->dynamic unit conversion: diag(1, bgi, 1, bgi, 1, bgi)
         *       acting on (x, px, y, py, t, pt);
         *   (2) solenoid-like focusing by phase advance theta (cf. Sol);
         *   (3) rotation by the same theta in the (x, y) plane;
         *   (4) dynamic->static unit conversion: diag(1, 1/bgf, 1, 1/bgf,
         *       1, 1/bgf) on the output;
         * plus a symplectic correction to t. Here
         *   theta = (alpha/ez) * log((pzf_ref - ptf_ref)/(pzi_ref - pti_ref))
         * with alpha = Bz/2, ez the normalized longitudinal electric field,
         * pzi_ref = sqrt(pti_ref^2 - 1), pzf_ref = sqrt(ptf_ref^2 - 1), and
         * (pti_ref, ptf_ref) = (initial, final) reference pt on this slice.
         *
         * This routine is called after the reference particle has been
         * advanced by one slice (see diagnostics/LinearMap.cpp), so
         * refpart.pt holds ptf_ref and pti_ref = ptf_ref + ez * slice_ds.
         *
         * The closed-form linear map was derived with SymPy by Jacobian
         * of the exact push at (x, px, y, py, t, pt) = 0. Writing
         * c = cos(theta0), s = sin(theta0) with theta0 the reference
         * phase advance above, the nonzero entries are
         *   R(1,1) = c^2,                R(1,2) = bgi * c*s / alpha,
         *   R(1,3) = c*s,                R(1,4) = bgi * s^2 / alpha,
         *   R(2,1) = -alpha * c*s / bgf, R(2,2) = (bgi/bgf) * c^2,
         *   R(2,3) = -alpha * s^2 / bgf, R(2,4) = (bgi/bgf) * c*s,
         *   R(3,1) = -c*s,               R(3,2) = -bgi * s^2 / alpha,
         *   R(3,3) = c^2,                R(3,4) = bgi * c*s / alpha,
         *   R(4,1) =  alpha * s^2 / bgf, R(4,2) = -(bgi/bgf) * c*s,
         *   R(4,3) = -alpha * c*s / bgf, R(4,4) = (bgi/bgf) * c^2,
         *   R(5,5) = 1,
         *   R(5,6) = (bgi/ez) * (ptf_ref/pzf_ref - pti_ref/pzi_ref),
         *   R(6,6) = bgi/bgf.
         * All other entries are zero. The quadratic-in-(px, py, x, y)
         * part of the t correction vanishes at linear order.
         *
         * References:
         *   - A. J. Dragt, Lie Methods for Nonlinear Dynamics with
         *     Applications to Accelerator Physics, Univ. of Maryland,
         *     2020 (MaryLie documentation), chapter on uniform-Ez
         *     acceleration with solenoidal focusing;
         *   - H. Wiedemann, Particle Accelerator Physics, 4th ed.,
         *     Springer 2015, chapters on electrostatic acceleration and
         *     solenoid focusing;
         *   - The Hamiltonian and nonlinear push follow Chad Mitchell's
         *     formulation already implemented in @see operator().
         *
         * Ideal behavior only: alignment and rotation errors and
         * feed-down from an off-axis reference orbit are ignored.
         *
         * @param[in] refpart reference particle after this slice's
         *                    reference advance (pt = ptf_ref)
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map (RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            // slice length and initial/final reference pt:
            // transport_map is called after the reference advance, so
            // refpart.pt holds the post-slice (final) pt.
            amrex::ParticleReal const slice_ds = m_ds / nslice();
            amrex::ParticleReal const ptf_ref = refpart.pt;
            amrex::ParticleReal const pti_ref = ptf_ref + m_ez * slice_ds;

            amrex::ParticleReal const pzf_ref = std::sqrt(powi<2>(ptf_ref) - 1.0_prt);
            amrex::ParticleReal const pzi_ref = std::sqrt(powi<2>(pti_ref) - 1.0_prt);
            amrex::ParticleReal const bgf = pzf_ref;  // |beta*gamma|_final
            amrex::ParticleReal const bgi = pzi_ref;  // |beta*gamma|_initial

            // Solenoid focusing constant alpha = Bz / 2. The "log"
            // factor is the integrated drift length in dynamic units and
            // is well-defined even when alpha = 0 (no solenoid).
            amrex::ParticleReal const alpha = 0.5_prt * m_bz;
            amrex::ParticleReal const drift_log =
                std::log((pzf_ref - ptf_ref) / (pzi_ref - pti_ref)) / m_ez;
            amrex::ParticleReal const rbg = bgi / bgf;

            Map6x6 R = Map6x6::Identity();

            // t couples to pt only: symplectic correction to the
            // time-of-flight, linearized at pt = 0.
            R(5,6) = bgi * (ptf_ref / pzf_ref - pti_ref / pzi_ref) / m_ez;
            // pt scales by the static<->dynamic unit change bgi/bgf.
            R(6,6) = rbg;

            if (alpha == 0.0_prt)
            {
                // No solenoid: pure accelerating "drift" with adiabatic
                // damping on the transverse momenta.
                R(1,2) = bgi * drift_log;
                R(2,2) = rbg;
                R(3,4) = bgi * drift_log;
                R(4,4) = rbg;
            }
            else
            {
                // Solenoid + acceleration: rotation * focusing, bracketed
                // by static<->dynamic unit conversions on (px, py).
                amrex::ParticleReal const theta0 = alpha * drift_log;
                auto const [s, c] = amrex::Math::sincos(theta0);

                amrex::ParticleReal const c2 = c * c;
                amrex::ParticleReal const s2 = s * s;
                amrex::ParticleReal const cs = c * s;

                R(1,1) =  c2;                R(1,2) =  bgi * cs / alpha;
                R(1,3) =  cs;                R(1,4) =  bgi * s2 / alpha;
                R(2,1) = -alpha * cs / bgf;  R(2,2) =  rbg * c2;
                R(2,3) = -alpha * s2 / bgf;  R(2,4) =  rbg * cs;
                R(3,1) = -cs;                R(3,2) = -bgi * s2 / alpha;
                R(3,3) =  c2;                R(3,4) =  bgi * cs / alpha;
                R(4,1) =  alpha * s2 / bgf;  R(4,2) = -rbg * cs;
                R(4,3) = -alpha * cs / bgf;  R(4,4) =  rbg * c2;
            }

            return R;
        }

        amrex::ParticleReal m_ez; //! electric field strength in 1/m
        amrex::ParticleReal m_bz; //! magnetic field strength in 1/m

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_slice_ds;  //! m_ds / nslice();
        amrex::ParticleReal m_ptf_ref, m_pti_ref, m_bgf, m_bgi;
        amrex::ParticleReal m_alpha, m_alpha_iez;
    };

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::ChrAcc)

#endif // IMPACTX_CHRACC_H
