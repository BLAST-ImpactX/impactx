/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl, Kyrre Sjobak
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_CHRPLASMALENS_H
#define IMPACTX_CHRPLASMALENS_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/pipeaperture.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"
#include "mixin/lineartransport.H"
#include "mixin/spintransport.H"
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
    struct ChrPlasmaLens
    : public mixin::Named,
      public mixin::BeamOptic<ChrPlasmaLens>,
      public mixin::LinearTransport<ChrPlasmaLens>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::SpinTransport,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "ChrPlasmaLens";
        using PType = ImpactXParticleContainer::ParticleType;

        /** An active cylindrically symmetric plasma lens with chromatic focusing
         *
         * The Hamiltonian is expanded through second order in the
         * transverse variables (x,px,y,py), with the exact pt dependence retained.
         *
         * @param ds Segment length in m.
         * @param k  plasma lens focusing strength in m^(-2)
         *           = (azimuthal field gradient in T/m) / (rigidity in T-m)
         *        OR azimuthal magnetic field gradient in T/m
         * @param unit  Unit specification
         *           unit = 0 focusing strength in m^(-2)
         *           unit = 1 focusing strength in T/m
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        ChrPlasmaLens (
            amrex::ParticleReal ds,
            amrex::ParticleReal k,
            int unit,
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
          m_k(k), m_unit(unit)
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

            // access reference particle values to find beta and gamma
            m_beta = refpart.beta();
            amrex::ParticleReal const gamma = refpart.gamma();

            // normalize focusing strength units to MAD-X convention if needed
            m_g = m_unit == 1 ? m_k / refpart.rigidity_Tm() : m_k;

            //For m_g=0 time propagation
            m_const1 = 1_prt / (2_prt * powi<3>(m_beta) * powi<2>(gamma));
        }

        /** This is a chrplasmalens functor, so that a variable of this type can be used like a chrplasmalens function.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index
         * @param refpart reference particle
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py,
            T_Real const & AMREX_RESTRICT pt,
            [[maybe_unused]] T_IdCpu const & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;
            using namespace std; // for cmath(float)

            // compute particle momentum deviation delta + 1
            T_Real const delta1 = sqrt(1_prt - 2_prt*pt/m_beta + powi<2>(pt));
            T_Real const delta = delta1 - 1_prt;

            // compute phase advance per unit length in s (in rad/m)
            // chromatic dependence on delta is included
            T_Real const omega = sqrt(abs(m_g)/delta1);

            // initialize output values
            T_Real xout = x;
            T_Real yout = y;
            T_Real tout = t;

            // initialize output values of momenta
            T_Real pxout = px;
            T_Real pyout = py;
            // T_Real const ptout = pt;

            // placeholder variables -- canonical position and momentum
            T_Real const q1 = x;
            T_Real const q2 = y;
            T_Real const p1 = px;
            T_Real const p2 = py;

            if (m_g != 0_prt)
            {
                // Focusing and defocusing differ only by the trig family
                // (sin/cos vs. sinh/cosh), the sign of the momentum kick, and
                // the sign of the q^2 w^2 term in the longitudinal update. The
                // m_g sign is uniform across particles and delta1 > 0 always,
                // so we select the trig factors branch-free. Only raw trig
                // values are selected, so the arithmetic order matches the
                // original branches byte-for-byte. The m_g == 0 drift case is
                // handled separately below to avoid the 1/(omega*delta1)
                // division.
                bool const focusing = m_g > 0_prt;
                auto const [sin_ods, cos_ods] = amrex::Math::sincos(omega*m_slice_ds);
                T_Real const sinh_ods = sinh(omega*m_slice_ds);
                T_Real const cosh_ods = cosh(omega*m_slice_ds);
                T_Real const c = focusing ? cos_ods : cosh_ods;
                T_Real const s = focusing ? sin_ods : sinh_ods;
                // momentum kick sign: -1 (focusing) / +1 (defocusing)
                T_Real const a = focusing ? -1_prt : 1_prt;

                // advance transverse position and momentum (both planes alike)
                xout = c * x + s / (omega * delta1) * px;
                pxout = a * omega * delta1 * s * x + c * px;

                yout = c * y + s / (omega * delta1) * py;
                pyout = a * omega * delta1 * s * y + c * py;

                // the corresponding symplectic update to t (both x and y)
                T_Real const term = pt + delta/m_beta;
                T_Real const t0 = t - term*m_slice_ds/delta1;

                // 2*omega*slice trig: sin/cos (focusing) vs. sinh/cosh (defocusing)
                T_Real const s2 = focusing ? sin(2_prt*m_slice_ds*omega) : sinh(2_prt*m_slice_ds*omega);
                T_Real const c2 = focusing ? cos(2_prt*m_slice_ds*omega) : cosh(2_prt*m_slice_ds*omega);
                // sign of the q^2 w^2 contribution: -1 (focusing) / +1 (defocusing)
                T_Real const sgn2 = focusing ? -1_prt : 1_prt;

                T_Real const w = omega*delta1;
                T_Real const term1 = -(powi<2>(p2) + sgn2 * powi<2>(q2) * powi<2>(w)) * s2;
                T_Real const term2 = -(powi<2>(p1) + sgn2 * powi<2>(q1) * powi<2>(w)) * s2;
                T_Real const term3 = -2_prt * q2 * p2 * w * c2;
                T_Real const term4 = -2_prt * q1 * p1 * w * c2;
                T_Real const term5 = 2_prt * omega * (q1*p1*delta1 + q2*p2*delta1
                    -(powi<2>(p1) + powi<2>(p2))*m_slice_ds - (powi<2>(q1) + powi<2>(q2)) * powi<2>(w)*m_slice_ds);
                tout = t0 + (-1_prt+m_beta*pt)
                          / (8_prt*m_beta * powi<3>(delta1)*omega)
                          * (term1+term2+term3+term4+term5);
            }
            else {
                xout = x + px*m_slice_ds / delta1;
                pxout = px;

                yout = y + py*m_slice_ds / delta1;
                pyout = py;

                // the corresponding symplectic update to t (zero strength = drift),
                // avoiding division by 0
                T_Real term = 2_prt * powi<2>(pt) + powi<2>(px) + powi<2>(py);
                term = 2_prt - 4_prt * m_beta * pt + powi<2>(m_beta) * term;
                term = -2_prt + powi<2>(refpart.gamma())*term;
                term = (-1_prt + m_beta * pt) * term;
                term = term * m_const1;
                tout = t - m_slice_ds * (1_prt / m_beta + term / powi<3>(delta1));
            }
            // advance longitudinal position and momentum


            // assign updated position & momenta
            x = xout;
            y = yout;
            t = tout;
            px = pxout;
            py = pyout;
            // pt = ptout;
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

            // assign intermediate parameter
            amrex::ParticleReal const step = slice_ds / std::sqrt(powi<2>(pt)-1.0_prt);

            // advance position and momentum (straight element)
            refpart.x = x + step*px;
            refpart.y = y + step*py;
            refpart.z = z + step*pz;
            refpart.t = t - step*pt;

            // advance integrated path length
            refpart.s = s + slice_ds;
        }


        /** This operator pushes the particle spin.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param sx particle spin x-component
         * @param sy particle spin y-component
         * @param sz particle spin z-component
         * @param idcpu particle global index
         * @param refpart reference particle
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void spin_and_phasespace_push (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py,
            T_Real const & AMREX_RESTRICT pt,
            T_Real & AMREX_RESTRICT sx,
            T_Real & AMREX_RESTRICT sy,
            T_Real & AMREX_RESTRICT sz,
            T_IdCpu const & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using namespace std; // for cmath(float)
            using amrex::Math::powi;

            // compute particle momentum deviation delta + 1
            T_Real const delta1 = sqrt(1_prt - 2_prt*pt/m_beta + powi<2>(pt));
            T_Real const inv_delta1 = 1_prt / delta1;

            // compute particle relativistic gamma
            T_Real const gamma = refpart.gamma() * (1_prt - m_beta*pt);
            T_Real const gyro_const = 1_prt + refpart.gyromagnetic_anomaly * gamma;

            // compute phase advance per unit length in s (in rad/m)
            // chromatic dependence on delta is included
            T_Real const omega = sqrt(abs(m_g)*inv_delta1);

            // compute trigonometric quantities
            auto const [sin_omega_ds, cos_omega_ds] = amrex::Math::sincos(omega*m_slice_ds);
            T_Real const sinh_omega_ds = sinh(omega*m_slice_ds);
            T_Real const cosh_omega_ds = cosh(omega*m_slice_ds);

            // A plasma lens focuses (or defocuses) both transverse planes
            // identically, so the x- and y-generators share the same factor and
            // differ only by an overall sign. The focusing/defocusing/drift
            // selection depends only on the (particle-uniform) sign of m_g, so
            // we select the trig factors branch-free. omega is real
            // (delta1 > 0 always); only raw trig values are selected so the
            // arithmetic order matches the original branches byte-for-byte. For
            // m_g == 0, omega == 0 makes (1-cos), sin, (cosh-1), sinh vanish, so
            // the generator is identically zero (drift) with no division by omega.
            bool const focusing = m_g > 0.0_prt;
            T_Real const fac = focusing ? (1_prt - cos_omega_ds) : (cosh_omega_ds - 1_prt);
            T_Real const sin_sel = focusing ? sin_omega_ds : sinh_omega_ds;
            // sign of lambdax: +gyro when focusing, -gyro when defocusing/drift
            T_Real const sgn = focusing ? gyro_const : -gyro_const;

            T_Real const lambdax =  sgn * ( py*inv_delta1*fac + y*omega*sin_sel );
            T_Real const lambday = -sgn * ( px*inv_delta1*fac + x*omega*sin_sel );
            T_Real const lambdaz = 0_prt;

            // push the spin vector using the generator just determined
            rotate_spin(lambdax,lambday,lambdaz,sx,sy,sz);

            // phase space push
            (*this)(x, y, t, px, py, pt, idcpu, refpart);

        }


        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
         *
         * @param[in] refpart reference particle
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map (RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = powi<2>(pt_ref) - 1_prt;

            // compute phase advance per unit length in s (in rad/m)
            amrex::ParticleReal const omega = std::sqrt(std::abs(m_k));

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();

            if (m_k > 0.0) {
                R(1,1) = std::cos(omega*slice_ds);
                R(1,2) = std::sin(omega*slice_ds)/omega;
                R(2,1) = -omega*std::sin(omega*slice_ds);
                R(2,2) = std::cos(omega*slice_ds);
                R(3,3) = std::cos(omega*slice_ds);
                R(3,4) = std::sin(omega*slice_ds)/omega;
                R(4,3) = -omega*std::sin(omega*slice_ds);
                R(4,4) = std::cos(omega*slice_ds);
                R(5,6) = slice_ds/betgam2;
            } else if (m_k < 0.0) {
                R(1,1) = std::cosh(omega*slice_ds);
                R(1,2) = std::sinh(omega*slice_ds)/omega;
                R(2,1) = omega*std::sinh(omega*slice_ds);
                R(2,2) = std::cosh(omega*slice_ds);
                R(3,3) = std::cosh(omega*slice_ds);
                R(3,4) = std::sinh(omega*slice_ds)/omega;
                R(4,3) = omega*std::sinh(omega*slice_ds);
                R(4,4) = std::cosh(omega*slice_ds);
                R(5,6) = slice_ds/betgam2;
            } else {
                R(1,2) = slice_ds;
                R(3,4) = slice_ds;
                R(5,6) = slice_ds / betgam2;
            }

            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(R);
        }

        amrex::ParticleReal m_k; //! focusing strength in 1/m^2 (or T/m)
        int m_unit; //! unit specification for focusing strength

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_slice_ds;  //! m_ds / nslice();
        amrex::ParticleReal m_beta;
        amrex::ParticleReal m_g;
        amrex::ParticleReal m_const1;
    };

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::ChrPlasmaLens)

#endif // IMPACTX_CHRPLASMALENS_H
