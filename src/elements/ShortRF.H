/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_SHORTRF_H
#define IMPACTX_SHORTRF_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"
#include "mixin/lineartransport.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>

#include <cmath>


namespace impactx::elements
{
    struct ShortRF
    : public mixin::Named,
      public mixin::BeamOptic<ShortRF>,
      public mixin::LinearTransport<ShortRF>,
      public mixin::Thin,
      public mixin::Alignment,
      public mixin::NoFinalize
    {
        static constexpr auto type = "ShortRF";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A short RF cavity element
         *
         * @param V Normalized RF voltage V = maximum energy gain/(m*c^2)
         * @param freq RF frequency in Hz
         * @param phase RF synchronous phase of the reference particle in deg
         *     phi = 0:  maximum energy gain (on-crest)
         *     phi = -90 deg:  zero-crossing for bunching
         *     phi = 90 deg:  zero-crossing for debunching
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param name a user defined and not necessarily unique name of the element
         */
        ShortRF (
            amrex::ParticleReal V,
            amrex::ParticleReal freq,
            amrex::ParticleReal phase,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            std::optional<std::string> name = std::nullopt
        )
        : Named(std::move(name)),
          Alignment(dx, dy, rotation_degree),
          m_V(V), m_freq(freq), m_phase(phase)
        {
        }

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

            Alignment::compute_constants(refpart);

            // Define parameters and intermediate constants
            using ablastr::constant::math::pi;
            using ablastr::constant::SI::c;
            m_k = 2.0_prt * pi / c * m_freq;
            m_phi = m_phase * pi / 180.0_prt;

            // access reference particle values (final, initial):
            amrex::ParticleReal const ptf_ref = refpart.pt;
            m_V_cos_phi = m_V * std::cos(m_phi);
            amrex::ParticleReal const pti_ref = ptf_ref + m_V_cos_phi;
            m_bgf = std::sqrt(std::pow(ptf_ref, 2) - 1.0_prt);
            m_bgi = std::sqrt(std::pow(pti_ref, 2) - 1.0_prt);
        }

        /** This is a shortrf functor, so that a variable of this type can be used like a
         *  shortrf function.
         *
         * The @see compute_constants method must be called before pushing particles through this operator.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index (unused)
         * @param refpart reference particle (unused)
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            amrex::ParticleReal & AMREX_RESTRICT x,
            amrex::ParticleReal & AMREX_RESTRICT y,
            amrex::ParticleReal & AMREX_RESTRICT t,
            amrex::ParticleReal & AMREX_RESTRICT px,
            amrex::ParticleReal & AMREX_RESTRICT py,
            amrex::ParticleReal & AMREX_RESTRICT pt,
            [[maybe_unused]] uint64_t & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {

            using namespace amrex::literals; // for _rt and _prt

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // initial conversion from static to dynamic units:
            px = px * m_bgi;
            py = py * m_bgi;
            pt = pt * m_bgi;

            // intialize output values
            amrex::ParticleReal xout = x;
            amrex::ParticleReal yout = y;
            amrex::ParticleReal tout = t;
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            // advance position and momentum in dynamic units
            // xout = x;
            pxout = px;

            // yout = y;
            pyout = py;

            // tout = t;
            ptout = pt - m_V * std::cos(m_k * t + m_phi) + m_V_cos_phi;

            // assign updated values
            x = xout;
            y = yout;
            t = tout;
            px = pxout;
            py = pyout;
            pt = ptout;

            // final conversion from dynamic to static units:
            px = px / m_bgf;
            py = py / m_bgf;
            pt = pt / m_bgf;

            // undo shift due to alignment errors of the element
            shift_out(x, y, px, py);
        }

        /** This pushes the reference particle. */
        using Thin::operator();

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // assign input reference particle values
            amrex::ParticleReal const x = refpart.x;
            amrex::ParticleReal const px = refpart.px;
            amrex::ParticleReal const y = refpart.y;
            amrex::ParticleReal const py = refpart.py;
            amrex::ParticleReal const z = refpart.z;
            amrex::ParticleReal const pz = refpart.pz;
            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;

            // Define parameters and intermediate constants
            using ablastr::constant::math::pi;
            amrex::ParticleReal const phi = m_phase*(pi/180.0_prt);

            // compute intial value of beta*gamma
            amrex::ParticleReal const bgi = std::sqrt(std::pow(pt, 2) - 1.0_prt);

            // advance pt
            refpart.pt = pt - m_V * std::cos(phi);

            // compute final value of beta*gamma
            amrex::ParticleReal const ptf = refpart.pt;
            amrex::ParticleReal const bgf = std::sqrt(std::pow(ptf, 2) - 1.0_prt);

            // advance position (x,y,z,t)
            refpart.x = x;
            refpart.y = y;
            refpart.z = z;
            refpart.t = t;

            // advance momentum (px,py,pz)
            refpart.px = px*bgf/bgi;
            refpart.py = py*bgf/bgi;
            refpart.pz = pz*bgf/bgi;

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
        transport_map ([[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // Define parameters and intermediate constants
            using ablastr::constant::math::pi;
            using ablastr::constant::SI::c;
            amrex::ParticleReal const k = (2.0_prt*pi/c)*m_freq;
            amrex::ParticleReal const phi = m_phase*(pi/180.0_prt);

            // access reference particle values (final, initial):
            amrex::ParticleReal const ptf_ref = refpart.pt;
            amrex::ParticleReal const pti_ref = ptf_ref + m_V * std::cos(phi);
            amrex::ParticleReal const bgf = std::sqrt(std::pow(ptf_ref, 2) - 1.0_prt);
            amrex::ParticleReal const bgi = std::sqrt(std::pow(pti_ref, 2) - 1.0_prt);

            // initialize transport matrix
            Map6x6 R = Map6x6::Identity();

            // This is a nonlinear element: the linearized map is returned here.
            // assign linear map matrix elements
            R(2,2) = bgi/bgf;
            R(4,4) = bgi/bgf;
            R(6,5) = k*m_V*std::sin(phi)/bgf;
            R(6,6) = bgi/bgf;

            return R;
        }

        amrex::ParticleReal m_V; //! normalized (max) RF voltage drop.
        amrex::ParticleReal m_freq; //! RF frequency in Hz.
        amrex::ParticleReal m_phase; //! reference RF phase in degrees.

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_k, m_phi, m_V_cos_phi, m_bgi, m_bgf;
    };

} // namespace impactx

#endif // IMPACTX_SHORTRF_H
