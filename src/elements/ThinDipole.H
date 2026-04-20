/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_THINDIPOLE_H
#define IMPACTX_THINDIPOLE_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"
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
    struct ThinDipole
    : public mixin::Named,
      public mixin::BeamOptic<ThinDipole>,
      public mixin::LinearTransport<ThinDipole>,
      public mixin::Thin,
      public mixin::Alignment,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "ThinDipole";
        using PType = ImpactXParticleContainer::ParticleType;

        static constexpr amrex::ParticleReal degree2rad = ablastr::constant::math::pi / 180.0;

        /** A general thin-kick dipole element with chromatic effects
         *
         * This model is equivalent to the model described in Section 3.1 of
         * G. Ripken, F. Schmidt, "A Symplectic Six-Dimensional Thin-Lens Formalism
         * for Tracking," CERN/SL/95-12 (AP), 1995.  It is intended to replicate
         * the thin-lens dipole model in MAD-X.
         *
         * @param theta - the total bending angle (degrees)
         * @param rc - the curvature radius (m)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param name a user defined and not necessarily unique name of the element
         */
        ThinDipole (
            amrex::ParticleReal theta,
            amrex::ParticleReal rc,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            std::optional<std::string> name = std::nullopt
        )
        : Named(std::move(name)),
          Alignment(dx, dy, rotation_degree),
          m_theta(theta * degree2rad),
          m_rc(rc)
        {
        }

        /** Reverse the element (negate bending angle) */
        void reverse () { m_theta = -m_theta; }

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

            // access reference particle to find relativistic beta
            m_beta = refpart.beta();

            // compute the effective (equivalent) arc length and curvature
            m_ds = m_theta * m_rc;
            m_kx = 1.0_prt / m_rc;
        }

        /** This is a multipole functor, so that a variable of this type can be used like a
         *  multipole function.
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
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py,
            T_Real & AMREX_RESTRICT pt,
            [[maybe_unused]] T_IdCpu & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;
            using namespace std; // for cmath(float)

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // access position data
            T_Real const xout = x;
            T_Real const yout = y;
            T_Real const tout = t;

            // initialize output values of momenta
            T_Real pxout = px;
            T_Real pyout = py;
            T_Real ptout = pt;

            // compute the function expressing dp/p in terms of pt (labeled f in Ripken etc.)
            T_Real const f = -1.0_prt + sqrt(1.0_prt - 2.0_prt * pt / m_beta + powi<2>(pt));
            T_Real const fprime = (1.0_prt - m_beta*pt) / (m_beta * (1.0_prt + f));

            // advance position and momentum
            x = xout;
            pxout = px - powi<2>(m_kx) * m_ds * xout + m_kx * m_ds * f; //eq (3.2b)

            y = yout;
            pyout = py;

            t = tout + m_kx * xout * m_ds * fprime; //eq (3.2e)
            ptout = pt;

            // assign updated momenta
            px = pxout;
            py = pyout;
            pt = ptout;

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
            // assign input reference particle values
            amrex::ParticleReal const px = refpart.px;
            amrex::ParticleReal const pz = refpart.pz;

            // calculate expensive terms once
            auto const [sin_theta, cos_theta] = amrex::Math::sincos(m_theta);

            // advance position and momentum (thin dipole)
            refpart.px = px*cos_theta - pz*sin_theta;
            refpart.pz = pz*cos_theta + px*sin_theta;
        }

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
         *
         * Linearization of the Ripken-Schmidt thin-dipole kick
         * (eqs. 3.2 of CERN/SL/95-12 (AP), 1995, already cited in the
         * class docstring) about the reference trajectory. With the
         * effective arc length ds = theta * rc and curvature kx = 1/rc,
         * the thin kick is
         *   dpx = - kx^2 * ds * x + kx * ds * f(pt)
         *   dt  = + kx * ds * x * f'(pt)
         * where f(pt) = -1 + sqrt(1 - 2*pt/beta + pt^2) satisfies
         * f(0) = 0, f'(0) = 1/beta. The Jacobian at the origin yields:
         *   R(2,1) = - kx^2 * ds     = - theta / rc
         *   R(2,6) = - kx * ds / beta = - theta / beta
         *   R(5,1) = + kx * ds / beta = + theta / beta
         * All other entries are identity.
         *
         * Cross-check: in the short-arc limit (theta << 1) this agrees
         * with Sbend::transport_map, specifically
         *   R(2,1) = -sin(theta)/rc  ~ -theta/rc,
         *   R(2,6) = -sin(theta)/bet ~ -theta/bet,
         *   R(5,1) = +sin(theta)/bet ~ +theta/bet,
         * after dropping the thick-lens drift contributions
         * (R(1,2) = rc*sin(theta), R(5,6) = ... etc.) that a thin
         * element does not carry.
         *
         * Sign convention: ImpactX uses pt = -Delta(gamma)/(beta0*gamma0)
         * (opposite sign to MAD-X's PT = +Delta(E)/(p0*c)); the signs of
         * the (2,6) and (5,1) entries above therefore match the sibling
         * Sbend element in ImpactX but differ by an overall sign from
         * the MAD-X thin-dipole formulas.
         *
         * Ideal behavior only: alignment and rotation errors and
         * feed-down from a nonzero reference orbit are ignored; the
         * linearization is taken at (x, y, t) = (0, 0, 0) and pt = 0.
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

            amrex::ParticleReal const beta = refpart.beta();
            amrex::ParticleReal const ds = m_theta * m_rc;     // effective arc length
            amrex::ParticleReal const kx = 1.0_prt / m_rc;     // curvature

            Map6x6 R = Map6x6::Identity();
            R(2,1) = -powi<2>(kx) * ds;      // geometric focusing
            R(2,6) = -kx * ds / beta;        // chromatic dispersion kick
            R(5,1) = +kx * ds / beta;        // path-length coupling to x

            return R;
        }

        amrex::ParticleReal m_theta;  //! dipole bending angle (rad)
        amrex::ParticleReal m_rc;     //! curvature radius (m)

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_beta;  //! beta of the reference particle
        amrex::ParticleReal m_ds;    //! effective (equivalent) arc length
        amrex::ParticleReal m_kx;    //! effective (equivalent) curvature
    };

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::ThinDipole)

#endif // IMPACTX_THINDIPOLE_H
