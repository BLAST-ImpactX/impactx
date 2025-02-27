/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PROT_H
#define IMPACTX_PROT_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"
#include "mixin/lineartransport.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"

#include <ablastr/constant.H>

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <stdexcept>


namespace impactx::elements
{
    struct PRot
    : public mixin::Named,
      public mixin::BeamOptic<PRot>,
      public mixin::LinearTransport<PRot>,
      public mixin::Thin,
      public mixin::NoFinalize
    {
        static constexpr auto type = "PRot";
        using PType = ImpactXParticleContainer::ParticleType;

        static constexpr amrex::ParticleReal degree2rad = ablastr::constant::math::pi / 180.0;

        /** An exact pole face rotation in the x-z plane, from a frame
         *  in which the reference orbit has angle phi_in with the z-axis,
         *  to a frame in which the reference orbit has angle phi_out
         *  with the z-axis.  The net rotation angle is phi_out - phi_in.
         *
         * @param phi_in Initial angle of reference trajectory w/r/t z (degrees)
         * @param phi_out Final angle of reference trajectory w/r/t/ z (degrees)
         * @param name a user defined and not necessarily unique name of the element
         */
        PRot (
            amrex::ParticleReal phi_in,
            amrex::ParticleReal phi_out,
            std::optional<std::string> name = std::nullopt
        )
        : Named(std::move(name)),
          m_phi_in(phi_in * degree2rad),
          m_phi_out(phi_out * degree2rad)
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

            // access reference particle values to find 1/beta
            m_inv_beta = 1.0_prt / refpart.beta();

            // store rotation angle and initial, final values of pz
            amrex::ParticleReal const theta = m_phi_out - m_phi_in;
            auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);
            m_sin_theta = sin_theta;
            m_cos_theta = cos_theta;
            auto const [sin_phi_in, cos_phi_in] = amrex::Math::sincos(m_phi_in);
            m_sin_phi_in = sin_phi_in;
            m_cos_phi_in = cos_phi_in;
        }

        /** This is a prot functor, so that a variable of this type can be used like a
         *  prot function.
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

            // intialize output values
            amrex::ParticleReal xout = x;
            amrex::ParticleReal yout = y;
            amrex::ParticleReal tout = t;
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal ptout = pt;

            amrex::ParticleReal const pz = std::sqrt(
                1.0_prt -
                2.0_prt * pt * m_inv_beta +
                std::pow(pt, 2) -
                std::pow(py, 2) -
                std::pow(px + m_sin_phi_in, 2)
            );
            amrex::ParticleReal const pzf = pz * m_cos_theta - (px + m_sin_phi_in) * m_sin_theta;

            // advance position and momentum
            xout = x * pz / pzf;
            pxout = px * m_cos_theta + (pz - m_cos_phi_in) * m_sin_theta;

            yout = y + py * x * m_sin_theta / pzf;
            pyout = py;

            tout = t - (pt - m_inv_beta) * x * m_sin_theta / pzf;
            ptout = pt;

            // assign updated values
            x = xout;
            y = yout;
            t = tout;
            px = pxout;
            py = pyout;
            pt = ptout;
        }

        /** This pushes the reference particle. */
        using Thin::operator();

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
            throw std::runtime_error(std::string(type) + ": Envelope tracking is not yet implemented!");
            return Map6x6::Identity();
        }

        amrex::ParticleReal m_phi_in; //! normalized (max) RF voltage drop.
        amrex::ParticleReal m_phi_out; //! RF wavenumber in 1/m.

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_inv_beta, m_sin_theta, m_cos_theta, m_sin_phi_in, m_cos_phi_in;
    };

} // namespace impactx

#endif // IMPACTX_PROT_H
