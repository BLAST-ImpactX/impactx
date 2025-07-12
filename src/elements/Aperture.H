/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_APERTURE_H
#define IMPACTX_APERTURE_H


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

#include <cmath>
#include <stdexcept>


namespace impactx::elements
{
    struct Aperture
    : public mixin::Named,
      public mixin::BeamOptic<Aperture>,
      public mixin::LinearTransport<Aperture>,
      public mixin::Thin,
      public mixin::Alignment,
      public mixin::NoFinalize
    {
        static constexpr auto type = "Aperture";
        using PType = ImpactXParticleContainer::ParticleType;

        // TODO: make AMREX_ENUM and simplify @see shape_name implementation with it
        enum Shape
        {
            rectangular,
            elliptical
        };
        enum Action
        {
            transmit,
            absorb
        };

        static std::string
        shape_name (Shape const & shape);

        static std::string
        action_name (Action const & action);

        /** A thin collimator element that applies a transverse aperture boundary.
         *  Particles outside the boundary are considered lost.
         *
         * @param shape aperture shape
         * @param action specify action of domain (transmit/absorb)
         * @param aperture_x horizontal half-aperture (m)
         * @param aperture_y vertical half-aperture (m)
         * @param repeat_x horizontal period for repeated masking, optional (m)
         * @param repeat_y vertical period for repeated masking, optional (m)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param name a user defined and not necessarily unique name of the element
         */
        Aperture (
            amrex::ParticleReal aperture_x,
            amrex::ParticleReal aperture_y,
            amrex::ParticleReal repeat_x,
            amrex::ParticleReal repeat_y,
            Shape shape,
            Action action,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            std::optional<std::string> name = std::nullopt
        )
        : Named(std::move(name)),
          Alignment(dx, dy, rotation_degree),
          m_shape(shape), m_action(action), m_repeat_x(repeat_x), m_repeat_y(repeat_y)
        {
            using namespace amrex::literals; // for _rt and _prt

            if (aperture_x > 0_prt) {
                m_inv_aperture_x = 1_prt / aperture_x;
            } else {
                throw std::runtime_error("Aperture: aperture_x must be > 0!");
            }

            if (aperture_y > 0_prt) {
                m_inv_aperture_y = 1_prt / aperture_y;
            } else {
                throw std::runtime_error("Aperture: aperture_y must be > 0!");
            }
        }

        /** Push all particles */
        using BeamOptic::operator();

        /** Compute and cache the constants for the push.
         *
         * In particular, used to pre-compute and cache variables that are
         * independent of the individually tracked particle.
         *
         * @param refpart reference particle (unused)
         */
        void compute_constants ([[maybe_unused]] RefPart const & refpart)
        {
            Alignment::compute_constants(refpart);
        }

        /** This is an aperture functor, so that a variable of this type can be used like an
         *  aperture function.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t (unused)
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t (unused)
         * @param idcpu particle global index
         * @param refpart reference particle (unused)
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            amrex::ParticleReal & AMREX_RESTRICT x,
            amrex::ParticleReal & AMREX_RESTRICT y,
            [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT t,
            amrex::ParticleReal & AMREX_RESTRICT px,
            amrex::ParticleReal & AMREX_RESTRICT py,
            [[maybe_unused]] amrex::ParticleReal & AMREX_RESTRICT pt,
            uint64_t & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // define intermediate variables
            amrex::ParticleReal u;
            amrex::ParticleReal v;
            amrex::ParticleReal dx = m_repeat_x * 0.5_prt;
            amrex::ParticleReal dy = m_repeat_y * 0.5_prt;

            // if the aperture is periodic, shift x,y coordinates to the fundamental domain
            u = (m_repeat_x==0.0) ? x : (std::fmod(std::abs(x)+dx,m_repeat_x)-dx);
            v = (m_repeat_y==0.0) ? y : (std::fmod(std::abs(y)+dy,m_repeat_y)-dy);

            // scale horizontal and vertical coordinates
            // performance optimization: we intentionally store the inverse of
            //   the aperture, so we can do a quick multiplication (instead of a
            //   costly division) here.
            u = u * m_inv_aperture_x;
            v = v * m_inv_aperture_y;

            // compare against the aperture boundary
            bool inside_aperture = true;

            switch (m_shape)
            {
                case Shape::rectangular:  // default
                    inside_aperture = (powi<2>(u) <= 1_prt && powi<2>(v) <= 1_prt);
                    break;

                case Shape::elliptical:
                    inside_aperture = (powi<2>(u) + powi<2>(v) <= 1_prt);
                    break;
            }

            // For transmit (default): invalidate if outside aperture
            // For absorb: invalidate if inside aperture
            if ((m_action == Action::transmit && !inside_aperture) ||
                (m_action == Action::absorb && inside_aperture)) {
                amrex::ParticleIDWrapper{idcpu}.make_invalid();
            }

            // undo shift due to alignment errors of the element
            shift_out(x, y, px, py);
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

        Shape m_shape; //! aperture type (rectangular, elliptical)
        Action m_action; //! action type (transmit, absorb)
        // performance optimization: we intentionally store the inverse of
        //   the aperture, so we can do a quick multiplication instead of a costly
        //   division in the comparison operation per particle
        amrex::ParticleReal m_inv_aperture_x; //! maximum horizontal coordinate
        amrex::ParticleReal m_inv_aperture_y; //! maximum vertical coordinate
        amrex::ParticleReal m_repeat_x; //! horizontal period for repeated masking
        amrex::ParticleReal m_repeat_y; //! vertical period for repeated masking

    };

} // namespace impactx

#endif // IMPACTX_APERTURE_H
