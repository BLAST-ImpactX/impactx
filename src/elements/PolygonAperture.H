/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Eric G. Stern, Chad Mitchell, Axel Huebl
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
#include <AMReX_SIMD.H>

#include <cmath>
#include <stdexcept>


namespace impactx::elements
{
    struct PolygonAperture
    : public mixin::Named,
      public mixin::BeamOptic<Aperture>,
      public mixin::LinearTransport<Aperture>,
      public mixin::Thin,
      public mixin::Alignment,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "PolygonAperture";
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

        /** A thin collimator element that applies a transverse aperture boundary defined
         *  by points of a polygon vertices.
         *  Particles outside the boundary are considered lost.
         *
         * @param action specify action of domain (transmit/absorb)
         * @param vertices_x array of coordinates of the horizontal vertex positions (m)
         * @param vertices_y array of coordinates of the vertical vertex positions (m)
         * @param min_radius2 a minimum radius**2 below which all particles are transmitted
         * @param repeat_x horizontal period for repeated masking, optional (m)
         * @param repeat_y vertical period for repeated masking, optional (m)
         * @param shift_odd_x for hexagonal/triangular mask patterns: horizontal shift of every 2nd (odd) vertical period by repeat_x / 2. Use alignment offsets dx,dy to move whole mask as needed.
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param name a user defined and not necessarily unique name of the element
         */
        using ParticleInt = std::int64_t;      // for indices/counters
        using ParticleIdT = std::uint64_t;     // for idcpu / global particle id

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
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            [[maybe_unused]] T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py,
            [[maybe_unused]] T_Real & AMREX_RESTRICT pt,
            T_IdCpu & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            namespace stdx = amrex::simd::stdx;
            using amrex::Math::powi;
            using std::fmod;
            using std::abs;
            using amrex::complex, amrex::arg;
            using VBool = decltype(x>0_prt);

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // define intermediate variables
            amrex::ParticleReal const dx = m_repeat_x * 0.5_prt;
            amrex::ParticleReal const dy = m_repeat_y * 0.5_prt;

            // shift every 2nd row by half of m_repeat_x
            T_Real sx = x;
            T_Real const sy = y;
            VBool const shift_x = m_shift_odd_x == false ? VBool{false} : fmod(floor((y+dy)/m_repeat_y), 2) != T_Real{0};
#ifdef AMREX_USE_SIMD
            amrex::simd::stdx::where(shift_x, sx) += dx;
#else
            sx += shift_x ? dx : 0_prt;
#endif

            // if the aperture is periodic, shift sx,sy coordinates to the fundamental domain
            T_Real u = (m_repeat_x == 0_prt) ? sx : (fmod(abs(sx)+dx, m_repeat_x)-dx);
            T_Real v = (m_repeat_y == 0_prt) ? sy : (fmod(abs(sy)+dy, m_repeat_y)-dy);

            // scale horizontal and vertical coordinates
            // performance optimization: we intentionally store the inverse of
            //   the aperture, so we can do a quick multiplication (instead of a
            //   costly division) here.
            u = u * m_inv_aperture_x;
            v = v * m_inv_aperture_y;

            // compare against the aperture boundary
            VBool inside_aperture;

            switch (m_shape)
            {
                case Shape::rectangular:  // default
                    inside_aperture = (powi<2>(u) <= 1_prt && powi<2>(v) <= 1_prt);
                    break;

                case Shape::elliptical:
                    inside_aperture = (powi<2>(u) + powi<2>(v) <= 1_prt);
                    break;

                default:
                    inside_aperture = VBool{true};
                    break;
            }

            // calculate whether point is inside polygon or not
            T_real r2 = powi<2>(u) + powi<2>(v);
            // are we inside the minimum radius?
            if (r2 < min_radius2) {
                inside_aperture = VBool(true); // yes!
            } else{
                // no, we have to do the polygon calculation
                T_real thetasum = 0.0;
                T_complex v1conj, v2;
                // vertices_x and vertices_y must be the same size
                for (ParticleInt i = 0; i <vertices_x.size()-1; ++i) {
                    v1conj = complex(u - vertices_x[i], -(v-vertices_y[i]));
                    v2 = complex(u - vertices_x[i+1], v-vertices_y[i+1]);
                    thetasum += amrex::arg(v2 * v1conj);
                }
                inside_aperture = abs(thetasum/(2*pi)) < 1.0e-12_prt;
            }

            // For transmit (default): invalidate if outside aperture
            // For absorb: invalidate if inside aperture
            auto const invalid_mask = m_action == Action::transmit ? !inside_aperture : inside_aperture;
            amrex::ParticleIDWrapper<T_IdCpu>{idcpu}.make_invalid(invalid_mask);

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
        bool m_shift_odd_x; //! for hexagonal/triangular mask patterns: horizontal shift of every 2nd (odd) vertical period by m_repeat_x / 2. Use Alignment offsets to move whole mask as needed.

    };

} // namespace impactx

#endif // IMPACTX_APERTURE_H
