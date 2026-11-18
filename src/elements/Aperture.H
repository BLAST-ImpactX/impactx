/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
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
#include "mixin/gpuelement.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"
#include "mixin/lineartransport.H"
#include "mixin/spintransport.H"
#include "mixin/nofinalize.H"
#include "mixin/repeatedmask.H"

#include <AMReX_Enum.H>
#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

namespace impactx::elements
{
    struct ApertureGPU
    : public mixin::BeamOptic<ApertureGPU>,
      public mixin::LinearTransport<ApertureGPU>,
      public mixin::Thin,
      public mixin::Alignment,
      public mixin::RepeatedMask,
      public mixin::SpinTransport,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "Aperture";
        using PType = ImpactXParticleContainer::ParticleType;

        AMREX_ENUM_IN_CLASS(Shape,
            rectangular,  //! rectangular aperture boundary
            elliptical,   //! elliptical aperture boundary
        );

        AMREX_ENUM_IN_CLASS(Action,
            transmit,  //! transmit particles inside the aperture
            absorb,    //! absorb particles inside the aperture
        );

        static constexpr Shape DEFAULT_shape = Shape::rectangular;
        static constexpr Action DEFAULT_action = Action::transmit;

        /** A thin collimator element that applies a transverse aperture boundary.
         *  Particles outside the boundary are considered lost.
         *
         *  A half-aperture of zero or less disables the boundary in that plane,
         *  the same convention as \see mixin::PipeAperture: each plane is bounded
         *  independently, so bounding a single plane gives a jaw (slit) collimator
         *  and bounding both an iris. With only ``aperture_y`` set, a particle is
         *  lost when ``|y| > aperture_y`` at any ``x``, and the rectangular and
         *  elliptical shapes degenerate to the same slab. Note that with
         *  ``action = absorb`` a disabled plane makes the absorbing domain
         *  unbounded in that direction: disabling both planes absorbs the whole
         *  beam, while for ``action = transmit`` the element has no effect on the beam.
         *
         * @param shape aperture shape
         * @param action specify action of domain (transmit/absorb)
         * @param aperture_x horizontal half-aperture (m), zero or less disables the horizontal boundary
         * @param aperture_y vertical half-aperture (m), zero or less disables the vertical boundary
         * @param repeat_x horizontal period for repeated masking, optional (m),
         *                 @see mixin::RepeatedMask
         * @param repeat_y vertical period for repeated masking, optional (m),
         *                 @see mixin::RepeatedMask
         * @param shift_odd_x for hexagonal/triangular mask patterns: horizontal shift of every
         *                    2nd (odd) vertical period by repeat_x / 2. Requires repeat_x and
         *                    repeat_y. Use alignment offsets dx,dy to move whole mask as needed.
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         */
        ApertureGPU (
            amrex::ParticleReal aperture_x,
            amrex::ParticleReal aperture_y,
            amrex::ParticleReal repeat_x = DEFAULT_repeat_x,
            amrex::ParticleReal repeat_y = DEFAULT_repeat_y,
            bool shift_odd_x = DEFAULT_shift_odd_x,
            Shape shape = DEFAULT_shape,
            Action action = DEFAULT_action,
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree
        )
        : Alignment(dx, dy, rotation_degree),
          RepeatedMask(repeat_x, repeat_y, shift_odd_x),
          m_shape(shape), m_action(action)
        {
            set_aperture_x(aperture_x);
            set_aperture_y(aperture_y);
        }

        /** Reverse the element (no-op: aperture is not invertible for absorbed particles) */
        void reverse () { /* no-op */ }

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
            T_Real const & AMREX_RESTRICT x,
            T_Real const & AMREX_RESTRICT y,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT t,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT px,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT py,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT pt,
            T_IdCpu & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            // if the aperture is periodic, translate the coordinates by whole periods
            // into the period cell around the origin
            T_Real u, v;
            to_period_cell(x, y, u, v);

            // compare against the aperture boundary
            auto inside_aperture = inside_shape(u, v);

            // overlapping shifted rows: the nearest copy of the nearest shifted row
            // can contain the position, too; rectangle and ellipse are convex and
            // symmetric about both axes, so no other copy can, @see mixin::RepeatedMask
            if (shifted_rows_overlap()) {
                T_Real un, vn;
                to_neighbour_row(u, v, un, vn);
                inside_aperture = inside_aperture || inside_shape(un, vn);
            }

            // For transmit (default): invalidate if outside aperture
            // For absorb: invalidate if inside aperture
            auto const invalid_mask = m_action == Action::transmit ? !inside_aperture : inside_aperture;
            amrex::ParticleIDWrapper<T_IdCpu>{idcpu}.make_invalid(invalid_mask);
        }

        /** Is a position inside the aperture shape centred on the origin?
         *
         * @param u horizontal position [m]
         * @param v vertical position [m]
         * @return true inside the boundary
         */
        template<typename T_Real>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        auto inside_shape (T_Real const & u, T_Real const & v) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;
            using VBool = decltype(u>0_prt);

            // scale horizontal and vertical coordinates
            // performance optimization: we intentionally store the inverse of
            //   the aperture, so we can do a quick multiplication (instead of a
            //   costly division) here.
            T_Real const su = u * m_inv_aperture_x;
            T_Real const sv = v * m_inv_aperture_y;

            VBool inside;
            switch (m_shape)
            {
                case Shape::rectangular:  // default
                    inside = (powi<2>(su) <= 1_prt && powi<2>(sv) <= 1_prt);
                    break;

                case Shape::elliptical:
                    inside = (powi<2>(su) + powi<2>(sv) <= 1_prt);
                    break;

                default:
                    inside = VBool{true};
                    break;
            }
            return inside;
        }

        /** This pushes the reference particle. */
        using Thin::operator();

        /** This operator pushes the particle spin.
         *
         * An aperture is field-free, so the spin is unaffected.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param sx particle spin x-component (unused)
         * @param sy particle spin y-component (unused)
         * @param sz particle spin z-component (unused)
         * @param idcpu particle global index
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void spin_and_phasespace_push (
            T_Real const & AMREX_RESTRICT x,
            T_Real const & AMREX_RESTRICT y,
            T_Real const & AMREX_RESTRICT t,
            T_Real const & AMREX_RESTRICT px,
            T_Real const & AMREX_RESTRICT py,
            T_Real const & AMREX_RESTRICT pt,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT sx,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT sy,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT sz,
            T_IdCpu & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            // spin is unaffected

            // phase space push
            (*this)(x, y, t, px, py, pt, idcpu, refpart);
        }

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
         *
         * An aperture (potentially) removes particles; it does not
         * deflect or drift the beam, so the linear transport map is the
         * identity.
         *
         * @param[in] refpart reference particle (unused)
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map ([[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            // apply the transverse rotation (roll) alignment error
            return rotate_aligned_map(Map6x6::Identity());
        }

        /** Maximum horizontal coordinate
         *
         * @return horizontal half-aperture size [m]; zero if the horizontal
         *         boundary is disabled
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal aperture_x () const
        {
            using namespace amrex::literals; // for _prt
            return m_inv_aperture_x > 0_prt ? 1_prt / m_inv_aperture_x : 0_prt;
        }

        /** Maximum vertical coordinate
         *
         * @return vertical half-aperture size [m]; zero if the vertical
         *         boundary is disabled
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal aperture_y () const
        {
            using namespace amrex::literals; // for _prt
            return m_inv_aperture_y > 0_prt ? 1_prt / m_inv_aperture_y : 0_prt;
        }

        /** Set the maximum horizontal coordinate
         *
         * A value of zero or less removes the horizontal boundary, the same
         * convention as \see mixin::PipeAperture::set_aperture_x.
         *
         * @param aperture_x horizontal half-aperture size [m]
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void set_aperture_x (amrex::ParticleReal aperture_x)
        {
            using namespace amrex::literals; // for _prt

            // performance optimization: we intentionally store the inverse of
            //   the aperture, @see m_inv_aperture_x
            m_inv_aperture_x = aperture_x > 0_prt ? 1_prt / aperture_x : 0_prt;
        }

        /** Set the maximum vertical coordinate
         *
         * A value of zero or less removes the vertical boundary, the same
         * convention as \see mixin::PipeAperture::set_aperture_y.
         *
         * @param aperture_y vertical half-aperture size [m]
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void set_aperture_y (amrex::ParticleReal aperture_y)
        {
            using namespace amrex::literals; // for _prt

            // performance optimization: we intentionally store the inverse of
            //   the aperture, @see m_inv_aperture_y
            m_inv_aperture_y = aperture_y > 0_prt ? 1_prt / aperture_y : 0_prt;
        }

        /** Can the shifted rows overlap?
         *
         * With shift_odd_x and aperture_y > repeat_y / 2, a position in the period cell can
         * also lie in the nearest copy of the nearest shifted row,
         * \see mixin::RepeatedMask::to_neighbour_row. A plane without a vertical boundary
         * is open across its whole period instead and does not reach into the next row.
         *
         * @return true if the shifted rows overlap
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool shifted_rows_overlap () const
        {
            using namespace amrex::literals; // for _prt

            // aperture_y > repeat_y / 2, with the stored inverse of aperture_y
            return m_shift_odd_x && m_inv_aperture_y > 0_prt
                && m_inv_aperture_y * m_repeat_y < 2_prt;
        }

        Shape m_shape; //! aperture type (rectangular, elliptical)
        Action m_action; //! action type (transmit, absorb)
        // performance optimization: we intentionally store the inverse of
        //   the aperture, so we can do a quick multiplication instead of a costly
        //   division in the comparison operation per particle.
        //   A stored zero scales the coordinate away and thus disables the
        //   boundary in that plane, @see set_aperture_x
        amrex::ParticleReal m_inv_aperture_x = 0; //! inverse of the horizontal half-aperture [1/m]
        amrex::ParticleReal m_inv_aperture_y = 0; //! inverse of the vertical half-aperture [1/m]
    };

    /** Aperture: GPU state plus common host-only metadata */
    using Aperture = Element<ApertureGPU>;

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::ApertureGPU)

#endif // IMPACTX_APERTURE_H
