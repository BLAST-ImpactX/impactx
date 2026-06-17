/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_ALIGNMENT_H
#define IMPACTX_ELEMENTS_MIXIN_ALIGNMENT_H

#include "particles/CovarianceMatrix.H"
#include "particles/ImpactXParticleContainer.H"

#include <ablastr/constant.H>

#include <AMReX_Math.H>
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_SmallMatrix.H>


namespace impactx::elements::mixin
{
    /** This is a helper class for lattice elements with horizontal/vertical alignment errors
     */
    struct Alignment
    {
        static constexpr amrex::ParticleReal degree2rad =
            ablastr::constant::math::pi_v<amrex::ParticleReal> / amrex::ParticleReal{180};

        static constexpr amrex::ParticleReal DEFAULT_dx = 0;
        static constexpr amrex::ParticleReal DEFAULT_dy = 0;
        static constexpr amrex::ParticleReal DEFAULT_rotation_degree = 0;

        /** A finite-length element
         *
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         */
        Alignment (
            amrex::ParticleReal dx,
            amrex::ParticleReal dy,
            amrex::ParticleReal rotation_degree
        )
        : m_dx(dx), m_dy(dy), m_rotation(rotation_degree * degree2rad)
        {
        }

        Alignment () = default;
        Alignment (Alignment const &) = default;
        Alignment& operator= (Alignment const &) = default;
        Alignment (Alignment&&) = default;
        Alignment& operator= (Alignment&& rhs) = default;

        ~Alignment () = default;

        /** Compute and cache the constants for the push.
         *
         * In particular, used to pre-compute and cache variables that are
         * independent of the individually tracked particle.
         *
         * @param refpart reference particle
         */
        void compute_constants ([[maybe_unused]] RefPart const & refpart)
        {
            auto const [sin_rotation, cos_rotation] = amrex::Math::sincos(m_rotation);
            m_sin_rotation = sin_rotation;
            m_cos_rotation = cos_rotation;
        }

        /** Shift the particle into the alignment error frame
         *
         * The @see compute_constants method must be called before pushing particles through this operator.
         *
         * @param[inout] x horizontal position relative to reference particle
         * @param[inout] y vertical position relative to reference particle
         * @param[inout] px horizontal momentum relative to reference particle
         * @param[inout] py vertical momentum relative to reference particle
         */
        template<typename T_Real>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void shift_in (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py
        ) const
        {
            // position
            T_Real const xc = x - m_dx;
            T_Real const yc = y - m_dy;
            x =  xc * m_cos_rotation + yc * m_sin_rotation;
            y = -xc * m_sin_rotation + yc * m_cos_rotation;

            // momentum
            T_Real const pxc = px;
            T_Real const pyc = py;
            px =  pxc * m_cos_rotation + pyc * m_sin_rotation;
            py = -pxc * m_sin_rotation + pyc * m_cos_rotation;
        }

        /** Rotate a transverse vector into the alignment error frame
         *
         * The @see compute_constants method must be called before pushing particles through this operator.
         *
         * @param[inout] sx horizontal spin/vector component
         * @param[inout] sy vertical spin/vector component
         */
        template<typename T_Real>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void rotate_in (
            T_Real & AMREX_RESTRICT sx,
            T_Real & AMREX_RESTRICT sy
        ) const
        {
            T_Real const sxc = sx;
            T_Real const syc = sy;
            sx =  sxc * m_cos_rotation + syc * m_sin_rotation;
            sy = -sxc * m_sin_rotation + syc * m_cos_rotation;
        }

        /** Shift the particle out of the alignment error frame
         *
         * The @see compute_constants method must be called before pushing particles through this operator.
         *
         * @param[inout] x horizontal position relative to reference particle
         * @param[inout] y vertical position relative to reference particle
         * @param[inout] px horizontal momentum relative to reference particle
         * @param[inout] py vertical momentum relative to reference particle
         */
        template<typename T_Real>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void shift_out (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py
        ) const
        {
            // position
            T_Real const xc = x;
            T_Real const yc = y;
            x = xc * m_cos_rotation - yc * m_sin_rotation;
            y = xc * m_sin_rotation + yc * m_cos_rotation;
            x += m_dx;
            y += m_dy;

            // momentum
            T_Real const pxc = px;
            T_Real const pyc = py;
            px = pxc * m_cos_rotation - pyc * m_sin_rotation;
            py = pxc * m_sin_rotation + pyc * m_cos_rotation;
        }

        /** Rotate a transverse vector out of the alignment error frame
         *
         * The @see compute_constants method must be called before pushing particles through this operator.
         *
         * @param[inout] sx horizontal spin/vector component
         * @param[inout] sy vertical spin/vector component
         */
        template<typename T_Real>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void rotate_out (
            T_Real & AMREX_RESTRICT sx,
            T_Real & AMREX_RESTRICT sy
        ) const
        {
            T_Real const sxc = sx;
            T_Real const syc = sy;
            sx = sxc * m_cos_rotation - syc * m_sin_rotation;
            sy = sxc * m_sin_rotation + syc * m_cos_rotation;
        }

        /** Horizontal translation error
         *
         * @return horizontal translation error in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal dx () const
        {
            return m_dx;
        }

        /** Vertical translation error
         *
         * @return vertical translation error in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal dy () const
        {
            return m_dy;
        }

        /** rotation error in the transverse plane
         *
         * @return rotation error in the transverse plane [degrees]
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal rotation () const
        {
            return m_rotation / degree2rad;
        }

        /** Whether the element has any (nonzero) translation or rotation error.
         *
         * Used to unswitch the alignment shift/rotate out of the per-particle push:
         * a perfectly aligned element skips the (otherwise identity) transforms.
         * Cheap to evaluate and called once per push launch (not per particle), so
         * it reads the current alignment errors directly instead of caching.
         *
         * @return true if dx, dy, or rotation is nonzero
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool misaligned () const
        {
            return (m_dx != amrex::ParticleReal(0)) ||
                   (m_dy != amrex::ParticleReal(0)) ||
                   (m_rotation != amrex::ParticleReal(0));
        }

        /** 6x6 linear map of the transverse rotation (roll) alignment error.
         *
         * Returns the lab-from-element rotation @c A = Rot(theta) acting on
         * the (x, y) and (px, py) sub-blocks and the identity on (t, pt). It is
         * the linear-map analogue of @see shift_in / @see shift_out: the
         * aligned (lab-frame) map of an element is @c A * R * A^T .
         *
         * sin/cos are evaluated from @c m_rotation directly rather than the
         * @see compute_constants cache, because the linear-optics path
         * (transfer_map, map_trace, envelope tracking) does not call
         * compute_constants().
         *
         * @return the 6x6 transverse-rotation map in (x, px, y, py, t, pt) order
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        rotation_map () const
        {
            auto const [sin_rotation, cos_rotation] = amrex::Math::sincos(m_rotation);

            Map6x6 A = Map6x6::Identity();
            A(1, 1) =  cos_rotation;  A(1, 3) = -sin_rotation;
            A(3, 1) =  sin_rotation;  A(3, 3) =  cos_rotation;
            A(2, 2) =  cos_rotation;  A(2, 4) = -sin_rotation;
            A(4, 2) =  sin_rotation;  A(4, 4) =  cos_rotation;
            return A;
        }

        /** Apply the transverse rotation (roll) alignment error to an element's
         *  intrinsic linear transport map.
         *
         * Elements call this on their own @c transport_map result just before
         * returning, so that linear-optics consumers (transfer_map, map_trace,
         * envelope tracking) see the same roll that the particle push applies
         * via @see shift_in / @see shift_out. A perfectly roll-aligned element
         * returns @c R unchanged.
         *
         * @param R the element's intrinsic (zero-rotation) 6x6 transport map
         * @return the lab-frame map @c A * R * A^T (see @see rotation_map)
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        rotate_aligned_map (Map6x6 const & R) const
        {
            if (m_rotation == amrex::ParticleReal(0)) { return R; }
            Map6x6 const A = rotation_map();
            return A * R * A.transpose();
        }

        amrex::ParticleReal m_dx = 0; //! horizontal translation error [m]
        amrex::ParticleReal m_dy = 0; //! vertical translation error [m]
        amrex::ParticleReal m_rotation = 0; //! rotation error in the transverse plane [rad]

    private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_sin_rotation;  //! std::sin(m_rotation)
        amrex::ParticleReal m_cos_rotation;  //! std::cos(m_rotation)
    };

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_ALIGNMENT_H
