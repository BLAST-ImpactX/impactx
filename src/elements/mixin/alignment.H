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


namespace impactx::elements::mixin
{
    /** This is a helper class for lattice elements with horizontal/vertical alignment errors
     */
    struct Alignment
    {
        static constexpr amrex::ParticleReal degree2rad = ablastr::constant::math::pi / 180.0;

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

        /** Shift a reference particle into the alignment error frame.
         *
         * This overload is stateless and computes the rotation factors
         * locally so it can be used from transport-map code paths that do
         * not call @see compute_constants.
         *
         * @param[in] refpart reference particle in lab coordinates
         * @returns reference particle in element-frame coordinates
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        RefPart
        shift_in (RefPart const & refpart) const
        {
            auto const [sin_rotation, cos_rotation] = amrex::Math::sincos(m_rotation);

            RefPart ref_local = refpart;

            // position
            amrex::ParticleReal const xc = refpart.x - m_dx;
            amrex::ParticleReal const yc = refpart.y - m_dy;
            ref_local.x =  xc * cos_rotation + yc * sin_rotation;
            ref_local.y = -xc * sin_rotation + yc * cos_rotation;

            // momentum
            amrex::ParticleReal const pxc = refpart.px;
            amrex::ParticleReal const pyc = refpart.py;
            ref_local.px =  pxc * cos_rotation + pyc * sin_rotation;
            ref_local.py = -pxc * sin_rotation + pyc * cos_rotation;

            return ref_local;
        }

        /** Shift a reference particle out of the alignment error frame.
         *
         * This overload is stateless and computes the rotation factors
         * locally so it can be used from transport-map code paths that do
         * not call @see compute_constants.
         *
         * @param[in] refpart reference particle in element-frame coordinates
         * @returns reference particle in lab coordinates
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        RefPart
        shift_out (RefPart const & refpart) const
        {
            auto const [sin_rotation, cos_rotation] = amrex::Math::sincos(m_rotation);

            RefPart ref_lab = refpart;

            // position
            ref_lab.x = refpart.x * cos_rotation - refpart.y * sin_rotation + m_dx;
            ref_lab.y = refpart.x * sin_rotation + refpart.y * cos_rotation + m_dy;

            // momentum
            ref_lab.px = refpart.px * cos_rotation - refpart.py * sin_rotation;
            ref_lab.py = refpart.px * sin_rotation + refpart.py * cos_rotation;

            return ref_lab;
        }

        /** Rotate an element-frame transport map into lab coordinates.
         *
         * The returned map is A^{-1} R A, where A is the linear lab->local
         * transverse rotation acting on (x, px, y, py). The longitudinal
         * coordinates (t, pt) are left untouched.
         *
         * @param[in] R_local transport map in the element frame
         * @returns transport map in lab coordinates
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        rotate_map (Map6x6 const & R_local) const
        {
            if (m_rotation == amrex::ParticleReal(0)) {
                return R_local;
            }

            auto const [sin_rotation, cos_rotation] = amrex::Math::sincos(m_rotation);

            Map6x6 A = Map6x6::Identity();
            A(1,1) = +cos_rotation;
            A(1,3) = +sin_rotation;
            A(2,2) = +cos_rotation;
            A(2,4) = +sin_rotation;
            A(3,1) = -sin_rotation;
            A(3,3) = +cos_rotation;
            A(4,2) = -sin_rotation;
            A(4,4) = +cos_rotation;

            Map6x6 Ainv = Map6x6::Identity();
            Ainv(1,1) = +cos_rotation;
            Ainv(1,3) = -sin_rotation;
            Ainv(2,2) = +cos_rotation;
            Ainv(2,4) = -sin_rotation;
            Ainv(3,1) = +sin_rotation;
            Ainv(3,3) = +cos_rotation;
            Ainv(4,2) = +sin_rotation;
            Ainv(4,4) = +cos_rotation;

            return Ainv * R_local * A;
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
