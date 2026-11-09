/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_REPEATED_MASK_H
#define IMPACTX_ELEMENTS_MIXIN_REPEATED_MASK_H

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>


namespace impactx::elements::mixin
{
    /** A transverse aperture mask that repeats periodically, e.g., a pepperpot
     *
     * The aperture shape is repeated with the horizontal and vertical periods
     * ``repeat_x`` and ``repeat_y``; a period of zero disables the repetition
     * in that plane. Optionally, every 2nd (odd) row of periods is shifted
     * horizontally by ``repeat_x / 2``, for hexagonal/triangular patterns.
     *
     * A particle is translated by whole periods into the period cell centred on
     * the origin, ``|x| < repeat_x / 2`` and ``|y| < repeat_y / 2``, @see to_period_cell,
     * and compared against the one copy of the shape there.
     * A shape that does not fit inside this cell would be cut at its edges, so an
     * element registers how far its shape extends, @see set_shape_extent, and both
     * its shape and the periods are validated against each other.
     */
    struct RepeatedMask
    {
        static constexpr amrex::ParticleReal DEFAULT_repeat_x = 0;
        static constexpr amrex::ParticleReal DEFAULT_repeat_y = 0;
        static constexpr bool DEFAULT_shift_odd_x = false;

        /** A periodically repeated aperture mask
         *
         * @param repeat_x horizontal period [m]; zero disables the repetition
         * @param repeat_y vertical period [m]; zero disables the repetition
         * @param shift_odd_x horizontal shift of every 2nd (odd) vertical period by
         *                    repeat_x / 2; requires both periods
         */
        RepeatedMask (
            amrex::ParticleReal repeat_x,
            amrex::ParticleReal repeat_y,
            bool shift_odd_x
        )
        {
            set_repeat(repeat_x, repeat_y, shift_odd_x);
        }

        /** Horizontal period
         *
         * @return horizontal period [m]; zero if the mask does not repeat horizontally
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal repeat_x () const { return m_repeat_x; }

        /** Vertical period
         *
         * @return vertical period [m]; zero if the mask does not repeat vertically
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal repeat_y () const { return m_repeat_y; }

        /** Is every 2nd (odd) vertical period shifted horizontally by repeat_x / 2? */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool shift_odd_x () const { return m_shift_odd_x; }

        /** Set the horizontal period
         *
         * @param repeat_x horizontal period [m]; zero disables the repetition
         */
        void set_repeat_x (amrex::ParticleReal repeat_x)
        {
            set_repeat(repeat_x, m_repeat_y, m_shift_odd_x);
        }

        /** Set the vertical period
         *
         * @param repeat_y vertical period [m]; zero disables the repetition
         */
        void set_repeat_y (amrex::ParticleReal repeat_y)
        {
            set_repeat(m_repeat_x, repeat_y, m_shift_odd_x);
        }

        /** Shift every 2nd (odd) vertical period horizontally by repeat_x / 2
         *
         * @param shift_odd_x enable the shift; requires both periods
         */
        void set_shift_odd_x (bool shift_odd_x)
        {
            set_repeat(m_repeat_x, m_repeat_y, shift_odd_x);
        }

        /** Set the periods and the odd-row shift at once
         *
         * All three are validated together before any is stored, so a rejected update
         * leaves the mask as it was.
         *
         * @param repeat_x horizontal period [m]; zero disables the repetition
         * @param repeat_y vertical period [m]; zero disables the repetition
         * @param shift_odd_x horizontal shift of every 2nd (odd) vertical period by
         *                    repeat_x / 2; requires both periods
         */
        void set_repeat (
            amrex::ParticleReal repeat_x,
            amrex::ParticleReal repeat_y,
            bool shift_odd_x
        )
        {
            using namespace amrex::literals; // for _prt

            require_valid(repeat_x, repeat_y, shift_odd_x, m_shape_extent_x, m_shape_extent_y);
            m_repeat_x = repeat_x;
            m_repeat_y = repeat_y;
            m_shift_odd_x = shift_odd_x;

            // performance optimization: we intentionally store the inverse of
            //   the periods, @see m_inv_repeat_x
            m_inv_repeat_x = repeat_x > 0_prt ? 1_prt / repeat_x : 0_prt;
            m_inv_repeat_y = repeat_y > 0_prt ? 1_prt / repeat_y : 0_prt;
        }

        /** Set how far the repeated shape extends from the origin
         *
         * The shape must lie inside the period cell centred on the origin in each
         * repeated plane, @see RepeatedMask. It is validated against the current
         * periods before it is stored.
         *
         * @param extent_x largest horizontal distance |x| of the shape from the origin [m]
         * @param extent_y largest vertical distance |y| of the shape from the origin [m]
         */
        void set_shape_extent (
            amrex::ParticleReal extent_x,
            amrex::ParticleReal extent_y
        )
        {
            require_valid(m_repeat_x, m_repeat_y, m_shift_odd_x, extent_x, extent_y);
            m_shape_extent_x = extent_x;
            m_shape_extent_y = extent_y;
        }

        /** Reject periods, odd-row shift and shape extent that do not fit together
         *
         * @param repeat_x horizontal period [m]
         * @param repeat_y vertical period [m]
         * @param shift_odd_x horizontal shift of every 2nd (odd) vertical period
         * @param extent_x largest horizontal distance |x| of the shape from the origin [m]
         * @param extent_y largest vertical distance |y| of the shape from the origin [m]
         */
        static void require_valid (
            amrex::ParticleReal repeat_x,
            amrex::ParticleReal repeat_y,
            bool shift_odd_x,
            amrex::ParticleReal extent_x,
            amrex::ParticleReal extent_y
        )
        {
            using namespace amrex::literals; // for _prt

            if (!(std::isfinite(repeat_x) && repeat_x >= 0_prt)) {
                throw std::invalid_argument("RepeatedMask: repeat_x must be finite and >= 0");
            }
            if (!(std::isfinite(repeat_y) && repeat_y >= 0_prt)) {
                throw std::invalid_argument("RepeatedMask: repeat_y must be finite and >= 0");
            }
            if (shift_odd_x && !(repeat_x > 0_prt && repeat_y > 0_prt)) {
                throw std::invalid_argument(
                    "RepeatedMask: shift_odd_x shifts every 2nd row of periods by repeat_x / 2, "
                    "so it requires repeat_x > 0 and repeat_y > 0");
            }
            require_inside_cell("x", extent_x, repeat_x);
            require_inside_cell("y", extent_y, repeat_y);
        }

        /** Reject a shape that does not fit inside the period cell in one plane
         *
         * @param plane "x" or "y", used in the error message
         * @param extent largest distance of the shape from the origin in this plane [m]
         * @param repeat period in this plane [m]; zero means no repetition
         */
        static void require_inside_cell (
            std::string const & plane,
            amrex::ParticleReal extent,
            amrex::ParticleReal repeat
        )
        {
            using namespace amrex::literals; // for _prt

            if (repeat > 0_prt && !(extent < 0.5_prt * repeat)) {
                std::ostringstream msg;
                msg << "RepeatedMask: the aperture shape reaches |" << plane << "| = " << extent
                    << " m, but a repeated shape must lie inside one period, |" << plane
                    << "| < repeat_" << plane << " / 2 = " << 0.5_prt * repeat << " m";
                throw std::invalid_argument(msg.str());
            }
        }

        /** Translate a position by whole periods into the period cell around the origin
         *
         * If shift_odd_x is set, every 2nd (odd) row of periods is shifted by
         * repeat_x / 2 first. A plane without repetition is passed through.
         *
         * @param[in] x horizontal position [m]
         * @param[in] y vertical position [m]
         * @param[out] u horizontal position in the period cell [m]
         * @param[out] v vertical position in the period cell [m]
         */
        template<typename T_Real>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void to_period_cell (
            T_Real const & AMREX_RESTRICT x,
            T_Real const & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT u,
            T_Real & AMREX_RESTRICT v
        ) const
        {
            using namespace amrex::literals; // for _prt
            using std::floor;
            using VBool = decltype(x > 0_prt);

            // the row of periods, needed for both the vertical fold and the odd-row shift
            T_Real const row = (m_repeat_y == 0_prt) ? T_Real(0_prt)
                : floor((y + 0.5_prt * m_repeat_y) * m_inv_repeat_y);
            v = (m_repeat_y == 0_prt) ? y : y - m_repeat_y * row;

            // shift every 2nd (odd) row by half of m_repeat_x
            // performance optimization: the parity of the integer-valued row without fmod
            T_Real sx = x;
            if (m_shift_odd_x) {
                VBool const odd_row = row != 2_prt * floor(0.5_prt * row);
#ifdef AMREX_USE_SIMD
                amrex::simd::stdx::where(odd_row, sx) += 0.5_prt * m_repeat_x;
#else
                sx += odd_row ? 0.5_prt * m_repeat_x : 0_prt;
#endif
            }

            u = (m_repeat_x == 0_prt) ? sx
                : sx - m_repeat_x * floor((sx + 0.5_prt * m_repeat_x) * m_inv_repeat_x);
        }

        //! horizontal period for repeated masking [m]
        amrex::ParticleReal m_repeat_x = DEFAULT_repeat_x;
        //! vertical period for repeated masking [m]
        amrex::ParticleReal m_repeat_y = DEFAULT_repeat_y;
        //! for hexagonal/triangular mask patterns: horizontal shift of every 2nd (odd)
        //! vertical period by m_repeat_x / 2. Use Alignment offsets to move whole mask as needed.
        bool m_shift_odd_x = DEFAULT_shift_odd_x;

        // performance optimization: we intentionally store the inverse of
        //   the periods, so we can do a quick multiplication instead of a costly
        //   division per particle, @see to_period_cell
        //! inverse of the horizontal period [1/m]; zero without repetition
        amrex::ParticleReal m_inv_repeat_x = 0;
        //! inverse of the vertical period [1/m]; zero without repetition
        amrex::ParticleReal m_inv_repeat_y = 0;

        //! largest horizontal distance |x| of the repeated shape from the origin [m];
        //! zero if the shape is not limited to one period
        amrex::ParticleReal m_shape_extent_x = 0;
        //! largest vertical distance |y| of the repeated shape from the origin [m];
        //! zero if the shape is not limited to one period
        amrex::ParticleReal m_shape_extent_y = 0;
    };

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_REPEATED_MASK_H
