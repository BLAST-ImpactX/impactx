/* Copyright 2022 The Regents of the University of California, through Lawrence
*            Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_COMPENSATED_REAL_H
#define IMPACTX_COMPENSATED_REAL_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>


namespace impactx
{
    /** A compensated floating-point type that automatically uses Kahan-Babuska summation.
     *
     * This class wraps a floating-point value and automatically applies the second-order
     * iterative Kahan-Babuska algorithm for all arithmetic operations to maintain high
     * precision when accumulating many small values.
     *
     * The algorithm follows Klein (2006) to provide high-precision arithmetic that avoids
     * floating-point precision errors when dealing with large numbers of small values.
     *
     * References:
     * - https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements
     * - Klein (2006). "A generalized Kahan–Babuška-Summation-Algorithm". in
     *   Computing. 76 (3–4). Springer-Verlag: 279–293. doi:10.1007/s00607-005-0139-x
     *
     * @tparam T Floating-point type (e.g., amrex::Real or amrex::ParticleReal)
     */
    template<typename T = amrex::Real>
    class CompensatedReal
    {
        T m_value = T(0.0);     /**< Main value */
        T m_cs = T(0.0);        /**< First-order compensation for lost low-order bits */
        T m_ccs = T(0.0);       /**< Second-order compensation for further lost bits */

    public:
        /** Default constructor */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        explicit CompensatedReal () = default;

        /** Constructor from value with proper compensation reset
         *
         * @param value Initial value
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        explicit CompensatedReal (T value) : m_value(value) {
            reset_compensation();
        }

        /** Copy constructor
         *
         * @param other Other CompensatedReal to copy from
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        explicit CompensatedReal (CompensatedReal const& other) = default;

        /** Assignment operator
         *
         * @param other Other CompensatedReal to assign from
         * @return Reference to this object
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        CompensatedReal& operator=(CompensatedReal const& other) = default;

        /** Assignment operator from value
         *
         * @param value Value to assign
         * @return Reference to this object
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        CompensatedReal& assign (T value) {
            m_value = value;
            reset_compensation();
            return *this;
        }
        // note: intentionally commented out to avoid usage like
        //   s = s + ds
        // CompensatedReal& operator=(T value) {

        /** Implicit conversion to the underlying type
         *
         * @return The compensated value
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        operator T() const {
            return value();
        }

        /** Get the final compensated value
         *
         * @return The compensated value (main + first-order + second-order compensation)
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T value () const {
            return m_value + m_cs + m_ccs;
        }

        /** Reset compensation variables to zero
         *
         * Call this when starting a new tracking sequence to ensure
         * compensation variables are properly initialized.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void reset_compensation () {
            m_cs = T(0.0);
            m_ccs = T(0.0);
        }

        /** Compensated addition operator
         *
         * Applies the second-order iterative Kahan-Babuška algorithm
         * to add a value with high precision.
         * Automatically promotes other floating-point types to T.
         *
         * @param rhs Value to add (any floating-point type)
         * @return Reference to this object
         */
        template<typename U>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        CompensatedReal& operator+=(U rhs) {
            T const other = static_cast<T>(rhs);

            // First-order Kahan-Babuška step
            T const t_sum = m_value + other;
            T c;
            if (amrex::Math::abs(m_value) >= amrex::Math::abs(other)) {
                c = (m_value - t_sum) + other;
            } else {
                c = (other - t_sum) + m_value;
            }
            m_value = t_sum;

            // Second-order compensation step
            T const t_comp = m_cs + c;
            T cc;
            if (amrex::Math::abs(m_cs) >= amrex::Math::abs(c)) {
                cc = (m_cs - t_comp) + c;
            } else {
                cc = (c - t_comp) + m_cs;
            }
            m_cs = t_comp;
            m_ccs += cc;

            return *this;
        }

        /** Compensated subtraction operator
         *
         * Applies compensated subtraction by adding the negative value.
         * Automatically promotes other floating-point types to T.
         *
         * @param rhs Value to subtract (any floating-point type)
         * @return Reference to this object
         */
        template<typename U>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        CompensatedReal& operator-=(U rhs) {
            return operator+=(-rhs);
        }

        /** Compensated addition with another CompensatedReal
         *
         * @param rhs Other CompensatedReal to add
         * @return Reference to this object
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        CompensatedReal& operator+=(CompensatedReal const& rhs) {
            // Add the compensated value
            T const rhs_value = rhs.value();
            return operator+=(rhs_value);
        }

        /** Compensated subtraction with another CompensatedReal
         *
         * @param rhs Other CompensatedReal to subtract
         * @return Reference to this object
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        CompensatedReal& operator-=(CompensatedReal const& rhs) {
            T const rhs_value = rhs.value();
            return operator+=(-rhs_value);
        }

        /** Addition operator
         *
         * @param rhs Value to add (any floating-point type)
         * @return New CompensatedReal with result
         */
        template<typename U>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        CompensatedReal operator+(U rhs) const {
            CompensatedReal result(*this);
            result += rhs;
            return CompensatedReal(result);
        }

        /** Subtraction operator
         *
         * @param rhs Value to subtract (any floating-point type)
         * @return New CompensatedReal with result
         */
        template<typename U>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        CompensatedReal operator-(U rhs) const {
            CompensatedReal result(*this);
            result -= rhs;
            return CompensatedReal(result);
        }

        /** Addition operator with CompensatedReal
         *
         * @param rhs Other CompensatedReal to add
         * @return New CompensatedReal with result
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        CompensatedReal operator+(CompensatedReal const& rhs) const {
            CompensatedReal result(*this);
            result += rhs;
            return CompensatedReal(result);
        }

        /** Subtraction operator with CompensatedReal
         *
         * @param rhs Other CompensatedReal to subtract
         * @return New CompensatedReal with result
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        CompensatedReal operator-(CompensatedReal const& rhs) const {
            CompensatedReal result(*this);
            result -= rhs;
            return CompensatedReal(result);
        }

        /** Equality comparison
         *
         * @param rhs Value to compare with
         * @return True if values are equal
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool operator==(T rhs) const {
            return value() == rhs;
        }

        /** Inequality comparison
         *
         * @param rhs Value to compare with
         * @return True if values are not equal
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool operator!=(T rhs) const {
            return value() != rhs;
        }

        /** Less than comparison
         *
         * @param rhs Value to compare with
         * @return True if this value is less than rhs
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool operator<(T rhs) const {
            return value() < rhs;
        }

        /** Greater than comparison
         *
         * @param rhs Value to compare with
         * @return True if this value is greater than rhs
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool operator>(T rhs) const {
            return value() > rhs;
        }

        /** Less than or equal comparison
         *
         * @param rhs Value to compare with
         * @return True if this value is less than or equal to rhs
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool operator<=(T rhs) const {
            return value() <= rhs;
        }

        /** Greater than or equal comparison
         *
         * @param rhs Value to compare with
         * @return True if this value is greater than or equal to rhs
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool operator>=(T rhs) const {
            return value() >= rhs;
        }

        /** Equality comparison with CompensatedReal
         *
         * @param rhs Other CompensatedReal to compare with
         * @return True if values are equal
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool operator==(CompensatedReal const& rhs) const {
            return value() == rhs.value();
        }

        /** Inequality comparison with CompensatedReal
         *
         * @param rhs Other CompensatedReal to compare with
         * @return True if values are not equal
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool operator!=(CompensatedReal const& rhs) const {
            return value() != rhs.value();
        }

        /** Less than comparison with CompensatedReal
         *
         * @param rhs Other CompensatedReal to compare with
         * @return True if this value is less than rhs
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool operator<(CompensatedReal const& rhs) const {
            return value() < rhs.value();
        }

        /** Greater than comparison with CompensatedReal
         *
         * @param rhs Other CompensatedReal to compare with
         * @return True if this value is greater than rhs
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool operator>(CompensatedReal const& rhs) const {
            return value() > rhs.value();
        }

        /** Less than or equal comparison with CompensatedReal
         *
         * @param rhs Other CompensatedReal to compare with
         * @return True if this value is less than or equal to rhs
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool operator<=(CompensatedReal const& rhs) const {
            return value() <= rhs.value();
        }

        /** Greater than or equal comparison with CompensatedReal
         *
         * @param rhs Other CompensatedReal to compare with
         * @return True if this value is greater than or equal to rhs
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool operator>=(CompensatedReal const& rhs) const {
            return value() >= rhs.value();
        }
    };

    /** Convenience typedef for typical usage with particle real type */
    using CompensatedParticleReal = CompensatedReal<amrex::ParticleReal>;

    // Free function operators for operand flexibility
    // These allow operations like: 2.5 + cr (where cr is CompensatedReal)
    // and ensure the result is always CompensatedReal, not the underlying type

    /** Free function addition operator
     *
     * @param lhs Left-hand side value
     * @param rhs Right-hand side CompensatedReal
     * @return CompensatedReal result
     */
    template<typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    CompensatedReal<T> operator+(T lhs, CompensatedReal<T> const& rhs) {
        return rhs + lhs;  // Delegate to member operator
    }

    /** Free function subtraction operator
     *
     * @param lhs Left-hand side value
     * @param rhs Right-hand side CompensatedReal
     * @return CompensatedReal result
     */
    template<typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    CompensatedReal<T> operator-(T lhs, CompensatedReal<T> const& rhs) {
        CompensatedReal<T> result(lhs);  // Create CompensatedReal from lhs
        result -= rhs;                   // Subtract rhs
        return CompensatedReal<T>(result);
    }

    /** Free function equality comparison
     *
     * @param lhs Left-hand side value
     * @param rhs Right-hand side CompensatedReal
     * @return True if values are equal
     */
    template<typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator==(T lhs, CompensatedReal<T> const& rhs) {
        return rhs == lhs;  // Delegate to member operator
    }

    /** Free function inequality comparison
     *
     * @param lhs Left-hand side value
     * @param rhs Right-hand side CompensatedReal
     * @return True if values are not equal
     */
    template<typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator!=(T lhs, CompensatedReal<T> const& rhs) {
        return rhs != lhs;
    }

    /** Free function less-than comparison
     *
     * @param lhs Left-hand side value
     * @param rhs Right-hand side CompensatedReal
     * @return True if lhs < rhs
     */
    template<typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator<(T lhs, CompensatedReal<T> const& rhs) {
        return rhs > lhs;  // Note: flipped comparison
    }

    /** Free function greater-than comparison
     *
     * @param lhs Left-hand side value
     * @param rhs Right-hand side CompensatedReal
     * @return True if lhs > rhs
     */
    template<typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator>(T lhs, CompensatedReal<T> const& rhs) {
        return rhs < lhs;
    }

    /** Free function less-than-or-equal comparison
     *
     * @param lhs Left-hand side value
     * @param rhs Right-hand side CompensatedReal
     * @return True if lhs <= rhs
     */
    template<typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator<=(T lhs, CompensatedReal<T> const& rhs) {
        return rhs >= lhs;
    }

    /** Free function greater-than-or-equal comparison
     *
     * @param lhs Left-hand side value
     * @param rhs Right-hand side CompensatedReal
     * @return True if lhs >= rhs
     */
    template<typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator>=(T lhs, CompensatedReal<T> const& rhs) {
        return rhs <= lhs;
    }

} // namespace impactx

#endif // IMPACTX_COMPENSATED_REAL_H
