/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_LINEAR_TRANSPORT_H
#define IMPACTX_ELEMENTS_MIXIN_LINEAR_TRANSPORT_H

#include "alignment.H"
#include "particles/CovarianceMatrix.H"

#include <ablastr/constant.H>

#include <AMReX_Math.H>
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_SmallMatrix.H>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>


namespace impactx::elements::mixin
{
    namespace detail
    {
        template <typename T_Element, typename = void>
        struct has_compute_constants : std::false_type {};

        template <typename T_Element>
        struct has_compute_constants<
            T_Element,
            std::void_t<decltype(std::declval<T_Element &>().compute_constants(std::declval<RefPart const &>()))>
        > : std::true_type {};

        template <typename T_Element, typename = void>
        struct has_particle_push : std::false_type {};

        template <typename T_Element>
        struct has_particle_push<
            T_Element,
            std::void_t<decltype(
                std::declval<T_Element &>()(
                    std::declval<amrex::ParticleReal &>(),
                    std::declval<amrex::ParticleReal &>(),
                    std::declval<amrex::ParticleReal &>(),
                    std::declval<amrex::ParticleReal &>(),
                    std::declval<amrex::ParticleReal &>(),
                    std::declval<amrex::ParticleReal &>(),
                    std::declval<std::uint64_t &>(),
                    std::declval<RefPart const &>()
                )
            )>
        > : std::true_type {};

        /** Push the centroid with the same particle-coordinate contract used in
         *  particle tracking, after the reference particle has been advanced.
         *
         *  This captures constant kicks without changing the shared reference
         *  particle push semantics across tracking modes.
         *
         * @param[in] element the current lattice element
         * @param[in,out] centroid the beam centroid relative to the reference
         * @param[in] ref the already-advanced reference particle
         */
        template <typename T_Element>
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        void
        push_centroid (
            T_Element const & element,
            PhaseSpaceVector & AMREX_RESTRICT centroid,
            RefPart const & AMREX_RESTRICT ref
        )
        {
            if constexpr (has_compute_constants<T_Element>::value && has_particle_push<T_Element>::value)
            {
                T_Element element_copy = element;
                element_copy.compute_constants(ref);

                auto & x = centroid(1);
                auto & px = centroid(2);
                auto & y = centroid(3);
                auto & py = centroid(4);
                auto & t = centroid(5);
                auto & pt = centroid(6);
                std::uint64_t idcpu = 0;

                element_copy(x, y, t, px, py, pt, idcpu, ref);
            }
        }

        /** Advance the reference particle and return the corresponding
         *  lab-frame linear transport map for one element slice.
         *
         *  Element-local contracts are without alignment errors (i.e., ideal):
         *    - @c operator()(RefPart&) advances the ideal reference particle
         *      without alignment errors.
         *    - @c transport_map(ref) returns the ideal/local linear map.
         *
         *  This helper centralizes @c Alignment handling for linear/envelope
         *  workflows by shifting the incoming reference particle into
         *  the element frame, applying the ideal reference push there,
         *  evaluating the ideal/local map at the updated local reference
         *  state, and finally transforming both the reference particle and
         *  linear map back to lab coordinates.
         *
         * @param[in] element the current lattice element
         * @param[in,out] ref the reference particle
         */
        template <typename T_Element>
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        advance_ref_and_transport_map (
            T_Element const & element,
            RefPart & AMREX_RESTRICT ref
        )
        {
            if constexpr (std::is_base_of_v<Alignment, T_Element>)
            {
                RefPart ref_local = element.shift_in(ref);
                element(ref_local);
                auto const R_local = element.transport_map(ref_local);
                ref = element.shift_out(ref_local);
                return element.rotate_map(R_local);
            }
            else
            {
                element(ref);
                return element.transport_map(ref);
            }
        }
    } // namespace detail

    /** Mixin for lattice elements that can be expressed as linear transport maps.
     *
     *  @tparam T_Element   the concrete element type (CRTP).
     *  @tparam Implemented whether the element provides a real linear
     *                      transport map. Defaults to @c true (the element
     *                      must define its own @c transport_map). Elements
     *                      that have not yet implemented a linear map
     *                      inherit @c LinearTransport<E, false>: the
     *                      partial specialization below provides a default
     *                      throwing @c transport_map so linear-optics
     *                      diagnostics and envelope tracking get a uniform,
     *                      self-documenting error. Such elements are still
     *                      free to define their own @c transport_map with a
     *                      more specific message; normal C++ name lookup
     *                      shadows the mixin's default.
     *
     *  Diagnostics use @c std::is_base_of_v<LinearTransport<T, true>, T>
     *  as the compile-time "has a real linear map" trait.
     */
    template <typename T_Element, bool Implemented = true>
    struct LinearTransport
    {
        static constexpr bool has_linear_transport = Implemented;

        /** Linear push of the covariance matrix through an element
         *
         * Advances the reference particle and pushes the covariance
         * matrix through the aligned lab-frame linear transport step.
         *
         * @param[in,out] cm covariance matrix
         * @param[in,out] ref reference particle seen by downstream elements
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        void
        operator() (
            Map6x6 & AMREX_RESTRICT cm,
            RefPart & AMREX_RESTRICT ref
        ) const
        {
            static_assert(
                std::is_base_of_v<LinearTransport, T_Element>,
                "LinearTransport can only be used as a mixin class!"
            );

            // small trick to force every derived class has to implement a method transport_map
            // (w/o using a purely virtual function)
            T_Element const & element = *static_cast<T_Element const*>(this);
            auto const R = detail::advance_ref_and_transport_map(element, ref);
            cm = R * cm * R.transpose();
        }

        /** Linear push of the covariance matrix and centroid through an
         *  element.
         *
         * Advances the reference particle first, matching particle tracking,
         * then propagates the centroid with the particle-coordinate contract
         * and the covariance matrix with the aligned linear transport map.
         *
         * @param[in,out] cm covariance matrix
         * @param[in,out] centroid beam centroid relative to the reference
         * @param[in,out] ref reference particle seen by downstream elements
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        void
        operator() (
            Map6x6 & AMREX_RESTRICT cm,
            PhaseSpaceVector & AMREX_RESTRICT centroid,
            RefPart & AMREX_RESTRICT ref
        ) const
        {
            static_assert(
                std::is_base_of_v<LinearTransport, T_Element>,
                "LinearTransport can only be used as a mixin class!"
            );

            T_Element const & element = *static_cast<T_Element const*>(this);
            auto const R = detail::advance_ref_and_transport_map(element, ref);
            detail::push_centroid(element, centroid, ref);
            cm = R * cm * R.transpose();
        }
    };

    /** Partial specialization for elements that have not yet implemented a
     *  linear transport map.
     *
     *  Provides:
     *    - a default @c transport_map that throws a uniform, self-documenting
     *      error (so neither the element nor its callers need custom error
     *      handling), and
     *    - the same @c operator() as the primary template, routed through
     *      the element's (possibly inherited) @c transport_map, which keeps
     *      downstream code (envelope tracking, @c lattice.transfer_map,
     *      @c map_trace) uniform across "implemented" and "not implemented"
     *      elements.
     *
     *  Elements that want a more specific error message can still provide
     *  their own @c transport_map, which shadows the inherited one via
     *  normal C++ name lookup.
     */
    template <typename T_Element>
    struct LinearTransport<T_Element, false>
    {
        static constexpr bool has_linear_transport = false;

        /** Uniform, self-documenting fallback for elements whose linear
         *  transport map is not yet implemented.
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map ([[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            throw std::runtime_error(
                std::string(T_Element::type)
                + ": Linear transport map is not yet implemented for this element."
            );
        }

        /** Covariance-matrix push; identical to the primary template. Routed
         *  through the element's @c transport_map so that any element-side
         *  override (custom error message) is picked up.
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        void
        operator() (
            [[maybe_unused]] Map6x6 & AMREX_RESTRICT cm,
            [[maybe_unused]] RefPart & AMREX_RESTRICT ref
        ) const
        {
            static_assert(
                std::is_base_of_v<LinearTransport, T_Element>,
                "LinearTransport can only be used as a mixin class!"
            );
            T_Element const & element = *static_cast<T_Element const*>(this);

            (void) element.transport_map(ref);  // expected to throw
        }

        /** Covariance-matrix and centroid push; identical throwing behavior as
         *  the primary overload above.
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        void
        operator() (
            [[maybe_unused]] Map6x6 & AMREX_RESTRICT cm,
            [[maybe_unused]] PhaseSpaceVector & AMREX_RESTRICT centroid,
            [[maybe_unused]] RefPart & AMREX_RESTRICT ref
        ) const
        {
            static_assert(
                std::is_base_of_v<LinearTransport, T_Element>,
                "LinearTransport can only be used as a mixin class!"
            );
            T_Element const & element = *static_cast<T_Element const*>(this);

            (void) element.transport_map(ref);  // expected to throw
        }
    };

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_LINEAR_TRANSPORT_H
