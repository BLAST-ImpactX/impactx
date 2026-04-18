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

#include "particles/CovarianceMatrix.H"

#include <ablastr/constant.H>

#include <AMReX_Math.H>
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_SmallMatrix.H>

#include <stdexcept>
#include <string>
#include <type_traits>


namespace impactx::elements::mixin
{
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
         * Expects that the reference particle was advanced first.
         *
         * @param[in,out] cm covariance matrix
         * @param[in] ref reference particle
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        void
        operator() (
            Map6x6 & AMREX_RESTRICT cm,
            RefPart const & AMREX_RESTRICT ref
        ) const
        {
            static_assert(
                std::is_base_of_v<LinearTransport, T_Element>,
                "LinearTransport can only be used as a mixin class!"
            );

            // small trick to force every derived class has to implement a method transport_map
            // (w/o using a purely virtual function)
            T_Element const & element = *static_cast<T_Element const*>(this);
            cm = element.transport_map(ref) * cm * element.transport_map(ref).transpose();
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
            Map6x6 & AMREX_RESTRICT cm,
            RefPart const & AMREX_RESTRICT ref
        ) const
        {
            static_assert(
                std::is_base_of_v<LinearTransport, T_Element>,
                "LinearTransport can only be used as a mixin class!"
            );

            T_Element const & element = *static_cast<T_Element const*>(this);
            cm = element.transport_map(ref) * cm * element.transport_map(ref).transpose();
        }
    };

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_LINEAR_TRANSPORT_H
