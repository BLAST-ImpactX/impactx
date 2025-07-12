/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
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


namespace impactx::elements::mixin
{
    /** This is a helper class for lattice elements that can be expressed as linear transport maps.
     */
    template<typename T_Element>
    struct LinearTransport
    {
        LinearTransport () = default;
        LinearTransport (LinearTransport const &) = default;
        LinearTransport& operator= (LinearTransport const &) = default;
        LinearTransport (LinearTransport&&) = default;
        LinearTransport& operator= (LinearTransport&& rhs) = default;

        ~LinearTransport () = default;

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

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_LINEAR_TRANSPORT_H
