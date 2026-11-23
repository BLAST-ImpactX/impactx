/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_HELPER_ARRAYVALIDATION_H
#define IMPACTX_ELEMENTS_HELPER_ARRAYVALIDATION_H

#include <AMReX_REAL.H>
#include <AMReX_TrackedVector.H>

#include <limits>
#include <stdexcept>
#include <string>
#include <vector>


namespace impactx::elements
{
    /** Reject a pair of element arrays the element cannot be evaluated with
     *
     * Several elements carry two arrays that only mean something together: cosine and
     * sine coefficients, normal and skew coefficients, horizontal and vertical vertices.
     * They must be equally long, and long enough for the field the element evaluates.
     *
     * @param first the first array
     * @param second the second array, which must have the same length
     * @param what element and array names, used in the error message
     * @param min_size the fewest entries this element can be evaluated with
     * @param max_size the most entries this element can be evaluated with
     */
    inline void
    require_array_pair (
        std::vector<amrex::ParticleReal> const & first,
        std::vector<amrex::ParticleReal> const & second,
        std::string const & what,
        int min_size,
        int max_size = std::numeric_limits<int>::max()
    )
    {
        if (first.size() != second.size())
        {
            throw std::invalid_argument(what + " arrays must have the same length!");
        }
        if (first.size() < static_cast<std::size_t>(min_size))
        {
            throw std::invalid_argument(
                what + " arrays need at least " + std::to_string(min_size) +
                (min_size == 1 ? " entry, got " : " entries, got ") +
                std::to_string(first.size()) + "!");
        }
        if (first.size() > static_cast<std::size_t>(max_size))
        {
            throw std::invalid_argument(
                what + " arrays hold at most " + std::to_string(max_size) +
                " entries, got " + std::to_string(first.size()) + "!");
        }
    }

    /** @see require_array_pair, for arrays the element already owns */
    inline void
    require_array_pair (
        amrex::Gpu::TrackedVector<amrex::ParticleReal> const & first,
        amrex::Gpu::TrackedVector<amrex::ParticleReal> const & second,
        std::string const & what,
        int min_size,
        int max_size = std::numeric_limits<int>::max()
    )
    {
        require_array_pair(first.host_const(), second.host_const(), what, min_size, max_size);
    }

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_HELPER_ARRAYVALIDATION_H
