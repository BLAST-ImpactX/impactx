/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_ACCESSORS_H
#define IMPACTX_ELEMENTS_MIXIN_ACCESSORS_H

#include "traits.H"

#include "elements/All.H"

#include <AMReX_BLassert.H>
#include <AMReX_REAL.H>

#include <string>
#include <type_traits>
#include <variant>


namespace impactx::elements
{
    /** @file
     *
     * Property access on the ``KnownElements`` variant.
     *
     * These free functions wrap the ``std::visit`` dispatch that would otherwise be spelled
     * out at every call site that only wants to read or write one element property. Getters
     * and setters are overloaded on arity and named after the property itself, mirroring the
     * element member functions: ``ds(element)`` reads, ``ds(element, 0.1)`` writes.
     *
     * Since ``KnownElements`` is a ``std::variant`` over types in this namespace, ADL finds
     * these unqualified; the explicit ``elements::`` spelling is preferred for clarity.
     *
     * Getters are total: every lattice element provides ``ds()``, ``nslice()``, ``name()``,
     * ``has_name()`` and a static ``type``. Setters are partial and assert on elements that
     * do not support them, @see mixin::has_settable_ds_v, mixin::has_settable_nslice_v and
     * mixin::is_named_v for the exact sets.
     *
     * @note This header includes ``elements/All.H`` and therefore sits *above* the element
     *       definitions. Do not include it from an element header or from another header in
     *       ``mixin/``. Those need the type-level traits instead, @see traits.H.
     */

    /** Return the segment length of an element
     *
     * @param element_variant the element to query
     * @return value in meters, zero for zero-length elements
     */
    inline amrex::ParticleReal
    ds (KnownElements const & element_variant)
    {
        return std::visit([](auto const & element) { return element.ds(); }, element_variant);
    }

    /** Set the segment length of an element
     *
     * @param[in,out] element_variant the element to modify
     * @param ds segment length in meters
     */
    inline void
    ds (KnownElements & element_variant, amrex::ParticleReal ds)
    {
        std::visit([ds](auto & element) {
            if constexpr (mixin::has_settable_ds_v<decltype(element)>) {
                element.m_ds = ds;
            } else {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
                    "elements::ds: this element has no settable segment length.");
            }
        }, element_variant);
    }

    /** Return the number of slices used for the application of collective effects
     *
     * @param element_variant the element to query
     * @return positive integer, one for elements that cannot be sliced
     */
    inline int
    nslice (KnownElements const & element_variant)
    {
        return std::visit([](auto const & element) { return element.nslice(); }, element_variant);
    }

    /** Set the number of slices used for the application of collective effects
     *
     * This subdivides the element's transport: the element then applies ``ds/nslice`` per
     * push. @see ScopedNslice to do this temporarily.
     *
     * @param[in,out] element_variant the element to modify
     * @param nslice number of slices, a positive integer
     */
    inline void
    nslice (KnownElements & element_variant, int nslice)
    {
        std::visit([nslice](auto & element) {
            if constexpr (mixin::has_settable_nslice_v<decltype(element)>) {
                element.m_nslice = nslice;
            } else {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
                    "elements::nslice: this element cannot be sliced.");
            }
        }, element_variant);
    }

    /** Return the length of a single slice of an element
     *
     * @param element_variant the element to query
     * @return ``ds/nslice`` in meters, zero for zero-length elements
     */
    inline amrex::ParticleReal
    slice_ds (KnownElements const & element_variant)
    {
        return std::visit([](auto const & element) {
            return element.ds() / element.nslice();
        }, element_variant);
    }

    /** Return the name of an element
     *
     * @param element_variant the element to query
     * @return the element name, or an empty string if it has none
     */
    inline std::string
    name (KnownElements const & element_variant)
    {
        return std::visit([](auto const & element) { return element.name(); }, element_variant);
    }

    /** Set the name of an element
     *
     * @param[in,out] element_variant the element to modify
     * @param name the new element name
     */
    inline void
    name (KnownElements & element_variant, std::string const & name)
    {
        std::visit([&name](auto & element) {
            if constexpr (mixin::is_named_v<decltype(element)>) {
                element.set_name(name);
            } else {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
                    "elements::name: this element has no settable name.");
            }
        }, element_variant);
    }

    /** Does this element carry a name?
     *
     * @param element_variant the element to query
     * @return true if the element has a name
     */
    inline bool
    has_name (KnownElements const & element_variant)
    {
        return std::visit([](auto const & element) { return element.has_name(); }, element_variant);
    }

    /** Return the element type, e.g. "quad"
     *
     * @param element_variant the element to query
     * @return the element's static type name
     */
    inline std::string
    type (KnownElements const & element_variant)
    {
        return std::visit([](auto const & element) {
            return std::string(std::decay_t<decltype(element)>::type);
        }, element_variant);
    }

    /** Does this element have a finite length? (@see mixin::Thick)
     *
     * @param element_variant the element to query
     * @return true for finite-length elements
     */
    inline bool
    is_thick (KnownElements const & element_variant)
    {
        return std::visit([](auto const & element) {
            return mixin::is_thick_v<decltype(element)>;
        }, element_variant);
    }

    /** Does this element have zero length? (@see mixin::Thin)
     *
     * @param element_variant the element to query
     * @return true for zero-length elements
     */
    inline bool
    is_thin (KnownElements const & element_variant)
    {
        return std::visit([](auto const & element) {
            return mixin::is_thin_v<decltype(element)>;
        }, element_variant);
    }

    /** Can this element be subdivided by overriding its slice count?
     *
     * True for the elements that derive their per-push transport length from their own
     * ``nslice()``, so that raising the slice count shortens one push. Note that this is
     * narrower than "has a slice count": @see mixin::can_slice_v for which elements are
     * excluded and why.
     *
     * @param element_variant the element to query
     * @return true if the element can be subdivided, @see ScopedNslice
     */
    inline bool
    can_slice (KnownElements const & element_variant)
    {
        return std::visit([](auto const & element) {
            return mixin::can_slice_v<decltype(element)>;
        }, element_variant);
    }

    /** Reverse an element in place
     *
     * After this call, pushing particles through the element reverses the effect of the
     * original element. This is an involution.
     *
     * @param[in,out] element_variant the element to reverse
     */
    inline void
    reverse (KnownElements & element_variant)
    {
        std::visit([](auto & element) { element.reverse(); }, element_variant);
    }

    /** Close and deallocate all data and handles of an element
     *
     * @param[in,out] element_variant the element to finalize
     */
    inline void
    finalize (KnownElements & element_variant)
    {
        std::visit([](auto & element) { element.finalize(); }, element_variant);
    }

    /** Temporarily override an element's slice count, restoring it on scope exit
     *
     * Subdividing an element is how a caller shortens its transport: with the slice count
     * temporarily doubled, each push advances ``ds/(2*nslice)`` instead of ``ds/nslice``.
     * Restoring on scope exit keeps the lattice consistent even if the pushes in between
     * throw, e.g. from a Programmable element's Python callback.
     *
     * The element must be sliceable, @see can_slice.
     */
    class ScopedNslice
    {
    public:
        /** Override an element's slice count for the lifetime of this object
         *
         * @param[in,out] element_variant the element to subdivide
         * @param nslice the temporary number of slices, a positive integer
         */
        ScopedNslice (KnownElements & element_variant, int nslice)
        : m_element_variant(element_variant),
          m_nslice_original(elements::nslice(element_variant))
        {
            elements::nslice(m_element_variant, nslice);
        }

        ~ScopedNslice ()
        {
            elements::nslice(m_element_variant, m_nslice_original);
        }

        ScopedNslice (ScopedNslice const &) = delete;
        ScopedNslice & operator= (ScopedNslice const &) = delete;
        ScopedNslice (ScopedNslice &&) = delete;
        ScopedNslice & operator= (ScopedNslice &&) = delete;

    private:
        KnownElements & m_element_variant; //! the subdivided element
        int m_nslice_original; //! slice count to restore on scope exit
    };

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_MIXIN_ACCESSORS_H
