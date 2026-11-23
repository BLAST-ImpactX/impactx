/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "Handle.H"

#include <memory>
#include <type_traits>
#include <utility>
#include <variant>


namespace impactx::elements
{
    ElementHandle make_handle (KnownElements const & element)
    {
        return std::visit(
            [](auto const & el) {
                return ElementHandle{std::make_shared<std::decay_t<decltype(el)>>(el)};
            },
            element
        );
    }

    ElementHandle make_handle (KnownElements && element)
    {
        return std::visit(
            [](auto && el) {
                return ElementHandle{
                    std::make_shared<std::decay_t<decltype(el)>>(std::move(el))};
            },
            std::move(element)
        );
    }

    ElementHandle copy_element (ElementHandle const & handle)
    {
        require_valid(handle, "copy");
        return std::visit(
            [](auto const & el) {
                return ElementHandle{std::make_shared<std::decay_t<decltype(*el)>>(*el)};
            },
            handle
        );
    }

} // namespace impactx::elements
