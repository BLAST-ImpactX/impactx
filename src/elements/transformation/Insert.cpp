/* Copyright 2022-2025 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "elements/transformation/Insert.H"

#include <stdexcept>
#include <type_traits>
#include <variant>


namespace impactx::elements::transformation
{
    std::list<elements::KnownElements>
    insert_element_every_ds (
        std::list<elements::KnownElements> list,
        amrex::ParticleReal ds,
        elements::KnownElements element
    )
    {
        // algorithm below is so far only implemented for thin elements to insert
        amrex::ParticleReal new_element_ds = 0.0;  // in meters
        std::visit([&new_element_ds](auto &&new_element)
        {
            new_element_ds = new_element.ds();
        }, element);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            new_element_ds == 0,
            "insert_element_ever_s: Only thin elements are supported."
        );

        std::list<elements::KnownElements> new_list;

        amrex::ParticleReal s = 0.0;  // in meters   // TODO: if we can avoid a global s, we can avoid wasting significant digits for long lattices
        amrex::ParticleReal s_next_insert = ds;  // in meters

        while (!list.empty())
        {
            // copy out front element
            elements::KnownElements cur_element_variant = list.front();
            list.pop_front();

            // check where the current element ends
            amrex::ParticleReal cur_s_out;  // in meters
            std::visit([&s, &cur_s_out](auto &&cur_element)
            {
                cur_s_out = s + cur_element.ds();
            }, cur_element_variant);

            // case 1: current element is thick and ends after next insert
            if (s_next_insert < cur_s_out)
            {
                amrex::ParticleReal const s_rel_insert = s_next_insert - s;

                // split element and shorten each part
                elements::KnownElements cur_element_leftover = cur_element_variant;
                std::visit([&s_rel_insert](auto &&cur_element)
                {
                    if constexpr(std::is_base_of_v<elements::mixin::Thin, std::decay_t<decltype(cur_element)>>)
                    {
                        throw std::runtime_error("insert_element_ever_s: Thin element cannot be split.");
                    }
                    else {
                        cur_element.m_ds = s_rel_insert;
                    }
                }, cur_element_variant);
                std::visit([&s_rel_insert](auto &&cur_element_left)
                {
                    if constexpr(std::is_base_of_v<elements::mixin::Thin, std::decay_t<decltype(cur_element_left)>>)
                    {
                        throw std::runtime_error("insert_element_ever_s: Thin element cannot be split.");
                    }
                    else {
                        cur_element_left.m_ds -= s_rel_insert;
                        cur_element_left.set_name(cur_element_left.name() + "_leftover");
                    }
                }, cur_element_leftover);

                // insert element in between
                new_list.push_back(cur_element_variant);
                new_list.push_back(element);

                // add leftover element to front of old list
                list.push_front(cur_element_leftover);

                s += s_rel_insert;
                s_next_insert += ds;
            }
            // case 2: current element ends exactly with next insert
            else if (s_next_insert == cur_s_out) {
                // copy current element
                new_list.push_back(cur_element_variant);
                // insert element
                new_list.push_back(element);

                s = cur_s_out;
                s_next_insert += ds;
            }
            // case 3: current element ends before next insert
            else {
                // thin element or element too thin to slice in ds
                new_list.push_back(cur_element_variant);

                s = cur_s_out;
                // unchanged: s_next_insert
            }
        }

        return new_list;
    }

} // namespace impactx::elements::transformation
