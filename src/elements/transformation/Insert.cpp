/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "elements/transformation/Insert.H"

#include "elements/mixin/accessors.H"

#include <stdexcept>


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
        double const new_element_ds = elements::ds(element);  // in meters
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            new_element_ds == 0,
            "insert_element_ever_s: Only thin elements are supported."
        );

        std::list<elements::KnownElements> new_list;

        double s = 0.0;  // in meters   // TODO: if we can avoid a global s, we can avoid wasting significant digits for long lattices
        double s_next_insert = ds;  // in meters

        while (!list.empty())
        {
            // copy out front element
            elements::KnownElements cur_element_variant = list.front();
            list.pop_front();

            // check where the current element ends
            double const cur_s_out = s + elements::ds(cur_element_variant);  // in meters

            // case 1: current element is thick and ends after next insert
            if (s_next_insert < cur_s_out)
            {
                double const s_rel_insert = s_next_insert - s;

                if (elements::is_thin(cur_element_variant))
                {
                    throw std::runtime_error("insert_element_ever_s: Thin element cannot be split.");
                }

                // split element and shorten each part
                elements::KnownElements cur_element_leftover = cur_element_variant;
                elements::ds(cur_element_variant, static_cast<amrex::ParticleReal>(s_rel_insert));
                elements::ds(
                    cur_element_leftover,
                    elements::ds(cur_element_leftover) - static_cast<amrex::ParticleReal>(s_rel_insert)
                );
                elements::name(cur_element_leftover, elements::name(cur_element_leftover) + "_leftover");

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
