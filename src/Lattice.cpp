/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "Lattice.H"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>


namespace impactx
{
    Lattice::Lattice (Lattice const & other)
      : m_elements(other.m_elements), m_generation(other.m_generation)
    {
    }

    Lattice & Lattice::operator= (Lattice const & other)
    {
        if (this != &other)
        {
            require_no_traversal("replace the lattice");
            m_elements = other.m_elements;
            ++m_generation;
        }
        return *this;
    }

    Lattice::Lattice (Lattice && other) noexcept
      : m_elements(std::move(other.m_elements)), m_generation(other.m_generation)
    {
    }

    Lattice & Lattice::operator= (Lattice && other) noexcept(false)
    {
        if (this != &other)
        {
            require_no_traversal("replace the lattice");
            m_elements = std::move(other.m_elements);
            ++m_generation;
        }
        return *this;
    }

    void Lattice::push_back (elements::ElementHandle handle)
    {
        elements::require_valid(handle, "insert");
        require_no_traversal("add an element");
        m_elements.push_back(std::move(handle));
        ++m_generation;
    }

    void Lattice::insert (size_type i, elements::ElementHandle handle)
    {
        elements::require_valid(handle, "insert");
        require_insertion_point(i);
        require_no_traversal("add an element");
        m_elements.insert(m_elements.begin() + static_cast<std::ptrdiff_t>(i), std::move(handle));
        ++m_generation;
    }

    void Lattice::replace (size_type i, elements::ElementHandle handle)
    {
        elements::require_valid(handle, "replace");
        require_element_at(i, "replace");
        require_no_traversal("replace an element");
        m_elements[i] = std::move(handle);
        ++m_generation;
    }

    void Lattice::erase (size_type i)
    {
        require_element_at(i, "remove");
        require_no_traversal("remove an element");
        m_elements.erase(m_elements.begin() + static_cast<std::ptrdiff_t>(i));
        ++m_generation;
    }

    void Lattice::pop_back ()
    {
        require_not_empty("remove an element");
        require_no_traversal("remove an element");
        m_elements.pop_back();
        ++m_generation;
    }

    void Lattice::clear ()
    {
        require_no_traversal("clear the lattice");
        m_elements.clear();
        ++m_generation;
    }

    void Lattice::assign (storage_type elements)
    {
        for (auto const & handle : elements)
        {
            elements::require_valid(handle, "replace the lattice");
        }
        require_no_traversal("replace the lattice");
        m_elements = std::move(elements);
        ++m_generation;
    }

    void Lattice::require_element_at (size_type i, char const * what) const
    {
        if (i >= m_elements.size())
        {
            throw std::out_of_range(
                std::string("Lattice: cannot ") + what + " position " + std::to_string(i) +
                ": the lattice holds " + std::to_string(m_elements.size()) + " elements.");
        }
    }

    void Lattice::require_insertion_point (size_type i) const
    {
        if (i > m_elements.size())
        {
            throw std::out_of_range(
                "Lattice: cannot insert at position " + std::to_string(i) +
                ": the lattice holds " + std::to_string(m_elements.size()) + " elements.");
        }
    }

    void Lattice::require_not_empty (char const * what) const
    {
        if (m_elements.empty())
        {
            throw std::out_of_range(
                std::string("Lattice: cannot ") + what + " an empty lattice.");
        }
    }

    void Lattice::require_no_traversal (char const * what) const
    {
        if (m_active_traversals > 0)
        {
            throw std::runtime_error(
                std::string("Lattice: cannot ") + what + " while tracking through it. "
                "Element parameters may be changed during tracking, but the sequence of "
                "elements may not.");
        }
    }

} // namespace impactx
