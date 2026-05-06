/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_DYNAMICDATA_H
#define IMPACTX_ELEMENTS_MIXIN_DYNAMICDATA_H

#include "impactx_export.H"

#include <AMReX_Gpu.H>

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>


namespace impactx::elements::mixin
{

    /** Registry for element data that outlives value-copied elements.
     *
     * Beamline elements are stored by value in std::variant and may be copied.
     * E.g., element data coefficient arrays live here in shared_ptr so all copies
     * share data.
     *
     * Example (in MyElement.H):
     * @code
     * struct MyCoefficients { ... };
     * using MyElementData = GPUDataRegistry<MyCoefficients>;
     * @endcode
     *
     * @warning This template holds its registry map and ID counter in exported
     *          static data members (`s_next_id`, `s_registry`) declared inside
     *          the class and defined out-of-line below. To guarantee a single
     *          shared registry across every translation unit and shared
     *          library, including downstream DSOs built with
     *          `-fvisibility=hidden` and/or LTO (e.g. the Python bindings
     *          module), you MUST add explicit instantiation control:
     *          - In the element header: `IMPACTX_GPUDATA_EXTERN(impactx::elements::MyElement)`
     *          - In the element .cpp:   `IMPACTX_GPUDATA_INSTANTIATE(impactx::elements::MyElement)`
     *
     * @tparam T element data struct with to_device() method
     */
    template<typename T>
    struct GPUDataRegistry
    {
        using DataType = T;  //!< The coefficient/data type stored in this registry

        // Storage for the registry's state.
        //
        // Declared here as plain `static` members and defined out-of-line below
        // (rather than `inline static ... = init;`) for two reasons:
        //
        // 1. NVCC 12.2 + C++20 evaluates the default-initializer of an
        //    `inline static` data member during implicit instantiation of the
        //    enclosing class template (any reference such as
        //    `DynamicData::get(id)` from an element header triggers it). With
        //    a `std::map<int, std::shared_ptr<T>>` value type that fails as
        //    "incomplete type is not allowed" on the brace-initializer. A
        //    declaration with no initializer sidesteps this entirely; the
        //    out-of-line definitions are only processed when an explicit
        //    instantiation actually requests them.
        // 2. Combined with `extern template struct` in element headers and
        //    `template struct` in the corresponding .cpp, this still provides
        //    exactly one strong definition per element type. With
        //    `IMPACTX_EXPORT` the resulting symbol has default visibility, so
        //    downstream shared libraries built with `-fvisibility=hidden`
        //    and/or LTO (e.g. the Python bindings module) resolve to that one
        //    definition at load time even if they re-instantiate the
        //    `emplace<Args...>` member template locally. The previous
        //    function-local-static design (PR #1368) produced a hidden
        //    duplicate per DSO via the inlined accessor body inside `emplace`
        //    and caused emplace/get to disagree across module boundaries.
        IMPACTX_EXPORT static int s_next_id;
        IMPACTX_EXPORT static std::map<int, std::shared_ptr<T>> s_registry;

        static int& next_id () { return s_next_id; }
        static std::map<int, std::shared_ptr<T>>& registry () { return s_registry; }

        /** Allocate unique ID for a new element instance. */
        static int allocate_id () { return next_id()++; }

        /** Retrieve data by element ID.
         *
         * @param[in] id Element ID assigned at construction
         * @return Shared pointer to coefficient data
         * @throws std::runtime_error if ID not found
         */
        static std::shared_ptr<T> const &
        get (int id)
        {
            auto& reg = registry();
            auto it = reg.find(id);
            if (it == reg.end() || !it->second)
                throw std::runtime_error(
                    "GPUDataRegistry::get failed for id=" + std::to_string(id));
            return it->second;
        }

        /** Store data for a new element.
         *
         * @param[in] id Element ID from allocate_id()
         * @param[in] data Shared pointer to coefficient data
         */
        static void store (int id, std::shared_ptr<T> data)
        {
            registry()[id] = std::move(data);
        }

        /** Construct and store data for a new element.
         *
         * @param[in] id Element ID from allocate_id()
         * @param[in] args Constructor arguments forwarded to T
         * @return Reference to the stored data
         */
        template<typename... Args>
        static T& emplace (int id, Args&&... args)
        {
            auto ptr = std::shared_ptr<T>(new T{std::forward<Args>(args)...});  // NOLINT(modernize-make-shared)
            registry()[id] = ptr;
            return *ptr;
        }

        /** Full teardown: clear all data including host. */
        static void clear ()
        {
            registry().clear();
            next_id() = 0;
        }
    };

    // Out-of-class template definitions for the static data members declared
    // inside `GPUDataRegistry`. These are deliberately not `inline static
    // ... = init;` inside the class: see the comment on the declarations for
    // the NVCC 12.2 + C++20 reason. Each element type gets exactly one strong
    // definition per program via `IMPACTX_GPUDATA_INSTANTIATE` in its .cpp.
    template<typename T>
    int GPUDataRegistry<T>::s_next_id = 0;

    template<typename T>
    std::map<int, std::shared_ptr<T>> GPUDataRegistry<T>::s_registry;

} // namespace impactx::elements::mixin


/** Suppress implicit instantiation of GPUDataRegistry for an element type.
 *
 * GPUDataRegistry holds its registry map and ID counter in exported static
 * data members (`s_next_id`, `s_registry`) declared inside the class and
 * defined out-of-line in this header. Without explicit instantiation control,
 * different translation units, especially downstream shared libraries built
 * with `-fvisibility=hidden` and/or LTO (e.g. the Python bindings module),
 * could materialize a separate copy of these statics, causing elements to be
 * registered in one registry instance but looked up in another (ODR violation).
 *
 * Use IMPACTX_GPUDATA_EXTERN in headers and IMPACTX_GPUDATA_INSTANTIATE in the
 * corresponding .cpp file to ensure a single shared registry across the program.
 *
 * @param ElementType The element class (must define `using DynamicData = ...`)
 */
#define IMPACTX_GPUDATA_EXTERN(ElementType) \
    extern template struct impactx::elements::mixin::GPUDataRegistry< \
        ElementType::DynamicData::DataType>;

/** Explicit instantiation of GPUDataRegistry for an element type.
 *
 * Provides the single definition of the registry's `s_next_id` / `s_registry`
 * static data members. Pair with IMPACTX_GPUDATA_EXTERN in the header.
 *
 * @param ElementType The element class (must define `using DynamicData = ...`)
 */
#define IMPACTX_GPUDATA_INSTANTIATE(ElementType) \
    template struct impactx::elements::mixin::GPUDataRegistry< \
        ElementType::DynamicData::DataType>;

#endif // IMPACTX_ELEMENTS_MIXIN_DYNAMICDATA_H
