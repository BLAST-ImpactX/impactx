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
     * @warning This template uses function-local static variables. To avoid ODR
     *          violations (separate registry instances per translation unit), you
     *          MUST add explicit instantiation control:
     *          - In the element header: `IMPACTX_GPUDATA_EXTERN(impactx::elements::MyElement)`
     *          - In the element .cpp:   `IMPACTX_GPUDATA_INSTANTIATE(impactx::elements::MyElement)`
     *
     * @tparam T element data struct with to_device() method
     */
    template<typename T>
    struct GPUDataRegistry
    {
        using DataType = T;  //!< The coefficient/data type stored in this registry

        IMPACTX_EXPORT static int& next_id () {
            static int id = 0;
            return id;
        }

        IMPACTX_EXPORT static std::map<int, std::shared_ptr<T>>& registry () {
            static std::map<int, std::shared_ptr<T>> reg;
            return reg;
        }

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

} // namespace impactx::elements::mixin


/** Suppress implicit instantiation of GPUDataRegistry for an element type.
 *
 * GPUDataRegistry uses function-local static variables for its registry map and
 * ID counter. Without explicit instantiation control, different translation units
 * may each get their own copy of these statics, causing elements to be registered
 * in one registry instance but looked up in another (ODR violation).
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
 * Provides the single definition of the registry's static variables.
 * Pair with IMPACTX_GPUDATA_EXTERN in the header.
 *
 * @param ElementType The element class (must define `using DynamicData = ...`)
 */
#define IMPACTX_GPUDATA_INSTANTIATE(ElementType) \
    template struct impactx::elements::mixin::GPUDataRegistry< \
        ElementType::DynamicData::DataType>;

#endif // IMPACTX_ELEMENTS_MIXIN_DYNAMICDATA_H
