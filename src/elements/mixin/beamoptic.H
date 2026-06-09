/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_BEAMOPTIC_H
#define IMPACTX_ELEMENTS_MIXIN_BEAMOPTIC_H

#include "spintransport.H"

#include "particles/ImpactXParticleContainer.H"
#include "particles/PushAll.H"

#include <AMReX_Extension.H> // for AMREX_RESTRICT
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <type_traits>


namespace impactx
{
    /** A generalized ParallelFor that dispatches to amrex::ParallelFor or amrex::ParallelForSIMD
     *  depending on whether the element type T_Element is vectorized.
     *
     * @tparam T_Element the element type
     * @tparam F the functor type
     * @param n the number of items to iterate over
     * @param f the functor to execute
     */
    template <typename T_Element, typename F>
    void ParallelFor (int n, F&& f) {
#ifdef AMREX_USE_SIMD
        if constexpr (amrex::simd::is_vectorized<T_Element>) {
            amrex::ParallelForSIMD<T_Element::simd_width>(n, std::forward<F>(f));
        } else
#endif
        {
            amrex::ParallelFor(n, std::forward<F>(f));
        }
    }
} // namespace impactx

namespace impactx::elements::mixin
{
namespace detail
{
    /** Load particle data from array pointers
     *
     * On GPU and CPU w/o SIMD, this dereferences a particle property at the
     * index position i.
     * On CPU with SIMD, this loads a SIMD register at the IndexType::width
     * SIMD-wide index position i.
     *
     * @tparam T data type (amrex::ParticleReal or uint64_t)
     * @tparam IndexType int or amrex::SIMDindex<SIMD_WIDTH, int>
     * @param ptr pointer to the array data
     * @param i index or SIMD index
     * @return a reference to the data (scalar) or a SIMD register (SIMD)
     */
    template <typename T, typename IndexType>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    decltype(auto) load_pdata (T* ptr, IndexType const i)
    {
        if constexpr (std::is_integral_v<IndexType>) {
            return ptr[i];
        } else {
#ifdef AMREX_USE_SIMD
            using namespace amrex::simd;
            using RealType = SIMDParticleReal<IndexType::width>;
            using IdCpuType = SIMDIdCpu<RealType>;
            using DataType = std::conditional_t<
                std::is_same_v<T, amrex::ParticleReal>,
                RealType,
                IdCpuType
            >;

            // initialize vector register
            // TODO stdx::vector_aligned needs alignment guarantees
            //      https://github.com/AMReX-Codes/amrex/issues/4592
            //      https://en.cppreference.com/w/cpp/experimental/simd/simd/copy_from
            DataType val;
            val.copy_from(&ptr[i.index], stdx::element_aligned);
            return val;

#else
            // error handling: we should never get here
            amrex::ignore_unused(ptr, i);
            amrex::Abort("SIMD index used but SIMD is not enabled");
            return ptr[0];
#endif
        }
    }

    /** Store particle data back to array pointers
     *
     * On GPU and CPU without SIMD, this does nothing because we already
     * modified the (global) RAM directly via pointer.
     *
     * On CPU with SIMD, this performs a conditional writeback of a SIMD register
     * to RAM (index in pointer array), but only if the argument was not passed
     * as const and thus was likely changed.
     *
     * Good optimizing compilers can eliminate writebacks of unchanged values
     * themselves, but we better help a little for robustness. Background:
     * https://github.com/AMReX-Codes/amrex/pull/4520#issuecomment-3064064215
     *
     * @tparam P_Method pointer to the push method (for is_nth_arg_non_const)
     * @tparam N the argument index (for is_nth_arg_non_const)
     * @tparam T data type
     * @tparam IndexType int or SIMD index
     * @tparam ValType the type of the value to store
     * @param val the value to store
     * @param ptr pointer to the SoA data
     * @param i index or SIMD index
     */
    template <auto P_Method, int N, typename T, typename IndexType, typename ValType>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void store_pdata (
        ValType const & AMREX_RESTRICT val,
        T * const AMREX_RESTRICT ptr,
        IndexType const i
    )
    {
#ifdef AMREX_USE_SIMD
        if constexpr (!std::is_integral_v<IndexType>) {
            if constexpr (amrex::simd::is_nth_arg_non_const(P_Method, N)) {
                // write back to memory
                // TODO stdx::vector_aligned needs alignment guarantees
                //      https://github.com/AMReX-Codes/amrex/issues/4592
                //      https://en.cppreference.com/w/cpp/experimental/simd/simd/copy_from
                val.copy_to(&ptr[i.index], amrex::simd::stdx::element_aligned);
            }
        }
#endif
        amrex::ignore_unused(val, ptr, i);
    }

    /** Push a single particle through an element
     *
     * Note: we usually would just write a C++ lambda below in ParallelFor. But, due to restrictions
     * in NVCC as of 11.5, we cannot write a lambda in a lambda as we also std::visit the element
     * types of our lattice elements list.
     *    error #3206-D: An extended __device__ lambda cannot be defined inside a generic lambda expression("operator()").
     * Thus, we fall back to writing a C++ functor here, instead of nesting two lambdas.
     *
     * Nvidia bug report: 3458976
     * Minimal demonstrator: https://cuda.godbolt.org/z/39e4q53Ye
     *
     * @tparam T_Element This can be a \see Drift, \see Quad, \see Sbend, etc.
     */
    template <typename T_Element, typename T_PT = amrex::ParticleReal>
    struct PushSingleParticleSpin
    {
        using PType = amrex::SoAParticle<RealSoA::nattribs, IntSoA::nattribs, T_PT>;
        using ParticleTileType = amrex::ParticleTile<amrex::SoAParticle<RealSoA::nattribs, IntSoA::nattribs, T_PT>,RealSoA::nattribs, IntSoA::nattribs>;

        /** Constructor taking in pointers to particle data
         *
         * @param element the beamline element to push through
         * @param part_x the array to the particle position (x)
         * @param part_y the array to the particle position (y)
         * @param part_t the array to the particle position (t)
         * @param part_px the array to the particle momentum (x)
         * @param part_py the array to the particle momentum (y)
         * @param part_pt the array to the particle momentum (t)
         * @param part_sx the array to the particle spin (x)
         * @param part_sy the array to the particle spin (y)
         * @param part_sz the array to the particle spin (z)
         * @param part_idcpu the array to the particle global index
         * @param ref_part the struct containing the reference particle
         */
        PushSingleParticleSpin (T_Element element,
                            T_PT* AMREX_RESTRICT part_x,
                            T_PT* AMREX_RESTRICT part_y,
                            T_PT* AMREX_RESTRICT part_t,
                            T_PT* AMREX_RESTRICT part_px,
                            T_PT* AMREX_RESTRICT part_py,
                            T_PT* AMREX_RESTRICT part_pt,
                            T_PT* AMREX_RESTRICT part_sx,
                            T_PT* AMREX_RESTRICT part_sy,
                            T_PT* AMREX_RESTRICT part_sz,
                            uint64_t* AMREX_RESTRICT part_idcpu,
                            RefPart ref_part)
            : m_element(std::move(element)),
              m_part_x(part_x), m_part_y(part_y), m_part_t(part_t),
              m_part_px(part_px), m_part_py(part_py), m_part_pt(part_pt),
              m_part_sx(part_sx), m_part_sy(part_sy), m_part_sz(part_sz),
              m_part_idcpu(part_idcpu),
              m_ref_part(std::move(ref_part))
        {
            // pre-compute and cache constants that are independent of the
            // individually tracked particle
            m_element.compute_constants(m_ref_part);
        }

        PushSingleParticleSpin () = delete;
        PushSingleParticleSpin (PushSingleParticleSpin const &) = default;
        PushSingleParticleSpin (PushSingleParticleSpin &&) = default;
        ~PushSingleParticleSpin () = default;

        /** Push a single particle through an element
         *
         * @param i particle index in the current box
         */
        template <typename IndexType>
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        void
        operator() (IndexType i) const
        {
            using namespace amrex::simd;

            // access SoA data
            // note: an optimizing compiler will eliminate loads of unused parameters
            decltype(auto) x = load_pdata(m_part_x, i);
            decltype(auto) y = load_pdata(m_part_y, i);
            decltype(auto) t = load_pdata(m_part_t, i);
            decltype(auto) px = load_pdata(m_part_px, i);
            decltype(auto) py = load_pdata(m_part_py, i);
            decltype(auto) pt = load_pdata(m_part_pt, i);
            decltype(auto) sx = load_pdata(m_part_sx, i);
            decltype(auto) sy = load_pdata(m_part_sy, i);
            decltype(auto) sz = load_pdata(m_part_sz, i);
            decltype(auto) idcpu = load_pdata(m_part_idcpu, i);

            // push spin & phase space through element
            m_element.spin_and_phasespace_push(x, y, t, px, py, pt, sx, sy, sz, idcpu, m_ref_part);

            // SIMD: write back to memory
#ifdef AMREX_USE_SIMD
            if constexpr (amrex::simd::is_vectorized<T_Element>)
            {
                using RealType = std::decay_t<decltype(x)>;
                using IdCpuType = std::decay_t<decltype(idcpu)>;
                constexpr auto P_Method = &T_Element::template spin_and_phasespace_push<RealType, IdCpuType>;

                store_pdata<P_Method, 0>(x, m_part_x, i);
                store_pdata<P_Method, 1>(y, m_part_y, i);
                store_pdata<P_Method, 2>(t, m_part_t, i);
                store_pdata<P_Method, 3>(px, m_part_px, i);
                store_pdata<P_Method, 4>(py, m_part_py, i);
                store_pdata<P_Method, 5>(pt, m_part_pt, i);
                store_pdata<P_Method, 6>(sx, m_part_sx, i);
                store_pdata<P_Method, 7>(sy, m_part_sy, i);
                store_pdata<P_Method, 8>(sz, m_part_sz, i);
                store_pdata<P_Method, 9>(idcpu, m_part_idcpu, i);
            }
#endif
        }

    private:
        T_Element m_element;
        T_PT* const AMREX_RESTRICT m_part_x;
        T_PT* const AMREX_RESTRICT m_part_y;
        T_PT* const AMREX_RESTRICT m_part_t;
        T_PT* const AMREX_RESTRICT m_part_px;
        T_PT* const AMREX_RESTRICT m_part_py;
        T_PT* const AMREX_RESTRICT m_part_pt;
        T_PT* const AMREX_RESTRICT m_part_sx;
        T_PT* const AMREX_RESTRICT m_part_sy;
        T_PT* const AMREX_RESTRICT m_part_sz;
        uint64_t* const AMREX_RESTRICT m_part_idcpu;
        RefPart m_ref_part;
    };

    /** Push a single particle through an element
     *
     * Note: we usually would just write a C++ lambda below in ParallelFor. But, due to restrictions
     * in NVCC as of 11.5, we cannot write a lambda in a lambda as we also std::visit the element
     * types of our lattice elements list.
     *    error #3206-D: An extended __device__ lambda cannot be defined inside a generic lambda expression("operator()").
     * Thus, we fall back to writing a C++ functor here, instead of nesting two lambdas.
     *
     * Nvidia bug report: 3458976
     * Minimal demonstrator: https://cuda.godbolt.org/z/39e4q53Ye
     *
     * @tparam T_Element This can be a \see Drift, \see Quad, \see Sbend, etc.
     */
    template <typename T_Element, typename T_PT = amrex::ParticleReal>
    struct PushSingleParticle
    {
        using PType = amrex::SoAParticle<RealSoA::nattribs, IntSoA::nattribs, T_PT>;
        using ParticleTileType = amrex::ParticleTile<amrex::SoAParticle<RealSoA::nattribs, IntSoA::nattribs, T_PT>,RealSoA::nattribs, IntSoA::nattribs>;

        /** Constructor taking in pointers to particle data
         *
         * @param element the beamline element to push through
         * @param part_x the array to the particle position (x)
         * @param part_y the array to the particle position (y)
         * @param part_t the array to the particle position (t)
         * @param part_px the array to the particle momentum (x)
         * @param part_py the array to the particle momentum (y)
         * @param part_pt the array to the particle momentum (t)
         * @param part_idcpu the array to the particle global index
         * @param ref_part the struct containing the reference particle
         */
        PushSingleParticle (T_Element element,
                            T_PT* AMREX_RESTRICT part_x,
                            T_PT* AMREX_RESTRICT part_y,
                            T_PT* AMREX_RESTRICT part_t,
                            T_PT* AMREX_RESTRICT part_px,
                            T_PT* AMREX_RESTRICT part_py,
                            T_PT* AMREX_RESTRICT part_pt,
                            uint64_t* AMREX_RESTRICT part_idcpu,
                            RefPart ref_part)
            : m_element(std::move(element)),
              m_part_x(part_x), m_part_y(part_y), m_part_t(part_t),
              m_part_px(part_px), m_part_py(part_py), m_part_pt(part_pt),
              m_part_idcpu(part_idcpu),
              m_ref_part(std::move(ref_part))
        {
            // pre-compute and cache constants that are independent of the
            // individually tracked particle
            m_element.compute_constants(m_ref_part);
        }

        PushSingleParticle () = delete;
        PushSingleParticle (PushSingleParticle const &) = default;
        PushSingleParticle (PushSingleParticle &&) = default;
        ~PushSingleParticle () = default;

        /** Push a single particle through an element
         *
         * @param i particle index in the current box
         */
        template <typename IndexType>
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        void
        operator() (IndexType i) const
        {
            using namespace amrex::simd;

            // access SoA data
            // note: an optimizing compiler will eliminate loads of unused parameters
            decltype(auto) x = load_pdata(m_part_x, i);
            decltype(auto) y = load_pdata(m_part_y, i);
            decltype(auto) t = load_pdata(m_part_t, i);
            decltype(auto) px = load_pdata(m_part_px, i);
            decltype(auto) py = load_pdata(m_part_py, i);
            decltype(auto) pt = load_pdata(m_part_pt, i);
            decltype(auto) idcpu = load_pdata(m_part_idcpu, i);

            // push through element
            m_element(x, y, t, px, py, pt, idcpu, m_ref_part);

            // write back to memory
#ifdef AMREX_USE_SIMD
            if constexpr (amrex::simd::is_vectorized<T_Element>)
            {
                using RealType = std::decay_t<decltype(x)>;
                using IdCpuType = std::decay_t<decltype(idcpu)>;
                constexpr auto P_Method = &T_Element::template operator()<RealType, IdCpuType>;

                store_pdata<P_Method, 0>(x, m_part_x, i);
                store_pdata<P_Method, 1>(y, m_part_y, i);
                store_pdata<P_Method, 2>(t, m_part_t, i);
                store_pdata<P_Method, 3>(px, m_part_px, i);
                store_pdata<P_Method, 4>(py, m_part_py, i);
                store_pdata<P_Method, 5>(pt, m_part_pt, i);
                store_pdata<P_Method, 6>(idcpu, m_part_idcpu, i);
            }
#endif
        }

    private:
        T_Element m_element;
        T_PT* const AMREX_RESTRICT m_part_x;
        T_PT* const AMREX_RESTRICT m_part_y;
        T_PT* const AMREX_RESTRICT m_part_t;
        T_PT* const AMREX_RESTRICT m_part_px;
        T_PT* const AMREX_RESTRICT m_part_py;
        T_PT* const AMREX_RESTRICT m_part_pt;
        uint64_t* const AMREX_RESTRICT m_part_idcpu;
        RefPart m_ref_part;
    };

    /** This pushes all particles on a particle iterator tile/box
     */
    template< typename T_Element, typename T_Iter >
    void push_all_particles (
        T_Iter & pti,
        RefPart & AMREX_RESTRICT ref_part,
        T_Element & element,
        bool spin /* = false */
    ) {
        const int np = pti.numParticles();

        // preparing access to particle data: SoA of Reals
        auto& soa = pti.GetStructOfArrays();
        auto& soa_real = soa.GetRealData();
        // particle real (SoA storage) type, e.g. float or double
        using T_PT = std::remove_pointer_t<std::decay_t<decltype(soa_real[RealSoA::x].dataPtr())>>;
        auto* const AMREX_RESTRICT part_x = soa_real[RealSoA::x].dataPtr();
        auto* const AMREX_RESTRICT part_y = soa_real[RealSoA::y].dataPtr();
        auto* const AMREX_RESTRICT part_t = soa_real[RealSoA::t].dataPtr();
        auto* const AMREX_RESTRICT part_px = soa_real[RealSoA::px].dataPtr();
        auto* const AMREX_RESTRICT part_py = soa_real[RealSoA::py].dataPtr();
        auto* const AMREX_RESTRICT part_pt = soa_real[RealSoA::pt].dataPtr();

        uint64_t* const AMREX_RESTRICT part_idcpu = soa.GetIdCPUData().dataPtr();

        if (spin) {
            if constexpr (std::is_base_of_v<mixin::SpinTransport, std::decay_t<T_Element>>) {
                auto* const AMREX_RESTRICT part_sx = soa_real[RealSoA::sx].dataPtr();
                auto* const AMREX_RESTRICT part_sy = soa_real[RealSoA::sy].dataPtr();
                auto* const AMREX_RESTRICT part_sz = soa_real[RealSoA::sz].dataPtr();

                detail::PushSingleParticleSpin<T_Element, T_PT> const pushSingleParticle(
                    element, part_x, part_y, part_t, part_px, part_py, part_pt, part_sx, part_sy, part_sz, part_idcpu, ref_part);

                // loop over beam particles in the box
                impactx::ParallelFor<T_Element>(np, pushSingleParticle);
            } else {
                throw std::runtime_error("Spin transport requested but element does not implement the `SpinTransport` interface class!");
            }
        }
        else {
            detail::PushSingleParticle<T_Element, T_PT> const pushSingleParticle(
                element, part_x, part_y, part_t, part_px, part_py, part_pt, part_idcpu, ref_part);

            // loop over beam particles in the box
            impactx::ParallelFor<T_Element>(np, pushSingleParticle);
        }
    }
} // namespace detail

    /** Mixin class for a regular beam optics lattice element.
     *
     * This element pushes particles relative to a reference particle.
     * All particles are pushed independently with the same logic.
     */
    template<typename T_Element>
    struct BeamOptic
    {
        /** Push first the reference particle, then all other particles
         *
         * @param[inout] pc container of the particles to push
         * @param[in] step global step for diagnostics
         * @param[in] period for periodic lattices, this is the current period (turn or cycle)
         */
        template<class T_PC>
        void operator() (
            T_PC & pc,
            int step,
            int period
        )
        {
            static_assert(
                std::is_base_of_v<BeamOptic, T_Element>,
                "BeamOptic can only be used as a mixin class!"
            );

            T_Element& element = *static_cast<T_Element*>(this);
            push_all(pc, element, step, period);
        }

        /** This pushes the particles on a particle iterator tile or box.
         *
         * Particles are relative to the reference particle.
         *
         * @param[in] pti particle iterator for a current tile or box.
         * @param[in] ref_part reference particle
         * @param[in] spin whether to perform spin tracking
         */
         template<class T_Iter>
         void operator() (
            T_Iter & pti,
            RefPart & AMREX_RESTRICT ref_part,
            bool spin
         )
         {
            static_assert(
                std::is_base_of_v<BeamOptic, T_Element>,
                "BeamOptic can only be used as a mixin class!"
            );

            T_Element& element = *static_cast<T_Element*>(this);
            detail::push_all_particles<T_Element>(pti, ref_part, element, spin);
         }
    };

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_BEAMOPTIC_H
