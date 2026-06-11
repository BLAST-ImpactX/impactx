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

#include "alignment.H"
#include "pipeaperture.H"
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
     * The `ForceWriteback` flag forces the writeback regardless of the push
     * method's argument constness. Use it for slots that may be modified by
     * something other than the named push method itself (e.g., a pre/post
     * step in the caller), where the push method's signature is not the
     * authoritative source of truth for whether the value changed.
     *
     * @tparam P_Method pointer to the push method (for is_nth_arg_non_const)
     * @tparam N the argument index (for is_nth_arg_non_const)
     * @tparam ForceWriteback if true, write back regardless of P_Method's argument constness
     * @tparam T data type
     * @tparam IndexType int or SIMD index
     * @tparam ValType the type of the value to store
     * @param val the value to store
     * @param ptr pointer to the SoA data
     * @param i index or SIMD index
     */
    template <auto P_Method, int N, bool ForceWriteback = false,
              typename T, typename IndexType, typename ValType>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void store_pdata (
        ValType const & AMREX_RESTRICT val,
        T * const AMREX_RESTRICT ptr,
        IndexType const i
    )
    {
#ifdef AMREX_USE_SIMD
        if constexpr (!std::is_integral_v<IndexType>) {
            if constexpr (ForceWriteback || amrex::simd::is_nth_arg_non_const(P_Method, N)) {
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
     * @tparam T_DoAlignment whether to apply or ignore element misalignment
     */
    template <typename T_Element, bool T_DoAlignment>
    struct PushSingleParticleSpin
    {
        using PType = ImpactXParticleContainer::ParticleType;
        using ParticleTileType = amrex::ParticleTile<amrex::SoAParticle<RealSoA::nattribs, IntSoA::nattribs>,RealSoA::nattribs, IntSoA::nattribs>;

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
                            amrex::ParticleReal* AMREX_RESTRICT part_x,
                            amrex::ParticleReal* AMREX_RESTRICT part_y,
                            amrex::ParticleReal* AMREX_RESTRICT part_t,
                            amrex::ParticleReal* AMREX_RESTRICT part_px,
                            amrex::ParticleReal* AMREX_RESTRICT part_py,
                            amrex::ParticleReal* AMREX_RESTRICT part_pt,
                            amrex::ParticleReal* AMREX_RESTRICT part_sx,
                            amrex::ParticleReal* AMREX_RESTRICT part_sy,
                            amrex::ParticleReal* AMREX_RESTRICT part_sz,
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

            // pre/post-push functionality declared by mixin classes
            constexpr bool has_alignment = std::is_base_of_v<mixin::Alignment, T_Element> && T_DoAlignment;
            constexpr bool has_aperture = std::is_base_of_v<mixin::PipeAperture, T_Element>;

            // shift into the alignment-error frame (compile-time gated)
            if constexpr (has_alignment) {
                m_element.shift_in(x, y, px, py);
                m_element.rotate_in(sx, sy);
            }

            // push spin & phase space through element
            m_element.spin_and_phasespace_push(x, y, t, px, py, pt, sx, sy, sz, idcpu, m_ref_part);

            // apply transverse aperture in the element-local frame (compile-time gated)
            if constexpr (has_aperture) {
                m_element.apply_aperture(x, y, idcpu);
            }

            // shift back to the lab frame (compile-time gated)
            if constexpr (has_alignment) {
                m_element.shift_out(x, y, px, py);
                m_element.rotate_out(sx, sy);
            }

            // SIMD: write back to memory
#ifdef AMREX_USE_SIMD
            if constexpr (amrex::simd::is_vectorized<T_Element>)
            {
                using RealType = std::decay_t<decltype(x)>;
                using IdCpuType = std::decay_t<decltype(idcpu)>;
                constexpr auto P_Method = &T_Element::template spin_and_phasespace_push<RealType, IdCpuType>;

                // shift_in/shift_out write x, y, px, py; rotate_in/rotate_out
                // write sx, sy; apply_aperture writes idcpu.
                // These wrapper-modified slots must be stored regardless of the push
                // method's argument constness.
                constexpr bool wb_align    = has_alignment;
                constexpr bool wb_aperture = has_aperture;

                store_pdata<P_Method, 0, wb_align>   (x,     m_part_x,     i);
                store_pdata<P_Method, 1, wb_align>   (y,     m_part_y,     i);
                store_pdata<P_Method, 2>             (t,     m_part_t,     i);
                store_pdata<P_Method, 3, wb_align>   (px,    m_part_px,    i);
                store_pdata<P_Method, 4, wb_align>   (py,    m_part_py,    i);
                store_pdata<P_Method, 5>             (pt,    m_part_pt,    i);
                store_pdata<P_Method, 6, wb_align>   (sx,    m_part_sx,    i);
                store_pdata<P_Method, 7, wb_align>   (sy,    m_part_sy,    i);
                store_pdata<P_Method, 8>             (sz,    m_part_sz,    i);
                store_pdata<P_Method, 9, wb_aperture>(idcpu, m_part_idcpu, i);
            }
#endif
        }

    private:
        T_Element m_element;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_x;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_y;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_t;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_px;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_py;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_pt;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_sx;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_sy;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_sz;
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
     * @tparam T_DoAlignment whether to apply or ignore element misalignment
     */
    template <typename T_Element, bool T_DoAlignment>
    struct PushSingleParticle
    {
        using PType = ImpactXParticleContainer::ParticleType;
        using ParticleTileType = amrex::ParticleTile<amrex::SoAParticle<RealSoA::nattribs, IntSoA::nattribs>,RealSoA::nattribs, IntSoA::nattribs>;

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
                            amrex::ParticleReal* AMREX_RESTRICT part_x,
                            amrex::ParticleReal* AMREX_RESTRICT part_y,
                            amrex::ParticleReal* AMREX_RESTRICT part_t,
                            amrex::ParticleReal* AMREX_RESTRICT part_px,
                            amrex::ParticleReal* AMREX_RESTRICT part_py,
                            amrex::ParticleReal* AMREX_RESTRICT part_pt,
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

            // pre/post-push functionality declared by mixin classes
            constexpr bool has_alignment = std::is_base_of_v<mixin::Alignment, T_Element> && T_DoAlignment;
            constexpr bool has_aperture = std::is_base_of_v<mixin::PipeAperture, T_Element>;

            // shift into the alignment-error frame (compile-time gated)
            if constexpr (has_alignment) {
                m_element.shift_in(x, y, px, py);
            }

            // push through element
            m_element(x, y, t, px, py, pt, idcpu, m_ref_part);

            // apply transverse aperture in the element-local frame (compile-time gated)
            if constexpr (has_aperture) {
                m_element.apply_aperture(x, y, idcpu);
            }

            // shift back to the lab frame (compile-time gated)
            if constexpr (has_alignment) {
                m_element.shift_out(x, y, px, py);
            }

            // write back to memory
#ifdef AMREX_USE_SIMD
            if constexpr (amrex::simd::is_vectorized<T_Element>)
            {
                using RealType = std::decay_t<decltype(x)>;
                using IdCpuType = std::decay_t<decltype(idcpu)>;
                constexpr auto P_Method = &T_Element::template operator()<RealType, IdCpuType>;

                // shift_in/shift_out write x, y, px, py; apply_aperture writes idcpu.
                // These wrapper-modified slots must be stored regardless of the push
                // method's argument constness.
                constexpr bool wb_align    = has_alignment;
                constexpr bool wb_aperture = has_aperture;

                store_pdata<P_Method, 0, wb_align>   (x,     m_part_x,     i);
                store_pdata<P_Method, 1, wb_align>   (y,     m_part_y,     i);
                store_pdata<P_Method, 2>             (t,     m_part_t,     i);
                store_pdata<P_Method, 3, wb_align>   (px,    m_part_px,    i);
                store_pdata<P_Method, 4, wb_align>   (py,    m_part_py,    i);
                store_pdata<P_Method, 5>             (pt,    m_part_pt,    i);
                store_pdata<P_Method, 6, wb_aperture>(idcpu, m_part_idcpu, i);
            }
#endif
        }

    private:
        T_Element m_element;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_x;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_y;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_t;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_px;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_py;
        amrex::ParticleReal* const AMREX_RESTRICT m_part_pt;
        uint64_t* const AMREX_RESTRICT m_part_idcpu;
        RefPart m_ref_part;
    };

    /** Run a callable with the element's runtime misalignment lifted to a compile-time
     *  bool, selecting specialized (branch-free) push paths for with and without alignment.
     */
    template< typename T_Element, typename F >
    void dispatch_misalignment (T_Element & element, F && f)
    {
        if constexpr (std::is_base_of_v<mixin::Alignment, std::decay_t<T_Element>>) {
            if (element.misaligned()) { f(std::true_type{}); }
            else                      { f(std::false_type{}); }
        } else {
            f(std::false_type{});
        }
    }

    /** This pushes all particles on a particle iterator tile/box
     */
    template< typename T_Element >
    void push_all_particles (
        ImpactXParticleContainer::iterator & pti,
        RefPart & AMREX_RESTRICT ref_part,
        T_Element & element,
        bool spin /* = false */
    ) {
        const int np = pti.numParticles();

        // preparing access to particle data: SoA of Reals
        auto& soa = pti.GetStructOfArrays();
        auto& soa_real = soa.GetRealData();
        amrex::ParticleReal* const AMREX_RESTRICT part_x = soa_real[RealSoA::x].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT part_y = soa_real[RealSoA::y].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT part_t = soa_real[RealSoA::t].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT part_px = soa_real[RealSoA::px].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT part_py = soa_real[RealSoA::py].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT part_pt = soa_real[RealSoA::pt].dataPtr();

        uint64_t* const AMREX_RESTRICT part_idcpu = soa.GetIdCPUData().dataPtr();

        if (spin) {
            if constexpr (std::is_base_of_v<mixin::SpinTransport, std::decay_t<T_Element>>) {
                amrex::ParticleReal* const AMREX_RESTRICT part_sx = soa_real[RealSoA::sx].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_sy = soa_real[RealSoA::sy].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_sz = soa_real[RealSoA::sz].dataPtr();

                // loop over beam particles in the box
                dispatch_misalignment(element, [&](auto do_alignment) {
                    detail::PushSingleParticleSpin<T_Element, decltype(do_alignment)::value> const pushSingleParticle(
                        element, part_x, part_y, part_t, part_px, part_py, part_pt, part_sx, part_sy, part_sz, part_idcpu, ref_part);
                    impactx::ParallelFor<T_Element>(np, pushSingleParticle);
                });
            } else {
                throw std::runtime_error("Spin transport requested but element does not implement the `SpinTransport` interface class!");
            }
        }
        else {
            // loop over beam particles in the box
            dispatch_misalignment(element, [&](auto do_alignment) {
                detail::PushSingleParticle<T_Element, decltype(do_alignment)::value> const pushSingleParticle(
                    element, part_x, part_y, part_t, part_px, part_py, part_pt, part_idcpu, ref_part);
                impactx::ParallelFor<T_Element>(np, pushSingleParticle);
            });
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
        void operator() (
            ImpactXParticleContainer & pc,
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
         void operator() (
            ImpactXParticleContainer::iterator & pti,
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
