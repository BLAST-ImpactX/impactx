/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
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
#include "traits.H"

#include "particles/ImpactXParticleContainer.H"
#include "particles/PushAll.H"

#include <AMReX_Extension.H> // for AMREX_RESTRICT
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <type_traits>
#include <utility>


namespace impactx
{
    /** A generalized ParallelFor that dispatches to amrex::ParallelFor or amrex::ParallelForSIMD
     *  depending on whether the element type T_Element is vectorized.
     *
     * @tparam T_Element the element type
     * @tparam UseSimd default-on if AMREX_USE_SIMD is set, allows disabling SIMD for performance corner cases
     * @tparam F the functor type
     * @param n the number of items to iterate over
     * @param f the functor to execute
     */
    template <typename T_Element, bool UseSimd = true, typename F>
    void ParallelFor (int n, F&& f) {
#ifdef AMREX_USE_SIMD
        if constexpr (UseSimd && amrex::simd::is_vectorized<T_Element>) {
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
    /** Per-push SIMD policy.
     *
     * Whether an element's phase space push (`operator()`) or spin
     * (`spin_and_phasespace_push`) push should be vectorized on CPU.
     *
     * By default a push follows the element's overall vectorization
     * (@see amrex::simd::is_vectorized). An element may override the choice
     * independently per push kind by declaring a `static constexpr bool`
     * member `simd_ps` and/or `simd_spin`. This allows, e.g., vectorizing
     * the heavy spin push while keeping a light no-spin push
     * scalar in double precision (where SIMD has a small overhead).
     */
    template <typename T, typename = void>
    struct push_simd_ps : std::bool_constant<amrex::simd::is_vectorized<T>> {};
    template <typename T>
    struct push_simd_ps<T, std::void_t<decltype(T::simd_ps)>>
        : std::bool_constant<T::simd_ps> {};

    template <typename T, typename = void>
    struct push_simd_spin : std::bool_constant<amrex::simd::is_vectorized<T>> {};
    template <typename T>
    struct push_simd_spin<T, std::void_t<decltype(T::simd_spin)>>
        : std::bool_constant<T::simd_spin> {};

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
        // a non-integral index means a SIMD lane-batch was loaded into registers
        constexpr bool is_simd = !std::is_integral_v<IndexType>;

        if constexpr (is_simd) {
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

    /** Store the particle index back to memory, if anything can have changed it
     *
     * Unlike the phase space slots, idcpu is only ever modified when a particle
     * is lost. If the push method itself does not write it, the only writer left
     * is \see PipeAperture::apply_aperture - which is a no-op unless an aperture
     * is actually set on this element. That is a runtime property, so guard the
     * write-back with it: an unset aperture then costs neither the store nor,
     * once the compiler sinks it into the branch, the load.
     *
     * @tparam P_Method the push method of the element, to inspect argument constness
     * @tparam N the argument position of idcpu in the push method
     * @tparam T_HasAperture whether the element carries the \see PipeAperture mixin
     * @tparam T_Element This can be a \see Drift, \see Quad, \see Sbend, etc.
     * @tparam IndexType int or SIMD index
     * @tparam ValType the type of the value to store
     * @param val the value to store
     * @param ptr pointer to the SoA data
     * @param i index or SIMD index
     * @param element the beamline element that was pushed through
     */
    template <auto P_Method, int N, bool T_HasAperture,
              typename T, typename IndexType, typename ValType, typename T_Element>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void store_idcpu (
        ValType const & AMREX_RESTRICT val,
        T * const AMREX_RESTRICT ptr,
        IndexType const i,
        [[maybe_unused]] T_Element const & AMREX_RESTRICT element
    )
    {
#ifdef AMREX_USE_SIMD
        // the push method writes idcpu itself: always store, as before
        if constexpr (amrex::simd::is_nth_arg_non_const(P_Method, N)) {
            store_pdata<P_Method, N, true>(val, ptr, i);
        } else if constexpr (T_HasAperture) {
            // only apply_aperture can have written idcpu, and only if it is set
            if (element.aperture_active()) {
                store_pdata<P_Method, N, true>(val, ptr, i);
            }
        }
#endif
        amrex::ignore_unused(val, ptr, i);
    }

    /** Positions of the particle slots in the arguments of a push
     *
     * The phase space push (`operator()`) and the spin push (`spin_and_phasespace_push`) of an
     * element take the real-valued slots in this order, followed by idcpu and the reference
     * particle.
     */
    struct Slot
    {
        static constexpr int x = 0;
        static constexpr int y = 1;
        static constexpr int t = 2;
        static constexpr int px = 3;
        static constexpr int py = 4;
        static constexpr int pt = 5;
        static constexpr int sx = 6;  //!< spin push only
        static constexpr int sy = 7;  //!< spin push only
        static constexpr int sz = 8;  //!< spin push only
    };

    /** Calls the phase space push of an element, @see push_only_reads */
    struct PhaseSpacePushCall
    {
        /** number of real-valued slots, which precede idcpu */
        static constexpr int n_real = Slot::pt + 1;

        template <typename T_Element, typename... Args>
        auto operator() (T_Element const & element, Args &&... args) const
            -> decltype(element(std::forward<Args>(args)...));
    };

    /** Calls the spin and phase space push of an element, @see push_only_reads */
    struct SpinPushCall
    {
        /** number of real-valued slots, which precede idcpu */
        static constexpr int n_real = Slot::sz + 1;

        template <typename T_Element, typename... Args>
        auto operator() (T_Element const & element, Args &&... args) const
            -> decltype(element.spin_and_phasespace_push(std::forward<Args>(args)...));
    };

    /** Implementation of @see push_only_reads */
    template <typename T_PushCall, typename T_Element, typename T_Real, typename T_IdCpu,
              int N0, int N1, int... I>
    AMREX_GPU_HOST_DEVICE
    constexpr bool push_only_reads_impl (std::integer_sequence<int, I...>)
    {
        return std::is_invocable_v<
            T_PushCall,
            T_Element const &,
            std::conditional_t<I == N0 || I == N1, T_Real, T_Real &>...,
            T_IdCpu &,
            RefPart const &
        >;
    }

    /** Whether the push of an element only reads a pair of particle slots
     *
     * A push declares the slots it changes as non-const references. This tests whether
     * the push accepts the slots at positions N0 and N1 as rvalues: a parameter that is a
     * non-const lvalue reference, also one deduced as `auto &`, does not bind to an rvalue,
     * while `T_Real const &` and a by-value parameter do. Unlike the member pointer that
     * @see store_pdata inspects, this also works for a push that is not a template.
     *
     * A push that does not accept the slots as rvalues counts as changing them. A push must
     * not take a slot as a forwarding reference (`auto &&`): it would accept an rvalue,
     * although it can change the slot.
     *
     * @tparam T_PushCall @see PhaseSpacePushCall or @see SpinPushCall
     * @tparam T_Element This can be a \see Drift, \see Quad, \see Sbend, etc.
     * @tparam T_Real amrex::ParticleReal or a SIMD register of it
     * @tparam T_IdCpu uint64_t or a SIMD register of it
     * @tparam N0 the position of the first slot of the pair, @see Slot
     * @tparam N1 the position of the second slot of the pair
     * @return true if the push only reads both slots
     */
    template <typename T_PushCall, typename T_Element, typename T_Real, typename T_IdCpu,
              int N0, int N1>
    AMREX_GPU_HOST_DEVICE
    constexpr bool push_only_reads ()
    {
        static_assert(N0 < T_PushCall::n_real && N1 < T_PushCall::n_real,
                      "push_only_reads: the slot pair is not an argument of this push");
        return push_only_reads_impl<T_PushCall, T_Element, T_Real, T_IdCpu, N0, N1>(
            std::make_integer_sequence<int, T_PushCall::n_real>{});
    }

    /** The alignment-error frame of an element, for one push of one particle
     *
     * The alignment transforms move pairs of slots into the element frame before the push
     * and back out after it: (x, y) is shifted and rotated, (px, py) and (sx, sy) are
     * rotated, @see Alignment.
     *
     * @tparam T_PushCall @see PhaseSpacePushCall or @see SpinPushCall
     * @tparam T_Element This can be a \see Drift, \see Quad, \see Sbend, etc.
     * @tparam T_Real amrex::ParticleReal or a SIMD register of it
     * @tparam T_IdCpu uint64_t or a SIMD register of it
     * @tparam T_Align whether the alignment transforms run
     */
    template <typename T_PushCall, typename T_Element, typename T_Real, typename T_IdCpu,
              bool T_Align>
    struct AlignmentFrame
    {
        /** A pair of slots in the alignment-error frame
         *
         * For a pair that the push only reads, moving it back out would reproduce the value
         * in memory, up to rounding. Such a pair is moved into local copies instead, which
         * are never written back: this skips the transform back out and the stores, and on
         * GPU the global memory accesses in between, @see push_only_reads.
         *
         * @tparam N0 the position of the first slot, @see Slot
         * @tparam N1 the position of the second slot
         */
        template <int N0, int N1>
        struct Pair
        {
            /** the push only reads the pair: it goes through the frame in local copies */
            static constexpr bool copy = T_Align &&
                push_only_reads<T_PushCall, T_Element, T_Real, T_IdCpu, N0, N1>();

            /** the transform back out changes the pair: store it, even if the push
             *  takes a slot of it as const (SIMD, @see store_pdata) */
            static constexpr bool writeback = T_Align && !copy;

            /** A slot of the pair, as the push sees it in the element frame
             *
             * @param val the slot: a reference to memory (GPU, CPU w/o SIMD) or a SIMD register
             * @return a copy of val if the pair is copied, otherwise val itself
             */
            AMREX_GPU_DEVICE AMREX_FORCE_INLINE
            static std::conditional_t<copy, T_Real, T_Real &> element_frame (T_Real & val)
            {
                return val;
            }

            /** Move the pair into the element frame
             *
             * @param element the element, which holds the alignment errors
             * @param a the first slot, as returned by element_frame
             * @param b the second slot, as returned by element_frame
             */
            AMREX_GPU_DEVICE AMREX_FORCE_INLINE
            static void move_in (T_Element const & element, T_Real & a, T_Real & b)
            {
                if constexpr (T_Align) {
                    if constexpr (N0 == Slot::x) { element.shift_in(a, b); }
                    else { element.rotate_in(a, b); }
                }
            }

            /** Move the pair back out of the element frame, if the push can change it
             *
             * @param element the element, which holds the alignment errors
             * @param a the first slot, as returned by element_frame
             * @param b the second slot, as returned by element_frame
             */
            AMREX_GPU_DEVICE AMREX_FORCE_INLINE
            static void move_out (T_Element const & element, T_Real & a, T_Real & b)
            {
                if constexpr (writeback) {
                    if constexpr (N0 == Slot::x) { element.shift_out(a, b); }
                    else { element.rotate_out(a, b); }
                }
            }
        };

        using XY = Pair<Slot::x, Slot::y>;
        using PxPy = Pair<Slot::px, Slot::py>;
        using SxSy = Pair<Slot::sx, Slot::sy>;  //!< spin push only
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
            // Note: the element arrives already prepared. Preparation happens exactly
            // once per element slice in impactx::elements::make_particle_pusher(), and
            // the result is copied into each tile's launch closure. This constructor
            // stores and executes; it must not prepare.
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

            using RealType = std::decay_t<decltype(x)>;
            using IdCpuType = std::decay_t<decltype(idcpu)>;

#ifdef AMREX_USE_SIMD
            // a non-integral index means a SIMD lane-batch was loaded into registers
            constexpr bool is_simd = !std::is_integral_v<IndexType>;
#endif

            // pre/post-push functionality declared by mixin classes
            constexpr bool has_alignment = mixin::has_alignment_v<T_Element> && T_DoAlignment;
            constexpr bool has_aperture = mixin::has_pipe_aperture_v<T_Element>;

            // slot pairs in the alignment-error frame, @see AlignmentFrame
            using PushCall = SpinPushCall;
            using Frame = AlignmentFrame<PushCall, T_Element, RealType, IdCpuType, has_alignment>;
            using XY = typename Frame::XY;
            using PxPy = typename Frame::PxPy;
            using SxSy = typename Frame::SxSy;

            decltype(auto) xe = XY::element_frame(x);
            decltype(auto) ye = XY::element_frame(y);
            decltype(auto) pxe = PxPy::element_frame(px);
            decltype(auto) pye = PxPy::element_frame(py);
            decltype(auto) sxe = SxSy::element_frame(sx);
            decltype(auto) sye = SxSy::element_frame(sy);

            // shift into the alignment-error frame (compile-time gated)
            XY::move_in(m_element, xe, ye);
            PxPy::move_in(m_element, pxe, pye);
            SxSy::move_in(m_element, sxe, sye);

            // push spin & phase space through element
            m_element.spin_and_phasespace_push(
                xe, ye, t, pxe, pye, pt, sxe, sye, sz, idcpu, m_ref_part);

            // apply transverse aperture in the element-local frame (compile-time gated)
            if constexpr (has_aperture) {
                m_element.apply_aperture(xe, ye, idcpu);
            }

            // shift back to the lab frame (compile-time gated): pairs the push can change
            XY::move_out(m_element, xe, ye);
            PxPy::move_out(m_element, pxe, pye);
            SxSy::move_out(m_element, sxe, sye);

            // SIMD: write back to memory
#ifdef AMREX_USE_SIMD
            if constexpr (is_simd)
            {
                constexpr auto P_Method = &T_Element::template spin_and_phasespace_push<RealType, IdCpuType>;

                // move_out writes the pairs that were not copied; apply_aperture writes idcpu.
                // These wrapper-modified slots must be stored regardless of the push
                // method's argument constness.
                store_pdata<P_Method, Slot::x,  XY::writeback>  (x,  m_part_x,  i);
                store_pdata<P_Method, Slot::y,  XY::writeback>  (y,  m_part_y,  i);
                store_pdata<P_Method, Slot::t>                  (t,  m_part_t,  i);
                store_pdata<P_Method, Slot::px, PxPy::writeback>(px, m_part_px, i);
                store_pdata<P_Method, Slot::py, PxPy::writeback>(py, m_part_py, i);
                store_pdata<P_Method, Slot::pt>                 (pt, m_part_pt, i);
                store_pdata<P_Method, Slot::sx, SxSy::writeback>(sx, m_part_sx, i);
                store_pdata<P_Method, Slot::sy, SxSy::writeback>(sy, m_part_sy, i);
                store_pdata<P_Method, Slot::sz>                 (sz, m_part_sz, i);
                store_idcpu<P_Method, PushCall::n_real, has_aperture>(
                    idcpu, m_part_idcpu, i, m_element);
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
            // Note: the element arrives already prepared, @see PushSingleParticleSpin.
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

            using RealType = std::decay_t<decltype(x)>;
            using IdCpuType = std::decay_t<decltype(idcpu)>;

#ifdef AMREX_USE_SIMD
            // a non-integral index means a SIMD lane-batch was loaded into registers
            constexpr bool is_simd = !std::is_integral_v<IndexType>;
#endif

            // pre/post-push functionality declared by mixin classes
            constexpr bool has_alignment = mixin::has_alignment_v<T_Element> && T_DoAlignment;
            constexpr bool has_aperture = mixin::has_pipe_aperture_v<T_Element>;

            // slot pairs in the alignment-error frame, @see AlignmentFrame
            using PushCall = PhaseSpacePushCall;
            using Frame = AlignmentFrame<PushCall, T_Element, RealType, IdCpuType, has_alignment>;
            using XY = typename Frame::XY;
            using PxPy = typename Frame::PxPy;

            decltype(auto) xe = XY::element_frame(x);
            decltype(auto) ye = XY::element_frame(y);
            decltype(auto) pxe = PxPy::element_frame(px);
            decltype(auto) pye = PxPy::element_frame(py);

            // shift into the alignment-error frame (compile-time gated)
            XY::move_in(m_element, xe, ye);
            PxPy::move_in(m_element, pxe, pye);

            // push through element
            m_element(xe, ye, t, pxe, pye, pt, idcpu, m_ref_part);

            // apply transverse aperture in the element-local frame (compile-time gated)
            if constexpr (has_aperture) {
                m_element.apply_aperture(xe, ye, idcpu);
            }

            // shift back to the lab frame (compile-time gated): pairs the push can change
            XY::move_out(m_element, xe, ye);
            PxPy::move_out(m_element, pxe, pye);

            // write back to memory
#ifdef AMREX_USE_SIMD
            if constexpr (is_simd)
            {
                constexpr auto P_Method = &T_Element::template operator()<RealType, IdCpuType>;

                // move_out writes the pairs that were not copied; apply_aperture writes idcpu.
                // These wrapper-modified slots must be stored regardless of the push
                // method's argument constness.
                store_pdata<P_Method, Slot::x,  XY::writeback>  (x,  m_part_x,  i);
                store_pdata<P_Method, Slot::y,  XY::writeback>  (y,  m_part_y,  i);
                store_pdata<P_Method, Slot::t>                  (t,  m_part_t,  i);
                store_pdata<P_Method, Slot::px, PxPy::writeback>(px, m_part_px, i);
                store_pdata<P_Method, Slot::py, PxPy::writeback>(py, m_part_py, i);
                store_pdata<P_Method, Slot::pt>                 (pt, m_part_pt, i);
                store_idcpu<P_Method, PushCall::n_real, has_aperture>(
                    idcpu, m_part_idcpu, i, m_element);
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
     *
     * With ImpactX_OPTIMIZE_ALIGNMENT (default), alignment-capable elements compile two
     * push variants (with/without the shift/rotate transforms) and pick the matching one
     * at runtime, so a perfectly aligned element skips the otherwise-identity transforms.
     * Without it, only the always-on (always shift/rotate) variant is compiled, which
     * halves the element-push compile time at the cost of running the transforms on every
     * particle even when the element is perfectly aligned.
     */
    template< typename T_Element, typename F >
    void dispatch_misalignment ([[maybe_unused]] T_Element const & element, F && f)
    {
        if constexpr (mixin::has_alignment_v<T_Element>) {
#ifdef ImpactX_OPTIMIZE_ALIGNMENT
            if (element.misaligned()) { f(std::true_type{}); }
            else                      { f(std::false_type{}); }
#else
            // always apply the alignment transforms (single instantiation)
            f(std::true_type{});
#endif
        } else {
            f(std::false_type{});
        }
    }

    /** This pushes all particles on a particle iterator tile/box
     *
     * @param[in,out] pti the particle iterator, pointing at the tile/box to push
     * @param[in] ref_part the reference particle, already advanced for this slice
     * @param[in] element an already prepared pusher,
     *                    @see impactx::elements::make_particle_pusher
     * @param[in] spin whether to perform spin tracking
     */
    template< typename T_Element >
    void push_all_particles (
        ImpactXParticleContainer::iterator & pti,
        RefPart & AMREX_RESTRICT ref_part,
        T_Element const & element,
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
            if constexpr (mixin::has_spin_transport_v<T_Element>) {
                amrex::ParticleReal* const AMREX_RESTRICT part_sx = soa_real[RealSoA::sx].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_sy = soa_real[RealSoA::sy].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_sz = soa_real[RealSoA::sz].dataPtr();

                // loop over beam particles in the box
                dispatch_misalignment(element, [&](auto do_alignment) {
                    detail::PushSingleParticleSpin<T_Element, decltype(do_alignment)::value> const pushSingleParticle(
                        element, part_x, part_y, part_t, part_px, part_py, part_pt, part_sx, part_sy, part_sz, part_idcpu, ref_part);
                    impactx::ParallelFor<T_Element, detail::push_simd_spin<T_Element>::value>(np, pushSingleParticle);
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
                impactx::ParallelFor<T_Element, detail::push_simd_ps<T_Element>::value>(np, pushSingleParticle);
            });
        }
    }
} // namespace detail

    /** Push beam particles on one tile, relative to the reference particle.
     *
     * This is deliberately separate from @see BeamOptic. A pusher -- the
     * prepared, resource-free thing that is copied into a launch closure -- needs exactly
     * this and nothing more. It must *not* be able to start a lattice push or advance the
     * reference particle: those are the canonical element's job, they run once per slice
     * on the host, and they may need resources the pusher does not own.
     *
     * Ordinary elements are their own pusher and get this through @see BeamOptic. The
     * resource-owning elements give it to their pusher alone.
     */
    template<typename T_Element>
    struct ParticlePush
    {
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
         ) const
         {
            static_assert(
                std::is_base_of_v<ParticlePush, T_Element>,
                "ParticlePush can only be used as a mixin class!"
            );

            // Note: *this must already be a prepared pusher here,
            // @see impactx::elements::make_particle_pusher and impactx::push_all.
            T_Element const & element = *static_cast<T_Element const *>(this);
            detail::push_all_particles<T_Element>(pti, ref_part, element, spin);
         }
    };

    /** Mixin class for a regular beam optics lattice element.
     *
     * This element pushes particles relative to a reference particle. All particles are
     * pushed independently with the same logic.
     *
     * On top of @see ParticlePush, which pushes the beam particles of one tile, this adds
     * the whole-container push: it advances the reference particle first and then every
     * beam particle, once per element slice. An element that owns resources a launch
     * closure must not carry gives @see ParticlePush to its pusher alone and
     * keeps this for itself.
     */
    template<typename T_Element>
    struct BeamOptic : ParticlePush<T_Element>
    {
        using ParticlePush<T_Element>::operator();

        /** Push the beam particles of one tile through the canonical element.
         *
         * @see ParticlePush does the same thing but assumes it is already prepared, which
         * is what the once-per-slice path hands it. The canonical element is not: its
         * constants belong to a reference particle and a slice, and on a fresh element
         * they have never been written at all. So prepare a copy here and push with that,
         * leaving the element itself alone.
         *
         * The tile loop in @see impactx::push_tiles reaches @see ParticlePush directly and
         * does not pay for this.
         *
         * @param[in] pti particle iterator for a current tile or box.
         * @param[in] ref_part reference particle
         * @param[in] spin whether to perform spin tracking
         */
        void operator() (
            ImpactXParticleContainer::iterator & pti,
            RefPart & AMREX_RESTRICT ref_part,
            bool spin
        ) const
        {
            // Prepare through @see elements::make_particle_pusher, as the whole-container
            // push does: what a prepared pusher *is* differs by element. Ordinary ones are
            // their own, while the elements that own array data return a separate,
            // resource-free type -- copying the element and pushing that would not even
            // compile for those, as their canonical type has no particle call operator.
            auto const pusher = make_particle_pusher(
                static_cast<T_Element const &>(*this), ref_part);
            using pusher_type = std::decay_t<decltype(pusher)>;

            // through ParticlePush, not the pusher's own call operator: for an ordinary
            // element the pusher is this very type, and calling it plainly would come
            // straight back here
            static_cast<ParticlePush<pusher_type> const &>(pusher)(pti, ref_part, spin);
        }

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
    };

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_BEAMOPTIC_H
