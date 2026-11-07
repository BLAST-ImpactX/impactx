/* Copyright 2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Marco Garten, Chad Mitchell, Yinjian Zhao, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */

#include "ReducedBeamCharacteristics.H"

#include "particles/ImpactXParticleContainer.H"
#include "particles/ReferenceParticle.H"
#include "particles/CovarianceMatrix.H"
#include "EmittanceInvariants.H"

#include <AMReX_BLProfiler.H>           // for TinyProfiler
#include <AMReX_Extension.H>            // for AMREX_FORCE_INLINE
#include <AMReX_GpuDevice.H>            // for dtoh_memcpy
#include <AMReX_GpuQualifiers.H>        // for AMREX_GPU_HOST_DEVICE
#include <AMReX_ParallelDescriptor.H>   // for ParallelDescriptor
#include <AMReX_ParticleReduceSIMD.H>   // for ParticleReduceSIMD
#include <AMReX_REAL.H>                 // for ParticleReal
#include <AMReX_Reduce.H>               // for ReduceOps
#include <AMReX_SIMD.H>                 // for simd::load_1d
#include <AMReX_SmallMatrix.H>          // for SmallMatrix
#include <AMReX_Tuple.H>                // for makeTuple
#include <AMReX_TypeList.H>             // for TypeMultiplier

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>


namespace impactx::diagnostics
{
namespace
{
    /** Square root of a quantity that is mathematically non-negative
     *
     * Variances and squared emittances are non-negative by construction, but the
     * finite-precision moment sums they are recovered from can put them marginally
     * below zero. std::sqrt of such a value raises FE_INVALID, which aborts the run
     * under amrex.fpe_trap_invalid, so clamp the (rounding-level) negative case to
     * zero instead.
     *
     * A NaN input is passed through rather than clamped: it means the beam moments
     * are already poisoned, which should stay visible in the output.
     *
     * @param x a quantity that is non-negative up to rounding
     * @returns sqrt(x), or zero if x is negative
     */
    AMREX_FORCE_INLINE
    amrex::ParticleReal
    sqrt_clamped (amrex::ParticleReal x)
    {
        using namespace amrex::literals; // for _prt

        return (x < 0.0_prt) ? 0.0_prt : std::sqrt(x);
    }

    //! constant coordinate shifts subtracted before forming the beam moments
    struct Shifts
    {
        amrex::ParticleReal x, y, t, px, py, pt, sx, sy, sz, w;
        //! slopes of the shear of the momenta, see momentum_slope
        amrex::ParticleReal slope_x, slope_y, slope_t;
    };

    /** Mean of a sample of particles
     *
     * If all values are the same, then this is exactly that value: the second moments of a
     * plane without any spread are 0.
     *
     * @param q a coordinate of a sample of particles
     * @param n size of the sample
     * @return the mean of the sample
     */
    amrex::ParticleReal
    sample_mean (amrex::ParticleReal const * q, int n)
    {
        double deviation = 0.0;
        for (int i = 1; i < n; ++i) {
            deviation += static_cast<double>(q[i]) - static_cast<double>(q[0]);
        }
        return static_cast<amrex::ParticleReal>(static_cast<double>(q[0]) + deviation / n);
    }

    /** Slope to shear the momentum of a strongly correlated phase plane
     *
     * The squared emittance <q^2><p^2> - <qp>^2 of a plane with a correlation r cancels, which
     * multiplies the relative error of its moments by 1/(1 - r^2). This is the case for a chirped
     * beam, or a beam far from its waist. In single precision, the moments of a beam with
     * 1 - r^2 close to 1e-7 do not hold its emittance any more.
     *
     * The shear p -> p - a*q is symplectic, so it does not change the emittance and the
     * eigenemittances, and with a slope a = <qp>/<q^2> the moments of the plane do not cancel.
     * Any constant slope is exact. So it is enough to estimate it from a few particles, before
     * the sums over all particles.
     *
     * @param q,p position and momentum of a sample of particles
     * @param n size of the sample
     * @return the slope, or 0 for a sample without a strong correlation: its plane does not
     *         cancel, or the sample is too narrow to tell the slope, e.g., a thin slice of a beam
     */
    double
    momentum_slope (amrex::ParticleReal const * q, amrex::ParticleReal const * p, int n)
    {
        // the correlation of 2 particles is always 1
        if (n < 3) { return 0.0; }

        double mean_q = 0.0, mean_p = 0.0;
        for (int i = 0; i < n; ++i) {
            mean_q += static_cast<double>(q[i]);
            mean_p += static_cast<double>(p[i]);
        }
        mean_q /= n;
        mean_p /= n;

        double qq = 0.0, pp = 0.0, qp = 0.0;
        for (int i = 0; i < n; ++i) {
            double const dq = static_cast<double>(q[i]) - mean_q;
            double const dp = static_cast<double>(p[i]) - mean_p;
            qq += dq*dq;
            pp += dp*dp;
            qp += dq*dp;
        }

        // r^2 > 0.9: the plane loses less than a digit below this
        return (qp*qp > 0.9*qq*pp) ? qp/qq : 0.0;
    }

    /** Per-particle beam-moments reduction kernel.
     *
     * Invoked by amrex::ParticleReduceSIMD as f(ptd, si). On CPU with SIMD
     * support (ImpactX_SIMD=ON), the vectorized main loop evaluates this for a
     * SIMD register of particles at a time; the scalar remainder loop, CPUs
     * without SIMD support, and GPUs evaluate it for one particle at a time.
     * amrex::simd::load_1d makes the same source cover all three cases.
     *
     * This is a functor with a templated operator() instead of a generic
     * lambda for CUDA (extended device lambda) portability.
     */
    struct BeamMomentsKernel
    {
        Shifts m_shift;

        template <typename PTD, typename SI>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        auto operator() (PTD const& ptd, SI const si) const noexcept
        {
            using amrex::simd::load_1d;

            // access SoA particle position, momentum, spin data and weighting
            auto const p_w  = load_1d(ptd.rdata(RealSoA::w), si);
            auto const p_x  = load_1d(ptd.rdata(RealSoA::x), si);
            auto const p_y  = load_1d(ptd.rdata(RealSoA::y), si);
            auto const p_t  = load_1d(ptd.rdata(RealSoA::t), si);
            auto const p_px = load_1d(ptd.rdata(RealSoA::px), si);
            auto const p_py = load_1d(ptd.rdata(RealSoA::py), si);
            auto const p_pt = load_1d(ptd.rdata(RealSoA::pt), si);
            auto const p_sx = load_1d(ptd.rdata(RealSoA::sx), si);
            auto const p_sy = load_1d(ptd.rdata(RealSoA::sy), si);
            auto const p_sz = load_1d(ptd.rdata(RealSoA::sz), si);

            // deviations from the shift: O(rms) rather than O(coordinate), which
            // keeps the (weighted) second moments below well-conditioned
            auto const dx  = p_x  - m_shift.x;
            auto const dy  = p_y  - m_shift.y;
            auto const dt  = p_t  - m_shift.t;
            // ... and the same for the momenta, which are also sheared, see momentum_slope
            auto const dpx = p_px - m_shift.px - m_shift.slope_x * dx;
            auto const dpy = p_py - m_shift.py - m_shift.slope_y * dy;
            auto const dpt = p_pt - m_shift.pt - m_shift.slope_t * dt;
            auto const dsx = p_sx - m_shift.sx;
            auto const dsy = p_sy - m_shift.sy;
            auto const dsz = p_sz - m_shift.sz;

            return amrex::makeTuple(
                // Sum(w), shifted: adding the same weight rounds in the same direction every time
                p_w - m_shift.w,
                // weighted first moments (shifted): x, y, t, px, py, pt, sx, sy, sz
                dx*p_w, dy*p_w, dt*p_w, dpx*p_w, dpy*p_w, dpt*p_w,
                dsx*p_w, dsy*p_w, dsz*p_w,
                // weighted second moments (shifted): diagonal x, y, t, px, py, pt
                dx*dx*p_w, dy*dy*p_w, dt*dt*p_w, dpx*dpx*p_w, dpy*dpy*p_w, dpt*dpt*p_w,
                // same-plane correlations: xpx, ypy, tpt
                dx*dpx*p_w, dy*dpy*p_w, dt*dpt*p_w,
                // dispersive correlations: xpt, pxpt, ypt, pypt
                dx*dpt*p_w, dpx*dpt*p_w, dy*dpt*p_w, dpy*dpt*p_w,
                // cross-plane correlations: xy, xpy, xt, pxy, pxpy, pxt, yt, pyt
                dx*dy*p_w, dx*dpy*p_w, dx*dt*p_w, dpx*dy*p_w, dpx*dpy*p_w, dpx*dt*p_w, dy*dt*p_w, dpy*dt*p_w,
                // spin second moments (diagonal): sx, sy, sz
                dsx*dsx*p_w, dsy*dsy*p_w, dsz*dsz*p_w,
                // min of x, y, t, px, py, pt
                p_x, p_y, p_t, p_px, p_py, p_pt,
                // max of x, y, t, px, py, pt
                p_x, p_y, p_t, p_px, p_py, p_pt
            );
        }
    };
} // anonymous namespace

    std::unordered_map<std::string, amrex::ParticleReal>
    reduced_beam_characteristics (ImpactXParticleContainer const & pc)
    {
        BL_PROFILE("impactx::diagnostics::reduced_beam_characteristics(pc)");

        using namespace amrex::literals; // for _prt

        auto const nan = std::numeric_limits<amrex::ParticleReal>::quiet_NaN();

        // preparing to access reference particle data: RefPart
        RefPart const ref_part = pc.GetRefParticle();
        // reference particle charge in C
        amrex::ParticleReal const q_C = ref_part.charge;
        /* A reference particle that carries no energy yet, e.g. before a Source
         * element loads it from file, has gamma = 0 and thus no valid beta.
         */
        bool const ref_is_initialized = ref_part.gamma() >= 1.0_prt;
        // reference particle relativistic beta*gamma
        amrex::ParticleReal const bg = ref_is_initialized ? ref_part.beta_gamma() : 0.0_prt;
        amrex::ParticleReal const bg2 = bg*bg;

        /* The reduced beam characteristics are computed in a single pass over the
         * particles from raw (weighted) power sums. For beams that are off-center
         * from the reference orbit, accumulating the sums relative to a shift near
         * the beam centroid keeps them at the O(rms^2) scale instead of O(offset^2).
         * The recovered central moments then carry a relative error ~eps*(offset/rms):
         * linear, as in a mean-subtracted two-pass reduction, rather than the
         * ~eps*(offset/rms)^2 of a raw <u^2> - <u>^2, which would lose accuracy in
         * single precision. The essential point is that we square the (small)
         * deviation from the shift, not the (large) coordinate.
         *
         * Any global constant is exact under the parallel-axis theorem. The mean of a
         * sample of particles needs only offset/rms (i.e. to land well within rms of the
         * centroid), not the true mean, and makes no assumption on the beam shape. A
         * sample instead of a single particle avoids amplifying the rounding errors of
         * the sums by 1 + (offset/rms)^2, if that particle is in the tail of the beam.
         * (In the offset >> rms limit the leading digits of the subtraction cancel
         * exactly per Sterbenz's lemma, but that is an asymptotic bonus, not the
         * mechanism.)
         *
         * The same goes for the shear of the momenta, see momentum_slope.
         */
        constexpr int comps[10] = {
            RealSoA::x, RealSoA::y, RealSoA::t,
            RealSoA::px, RealSoA::py, RealSoA::pt,
            RealSoA::sx, RealSoA::sy, RealSoA::sz,
            RealSoA::w
        };
        // the 10 shifts of comps, then the slopes in x, y, t
        std::array<amrex::ParticleReal, 13> shift{};
        bool has_particles = false;
        {
            // Sample the first particles held locally, then agree on a single global
            // shift taken from the lowest rank that actually owns a particle. Any
            // global constant is exact under the parallel-axis theorem, so a rank
            // with no particles (e.g. one MPI rank of many, or a container that is
            // momentarily empty right after a restart) simply contributes nothing.
            // For a globally empty beam the zero shift is kept.
            bool found = false;
            for (int lev = 0; lev <= pc.finestLevel() && !found; ++lev) {
                for (auto const & kv : pc.GetParticles(lev)) {
                    auto const & ptile = kv.second;
                    if (ptile.numParticles() > 0) {
                        auto const & soa = ptile.GetStructOfArrays().GetRealData();

                        constexpr int max_sample = 32;
                        int const num_sample = static_cast<int>(
                            std::min<amrex::Long>(ptile.numParticles(), max_sample));
                        std::array<std::array<amrex::ParticleReal, max_sample>, 10> sample{};
                        for (int c = 0; c < 10; ++c) {
                            amrex::Gpu::dtoh_memcpy(
                                sample[c].data(), soa[comps[c]].dataPtr(),
                                num_sample * sizeof(amrex::ParticleReal));
                            shift[c] = sample_mean(sample[c].data(), num_sample);
                        }
                        for (int plane = 0; plane < 3; ++plane) {
                            shift[10 + plane] = static_cast<amrex::ParticleReal>(momentum_slope(
                                sample[plane].data(), sample[3 + plane].data(), num_sample));
                        }
                        found = true;
                        break;
                    }
                }
            }
            int src_rank = found ? amrex::ParallelDescriptor::MyProc()
                                 : amrex::ParallelDescriptor::NProcs();
            amrex::ParallelAllReduce::Min(src_rank, amrex::ParallelDescriptor::Communicator());
            // Does the beam hold any particle at all:
            has_particles = src_rank < amrex::ParallelDescriptor::NProcs();
            if (has_particles) {
                amrex::ParallelDescriptor::Bcast(shift.data(), shift.size(), src_rank);
            }
        }
        amrex::ParticleReal const shift_x  = shift[0];
        amrex::ParticleReal const shift_y  = shift[1];
        amrex::ParticleReal const shift_t  = shift[2];
        amrex::ParticleReal const shift_px = shift[3];
        amrex::ParticleReal const shift_py = shift[4];
        amrex::ParticleReal const shift_pt = shift[5];
        amrex::ParticleReal const shift_sx = shift[6];
        amrex::ParticleReal const shift_sy = shift[7];
        amrex::ParticleReal const shift_sz = shift[8];
        amrex::ParticleReal const shift_w  = shift[9];
        amrex::ParticleReal const slope_x  = shift[10];
        amrex::ParticleReal const slope_y  = shift[11];
        amrex::ParticleReal const slope_t  = shift[12];

        /* The variables below need to be static to work around an MSVC bug
         * https://stackoverflow.com/questions/55136414/constexpr-variable-captured-inside-lambda-loses-its-constexpr-ness
         */
        // numbers of same-type reduction operations, fused into a single pass:
        //   Sum(w), 9 weighted first moments, 24 weighted (shifted) second moments,
        //   and the min/max of the 6 phase-space coordinates.
        static constexpr std::size_t num_sum = 34;
        static constexpr std::size_t num_min = 6;
        static constexpr std::size_t num_max = 6;

        amrex::TypeMultiplier<amrex::ReduceOps,
            amrex::ReduceOpSum[num_sum],  // Sum(w) + first and second moments, all shifted
            amrex::ReduceOpMin[num_min],  // min of x, y, t, px, py, pt
            amrex::ReduceOpMax[num_max]   // max of x, y, t, px, py, pt
        > reduce_ops;
        using ReducedDataT = amrex::TypeMultiplier<amrex::ReduceData, amrex::ParticleReal[num_sum + num_min + num_max]>;

        /* Fused single-pass reduction over all particles. With ImpactX_SIMD=ON,
         * the sums accumulate in per-lane partial sums that are folded at the
         * end, which reassociates the floating-point additions relative to the
         * scalar evaluation order (rounding-level differences only); min/max
         * are exact either way.
         */
        BeamMomentsKernel const beam_moments{
            Shifts{shift_x, shift_y, shift_t, shift_px, shift_py, shift_pt,
                   shift_sx, shift_sy, shift_sz, shift_w, slope_x, slope_y, slope_t}};
        auto r = amrex::ParticleReduceSIMD<ReducedDataT>(pc, beam_moments, reduce_ops);

        // extract this rank's partial sums, minima and maxima
        std::vector<amrex::ParticleReal> values_sum(num_sum);
        amrex::constexpr_for<0, num_sum> ([&](auto i) {
            values_sum[i] = amrex::get<i>(r);
        });
        std::vector<amrex::ParticleReal> values_min(num_min);
        amrex::constexpr_for<0, num_min> ([&](auto i) {
            constexpr std::size_t idx = i + num_sum;
            values_min[i] = amrex::get<idx>(r);
        });
        std::vector<amrex::ParticleReal> values_max(num_max);
        amrex::constexpr_for<0, num_max> ([&](auto i) {
            constexpr std::size_t idx = i + num_sum + num_min;
            values_max[i] = amrex::get<idx>(r);
        });

        // undo the shift of Sum(w) on this rank
        amrex::Long num_particles_local = 0;
        for (int lev = 0; lev <= pc.finestLevel(); ++lev) {
            for (auto const & kv : pc.GetParticles(lev)) {
                num_particles_local += kv.second.numParticles();
            }
        }
        values_sum[0] += static_cast<amrex::ParticleReal>(
            static_cast<double>(shift_w) * static_cast<double>(num_particles_local));

        // reduce across MPI ranks (allreduce)
        amrex::ParallelAllReduce::Sum(
            values_sum.data(), values_sum.size(), amrex::ParallelDescriptor::Communicator());
        amrex::ParallelAllReduce::Min(
            values_min.data(), values_min.size(), amrex::ParallelDescriptor::Communicator());
        amrex::ParallelAllReduce::Max(
            values_max.data(), values_max.size(), amrex::ParallelDescriptor::Communicator());

        // Recover the beam moments from the raw (shifted, sheared) power sums via the
        // parallel-axis theorem. The shift keeps every sum at the O(rms^2) scale,
        // so the central moments below stay well-conditioned even in single precision.
        amrex::ParticleReal const w_sum = values_sum[0];

        /* Every moment below is normalized by the beam's total weight, which is zero
         * for a beam that holds no particles at all - e.g. the initial diagnostic of
         * a simulation whose beam a Source element loads further down the lattice -
         * and also for one of zero charge, as used for massless "test" particles.
         * Both would turn each normalization into invalid divisions by zero.
         * Normalize by one instead, which carries the
         * (all-zero) sums through the derived quantities without an exception, and
         * report the undefined results as NaN at the end of this function.
         *
         * The comparison is an equality on purpose: it is quiet for a NaN w_sum,
         * so a beam poisoned by a non-finite particle coordinate is not silently
         * absorbed here, but keeps tripping FPE diagnostics when enabled.
         */
        bool const has_weight = !(w_sum == 0.0_prt);
        amrex::ParticleReal const w_norm = has_weight ? w_sum : 1.0_prt;

        // shifted means: <u - shift_u>, and sheared for the momenta
        amrex::ParticleReal const dmean_x  = values_sum[1] / w_norm;
        amrex::ParticleReal const dmean_y  = values_sum[2] / w_norm;
        amrex::ParticleReal const dmean_t  = values_sum[3] / w_norm;
        amrex::ParticleReal const dmean_px = values_sum[4] / w_norm;
        amrex::ParticleReal const dmean_py = values_sum[5] / w_norm;
        amrex::ParticleReal const dmean_pt = values_sum[6] / w_norm;
        amrex::ParticleReal const dmean_sx = values_sum[7] / w_norm;
        amrex::ParticleReal const dmean_sy = values_sum[8] / w_norm;
        amrex::ParticleReal const dmean_sz = values_sum[9] / w_norm;
        // means
        amrex::ParticleReal const mean_x  = shift_x  + dmean_x;
        amrex::ParticleReal const mean_y  = shift_y  + dmean_y;
        amrex::ParticleReal const mean_t  = shift_t  + dmean_t;
        amrex::ParticleReal const mean_px = shift_px + dmean_px + slope_x * dmean_x;
        amrex::ParticleReal const mean_py = shift_py + dmean_py + slope_y * dmean_y;
        amrex::ParticleReal const mean_pt = shift_pt + dmean_pt + slope_t * dmean_t;
        amrex::ParticleReal const mean_sx = shift_sx + dmean_sx;
        amrex::ParticleReal const mean_sy = shift_sy + dmean_sy;
        amrex::ParticleReal const mean_sz = shift_sz + dmean_sz;
        // minimum values
        amrex::ParticleReal const min_x = values_min.at(0);
        amrex::ParticleReal const min_y = values_min.at(1);
        amrex::ParticleReal const min_t = values_min.at(2);
        amrex::ParticleReal const min_px = values_min.at(3);
        amrex::ParticleReal const min_py = values_min.at(4);
        amrex::ParticleReal const min_pt = values_min.at(5);
        // maximum values
        amrex::ParticleReal const max_x = values_max.at(0);
        amrex::ParticleReal const max_y = values_max.at(1);
        amrex::ParticleReal const max_t = values_max.at(2);
        amrex::ParticleReal const max_px = values_max.at(3);
        amrex::ParticleReal const max_py = values_max.at(4);
        amrex::ParticleReal const max_pt = values_max.at(5);
        // mean square and correlation values (central moments via parallel-axis theorem)
        // of the positions and the sheared momenta (s_)
        amrex::ParticleReal const x_ms     = values_sum[10] / w_norm - dmean_x  * dmean_x;
        amrex::ParticleReal const y_ms     = values_sum[11] / w_norm - dmean_y  * dmean_y;
        amrex::ParticleReal const t_ms     = values_sum[12] / w_norm - dmean_t  * dmean_t;
        amrex::ParticleReal const s_px_ms  = values_sum[13] / w_norm - dmean_px * dmean_px;
        amrex::ParticleReal const s_py_ms  = values_sum[14] / w_norm - dmean_py * dmean_py;
        amrex::ParticleReal const s_pt_ms  = values_sum[15] / w_norm - dmean_pt * dmean_pt;
        amrex::ParticleReal const s_xpx    = values_sum[16] / w_norm - dmean_x  * dmean_px;
        amrex::ParticleReal const s_ypy    = values_sum[17] / w_norm - dmean_y  * dmean_py;
        amrex::ParticleReal const s_tpt    = values_sum[18] / w_norm - dmean_t  * dmean_pt;
        amrex::ParticleReal const s_xpt    = values_sum[19] / w_norm - dmean_x  * dmean_pt;
        amrex::ParticleReal const s_pxpt   = values_sum[20] / w_norm - dmean_px * dmean_pt;
        amrex::ParticleReal const s_ypt    = values_sum[21] / w_norm - dmean_y  * dmean_pt;
        amrex::ParticleReal const s_pypt   = values_sum[22] / w_norm - dmean_py * dmean_pt;
        amrex::ParticleReal const xy       = values_sum[23] / w_norm - dmean_x  * dmean_y;
        amrex::ParticleReal const s_xpy    = values_sum[24] / w_norm - dmean_x  * dmean_py;
        amrex::ParticleReal const xt       = values_sum[25] / w_norm - dmean_x  * dmean_t;
        amrex::ParticleReal const s_pxy    = values_sum[26] / w_norm - dmean_px * dmean_y;
        amrex::ParticleReal const s_pxpy   = values_sum[27] / w_norm - dmean_px * dmean_py;
        amrex::ParticleReal const s_pxt    = values_sum[28] / w_norm - dmean_px * dmean_t;
        amrex::ParticleReal const yt       = values_sum[29] / w_norm - dmean_y  * dmean_t;
        amrex::ParticleReal const s_pyt    = values_sum[30] / w_norm - dmean_py * dmean_t;
        amrex::ParticleReal const sx_ms    = values_sum[31] / w_norm - dmean_sx * dmean_sx;
        amrex::ParticleReal const sy_ms    = values_sum[32] / w_norm - dmean_sy * dmean_sy;
        amrex::ParticleReal const sz_ms    = values_sum[33] / w_norm - dmean_sz * dmean_sz;
        // the same for the momenta, undoing the shear: p = s_p + slope*q
        amrex::ParticleReal const px_ms = s_px_ms + slope_x * (2.0_prt * s_xpx + slope_x * x_ms);
        amrex::ParticleReal const py_ms = s_py_ms + slope_y * (2.0_prt * s_ypy + slope_y * y_ms);
        amrex::ParticleReal const pt_ms = s_pt_ms + slope_t * (2.0_prt * s_tpt + slope_t * t_ms);
        amrex::ParticleReal const tpt   = s_tpt + slope_t * t_ms;
        amrex::ParticleReal const xpt   = s_xpt + slope_t * xt;
        amrex::ParticleReal const ypt   = s_ypt + slope_t * yt;
        amrex::ParticleReal const pxpt  = s_pxpt + slope_x * xpt + slope_t * s_pxt;
        amrex::ParticleReal const pypt  = s_pypt + slope_y * ypt + slope_t * s_pyt;
        // beam charge
        amrex::ParticleReal const charge = q_C * w_sum;
        // standard deviations of positions
        amrex::ParticleReal const sigma_x = sqrt_clamped(x_ms);
        amrex::ParticleReal const sigma_y = sqrt_clamped(y_ms);
        amrex::ParticleReal const sigma_t = sqrt_clamped(t_ms);
        // standard deviations of momenta
        amrex::ParticleReal const sigma_px = sqrt_clamped(px_ms);
        amrex::ParticleReal const sigma_py = sqrt_clamped(py_ms);
        amrex::ParticleReal const sigma_pt = sqrt_clamped(pt_ms);
        // standard deviations of spin
        amrex::ParticleReal const sigma_sx = sqrt_clamped(sx_ms);
        amrex::ParticleReal const sigma_sy = sqrt_clamped(sy_ms);
        amrex::ParticleReal const sigma_sz = sqrt_clamped(sz_ms);
        // RMS emittances: a shear does not change them, and does not cancel
        amrex::ParticleReal const e2_x = x_ms*s_px_ms-s_xpx*s_xpx;
        amrex::ParticleReal const e2_y = y_ms*s_py_ms-s_ypy*s_ypy;
        amrex::ParticleReal const e2_t = t_ms*s_pt_ms-s_tpt*s_tpt;
        amrex::ParticleReal const emittance_x = sqrt_clamped(e2_x);
        amrex::ParticleReal const emittance_y = sqrt_clamped(e2_y);
        amrex::ParticleReal const emittance_t = sqrt_clamped(e2_t);
        // Dispersion and dispersive beam moments
        amrex::ParticleReal const dispersion_x = ((pt_ms > 0.0) ? (- xpt / pt_ms) : 0.0_prt);
        amrex::ParticleReal const dispersion_px = ((pt_ms > 0.0) ? (- pxpt / pt_ms) : 0.0_prt);
        amrex::ParticleReal const dispersion_y = ((pt_ms > 0.0) ? (- ypt / pt_ms) : 0.0_prt);
        amrex::ParticleReal const dispersion_py = ((pt_ms > 0.0) ? (- pypt / pt_ms) : 0.0_prt);
        // the correlation of the sheared momenta with pt, as in dispersion_px, dispersion_py
        amrex::ParticleReal const s_dispersion_px =
            (pt_ms > 0.0) ? (- (s_pxpt + slope_t * s_pxt) / pt_ms) : 0.0_prt;
        amrex::ParticleReal const s_dispersion_py =
            (pt_ms > 0.0) ? (- (s_pypt + slope_t * s_pyt) / pt_ms) : 0.0_prt;
        amrex::ParticleReal const x_msd = x_ms - pt_ms*dispersion_x*dispersion_x;
        amrex::ParticleReal const s_px_msd = s_px_ms - pt_ms*s_dispersion_px*s_dispersion_px;
        amrex::ParticleReal const s_xpx_d = s_xpx - pt_ms*dispersion_x*s_dispersion_px;
        amrex::ParticleReal const xpx_d = s_xpx_d + slope_x * x_msd;
        amrex::ParticleReal const emittance_xd = sqrt_clamped(x_msd*s_px_msd-s_xpx_d*s_xpx_d);
        amrex::ParticleReal const y_msd = y_ms - pt_ms*dispersion_y*dispersion_y;
        amrex::ParticleReal const s_py_msd = s_py_ms - pt_ms*s_dispersion_py*s_dispersion_py;
        amrex::ParticleReal const s_ypy_d = s_ypy - pt_ms*dispersion_y*s_dispersion_py;
        amrex::ParticleReal const ypy_d = s_ypy_d + slope_y * y_msd;
        amrex::ParticleReal const emittance_yd = sqrt_clamped(y_msd*s_py_msd-s_ypy_d*s_ypy_d);
        /* Courant-Snyder (Twiss) beta-function and alpha
         *
         * Both are ratios to an rms emittance, and are undefined where that emittance
         * vanishes: a single particle, a cold plane, or a longitudinal plane without
         * energy spread.
         */
        amrex::ParticleReal const beta_x = (emittance_xd > 0.0) ? x_msd / emittance_xd : nan;
        amrex::ParticleReal const beta_y = (emittance_yd > 0.0) ? y_msd / emittance_yd : nan;
        amrex::ParticleReal const beta_t = (emittance_t > 0.0) ? t_ms / emittance_t : nan;
        amrex::ParticleReal const alpha_x = (emittance_xd > 0.0) ? - xpx_d / emittance_xd : nan;
        amrex::ParticleReal const alpha_y = (emittance_yd > 0.0) ? - ypy_d / emittance_yd : nan;
        amrex::ParticleReal const alpha_t = (emittance_t > 0.0) ? - tpt / emittance_t : nan;

        // Calculate normalized emittances
        amrex::ParticleReal emittance_xn = emittance_x * bg;
        amrex::ParticleReal emittance_yn = emittance_y * bg;
        amrex::ParticleReal emittance_tn = emittance_t * bg;

        // Determine whether to calculate eigenemittances, and initialize
        amrex::ParmParse pp_diag("diag");
        bool compute_eigenemittances = false;
        pp_diag.queryAdd("eigenemittances", compute_eigenemittances);
        amrex::ParticleReal emittance_1 = emittance_xn;
        amrex::ParticleReal emittance_2 = emittance_yn;
        amrex::ParticleReal emittance_3 = emittance_tn;

        if (compute_eigenemittances) {
           // Store the covariance matrix in dynamical variables, with the sheared momenta:
           // a shear does not change the eigenemittances, and does not cancel
           amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1> Sigma;
           Sigma(1,1) = x_ms;
           Sigma(1,2) = s_xpx * bg;
           Sigma(1,3) = xy;
           Sigma(1,4) = s_xpy * bg;
           Sigma(1,5) = xt;
           Sigma(1,6) = s_xpt * bg;
           Sigma(2,1) = s_xpx * bg;
           Sigma(2,2) = s_px_ms * bg2;
           Sigma(2,3) = s_pxy * bg;
           Sigma(2,4) = s_pxpy * bg2;
           Sigma(2,5) = s_pxt * bg;
           Sigma(2,6) = s_pxpt * bg2;
           Sigma(3,1) = xy;
           Sigma(3,2) = s_pxy * bg;
           Sigma(3,3) = y_ms;
           Sigma(3,4) = s_ypy * bg;
           Sigma(3,5) = yt;
           Sigma(3,6) = s_ypt * bg;
           Sigma(4,1) = s_xpy * bg;
           Sigma(4,2) = s_pxpy * bg2;
           Sigma(4,3) = s_ypy * bg;
           Sigma(4,4) = s_py_ms * bg2;
           Sigma(4,5) = s_pyt * bg;
           Sigma(4,6) = s_pypt * bg2;
           Sigma(5,1) = xt;
           Sigma(5,2) = s_pxt * bg;
           Sigma(5,3) = yt;
           Sigma(5,4) = s_pyt * bg;
           Sigma(5,5) = t_ms;
           Sigma(5,6) = s_tpt * bg;
           Sigma(6,1) = s_xpt * bg;
           Sigma(6,2) = s_pxpt * bg2;
           Sigma(6,3) = s_ypt * bg;
           Sigma(6,4) = s_pypt * bg2;
           Sigma(6,5) = s_tpt * bg;
           Sigma(6,6) = s_pt_ms * bg2;
           // Calculate eigenemittances
           std::tuple<amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal> emittances = Eigenemittances(Sigma);
           emittance_1 = std::get<0>(emittances);
           emittance_2 = std::get<1>(emittances);
           emittance_3 = std::get<2>(emittances);
        }

        std::unordered_map<std::string, amrex::ParticleReal> data;
        data["mean_x"] = mean_x;
        data["min_x"] = min_x;
        data["max_x"] = max_x;
        data["mean_y"] = mean_y;
        data["min_y"] = min_y;
        data["max_y"] = max_y;
        data["mean_t"] = mean_t;
        data["min_t"] = min_t;
        data["max_t"] = max_t;
        data["sigma_x"] = sigma_x;
        data["sigma_y"] = sigma_y;
        data["sigma_t"] = sigma_t;
        data["mean_px"] = mean_px;
        data["min_px"] = min_px;
        data["max_px"] = max_px;
        data["mean_py"] = mean_py;
        data["min_py"] = min_py;
        data["max_py"] = max_py;
        data["mean_pt"] = mean_pt;
        data["min_pt"] = min_pt;
        data["max_pt"] = max_pt;
        data["sigma_px"] = sigma_px;
        data["sigma_py"] = sigma_py;
        data["sigma_pt"] = sigma_pt;
        // start deprecated attributes
        data["x_mean"] = mean_x;
        data["x_min"] = min_x;
        data["x_max"] = max_x;
        data["y_mean"] = mean_y;
        data["y_min"] = min_y;
        data["y_max"] = max_y;
        data["t_mean"] = mean_t;
        data["t_min"] = min_t;
        data["t_max"] = max_t;
        data["sig_x"] = sigma_x;
        data["sig_y"] = sigma_y;
        data["sig_t"] = sigma_t;
        data["px_mean"] = mean_px;
        data["px_min"] = min_px;
        data["px_max"] = max_px;
        data["py_mean"] = mean_py;
        data["py_min"] = min_py;
        data["py_max"] = max_py;
        data["pt_mean"] = mean_pt;
        data["pt_min"] = min_pt;
        data["pt_max"] = max_pt;
        data["sig_px"] = sigma_px;
        data["sig_py"] = sigma_py;
        data["sig_pt"] = sigma_pt;
        // end deprecated attributes
        data["emittance_x"] = emittance_x;
        data["emittance_y"] = emittance_y;
        data["emittance_t"] = emittance_t;
        data["alpha_x"] = alpha_x;
        data["alpha_y"] = alpha_y;
        data["alpha_t"] = alpha_t;
        data["beta_x"] = beta_x;
        data["beta_y"] = beta_y;
        data["beta_t"] = beta_t;
        data["dispersion_x"] = dispersion_x;
        data["dispersion_px"] = dispersion_px;
        data["dispersion_y"] = dispersion_y;
        data["dispersion_py"] = dispersion_py;
        data["emittance_xn"] = emittance_xn;
        data["emittance_yn"] = emittance_yn;
        data["emittance_tn"] = emittance_tn;
        if (compute_eigenemittances) {
           data["emittance_1"] = emittance_1;
           data["emittance_2"] = emittance_2;
           data["emittance_3"] = emittance_3;
        }
        data["charge_C"] = charge;
        data["mean_sx"] = mean_sx;
        data["mean_sy"] = mean_sy;
        data["mean_sz"] = mean_sz;
        data["sigma_sx"] = sigma_sx;
        data["sigma_sy"] = sigma_sy;
        data["sigma_sz"] = sigma_sz;

        /* Without a total weight to normalize by, or without a reference particle to
         * measure against, none of the moments above carries information: they are
         * the zeros that the substitute normalization produced. Report them as
         * "undefined" so that they cannot be mistaken for a measurement of a real,
         * if degenerate, beam.
         *
         * Two kinds of entry survive. The charge is exactly zero, as reduced. And the
         * minima and maxima are unweighted: they stay meaningful for a beam of zero
         * charge, and are undefined only without any particle to take them over - in
         * which case they hold the identity elements of their reductions.
         */
        if (!has_weight || !ref_is_initialized)
        {
            auto const is_defined = [has_particles](std::string const & moment)
            {
                if (moment == "charge_C") { return true; }
                return has_particles &&
                       (moment.starts_with("min_") || moment.starts_with("max_") ||
                        moment.ends_with("_min") || moment.ends_with("_max"));
            };

            for (auto & moment : data) {
                if (!is_defined(moment.first)) { moment.second = nan; }
            }
        }

        return data;
    }

    std::unordered_map<std::string, amrex::ParticleReal>
    reduced_beam_characteristics (Map6x6 const & cm, RefPart const & ref_part)
    {
        BL_PROFILE("impactx::diagnostics::reduced_beam_characteristics(cm)");

        using namespace amrex::literals; // for _prt

        auto const nan = std::numeric_limits<amrex::ParticleReal>::quiet_NaN();

        // reference particle relativistic beta*gamma
        amrex::ParticleReal const bg = ref_part.beta_gamma();
        amrex::ParticleReal const bg2 = bg*bg;

       // mean square and correlation values
        amrex::ParticleReal const x_ms   = cm(1,1);
        amrex::ParticleReal const y_ms   = cm(3,3);
        amrex::ParticleReal const t_ms   = cm(5,5);
        amrex::ParticleReal const px_ms  = cm(2,2);
        amrex::ParticleReal const py_ms  = cm(4,4);
        amrex::ParticleReal const pt_ms  = cm(6,6);
        amrex::ParticleReal const xpx    = cm(1,2);
        amrex::ParticleReal const ypy    = cm(3,4);
        amrex::ParticleReal const tpt    = cm(5,6);
        amrex::ParticleReal const xpt    = cm(1,6);
        amrex::ParticleReal const pxpt   = cm(2,6);
        amrex::ParticleReal const ypt    = cm(3,6);
        amrex::ParticleReal const pypt   = cm(4,6);
        amrex::ParticleReal const xy     = cm(1,3);
        amrex::ParticleReal const xpy    = cm(1,4);
        amrex::ParticleReal const xt     = cm(1,5);
        amrex::ParticleReal const pxy    = cm(2,3);
        amrex::ParticleReal const pxpy   = cm(2,4);
        amrex::ParticleReal const pxt    = cm(2,5);
        amrex::ParticleReal const yt     = cm(3,5);
        amrex::ParticleReal const pyt    = cm(4,5);
        // standard deviations of positions
        amrex::ParticleReal const sig_x = sqrt_clamped(x_ms);
        amrex::ParticleReal const sig_y = sqrt_clamped(y_ms);
        amrex::ParticleReal const sig_t = sqrt_clamped(t_ms);
        // standard deviations of momenta
        amrex::ParticleReal const sig_px = sqrt_clamped(px_ms);
        amrex::ParticleReal const sig_py = sqrt_clamped(py_ms);
        amrex::ParticleReal const sig_pt = sqrt_clamped(pt_ms);
        // RMS emittances
        amrex::ParticleReal const e2_x = x_ms*px_ms-xpx*xpx;
        amrex::ParticleReal const e2_y = y_ms*py_ms-ypy*ypy;
        amrex::ParticleReal const e2_t = t_ms*pt_ms-tpt*tpt;
        amrex::ParticleReal const emittance_x = sqrt_clamped(e2_x);
        amrex::ParticleReal const emittance_y = sqrt_clamped(e2_y);
        amrex::ParticleReal const emittance_t = sqrt_clamped(e2_t);
        // Dispersion and dispersive beam moments
        amrex::ParticleReal const dispersion_x = ((pt_ms > 0.0) ? (- xpt / pt_ms) : 0.0_prt);
        amrex::ParticleReal const dispersion_px = ((pt_ms > 0.0) ? (- pxpt / pt_ms) : 0.0_prt);
        amrex::ParticleReal const dispersion_y = ((pt_ms > 0.0) ? (- ypt / pt_ms) : 0.0_prt);
        amrex::ParticleReal const dispersion_py = ((pt_ms > 0.0) ? (- pypt / pt_ms) : 0.0_prt);
        amrex::ParticleReal const x_msd = x_ms - pt_ms*dispersion_x*dispersion_x;
        amrex::ParticleReal const px_msd = px_ms - pt_ms*dispersion_px*dispersion_px;
        amrex::ParticleReal const xpx_d = xpx - pt_ms*dispersion_x*dispersion_px;
        amrex::ParticleReal const emittance_xd = sqrt_clamped(x_msd*px_msd-xpx_d*xpx_d);
        amrex::ParticleReal const y_msd = y_ms - pt_ms*dispersion_y*dispersion_y;
        amrex::ParticleReal const py_msd = py_ms - pt_ms*dispersion_py*dispersion_py;
        amrex::ParticleReal const ypy_d = ypy - pt_ms*dispersion_y*dispersion_py;
        amrex::ParticleReal const emittance_yd = sqrt_clamped(y_msd*py_msd-ypy_d*ypy_d);
        /* Courant-Snyder (Twiss) beta-function and alpha
         *
         * Both are ratios to an rms emittance, and are undefined where that emittance
         * vanishes: a single particle, a cold plane, or a longitudinal plane without
         * energy spread.
         */
        amrex::ParticleReal const beta_x = (emittance_xd > 0.0) ? x_msd / emittance_xd : nan;
        amrex::ParticleReal const beta_y = (emittance_yd > 0.0) ? y_msd / emittance_yd : nan;
        amrex::ParticleReal const beta_t = (emittance_t > 0.0) ? t_ms / emittance_t : nan;
        amrex::ParticleReal const alpha_x = (emittance_xd > 0.0) ? - xpx_d / emittance_xd : nan;
        amrex::ParticleReal const alpha_y = (emittance_yd > 0.0) ? - ypy_d / emittance_yd : nan;
        amrex::ParticleReal const alpha_t = (emittance_t > 0.0) ? - tpt / emittance_t : nan;

        // Calculate normalized emittances
        amrex::ParticleReal emittance_xn = emittance_x * bg;
        amrex::ParticleReal emittance_yn = emittance_y * bg;
        amrex::ParticleReal emittance_tn = emittance_t * bg;

        // Determine whether to calculate eigenemittances, and initialize
        amrex::ParmParse pp_diag("diag");
        bool compute_eigenemittances = false;
        pp_diag.queryAdd("eigenemittances", compute_eigenemittances);
        amrex::ParticleReal emittance_1 = emittance_xn;
        amrex::ParticleReal emittance_2 = emittance_yn;
        amrex::ParticleReal emittance_3 = emittance_tn;

        if (compute_eigenemittances) {
           // Store the covariance matrix in dynamical variables:
           amrex::SmallMatrix<amrex::ParticleReal, 6, 6, amrex::Order::F, 1> Sigma;
           Sigma(1,1) = x_ms;
           Sigma(1,2) = xpx * bg;
           Sigma(1,3) = xy;
           Sigma(1,4) = xpy * bg;
           Sigma(1,5) = xt;
           Sigma(1,6) = xpt * bg;
           Sigma(2,1) = xpx * bg;
           Sigma(2,2) = px_ms * bg2;
           Sigma(2,3) = pxy * bg;
           Sigma(2,4) = pxpy * bg2;
           Sigma(2,5) = pxt * bg;
           Sigma(2,6) = pxpt * bg2;
           Sigma(3,1) = xy;
           Sigma(3,2) = pxy * bg;
           Sigma(3,3) = y_ms;
           Sigma(3,4) = ypy * bg;
           Sigma(3,5) = yt;
           Sigma(3,6) = ypt * bg;
           Sigma(4,1) = xpy * bg;
           Sigma(4,2) = pxpy * bg2;
           Sigma(4,3) = ypy * bg;
           Sigma(4,4) = py_ms * bg2;
           Sigma(4,5) = pyt * bg;
           Sigma(4,6) = pypt * bg2;
           Sigma(5,1) = xt;
           Sigma(5,2) = pxt * bg;
           Sigma(5,3) = yt;
           Sigma(5,4) = pyt * bg;
           Sigma(5,5) = t_ms;
           Sigma(5,6) = tpt * bg;
           Sigma(6,1) = xpt * bg;
           Sigma(6,2) = pxpt * bg2;
           Sigma(6,3) = ypt * bg;
           Sigma(6,4) = pypt * bg2;
           Sigma(6,5) = tpt * bg;
           Sigma(6,6) = pt_ms * bg2;
           // Calculate eigenemittances
           std::tuple<amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal> emittances = Eigenemittances(Sigma);
           emittance_1 = std::get<0>(emittances);
           emittance_2 = std::get<1>(emittances);
           emittance_3 = std::get<2>(emittances);
        }

        std::unordered_map<std::string, amrex::ParticleReal> data;
        data["mean_x"] = 0.0_prt;
        data["min_x"] = nan;
        data["max_x"] = nan;
        data["mean_y"] = 0.0_prt;
        data["min_y"] = nan;
        data["max_y"] = nan;
        data["mean_t"] = 0.0_prt;
        data["min_t"] = nan;
        data["max_t"] = nan;
        data["sigma_x"] = sig_x;
        data["sigma_y"] = sig_y;
        data["sigma_t"] = sig_t;
        data["mean_px"] = 0.0_prt;
        data["min_px"] = nan;
        data["max_px"] = nan;
        data["mean_py"] = 0.0_prt;
        data["min_py"] = nan;
        data["max_py"] = nan;
        data["mean_pt"] = 0.0_prt;
        data["min_pt"] = nan;
        data["max_pt"] = nan;
        data["sigma_px"] = sig_px;
        data["sigma_py"] = sig_py;
        data["sigma_pt"] = sig_pt;
        // start deprecated attributes
        data["x_mean"] = 0.0_prt;
        data["x_min"] = nan;
        data["x_max"] = nan;
        data["y_mean"] = 0.0_prt;
        data["y_min"] = nan;
        data["y_max"] = nan;
        data["t_mean"] = 0.0_prt;
        data["t_min"] = nan;
        data["t_max"] = nan;
        data["sig_x"] = sig_x;
        data["sig_y"] = sig_y;
        data["sig_t"] = sig_t;
        data["px_mean"] = 0.0_prt;
        data["px_min"] = nan;
        data["px_max"] = nan;
        data["py_mean"] = 0.0_prt;
        data["py_min"] = nan;
        data["py_max"] = nan;
        data["pt_mean"] = 0.0_prt;
        data["pt_min"] = nan;
        data["pt_max"] = nan;
        data["sig_px"] = sig_px;
        data["sig_py"] = sig_py;
        data["sig_pt"] = sig_pt;
        // end deprecated attributes
        data["emittance_x"] = emittance_x;
        data["emittance_y"] = emittance_y;
        data["emittance_t"] = emittance_t;
        data["alpha_x"] = alpha_x;
        data["alpha_y"] = alpha_y;
        data["alpha_t"] = alpha_t;
        data["beta_x"] = beta_x;
        data["beta_y"] = beta_y;
        data["beta_t"] = beta_t;
        data["dispersion_x"] = dispersion_x;
        data["dispersion_px"] = dispersion_px;
        data["dispersion_y"] = dispersion_y;
        data["dispersion_py"] = dispersion_py;
        data["emittance_xn"] = emittance_xn;
        data["emittance_yn"] = emittance_yn;
        data["emittance_tn"] = emittance_tn;
        if (compute_eigenemittances) {
           data["emittance_1"] = emittance_1;
           data["emittance_2"] = emittance_2;
           data["emittance_3"] = emittance_3;
        }
        data["charge_C"] = nan;  // TODO: with space charge
        data["mean_sx"] = nan;
        data["mean_sy"] = nan;
        data["mean_sz"] = nan;
        data["sigma_sx"] = nan;
        data["sigma_sy"] = nan;
        data["sigma_sz"] = nan;

        return data;
    }

} // namespace impactx::diagnostics
