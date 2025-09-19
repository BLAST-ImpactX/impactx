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
#include <AMReX_GpuQualifiers.H>        // for AMREX_GPU_DEVICE
#include <AMReX_ParallelDescriptor.H>   // for ParallelDescriptor
#include <AMReX_ParticleReduce.H>       // for ParticleReduce
#include <AMReX_REAL.H>                 // for ParticleReal
#include <AMReX_Reduce.H>               // for ReduceOps
#include <AMReX_SmallMatrix.H>          // for SmallMatrix
#include <AMReX_TypeList.H>             // for TypeMultiplier

#include <limits>


namespace impactx::diagnostics
{
namespace detail
{
    /** prepare reduction operations for calculation of mean and min/max values in 6D phase space */
    template<std::size_t num_red_ops_1_sum, std::size_t num_red_ops_1_min, std::size_t num_red_ops_1_max>
    constexpr auto
    get_reduce_ops1 ()
    {
        if constexpr (num_red_ops_1_min + num_red_ops_1_max > 0) {
            return amrex::TypeMultiplier<amrex::ReduceOps,
                amrex::ReduceOpSum[num_red_ops_1_sum],  // preparing mean values for w, x, y, t, px, py, pt
                amrex::ReduceOpMin[num_red_ops_1_min],  // preparing min values for x, y, t, px, py, pt
                amrex::ReduceOpMax[num_red_ops_1_max]   // preparing max values for x, y, t, px, py, pt
            >{};
        } else {
            return amrex::TypeMultiplier<amrex::ReduceOps,
                amrex::ReduceOpSum[num_red_ops_1_sum]  // preparing mean values for w, x, y, t, px, py, pt
            >{};
        }
    }
} // namespace detail

    template<MomentSelection selection>
    std::unordered_map<std::string, amrex::ParticleReal>
    reduced_beam_characteristics (ImpactXParticleContainer const & pc)
    {
        BL_PROFILE("impactx::diagnostics::reduced_beam_characteristics(pc)");

        // preparing to access reference particle data: RefPart
        RefPart const ref_part = pc.GetRefParticle();
        // reference particle charge in C
        amrex::ParticleReal const q_C = ref_part.charge;
        // reference particle relativistic beta*gamma
        amrex::ParticleReal const bg = ref_part.beta_gamma();
        amrex::ParticleReal const bg2 = bg*bg;

        // preparing access to particle data: SoA
        using PType = typename ImpactXParticleContainer::SuperParticleType;

        /* The variables below need to be static to work around an MSVC bug
         * https://stackoverflow.com/questions/55136414/constexpr-variable-captured-inside-lambda-loses-its-constexpr-ness
         */
        // numbers of same type reduction operations in first concurrent batch
        static constexpr std::size_t num_red_ops_1_sum = 7;  // summation
        static constexpr std::size_t num_red_ops_1_min = (selection & MomentSelection::MinMax) != MomentSelection::None ? 6 : 0;  // minimum
        static constexpr std::size_t num_red_ops_1_max = (selection & MomentSelection::MinMax) != MomentSelection::None ? 6 : 0;  // maximum

        // prepare reduction operations for calculation of mean and min/max values in 6D phase space
        auto reduce_ops_1 = detail::get_reduce_ops1<num_red_ops_1_sum, num_red_ops_1_min, num_red_ops_1_max>();
        using ReducedDataT1 = amrex::TypeMultiplier<amrex::ReduceData, amrex::ParticleReal[num_red_ops_1_sum + num_red_ops_1_min + num_red_ops_1_max]>;

        auto r1 = amrex::ParticleReduce<ReducedDataT1>(
            pc,
            [=] AMREX_GPU_DEVICE(const PType& p) noexcept -> ReducedDataT1::Type
            {
                // access particle position data
                const amrex::ParticleReal p_x = p.rdata(RealSoA::x);
                const amrex::ParticleReal p_y = p.rdata(RealSoA::y);
                const amrex::ParticleReal p_t = p.rdata(RealSoA::t);

                // access SoA particle momentum data and weighting
                const amrex::ParticleReal p_w = p.rdata(RealSoA::w);
                const amrex::ParticleReal p_px = p.rdata(RealSoA::px);
                const amrex::ParticleReal p_py = p.rdata(RealSoA::py);
                const amrex::ParticleReal p_pt = p.rdata(RealSoA::pt);

                // prepare mean position values
                const amrex::ParticleReal p_x_mean = p_x * p_w;
                const amrex::ParticleReal p_y_mean = p_y * p_w;
                const amrex::ParticleReal p_t_mean = p_t * p_w;

                const amrex::ParticleReal p_px_mean = p_px * p_w;
                const amrex::ParticleReal p_py_mean = p_py * p_w;
                const amrex::ParticleReal p_pt_mean = p_pt * p_w;

                if constexpr ((selection & MomentSelection::MinMax) != MomentSelection::None) {
                    return {p_w,
                            p_x_mean, p_y_mean, p_t_mean,
                            p_px_mean, p_py_mean, p_pt_mean,
                            p_x, p_y, p_t, p_px, p_py, p_pt,
                            p_x, p_y, p_t, p_px, p_py, p_pt};
                } else {
                    return {p_w,
                            p_x_mean, p_y_mean, p_t_mean,
                            p_px_mean, p_py_mean, p_pt_mean};
                }
            },
            reduce_ops_1
        );

        std::vector<amrex::ParticleReal> values_per_rank_1st(num_red_ops_1_sum);

        /* contains in this order:
         * w, x_mean, y_mean, t_mean
         * px_mean, py_mean, pt_mean
         */
        amrex::constexpr_for<0, num_red_ops_1_sum> ([&](auto i) {
            values_per_rank_1st[i] = amrex::get<i>(r1);
        });

        // reduced sum over mpi ranks (allreduce)
        amrex::ParallelAllReduce::Sum(
            values_per_rank_1st.data(),
            values_per_rank_1st.size(),
            amrex::ParallelDescriptor::Communicator()
        );

        amrex::ParticleReal const w_sum   = values_per_rank_1st.at(0);
        amrex::ParticleReal const x_mean  = values_per_rank_1st.at(1) /= w_sum;
        amrex::ParticleReal const y_mean  = values_per_rank_1st.at(2) /= w_sum;
        amrex::ParticleReal const t_mean  = values_per_rank_1st.at(3) /= w_sum;
        amrex::ParticleReal const px_mean = values_per_rank_1st.at(4) /= w_sum;
        amrex::ParticleReal const py_mean = values_per_rank_1st.at(5) /= w_sum;
        amrex::ParticleReal const pt_mean = values_per_rank_1st.at(6) /= w_sum;

        std::vector<amrex::ParticleReal> values_per_rank_min(num_red_ops_1_min);

        /* contains in this order:
         * x_min, y_min, t_min
         * px_min, py_min, pt_min
         */
        amrex::constexpr_for<0, num_red_ops_1_min> ([&](auto i) {
            constexpr std::size_t idx = i + num_red_ops_1_sum;
            values_per_rank_min[i] = amrex::get<idx>(r1);
        });

        std::vector<amrex::ParticleReal> values_per_rank_max(num_red_ops_1_max);

        /* contains in this order:
         * x_max, y_max, t_max
         * px_max, py_max, pt_max
         */
        amrex::constexpr_for<0, num_red_ops_1_max> ([&](auto i) {
            constexpr std::size_t idx = i + num_red_ops_1_sum + num_red_ops_1_min;
            values_per_rank_max[i] = amrex::get<idx>(r1);
        });

        // reduced sum over mpi ranks (allreduce)
        amrex::ParallelAllReduce::Min(
                values_per_rank_min.data(),
                values_per_rank_min.size(),
                amrex::ParallelDescriptor::Communicator()
        );

        // reduced sum over mpi ranks (allreduce)
        amrex::ParallelAllReduce::Max(
                values_per_rank_max.data(),
                values_per_rank_max.size(),
                amrex::ParallelDescriptor::Communicator()
        );

        /* The variable below needs to be static to work around an MSVC bug
         * https://stackoverflow.com/questions/55136414/constexpr-variable-captured-inside-lambda-loses-its-constexpr-ness
         */
        // number of reduction operations in second concurrent batch
        static constexpr std::size_t num_red_ops_2 =
            ((selection & MomentSelection::Moment) != MomentSelection::None ? 6 : 0) +      // x_ms, y_ms, t_ms, px_ms, py_ms, pt_ms
            ((selection & MomentSelection::Covariance) != MomentSelection::None ? 3 : 0) +  // xpx, ypy, tpt
            ((selection & MomentSelection::Dispersion) != MomentSelection::None ? 4 : 0) +  // xpt, pxpt, ypt, pypt
            ((selection & MomentSelection::Covariance) != MomentSelection::None ? 8 : 0) +  // xy, xpy, xt, pxy, pxpy, pxt, yt, pyt
            ((selection & MomentSelection::Moment) != MomentSelection::None ? 1 : 0);       // charge

        // prepare reduction operations for calculation of mean square and correlation values
        amrex::TypeMultiplier<amrex::ReduceOps, amrex::ReduceOpSum[num_red_ops_2]> reduce_ops_2;
        using ReducedDataT2 = amrex::TypeMultiplier<amrex::ReduceData, amrex::ParticleReal[num_red_ops_2]>;

        auto r2 = amrex::ParticleReduce<ReducedDataT2>(
                pc,
                [=] AMREX_GPU_DEVICE(const PType& p) noexcept
            -> ReducedDataT2::Type
            {
                // access SoA particle momentum data and weighting
                const amrex::ParticleReal p_w = p.rdata(RealSoA::w);
                const amrex::ParticleReal p_px = p.rdata(RealSoA::px);
                const amrex::ParticleReal p_py = p.rdata(RealSoA::py);
                const amrex::ParticleReal p_pt = p.rdata(RealSoA::pt);
                // access position data
                const amrex::ParticleReal p_x = p.rdata(RealSoA::x);
                const amrex::ParticleReal p_y = p.rdata(RealSoA::y);
                const amrex::ParticleReal p_t = p.rdata(RealSoA::t);

                // Initialize return values
                amrex::ParticleReal x_ms_val = 0.0, y_ms_val = 0.0, t_ms_val = 0.0;
                amrex::ParticleReal px_ms_val = 0.0, py_ms_val = 0.0, pt_ms_val = 0.0;
                amrex::ParticleReal xpx_val = 0.0, ypy_val = 0.0, tpt_val = 0.0;
                amrex::ParticleReal xpt_val = 0.0, pxpt_val = 0.0, ypt_val = 0.0, pypt_val = 0.0;
                amrex::ParticleReal xy_val = 0.0, xpy_val = 0.0, xt_val = 0.0, pxy_val = 0.0;
                amrex::ParticleReal pxpy_val = 0.0, pxt_val = 0.0, yt_val = 0.0, pyt_val = 0.0;
                amrex::ParticleReal charge_val = 0.0;

                // Moment calculations (mean square values)
                if constexpr ((selection & MomentSelection::Moment) != MomentSelection::None) {
                // prepare mean square for positions
                    x_ms_val = (p_x-x_mean)*(p_x-x_mean)*p_w;
                    y_ms_val = (p_y-y_mean)*(p_y-y_mean)*p_w;
                    t_ms_val = (p_t-t_mean)*(p_t-t_mean)*p_w;
                // prepare mean square for momenta
                    px_ms_val = (p_px-px_mean)*(p_px-px_mean)*p_w;
                    py_ms_val = (p_py-py_mean)*(p_py-py_mean)*p_w;
                    pt_ms_val = (p_pt-pt_mean)*(p_pt-pt_mean)*p_w;
                }

                // Covariance calculations (position-momentum correlations)
                if constexpr ((selection & MomentSelection::Covariance) != MomentSelection::None) {
                    xpx_val = (p_x-x_mean)*(p_px-px_mean)*p_w;
                    ypy_val = (p_y-y_mean)*(p_py-py_mean)*p_w;
                    tpt_val = (p_t-t_mean)*(p_pt-pt_mean)*p_w;
                }

                // Dispersion calculations
                if constexpr ((selection & MomentSelection::Dispersion) != MomentSelection::None) {
                    xpt_val = (p_x-x_mean)*(p_pt-pt_mean)*p_w;
                    pxpt_val = (p_px-px_mean)*(p_pt-pt_mean)*p_w;
                    ypt_val = (p_y-y_mean)*(p_pt-pt_mean)*p_w;
                    pypt_val = (p_py-py_mean)*(p_pt-pt_mean)*p_w;
                }

                // Additional cross-plane correlations
                if constexpr ((selection & MomentSelection::Covariance) != MomentSelection::None) {
                    xy_val = (p_x-x_mean)*(p_y-y_mean)*p_w;
                    xpy_val = (p_x-x_mean)*(p_py-py_mean)*p_w;
                    xt_val = (p_x-x_mean)*(p_t-t_mean)*p_w;
                    pxy_val = (p_px-px_mean)*(p_y-y_mean)*p_w;
                    pxpy_val = (p_px-px_mean)*(p_py-py_mean)*p_w;
                    pxt_val = (p_px-px_mean)*(p_t-t_mean)*p_w;
                    yt_val = (p_y-y_mean)*(p_t-t_mean)*p_w;
                    pyt_val = (p_py-py_mean)*(p_t-t_mean)*p_w;
                }

                // Charge calculation
                if constexpr ((selection & MomentSelection::Moment) != MomentSelection::None) {
                    charge_val = q_C*p_w;
                }

                // Build return tuple using conditional assignment - much more readable!
                // Define all possible values in a structured way
                constexpr bool has_moment = (selection & MomentSelection::Moment) != MomentSelection::None;
                constexpr bool has_covariance = (selection & MomentSelection::Covariance) != MomentSelection::None;
                constexpr bool has_dispersion = (selection & MomentSelection::Dispersion) != MomentSelection::None;

                // Helper function to build the appropriate tuple based on selection
                if constexpr (has_moment && has_covariance && has_dispersion) {
                    // All flags: 22 values
                    return {x_ms_val, y_ms_val, t_ms_val, px_ms_val, py_ms_val, pt_ms_val,
                            xpx_val, ypy_val, tpt_val,
                            xpt_val, pxpt_val, ypt_val, pypt_val,
                            xy_val, xpy_val, xt_val, pxy_val, pxpy_val, pxt_val, yt_val, pyt_val,
                            charge_val};
                } else if constexpr (has_moment && has_covariance) {
                    // Moment + Covariance: 17 values
                    return {x_ms_val, y_ms_val, t_ms_val, px_ms_val, py_ms_val, pt_ms_val,
                            xpx_val, ypy_val, tpt_val,
                            xy_val, xpy_val, xt_val, pxy_val, pxpy_val, pxt_val, yt_val, pyt_val,
                            charge_val};
                } else if constexpr (has_moment && has_dispersion) {
                    // Moment + Dispersion: 10 values
                    return {x_ms_val, y_ms_val, t_ms_val, px_ms_val, py_ms_val, pt_ms_val,
                            xpt_val, pxpt_val, ypt_val, pypt_val,
                            charge_val};
                } else if constexpr (has_moment) {
                    // Moment only: 7 values
                    return {x_ms_val, y_ms_val, t_ms_val, px_ms_val, py_ms_val, pt_ms_val,
                            charge_val};
                } else if constexpr (has_covariance && has_dispersion) {
                    // Covariance + Dispersion: 15 values
                    return {xpx_val, ypy_val, tpt_val,
                            xpt_val, pxpt_val, ypt_val, pypt_val,
                            xy_val, xpy_val, xt_val, pxy_val, pxpy_val, pxt_val, yt_val, pyt_val};
                } else if constexpr (has_covariance) {
                    // Covariance only: 11 values
                    return {xpx_val, ypy_val, tpt_val,
                            xy_val, xpy_val, xt_val, pxy_val, pxpy_val, pxt_val, yt_val, pyt_val};
                } else if constexpr (has_dispersion) {
                    // Dispersion only: 4 values
                    return {xpt_val, pxpt_val, ypt_val, pypt_val};
                } else {
                    // No flags: 0 values (should not happen)
                    return {};
                }
            },
                reduce_ops_2
        );

        // Extract values conditionally based on selection - much more readable!
        // Define flag checks once at the top (similar to num_red_ops_2 calculation)
        constexpr bool has_moment = (selection & MomentSelection::Moment) != MomentSelection::None;
        constexpr bool has_covariance = (selection & MomentSelection::Covariance) != MomentSelection::None;
        constexpr bool has_dispersion = (selection & MomentSelection::Dispersion) != MomentSelection::None;

        // Helper function to get the correct index based on what's enabled
        constexpr auto get_index = [](auto base, auto offset) -> std::size_t {
            return base + offset;
        };

        // Calculate base indices for each section based on what's enabled
        constexpr std::size_t moment_base = 0;
        constexpr std::size_t covar_base = has_moment ? 6 : 0;
        constexpr std::size_t disp_base = has_moment ? (has_covariance ? 9 : 6) : (has_covariance ? 3 : 0);
        constexpr std::size_t cross_base = has_moment ? (has_covariance ? (has_dispersion ? 13 : 9) : 9) :
                                                       (has_covariance ? (has_dispersion ? 7 : 3) : 3);
        constexpr std::size_t charge_idx = has_moment ? (has_covariance ? (has_dispersion ? 21 : 17) :
                                                        (has_dispersion ? 10 : 6)) : 0;

        // Moment calculations (mean square values)
        amrex::ParticleReal x_ms = 0.0, y_ms = 0.0, t_ms = 0.0;
        amrex::ParticleReal px_ms = 0.0, py_ms = 0.0, pt_ms = 0.0;
        if constexpr (has_moment) {
            x_ms = amrex::get<get_index(moment_base, 0)>(r2) /= w_sum;
            y_ms = amrex::get<get_index(moment_base, 1)>(r2) /= w_sum;
            t_ms = amrex::get<get_index(moment_base, 2)>(r2) /= w_sum;
            px_ms = amrex::get<get_index(moment_base, 3)>(r2) /= w_sum;
            py_ms = amrex::get<get_index(moment_base, 4)>(r2) /= w_sum;
            pt_ms = amrex::get<get_index(moment_base, 5)>(r2) /= w_sum;
        }

        // Covariance calculations (position-momentum correlations)
        amrex::ParticleReal xpx = 0.0, ypy = 0.0, tpt = 0.0;
        if constexpr (has_covariance) {
            xpx = amrex::get<get_index(covar_base, 0)>(r2) /= w_sum;
            ypy = amrex::get<get_index(covar_base, 1)>(r2) /= w_sum;
            tpt = amrex::get<get_index(covar_base, 2)>(r2) /= w_sum;
        }

        // Dispersion calculations
        amrex::ParticleReal xpt = 0.0, pxpt = 0.0, ypt = 0.0, pypt = 0.0;
        if constexpr (has_dispersion) {
            xpt = amrex::get<get_index(disp_base, 0)>(r2) /= w_sum;
            pxpt = amrex::get<get_index(disp_base, 1)>(r2) /= w_sum;
            ypt = amrex::get<get_index(disp_base, 2)>(r2) /= w_sum;
            pypt = amrex::get<get_index(disp_base, 3)>(r2) /= w_sum;
        }

        // Additional cross-plane correlations
        amrex::ParticleReal xy = 0.0, xpy = 0.0, xt = 0.0, pxy = 0.0, pxpy = 0.0, pxt = 0.0, yt = 0.0, pyt = 0.0;
        if constexpr (has_covariance) {
            xy = amrex::get<get_index(cross_base, 0)>(r2) /= w_sum;
            xpy = amrex::get<get_index(cross_base, 1)>(r2) /= w_sum;
            xt = amrex::get<get_index(cross_base, 2)>(r2) /= w_sum;
            pxy = amrex::get<get_index(cross_base, 3)>(r2) /= w_sum;
            pxpy = amrex::get<get_index(cross_base, 4)>(r2) /= w_sum;
            pxt = amrex::get<get_index(cross_base, 5)>(r2) /= w_sum;
            yt = amrex::get<get_index(cross_base, 6)>(r2) /= w_sum;
            pyt = amrex::get<get_index(cross_base, 7)>(r2) /= w_sum;
        }

        // Charge calculation
        amrex::ParticleReal charge = 0.0;
        if constexpr (has_moment) {
            charge = amrex::get<charge_idx>(r2);
        }

        [[maybe_unused]] amrex::ParticleReal x_min, y_min, t_min,
                                             px_min, py_min, pt_min,
                                             x_max, y_max, t_max,
                                             px_max, py_max, pt_max;
        if constexpr ((selection & MomentSelection::MinMax) != MomentSelection::None) {
            // minimum values
            x_min = values_per_rank_min.at(0);
            y_min = values_per_rank_min.at(1);
            t_min = values_per_rank_min.at(2);
            px_min = values_per_rank_min.at(3);
            py_min = values_per_rank_min.at(4);
            pt_min = values_per_rank_min.at(5);
            // maximum values
            x_max = values_per_rank_max.at(0);
            y_max = values_per_rank_max.at(1);
            t_max = values_per_rank_max.at(2);
            px_max = values_per_rank_max.at(3);
            py_max = values_per_rank_max.at(4);
            pt_max = values_per_rank_max.at(5);
        }
        // standard deviations of positions
        amrex::ParticleReal sig_x = 0.0, sig_y = 0.0, sig_t = 0.0;
        if constexpr ((selection & MomentSelection::Moment) != MomentSelection::None) {
            sig_x = std::sqrt(x_ms);
            sig_y = std::sqrt(y_ms);
            sig_t = std::sqrt(t_ms);
        }
        // standard deviations of momenta
        amrex::ParticleReal sig_px = 0.0, sig_py = 0.0, sig_pt = 0.0;
        if constexpr ((selection & MomentSelection::Moment) != MomentSelection::None) {
            sig_px = std::sqrt(px_ms);
            sig_py = std::sqrt(py_ms);
            sig_pt = std::sqrt(pt_ms);
        }
        // RMS emittances
        amrex::ParticleReal emittance_x = 0.0, emittance_y = 0.0, emittance_t = 0.0;
        if constexpr ((selection & MomentSelection::Covariance) != MomentSelection::None) {
        amrex::ParticleReal const e2_x = x_ms*px_ms-xpx*xpx;
        amrex::ParticleReal const e2_y = y_ms*py_ms-ypy*ypy;
        amrex::ParticleReal const e2_t = t_ms*pt_ms-tpt*tpt;
            emittance_x = (e2_x > 0.0)? std::sqrt(e2_x) : 0.0;
            emittance_y = (e2_y > 0.0)? std::sqrt(e2_y) : 0.0;
            emittance_t = (e2_t > 0.0)? std::sqrt(e2_t) : 0.0;
        }
        // Dispersion and dispersive beam moments
        amrex::ParticleReal dispersion_x = 0.0, dispersion_px = 0.0, dispersion_y = 0.0, dispersion_py = 0.0;
        amrex::ParticleReal emittance_xd = 0.0, emittance_yd = 0.0;
        if constexpr ((selection & MomentSelection::Dispersion) != MomentSelection::None) {
            dispersion_x = ((pt_ms > 0.0) ? (- xpt / pt_ms) : 0.0);
            dispersion_px = ((pt_ms > 0.0) ? (- pxpt / pt_ms) : 0.0);
            dispersion_y = ((pt_ms > 0.0) ? (- ypt / pt_ms) : 0.0);
            dispersion_py = ((pt_ms > 0.0) ? (- pypt / pt_ms) : 0.0);
        }
        if constexpr (((selection & MomentSelection::Dispersion) != MomentSelection::None) && ((selection & MomentSelection::Covariance) != MomentSelection::None)) {
        amrex::ParticleReal const x_msd = x_ms - pt_ms*dispersion_x*dispersion_x;
        amrex::ParticleReal const px_msd = px_ms - pt_ms*dispersion_px*dispersion_px;
        amrex::ParticleReal const xpx_d = xpx - pt_ms*dispersion_x*dispersion_px;
            emittance_xd = std::sqrt(x_msd*px_msd-xpx_d*xpx_d);
        amrex::ParticleReal const y_msd = y_ms - pt_ms*dispersion_y*dispersion_y;
        amrex::ParticleReal const py_msd = py_ms - pt_ms*dispersion_py*dispersion_py;
        amrex::ParticleReal const ypy_d = ypy - pt_ms*dispersion_y*dispersion_py;
            emittance_yd = std::sqrt(y_msd*py_msd-ypy_d*ypy_d);
        }
        // Courant-Snyder (Twiss) beta-function and alpha
        amrex::ParticleReal beta_x = 0.0, beta_y = 0.0, beta_t = 0.0;
        amrex::ParticleReal alpha_x = 0.0, alpha_y = 0.0, alpha_t = 0.0;
        if constexpr ((selection & MomentSelection::Covariance) != MomentSelection::None) {
            if constexpr ((selection & MomentSelection::Dispersion) != MomentSelection::None) {
                beta_x = (emittance_xd > 0.0) ? (x_ms - pt_ms*dispersion_x*dispersion_x) / emittance_xd : 0.0;
                beta_y = (emittance_yd > 0.0) ? (y_ms - pt_ms*dispersion_y*dispersion_y) / emittance_yd : 0.0;
                alpha_x = (emittance_xd > 0.0) ? - (xpx - pt_ms*dispersion_x*dispersion_px) / emittance_xd : 0.0;
                alpha_y = (emittance_yd > 0.0) ? - (ypy - pt_ms*dispersion_y*dispersion_py) / emittance_yd : 0.0;
            }
            beta_t = (emittance_t > 0.0) ? t_ms / emittance_t : 0.0;
            alpha_t = (emittance_t > 0.0) ? - tpt / emittance_t : 0.0;
        }

        // Calculate normalized emittances
        amrex::ParticleReal emittance_xn = 0.0, emittance_yn = 0.0, emittance_tn = 0.0;
        if constexpr ((selection & MomentSelection::Covariance) != MomentSelection::None) {
            emittance_xn = emittance_x * bg;
            emittance_yn = emittance_y * bg;
            emittance_tn = emittance_t * bg;
        }

        // Calculate eigenemittances conditionally based on selection
        amrex::ParticleReal emittance_1 = emittance_xn;
        amrex::ParticleReal emittance_2 = emittance_yn;
        amrex::ParticleReal emittance_3 = emittance_tn;

        if constexpr ((selection & MomentSelection::Eigenemittance) != MomentSelection::None) {
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

        // Always add mean values (they're needed for basic calculations)
        data["x_mean"] = x_mean;
        data["y_mean"] = y_mean;
        data["t_mean"] = t_mean;
        data["px_mean"] = px_mean;
        data["py_mean"] = py_mean;
        data["pt_mean"] = pt_mean;

        // Add min/max values conditionally
        if constexpr ((selection & MomentSelection::MinMax) != MomentSelection::None) {
            data["x_min"] = x_min;
            data["x_max"] = x_max;
            data["y_min"] = y_min;
            data["y_max"] = y_max;
            data["t_min"] = t_min;
            data["t_max"] = t_max;
            data["px_min"] = px_min;
            data["px_max"] = px_max;
            data["py_min"] = py_min;
            data["py_max"] = py_max;
            data["pt_min"] = pt_min;
            data["pt_max"] = pt_max;
        }

        // Add moment-based values conditionally
        if constexpr ((selection & MomentSelection::Moment) != MomentSelection::None) {
            data["sig_x"] = sig_x;
            data["sig_y"] = sig_y;
            data["sig_t"] = sig_t;
        data["sig_px"] = sig_px;
        data["sig_py"] = sig_py;
        data["sig_pt"] = sig_pt;
            data["charge_C"] = charge;
        }

        // Add covariance-based values conditionally
        if constexpr ((selection & MomentSelection::Covariance) != MomentSelection::None) {
        data["emittance_x"] = emittance_x;
        data["emittance_y"] = emittance_y;
        data["emittance_t"] = emittance_t;
        data["alpha_x"] = alpha_x;
        data["alpha_y"] = alpha_y;
        data["alpha_t"] = alpha_t;
        data["beta_x"] = beta_x;
        data["beta_y"] = beta_y;
        data["beta_t"] = beta_t;
            data["emittance_xn"] = emittance_xn;
            data["emittance_yn"] = emittance_yn;
            data["emittance_tn"] = emittance_tn;
        }

        // Add dispersion-based values conditionally
        if constexpr ((selection & MomentSelection::Dispersion) != MomentSelection::None) {
        data["dispersion_x"] = dispersion_x;
        data["dispersion_px"] = dispersion_px;
        data["dispersion_y"] = dispersion_y;
        data["dispersion_py"] = dispersion_py;
        }

        // Add eigenemittance values conditionally
        if constexpr ((selection & MomentSelection::Eigenemittance) != MomentSelection::None) {
           data["emittance_1"] = emittance_1;
           data["emittance_2"] = emittance_2;
           data["emittance_3"] = emittance_3;
        }

        return data;
    }

    std::unordered_map<std::string, amrex::ParticleReal>
    reduced_beam_characteristics (Map6x6 const & cm, RefPart const & ref_part)
    {
        BL_PROFILE("impactx::diagnostics::reduced_beam_characteristics(cm)");

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
        amrex::ParticleReal const sig_x = std::sqrt(x_ms);
        amrex::ParticleReal const sig_y = std::sqrt(y_ms);
        amrex::ParticleReal const sig_t = std::sqrt(t_ms);
        // standard deviations of momenta
        amrex::ParticleReal const sig_px = std::sqrt(px_ms);
        amrex::ParticleReal const sig_py = std::sqrt(py_ms);
        amrex::ParticleReal const sig_pt = std::sqrt(pt_ms);
        // RMS emittances
        amrex::ParticleReal const e2_x = x_ms*px_ms-xpx*xpx;
        amrex::ParticleReal const e2_y = y_ms*py_ms-ypy*ypy;
        amrex::ParticleReal const e2_t = t_ms*pt_ms-tpt*tpt;
        amrex::ParticleReal const emittance_x = (e2_x > 0.0)? std::sqrt(e2_x) : 0.0;
        amrex::ParticleReal const emittance_y = (e2_y > 0.0)? std::sqrt(e2_y) : 0.0;
        amrex::ParticleReal const emittance_t = (e2_t > 0.0)? std::sqrt(e2_t) : 0.0;
        // Dispersion and dispersive beam moments
        amrex::ParticleReal const dispersion_x = ((pt_ms > 0.0) ? (- xpt / pt_ms) : 0.0);
        amrex::ParticleReal const dispersion_px = ((pt_ms > 0.0) ? (- pxpt / pt_ms) : 0.0);
        amrex::ParticleReal const dispersion_y = ((pt_ms > 0.0) ? (- ypt / pt_ms) : 0.0);
        amrex::ParticleReal const dispersion_py = ((pt_ms > 0.0) ? (- pypt / pt_ms) : 0.0);
        amrex::ParticleReal const x_msd = x_ms - pt_ms*dispersion_x*dispersion_x;
        amrex::ParticleReal const px_msd = px_ms - pt_ms*dispersion_px*dispersion_px;
        amrex::ParticleReal const xpx_d = xpx - pt_ms*dispersion_x*dispersion_px;
        amrex::ParticleReal const emittance_xd = std::sqrt(x_msd*px_msd-xpx_d*xpx_d);
        amrex::ParticleReal const y_msd = y_ms - pt_ms*dispersion_y*dispersion_y;
        amrex::ParticleReal const py_msd = py_ms - pt_ms*dispersion_py*dispersion_py;
        amrex::ParticleReal const ypy_d = ypy - pt_ms*dispersion_y*dispersion_py;
        amrex::ParticleReal const emittance_yd = std::sqrt(y_msd*py_msd-ypy_d*ypy_d);
        // Courant-Snyder (Twiss) beta-function
        amrex::ParticleReal const beta_x = x_msd / emittance_xd;
        amrex::ParticleReal const beta_y = y_msd / emittance_yd;
        amrex::ParticleReal const beta_t = t_ms / emittance_t;
        // Courant-Snyder (Twiss) alpha
        amrex::ParticleReal const alpha_x = - xpx_d / emittance_xd;
        amrex::ParticleReal const alpha_y = - ypy_d / emittance_yd;
        amrex::ParticleReal const alpha_t = - tpt / emittance_t;

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

        auto const nan = std::numeric_limits<amrex::ParticleReal>::quiet_NaN();

        std::unordered_map<std::string, amrex::ParticleReal> data;
        data["x_mean"] = nan;
        data["x_min"] = nan;
        data["x_max"] = nan;
        data["y_mean"] = nan;
        data["y_min"] = nan;
        data["y_max"] = nan;
        data["t_mean"] = nan;
        data["t_min"] = nan;
        data["t_max"] = nan;
        data["sig_x"] = sig_x;
        data["sig_y"] = sig_y;
        data["sig_t"] = sig_t;
        data["px_mean"] = nan;
        data["px_min"] = nan;
        data["px_max"] = nan;
        data["py_mean"] = nan;
        data["py_min"] = nan;
        data["py_max"] = nan;
        data["pt_mean"] = nan;
        data["pt_min"] = nan;
        data["pt_max"] = nan;
        data["sig_px"] = sig_px;
        data["sig_py"] = sig_py;
        data["sig_pt"] = sig_pt;
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

        return data;
    }

    namespace detail {
        // Elegant dispatcher using a constexpr array of valid combinations
        // This only instantiates templates for valid flag combinations
        template<auto Selection>
        auto dispatch_to_template(ImpactXParticleContainer const & pc) {
            return reduced_beam_characteristics<static_cast<MomentSelection>(Selection)>(pc);
        }

        // All valid combinations (always including Moment flag)
        constexpr std::array<std::underlying_type_t<MomentSelection>, 16> valid_combinations = {{
            0x01, // Moment
            0x03, // Moment|MinMax
            0x05, // Moment|Dispersion
            0x07, // Moment|MinMax|Dispersion
            0x09, // Moment|Covariance
            0x0B, // Moment|MinMax|Covariance
            0x0D, // Moment|Dispersion|Covariance
            0x0F, // Moment|MinMax|Dispersion|Covariance
            0x11, // Moment|Eigenemittance
            0x13, // Moment|MinMax|Eigenemittance
            0x15, // Moment|Dispersion|Eigenemittance
            0x17, // Moment|MinMax|Dispersion|Eigenemittance
            0x19, // Moment|Covariance|Eigenemittance
            0x1B, // Moment|MinMax|Covariance|Eigenemittance
            0x1D, // Moment|Dispersion|Covariance|Eigenemittance
            0x1F  // All flags
        }};

        // Recursive template dispatcher that only tries valid combinations
        template<std::size_t Index = 0>
        auto try_valid_dispatch(ImpactXParticleContainer const & pc, std::underlying_type_t<MomentSelection> selection) {
            if constexpr (Index < valid_combinations.size()) {
                constexpr auto test_value = valid_combinations[Index];
                if (selection == test_value) {
                    return dispatch_to_template<test_value>(pc);
                } else {
                    return try_valid_dispatch<Index + 1>(pc, selection);
                }
            } else {
                // Fallback - should never reach here with valid input
                return dispatch_to_template<0x1F>(pc);
            }
        }

        // Helper to normalize selection: always include Moment, mask out invalid bits
        constexpr auto normalize_selection(MomentSelection selection) {
            constexpr auto all_valid_flags = static_cast<std::underlying_type_t<MomentSelection>>(
                MomentSelection::Moment | MomentSelection::MinMax | MomentSelection::Dispersion |
                MomentSelection::Covariance | MomentSelection::Eigenemittance
            );
            auto normalized = static_cast<std::underlying_type_t<MomentSelection>>(selection) & all_valid_flags;
            // Always include Moment flag
            normalized |= static_cast<std::underlying_type_t<MomentSelection>>(MomentSelection::Moment);
            return normalized;
        }
    }

    std::unordered_map<std::string, amrex::ParticleReal>
    reduced_beam_characteristics (ImpactXParticleContainer const & pc, MomentSelection selection)
    {
        // Normalize selection and dispatch using valid combinations only
        auto normalized = detail::normalize_selection(selection);
        return detail::try_valid_dispatch(pc, normalized);
    }

} // namespace impactx::diagnostics
