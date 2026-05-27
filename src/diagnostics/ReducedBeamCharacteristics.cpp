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
        static constexpr std::size_t num_red_ops_1_sum = 10;  // summation
        static constexpr std::size_t num_red_ops_1_min = 6;  // minimum
        static constexpr std::size_t num_red_ops_1_max = 6;  // maximum

        // prepare reduction operations for calculation of mean and min/max values in 6D phase space
        amrex::TypeMultiplier<amrex::ReduceOps,
            amrex::ReduceOpSum[num_red_ops_1_sum],  // preparing mean values for w, x, y, t, px, py, pt, sx, sy, sz
            amrex::ReduceOpMin[num_red_ops_1_min],  // preparing min values for x, y, t, px, py, pt
            amrex::ReduceOpMax[num_red_ops_1_max]   // preparing max values for x, y, t, px, py, pt
        > reduce_ops_1;
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

                // access particle spin data
                const amrex::ParticleReal p_sx = p.rdata(RealSoA::sx);
                const amrex::ParticleReal p_sy = p.rdata(RealSoA::sy);
                const amrex::ParticleReal p_sz = p.rdata(RealSoA::sz);

                // prepare mean position values
                const amrex::ParticleReal p_x_mean = p_x * p_w;
                const amrex::ParticleReal p_y_mean = p_y * p_w;
                const amrex::ParticleReal p_t_mean = p_t * p_w;

                const amrex::ParticleReal p_px_mean = p_px * p_w;
                const amrex::ParticleReal p_py_mean = p_py * p_w;
                const amrex::ParticleReal p_pt_mean = p_pt * p_w;

                const amrex::ParticleReal p_sx_mean = p_sx * p_w;
                const amrex::ParticleReal p_sy_mean = p_sy * p_w;
                const amrex::ParticleReal p_sz_mean = p_sz * p_w;

                return {p_w,
                        p_x_mean, p_y_mean, p_t_mean,
                        p_px_mean, p_py_mean, p_pt_mean,
                        p_sx_mean, p_sy_mean, p_sz_mean,
                        p_x, p_y, p_t, p_px, p_py, p_pt,
                        p_x, p_y, p_t, p_px, p_py, p_pt};
            },
            reduce_ops_1
        );

        std::vector<amrex::ParticleReal> values_per_rank_1st(num_red_ops_1_sum);

        /* contains in this order:
         * w, x_mean, y_mean, t_mean
         * px_mean, py_mean, pt_mean
         * sx_mean, sy_mean, sz_mean
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
        amrex::ParticleReal const mean_x  = values_per_rank_1st.at(1) /= w_sum;
        amrex::ParticleReal const mean_y  = values_per_rank_1st.at(2) /= w_sum;
        amrex::ParticleReal const mean_t  = values_per_rank_1st.at(3) /= w_sum;
        amrex::ParticleReal const mean_px = values_per_rank_1st.at(4) /= w_sum;
        amrex::ParticleReal const mean_py = values_per_rank_1st.at(5) /= w_sum;
        amrex::ParticleReal const mean_pt = values_per_rank_1st.at(6) /= w_sum;
        amrex::ParticleReal const mean_sx = values_per_rank_1st.at(7) /= w_sum;
        amrex::ParticleReal const mean_sy = values_per_rank_1st.at(8) /= w_sum;
        amrex::ParticleReal const mean_sz = values_per_rank_1st.at(9) /= w_sum;

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
        static constexpr std::size_t num_red_ops_2 = 25;
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
                // access spin data
                const amrex::ParticleReal p_sx = p.rdata(RealSoA::sx);
                const amrex::ParticleReal p_sy = p.rdata(RealSoA::sy);
                const amrex::ParticleReal p_sz = p.rdata(RealSoA::sz);
                // prepare mean square for positions
                const amrex::ParticleReal p_x_ms = (p_x-mean_x)*(p_x-mean_x)*p_w;
                const amrex::ParticleReal p_y_ms = (p_y-mean_y)*(p_y-mean_y)*p_w;
                const amrex::ParticleReal p_t_ms = (p_t-mean_t)*(p_t-mean_t)*p_w;
                // prepare mean square for momenta
                const amrex::ParticleReal p_px_ms = (p_px-mean_px)*(p_px-mean_px)*p_w;
                const amrex::ParticleReal p_py_ms = (p_py-mean_py)*(p_py-mean_py)*p_w;
                const amrex::ParticleReal p_pt_ms = (p_pt-mean_pt)*(p_pt-mean_pt)*p_w;
                // prepare mean square for spin
                const amrex::ParticleReal p_sx_ms = (p_sx-mean_sx)*(p_sx-mean_sx)*p_w;
                const amrex::ParticleReal p_sy_ms = (p_sy-mean_sy)*(p_sy-mean_sy)*p_w;
                const amrex::ParticleReal p_sz_ms = (p_sz-mean_sz)*(p_sz-mean_sz)*p_w;
                // prepare position-momentum correlations
                const amrex::ParticleReal p_xpx = (p_x-mean_x)*(p_px-mean_px)*p_w;
                const amrex::ParticleReal p_ypy = (p_y-mean_y)*(p_py-mean_py)*p_w;
                const amrex::ParticleReal p_tpt = (p_t-mean_t)*(p_pt-mean_pt)*p_w;
                // prepare correlations for dispersion (4 required)
                const amrex::ParticleReal p_xpt = (p_x-mean_x)*(p_pt-mean_pt)*p_w;
                const amrex::ParticleReal p_pxpt = (p_px-mean_px)*(p_pt-mean_pt)*p_w;
                const amrex::ParticleReal p_ypt = (p_y-mean_y)*(p_pt-mean_pt)*p_w;
                const amrex::ParticleReal p_pypt = (p_py-mean_py)*(p_pt-mean_pt)*p_w;
                // prepare additional cross-plane correlations (8 required)
                const amrex::ParticleReal p_xy = (p_x-mean_x)*(p_y-mean_y)*p_w;
                const amrex::ParticleReal p_xpy = (p_x-mean_x)*(p_py-mean_py)*p_w;
                const amrex::ParticleReal p_xt = (p_x-mean_x)*(p_t-mean_t)*p_w;
                const amrex::ParticleReal p_pxy = (p_px-mean_px)*(p_y-mean_y)*p_w;
                const amrex::ParticleReal p_pxpy = (p_px-mean_px)*(p_py-mean_py)*p_w;
                const amrex::ParticleReal p_pxt = (p_px-mean_px)*(p_t-mean_t)*p_w;
                const amrex::ParticleReal p_yt = (p_y-mean_y)*(p_t-mean_t)*p_w;
                const amrex::ParticleReal p_pyt = (p_py-mean_py)*(p_t-mean_t)*p_w;

                const amrex::ParticleReal p_charge = q_C*p_w;

                return {p_x_ms, p_y_ms, p_t_ms,
                        p_px_ms, p_py_ms, p_pt_ms,
                        p_xpx, p_ypy, p_tpt,
                        p_xpt, p_pxpt, p_ypt, p_pypt,
                        p_xy, p_xpy, p_xt, p_pxy, p_pxpy, p_pxt, p_yt, p_pyt,
                        p_charge,
                        p_sx_ms, p_sy_ms, p_sz_ms};
            },
                reduce_ops_2
        );

        std::vector<amrex::ParticleReal> values_per_rank_2nd(num_red_ops_2);

        /* contains in this order:
         * x_ms, y_ms, t_ms
         * px_ms, py_ms, pt_ms,
         * xpx, ypy, tpt,
         * p_xpt, p_pxpt, p_ypt, p_pypt,
         * p_xy, p_xpy, p_xt, p_pxy, p_pxpy, p_pxt, p_yt, p_pyt,
         * charge,
         * sx_ms, sy_ms, sz_ms
         */
        amrex::constexpr_for<0, num_red_ops_2> ([&](auto i) {
            values_per_rank_2nd[i] = amrex::get<i>(r2);
        });

        // reduced sum over mpi ranks (allreduce)
        amrex::ParallelAllReduce::Sum(
            values_per_rank_2nd.data(),
            values_per_rank_2nd.size(),
            amrex::ParallelDescriptor::Communicator()
        );

        // minimum values
        amrex::ParticleReal const min_x = values_per_rank_min.at(0);
        amrex::ParticleReal const min_y = values_per_rank_min.at(1);
        amrex::ParticleReal const min_t = values_per_rank_min.at(2);
        amrex::ParticleReal const min_px = values_per_rank_min.at(3);
        amrex::ParticleReal const min_py = values_per_rank_min.at(4);
        amrex::ParticleReal const min_pt = values_per_rank_min.at(5);
        // maximum values
        amrex::ParticleReal const max_x = values_per_rank_max.at(0);
        amrex::ParticleReal const max_y = values_per_rank_max.at(1);
        amrex::ParticleReal const max_t = values_per_rank_max.at(2);
        amrex::ParticleReal const max_px = values_per_rank_max.at(3);
        amrex::ParticleReal const max_py = values_per_rank_max.at(4);
        amrex::ParticleReal const max_pt = values_per_rank_max.at(5);
        // mean square and correlation values
        amrex::ParticleReal const x_ms   = values_per_rank_2nd.at(0) /= w_sum;
        amrex::ParticleReal const y_ms   = values_per_rank_2nd.at(1) /= w_sum;
        amrex::ParticleReal const t_ms   = values_per_rank_2nd.at(2) /= w_sum;
        amrex::ParticleReal const px_ms  = values_per_rank_2nd.at(3) /= w_sum;
        amrex::ParticleReal const py_ms  = values_per_rank_2nd.at(4) /= w_sum;
        amrex::ParticleReal const pt_ms  = values_per_rank_2nd.at(5) /= w_sum;
        amrex::ParticleReal const xpx    = values_per_rank_2nd.at(6) /= w_sum;
        amrex::ParticleReal const ypy    = values_per_rank_2nd.at(7) /= w_sum;
        amrex::ParticleReal const tpt    = values_per_rank_2nd.at(8) /= w_sum;
        amrex::ParticleReal const xpt    = values_per_rank_2nd.at(9) /= w_sum;
        amrex::ParticleReal const pxpt   = values_per_rank_2nd.at(10) /= w_sum;
        amrex::ParticleReal const ypt    = values_per_rank_2nd.at(11) /= w_sum;
        amrex::ParticleReal const pypt   = values_per_rank_2nd.at(12) /= w_sum;
        amrex::ParticleReal const xy     = values_per_rank_2nd.at(13) /= w_sum;
        amrex::ParticleReal const xpy    = values_per_rank_2nd.at(14) /= w_sum;
        amrex::ParticleReal const xt     = values_per_rank_2nd.at(15) /= w_sum;
        amrex::ParticleReal const pxy    = values_per_rank_2nd.at(16) /= w_sum;
        amrex::ParticleReal const pxpy   = values_per_rank_2nd.at(17) /= w_sum;
        amrex::ParticleReal const pxt    = values_per_rank_2nd.at(18) /= w_sum;
        amrex::ParticleReal const yt     = values_per_rank_2nd.at(19) /= w_sum;
        amrex::ParticleReal const pyt    = values_per_rank_2nd.at(20) /= w_sum;
        amrex::ParticleReal const charge = values_per_rank_2nd.at(21);
        amrex::ParticleReal const sx_ms  = values_per_rank_2nd.at(22) /= w_sum;
        amrex::ParticleReal const sy_ms  = values_per_rank_2nd.at(23) /= w_sum;
        amrex::ParticleReal const sz_ms  = values_per_rank_2nd.at(24) /= w_sum;
        // standard deviations of positions
        amrex::ParticleReal const sigma_x = std::sqrt(x_ms);
        amrex::ParticleReal const sigma_y = std::sqrt(y_ms);
        amrex::ParticleReal const sigma_t = std::sqrt(t_ms);
        // standard deviations of momenta
        amrex::ParticleReal const sigma_px = std::sqrt(px_ms);
        amrex::ParticleReal const sigma_py = std::sqrt(py_ms);
        amrex::ParticleReal const sigma_pt = std::sqrt(pt_ms);
        // standard deviations of spin
        amrex::ParticleReal const sigma_sx = std::sqrt(sx_ms);
        amrex::ParticleReal const sigma_sy = std::sqrt(sy_ms);
        amrex::ParticleReal const sigma_sz = std::sqrt(sz_ms);
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
        data["beta_gamma"] = bg;
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
        data["mean_x"] = nan;
        data["min_x"] = nan;
        data["max_x"] = nan;
        data["mean_y"] = nan;
        data["min_y"] = nan;
        data["max_y"] = nan;
        data["mean_t"] = nan;
        data["min_t"] = nan;
        data["max_t"] = nan;
        data["sigma_x"] = sig_x;
        data["sigma_y"] = sig_y;
        data["sigma_t"] = sig_t;
        data["mean_px"] = nan;
        data["min_px"] = nan;
        data["max_px"] = nan;
        data["mean_py"] = nan;
        data["min_py"] = nan;
        data["max_py"] = nan;
        data["mean_pt"] = nan;
        data["min_pt"] = nan;
        data["max_pt"] = nan;
        data["sigma_px"] = sig_px;
        data["sigma_py"] = sig_py;
        data["sigma_pt"] = sig_pt;
        // start deprecated attributes
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
        data["beta_gamma"] = bg;
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
