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

#include "BeamMomentsSelection.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/ReferenceParticle.H"
#include "particles/CovarianceMatrix.H"
#include "EmittanceInvariants.H"

#include <AMReX_Array.H>                // for GpuArray
#include <AMReX_BLProfiler.H>           // for TinyProfiler
#include <AMReX_GpuDevice.H>            // for dtoh_memcpy
#include <AMReX_GpuQualifiers.H>        // for AMREX_GPU_DEVICE
#include <AMReX_ParallelDescriptor.H>   // for ParallelDescriptor
#include <AMReX_ParticleReduce.H>       // for ParticleReduce
#include <AMReX_REAL.H>                 // for ParticleReal
#include <AMReX_Reduce.H>               // for ReduceOps
#include <AMReX_SmallMatrix.H>          // for SmallMatrix
#include <AMReX_TypeList.H>             // for TypeMultiplier

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>


namespace impactx::diagnostics
{
namespace
{
    /** Central (weighted) beam moments, first moments, per-coordinate extremes
     *  and beam charge, gathered in one struct so that the derived-quantity
     *  recovery (sigmas, emittances, Twiss, dispersion, eigenemittances) and the
     *  output-map assembly are shared between the particle-based and the
     *  covariance-matrix-based overloads below (single source of truth for the
     *  parallel-axis recovery math).
     */
    struct RawMoments
    {
        // second central moments (phase space)
        amrex::ParticleReal x_ms, y_ms, t_ms, px_ms, py_ms, pt_ms;
        amrex::ParticleReal xpx, ypy, tpt;
        amrex::ParticleReal xpt, pxpt, ypt, pypt;
        amrex::ParticleReal xy, xpy, xt, pxy, pxpy, pxt, yt, pyt;
        // spin second central moments
        amrex::ParticleReal sx_ms, sy_ms, sz_ms;
        // first moments (means)
        amrex::ParticleReal mean_x, mean_y, mean_t, mean_px, mean_py, mean_pt;
        amrex::ParticleReal mean_sx, mean_sy, mean_sz;
        // per-coordinate extremes
        amrex::ParticleReal min_x, min_y, min_t, min_px, min_py, min_pt;
        amrex::ParticleReal max_x, max_y, max_t, max_px, max_py, max_pt;
        // beam charge [C]
        amrex::ParticleReal charge;
    };

    //! minimal requirement to compute an output key: the smallest profile plus
    //! whether the spin / min-max / eigenemittance reductions are also needed
    struct KeyReq { MomentsProfile profile; bool spin; bool minmax; bool eigen; };

    /** Canonicalize a (possibly deprecated) output name. Deprecated aliases
     *  (e.g. x_mean, sig_x, pt_min) map to their canonical spelling; any other
     *  name is returned unchanged.
     */
    std::string canonicalize (std::string const & key)
    {
        static const std::unordered_map<std::string, std::string> aliases = {
            {"x_mean", "mean_x"}, {"y_mean", "mean_y"}, {"t_mean", "mean_t"},
            {"px_mean", "mean_px"}, {"py_mean", "mean_py"}, {"pt_mean", "mean_pt"},
            {"x_min", "min_x"}, {"y_min", "min_y"}, {"t_min", "min_t"},
            {"px_min", "min_px"}, {"py_min", "min_py"}, {"pt_min", "min_pt"},
            {"x_max", "max_x"}, {"y_max", "max_y"}, {"t_max", "max_t"},
            {"px_max", "max_px"}, {"py_max", "max_py"}, {"pt_max", "max_pt"},
            {"sig_x", "sigma_x"}, {"sig_y", "sigma_y"}, {"sig_t", "sigma_t"},
            {"sig_px", "sigma_px"}, {"sig_py", "sigma_py"}, {"sig_pt", "sigma_pt"}
        };
        auto const it = aliases.find(key);
        return (it == aliases.end()) ? key : it->second;
    }

    /** Requirement to compute a canonical output key. Throws std::runtime_error
     *  for an unknown key (already canonicalized).
     */
    KeyReq key_requirement (std::string const & key)
    {
        using P = MomentsProfile;
        static const std::unordered_map<std::string, KeyReq> table = {
            // Positions: means/sigmas of x, y, t and the beam charge
            {"mean_x",   {P::Positions, false, false, false}},
            {"mean_y",   {P::Positions, false, false, false}},
            {"mean_t",   {P::Positions, false, false, false}},
            {"sigma_x",  {P::Positions, false, false, false}},
            {"sigma_y",  {P::Positions, false, false, false}},
            {"sigma_t",  {P::Positions, false, false, false}},
            {"charge_C", {P::Positions, false, false, false}},
            // Sizes: means/sigmas of px, py, pt
            {"mean_px",  {P::Sizes, false, false, false}},
            {"mean_py",  {P::Sizes, false, false, false}},
            {"mean_pt",  {P::Sizes, false, false, false}},
            {"sigma_px", {P::Sizes, false, false, false}},
            {"sigma_py", {P::Sizes, false, false, false}},
            {"sigma_pt", {P::Sizes, false, false, false}},
            // per-coordinate min/max (min/max of momenta needs the momenta loaded)
            {"min_x",  {P::Positions, false, true, false}},
            {"max_x",  {P::Positions, false, true, false}},
            {"min_y",  {P::Positions, false, true, false}},
            {"max_y",  {P::Positions, false, true, false}},
            {"min_t",  {P::Positions, false, true, false}},
            {"max_t",  {P::Positions, false, true, false}},
            {"min_px", {P::Sizes, false, true, false}},
            {"max_px", {P::Sizes, false, true, false}},
            {"min_py", {P::Sizes, false, true, false}},
            {"max_py", {P::Sizes, false, true, false}},
            {"min_pt", {P::Sizes, false, true, false}},
            {"max_pt", {P::Sizes, false, true, false}},
            // Twiss: emittances, dispersion and (dispersion-corrected) Twiss
            {"emittance_x",   {P::Twiss, false, false, false}},
            {"emittance_y",   {P::Twiss, false, false, false}},
            {"emittance_t",   {P::Twiss, false, false, false}},
            {"emittance_xn",  {P::Twiss, false, false, false}},
            {"emittance_yn",  {P::Twiss, false, false, false}},
            {"emittance_tn",  {P::Twiss, false, false, false}},
            {"alpha_x",       {P::Twiss, false, false, false}},
            {"alpha_y",       {P::Twiss, false, false, false}},
            {"alpha_t",       {P::Twiss, false, false, false}},
            {"beta_x",        {P::Twiss, false, false, false}},
            {"beta_y",        {P::Twiss, false, false, false}},
            {"beta_t",        {P::Twiss, false, false, false}},
            {"dispersion_x",  {P::Twiss, false, false, false}},
            {"dispersion_px", {P::Twiss, false, false, false}},
            {"dispersion_y",  {P::Twiss, false, false, false}},
            {"dispersion_py", {P::Twiss, false, false, false}},
            // eigenemittances need the full 6x6 covariance
            {"emittance_1", {P::Full, false, false, true}},
            {"emittance_2", {P::Full, false, false, true}},
            {"emittance_3", {P::Full, false, false, true}},
            // spin first and second moments
            {"mean_sx",  {P::Positions, true, false, false}},
            {"mean_sy",  {P::Positions, true, false, false}},
            {"mean_sz",  {P::Positions, true, false, false}},
            {"sigma_sx", {P::Positions, true, false, false}},
            {"sigma_sy", {P::Positions, true, false, false}},
            {"sigma_sz", {P::Positions, true, false, false}}
        };
        auto const it = table.find(key);
        if (it == table.end())
        {
            throw std::runtime_error(
                "reduced beam characteristics: unknown moment name '" + key + "'");
        }
        return it->second;
    }

    /** Whether a key's requirement is satisfied by a selection. */
    bool key_covered (std::string const & key, MomentsSelection const & sel)
    {
        KeyReq const req = key_requirement(canonicalize(key));
        return static_cast<int>(req.profile) <= static_cast<int>(sel.profile)
            && (!req.spin   || sel.spin)
            && (!req.minmax || sel.minmax)
            && (!req.eigen  || sel.eigen);
    }

    /** Recover the derived beam characteristics from the (central) moments and
     *  assemble the output map. The arithmetic below is unchanged from the
     *  original inlined recovery; it now lives in exactly one place.
     */
    std::unordered_map<std::string, amrex::ParticleReal>
    derive_and_assemble (RawMoments const & m,
                         amrex::ParticleReal const bg,
                         amrex::ParticleReal const bg2,
                         MomentsSelection const & sel)
    {
        using namespace amrex::literals; // for _prt

        // unpack the (central) moments recovered by the caller
        amrex::ParticleReal const x_ms   = m.x_ms;
        amrex::ParticleReal const y_ms   = m.y_ms;
        amrex::ParticleReal const t_ms   = m.t_ms;
        amrex::ParticleReal const px_ms  = m.px_ms;
        amrex::ParticleReal const py_ms  = m.py_ms;
        amrex::ParticleReal const pt_ms  = m.pt_ms;
        amrex::ParticleReal const xpx    = m.xpx;
        amrex::ParticleReal const ypy    = m.ypy;
        amrex::ParticleReal const tpt    = m.tpt;
        amrex::ParticleReal const xpt    = m.xpt;
        amrex::ParticleReal const pxpt   = m.pxpt;
        amrex::ParticleReal const ypt    = m.ypt;
        amrex::ParticleReal const pypt   = m.pypt;
        amrex::ParticleReal const xy     = m.xy;
        amrex::ParticleReal const xpy    = m.xpy;
        amrex::ParticleReal const xt     = m.xt;
        amrex::ParticleReal const pxy    = m.pxy;
        amrex::ParticleReal const pxpy   = m.pxpy;
        amrex::ParticleReal const pxt    = m.pxt;
        amrex::ParticleReal const yt     = m.yt;
        amrex::ParticleReal const pyt    = m.pyt;
        amrex::ParticleReal const sx_ms  = m.sx_ms;
        amrex::ParticleReal const sy_ms  = m.sy_ms;
        amrex::ParticleReal const sz_ms  = m.sz_ms;
        amrex::ParticleReal const mean_x  = m.mean_x;
        amrex::ParticleReal const mean_y  = m.mean_y;
        amrex::ParticleReal const mean_t  = m.mean_t;
        amrex::ParticleReal const mean_px = m.mean_px;
        amrex::ParticleReal const mean_py = m.mean_py;
        amrex::ParticleReal const mean_pt = m.mean_pt;
        amrex::ParticleReal const mean_sx = m.mean_sx;
        amrex::ParticleReal const mean_sy = m.mean_sy;
        amrex::ParticleReal const mean_sz = m.mean_sz;
        amrex::ParticleReal const min_x  = m.min_x;
        amrex::ParticleReal const min_y  = m.min_y;
        amrex::ParticleReal const min_t  = m.min_t;
        amrex::ParticleReal const min_px = m.min_px;
        amrex::ParticleReal const min_py = m.min_py;
        amrex::ParticleReal const min_pt = m.min_pt;
        amrex::ParticleReal const max_x  = m.max_x;
        amrex::ParticleReal const max_y  = m.max_y;
        amrex::ParticleReal const max_t  = m.max_t;
        amrex::ParticleReal const max_px = m.max_px;
        amrex::ParticleReal const max_py = m.max_py;
        amrex::ParticleReal const max_pt = m.max_pt;
        amrex::ParticleReal const charge = m.charge;

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
        amrex::ParticleReal const emittance_x = (e2_x > 0.0)? std::sqrt(e2_x) : 0.0_prt;
        amrex::ParticleReal const emittance_y = (e2_y > 0.0)? std::sqrt(e2_y) : 0.0_prt;
        amrex::ParticleReal const emittance_t = (e2_t > 0.0)? std::sqrt(e2_t) : 0.0_prt;
        // Dispersion and dispersive beam moments
        amrex::ParticleReal const dispersion_x = ((pt_ms > 0.0) ? (- xpt / pt_ms) : 0.0_prt);
        amrex::ParticleReal const dispersion_px = ((pt_ms > 0.0) ? (- pxpt / pt_ms) : 0.0_prt);
        amrex::ParticleReal const dispersion_y = ((pt_ms > 0.0) ? (- ypt / pt_ms) : 0.0_prt);
        amrex::ParticleReal const dispersion_py = ((pt_ms > 0.0) ? (- pypt / pt_ms) : 0.0_prt);
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
        bool const compute_eigenemittances = sel.eigen;
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

        // Filter to the requested selection: a non-empty key list keeps exactly
        // those keys (honoring the requested spelling); an empty list keeps every
        // key the profile and flags cover (used by the "all" and default
        // selections).
        if (sel.keys.empty())
        {
            for (auto it = data.begin(); it != data.end(); )
            {
                if (key_covered(it->first, sel)) { ++it; }
                else { it = data.erase(it); }
            }
        }
        else
        {
            std::unordered_map<std::string, amrex::ParticleReal> filtered;
            for (auto const & k : sel.keys)
            {
                auto const it = data.find(k);
                if (it != data.end()) { filtered.emplace(k, it->second); }
            }
            data = std::move(filtered);
        }

        return data;
    }

    // -----------------------------------------------------------------------
    // Selective single-pass reduction. One templated kernel body generates the
    // reduction for every profile; a profile is the ordered list of weighted
    // power sums (and the coordinates to min/max) needed to recover a requested
    // subset of moments. Lighter profiles read fewer SoA arrays and reduce fewer
    // values. The full profile reproduces the previous reduction bit-for-bit.
    // -----------------------------------------------------------------------

    //! phase-space and spin coordinate identifiers; `none` is a
    //! multiplicative-identity placeholder used by the generic slot formula
    enum class Coord : int { x = 0, y, t, px, py, pt, sx, sy, sz, none };

    // (MomentsProfile, the reduction-profile enum, is declared in
    //  BeamMomentsSelection.H so it can be shared with MomentsSelection.)

    //! one reduced (weighted) sum slot: sum over particles of w * dev(a) * dev(b),
    //! with dev(u) = u - shift_u and dev(none) = 1. Hence {none,none} = Sum(w),
    //! {a,none} = first moment of a, {a,b} = second moment of a and b.
    struct SumSpec { Coord a; Coord b; };

    //! Full has 1 + 6 first + 6 diagonal + 3 same-plane + 4 dispersive + 8
    //! cross-plane = 28 phase-space sums; the spin toggle adds 3 + 3 = 6.
    static constexpr int max_sums = 34;
    static constexpr int max_minmax = 6;

    struct ProfileDesc
    {
        SumSpec sums[max_sums] {};
        int n_sums = 0;
        Coord minmax[max_minmax] {};
        int n_minmax = 0;
    };

    /** Compile-time description of a reduction profile: the ordered list of
     *  weighted-sum slots and the coordinates to reduce min/max over.
     */
    constexpr ProfileDesc make_desc (MomentsProfile const profile, bool const spin, bool const minmax)
    {
        ProfileDesc d {};
        auto add_sum = [&d] (Coord const a, Coord const b)
        {
            d.sums[d.n_sums] = SumSpec{a, b};
            ++d.n_sums;
        };

        // Sum(w)
        add_sum(Coord::none, Coord::none);

        // positions: first moments, then diagonal second moments
        Coord const pos[3] = {Coord::x, Coord::y, Coord::t};
        for (Coord const c : pos) { add_sum(c, Coord::none); }
        for (Coord const c : pos) { add_sum(c, c); }

        // momenta (Sizes and up)
        if (profile >= MomentsProfile::Sizes)
        {
            Coord const mom[3] = {Coord::px, Coord::py, Coord::pt};
            for (Coord const c : mom) { add_sum(c, Coord::none); }
            for (Coord const c : mom) { add_sum(c, c); }
        }

        // correlations for emittance, dispersion and (dispersion-corrected) Twiss
        if (profile >= MomentsProfile::Twiss)
        {
            // same-plane
            add_sum(Coord::x, Coord::px);
            add_sum(Coord::y, Coord::py);
            add_sum(Coord::t, Coord::pt);
            // dispersive
            add_sum(Coord::x, Coord::pt);
            add_sum(Coord::px, Coord::pt);
            add_sum(Coord::y, Coord::pt);
            add_sum(Coord::py, Coord::pt);
        }

        // cross-plane correlations (only eigenemittances consume these)
        if (profile >= MomentsProfile::Full)
        {
            add_sum(Coord::x, Coord::y);
            add_sum(Coord::x, Coord::py);
            add_sum(Coord::x, Coord::t);
            add_sum(Coord::px, Coord::y);
            add_sum(Coord::px, Coord::py);
            add_sum(Coord::px, Coord::t);
            add_sum(Coord::y, Coord::t);
            add_sum(Coord::py, Coord::t);
        }

        // spin first and diagonal second moments
        if (spin)
        {
            add_sum(Coord::sx, Coord::none);
            add_sum(Coord::sy, Coord::none);
            add_sum(Coord::sz, Coord::none);
            add_sum(Coord::sx, Coord::sx);
            add_sum(Coord::sy, Coord::sy);
            add_sum(Coord::sz, Coord::sz);
        }

        // min/max over the coordinates the profile already loads
        if (minmax)
        {
            d.minmax[d.n_minmax++] = Coord::x;
            d.minmax[d.n_minmax++] = Coord::y;
            d.minmax[d.n_minmax++] = Coord::t;
            if (profile >= MomentsProfile::Sizes)
            {
                d.minmax[d.n_minmax++] = Coord::px;
                d.minmax[d.n_minmax++] = Coord::py;
                d.minmax[d.n_minmax++] = Coord::pt;
            }
        }

        return d;
    }

    /** SoA component index for a coordinate (compile-time). */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr int soa_index (Coord const c)
    {
        switch (c)
        {
            case Coord::x:  return RealSoA::x;
            case Coord::y:  return RealSoA::y;
            case Coord::t:  return RealSoA::t;
            case Coord::px: return RealSoA::px;
            case Coord::py: return RealSoA::py;
            case Coord::pt: return RealSoA::pt;
            case Coord::sx: return RealSoA::sx;
            case Coord::sy: return RealSoA::sy;
            case Coord::sz: return RealSoA::sz;
            case Coord::none: return RealSoA::x;  // never dereferenced (guarded by deviation)
        }
        return RealSoA::x;
    }

    /** Shifted deviation dev(C) = (C - shift_C), or the exact constant 1 for
     *  Coord::none so that a first-moment slot {a,none} is dev(a)*1*w and Sum(w)
     *  is 1*1*w. Multiplying by 1 leaves the value bit-identical.
     */
    template <Coord C, typename PType>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::ParticleReal deviation (
        PType const & p,
        amrex::GpuArray<amrex::ParticleReal, 9> const & shift)
    {
        if constexpr (C == Coord::none)
        {
            return amrex::ParticleReal(1);
        }
        else
        {
            return p.rdata(soa_index(C)) - shift[static_cast<int>(C)];
        }
    }

    /** Fill the weighted power-sum slots (the first desc.n_sums tuple entries)
     *  for a profile. Shared by the with- and without-min/max reduction paths so
     *  the generic slot formula w*dev(a)*dev(b) lives in exactly one place.
     */
    template <MomentsProfile P, bool Spin, typename TupleT, typename PType>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void fill_sum_slots (
        TupleT & out,
        PType const & p,
        amrex::GpuArray<amrex::ParticleReal, 9> const & shift)
    {
        static constexpr ProfileDesc desc = make_desc(P, Spin, false);
        amrex::ParticleReal const p_w = p.rdata(RealSoA::w);
        amrex::constexpr_for<0, desc.n_sums>([&] (auto i)
        {
            constexpr SumSpec s = desc.sums[i];
            amrex::get<i>(out) =
                deviation<s.a>(p, shift) * deviation<s.b>(p, shift) * p_w;
        });
    }

    /** Reduce a single pass over the particles for the given profile and recover
     *  the central moments into a RawMoments. The full profile
     *  (Profile::Full, spin, minmax) reproduces the previous hand-written
     *  reduction bit-for-bit: every slot is w*dev(a)*dev(b) in the same
     *  multiplication order, and the parallel-axis recovery keeps the same order
     *  of operations, keyed off the SumSpec instead of hard-coded indices.
     */
    template <MomentsProfile P, bool Spin, bool MinMax>
    RawMoments reduce_and_recover (
        ImpactXParticleContainer const & pc,
        std::array<amrex::ParticleReal, 9> const & shift,
        amrex::ParticleReal const q_C)
    {
        using namespace amrex::literals; // for _prt
        using PType = typename ImpactXParticleContainer::SuperParticleType;

        static constexpr ProfileDesc desc = make_desc(P, Spin, MinMax);
        static constexpr int n_sum = desc.n_sums;
        static constexpr int n_mm  = desc.n_minmax;

        amrex::GpuArray<amrex::ParticleReal, 9> shift_d {};
        for (int i = 0; i < 9; ++i) { shift_d[i] = shift[i]; }

        std::array<amrex::ParticleReal, n_sum> values_sum {};
        [[maybe_unused]] std::array<amrex::ParticleReal, (MinMax ? max_minmax : 1)> values_min {};
        [[maybe_unused]] std::array<amrex::ParticleReal, (MinMax ? max_minmax : 1)> values_max {};

        // A zero-count ReduceOpMin/Max would decay to a bogus pointer op, so the
        // min/max operations are only ever named on the MinMax path.
        if constexpr (MinMax)
        {
            amrex::TypeMultiplier<amrex::ReduceOps,
                amrex::ReduceOpSum[n_sum],
                amrex::ReduceOpMin[n_mm],
                amrex::ReduceOpMax[n_mm]
            > reduce_ops;
            using ReducedDataT = amrex::TypeMultiplier<amrex::ReduceData,
                amrex::ParticleReal[n_sum + 2 * n_mm]>;

            auto r = amrex::ParticleReduce<ReducedDataT>(
                pc,
                [=] AMREX_GPU_DEVICE (PType const & p) noexcept -> typename ReducedDataT::Type
                {
                    typename ReducedDataT::Type out;
                    fill_sum_slots<P, Spin>(out, p, shift_d);
                    amrex::constexpr_for<0, n_mm>([&] (auto j)
                    {
                        constexpr Coord c = desc.minmax[j];
                        amrex::ParticleReal const v = p.rdata(soa_index(c));
                        amrex::get<n_sum + j>(out) = v;
                        amrex::get<n_sum + n_mm + j>(out) = v;
                    });
                    return out;
                },
                reduce_ops
            );

            amrex::constexpr_for<0, n_sum>([&] (auto i) { values_sum[i] = amrex::get<i>(r); });
            amrex::constexpr_for<0, n_mm>([&] (auto j)
            {
                values_min[j] = amrex::get<n_sum + j>(r);
                values_max[j] = amrex::get<n_sum + n_mm + j>(r);
            });
        }
        else
        {
            amrex::TypeMultiplier<amrex::ReduceOps,
                amrex::ReduceOpSum[n_sum]
            > reduce_ops;
            using ReducedDataT = amrex::TypeMultiplier<amrex::ReduceData,
                amrex::ParticleReal[n_sum]>;

            auto r = amrex::ParticleReduce<ReducedDataT>(
                pc,
                [=] AMREX_GPU_DEVICE (PType const & p) noexcept -> typename ReducedDataT::Type
                {
                    typename ReducedDataT::Type out;
                    fill_sum_slots<P, Spin>(out, p, shift_d);
                    return out;
                },
                reduce_ops
            );

            amrex::constexpr_for<0, n_sum>([&] (auto i) { values_sum[i] = amrex::get<i>(r); });
        }

        // reduce across MPI ranks (allreduce)
        amrex::ParallelAllReduce::Sum(
            values_sum.data(), n_sum, amrex::ParallelDescriptor::Communicator());
        if constexpr (MinMax)
        {
            amrex::ParallelAllReduce::Min(
                values_min.data(), n_mm, amrex::ParallelDescriptor::Communicator());
            amrex::ParallelAllReduce::Max(
                values_max.data(), n_mm, amrex::ParallelDescriptor::Communicator());
        }

        // Recover central moments via the parallel-axis theorem. Each slot is an
        // independent sum, so slot order does not affect the results; the
        // recovery is keyed off the SumSpec, matching the previous arithmetic.
        amrex::ParticleReal w_sum = 0.0_prt;
        amrex::constexpr_for<0, n_sum>([&] (auto i)
        {
            constexpr SumSpec s = desc.sums[i];
            if constexpr (s.a == Coord::none && s.b == Coord::none) { w_sum = values_sum[i]; }
        });
        amrex::ParticleReal dmean[9] = {};
        amrex::constexpr_for<0, n_sum>([&] (auto i)
        {
            constexpr SumSpec s = desc.sums[i];
            if constexpr (s.a != Coord::none && s.b == Coord::none)
            {
                dmean[static_cast<int>(s.a)] = values_sum[i] / w_sum;
            }
        });
        amrex::ParticleReal cov[9][9] = {};
        amrex::constexpr_for<0, n_sum>([&] (auto i)
        {
            constexpr SumSpec s = desc.sums[i];
            if constexpr (s.a != Coord::none && s.b != Coord::none)
            {
                cov[static_cast<int>(s.a)][static_cast<int>(s.b)] =
                    values_sum[i] / w_sum
                    - dmean[static_cast<int>(s.a)] * dmean[static_cast<int>(s.b)];
            }
        });

        auto const nan = std::numeric_limits<amrex::ParticleReal>::quiet_NaN();
        amrex::ParticleReal cmin[6] = {nan, nan, nan, nan, nan, nan};
        amrex::ParticleReal cmax[6] = {nan, nan, nan, nan, nan, nan};
        if constexpr (MinMax)
        {
            amrex::constexpr_for<0, n_mm>([&] (auto j)
            {
                constexpr int c = static_cast<int>(desc.minmax[j]);
                cmin[c] = values_min[j];
                cmax[c] = values_max[j];
            });
        }

        // map into the shared RawMoments layout
        constexpr int X  = static_cast<int>(Coord::x);
        constexpr int Y  = static_cast<int>(Coord::y);
        constexpr int T  = static_cast<int>(Coord::t);
        constexpr int PX = static_cast<int>(Coord::px);
        constexpr int PY = static_cast<int>(Coord::py);
        constexpr int PT = static_cast<int>(Coord::pt);
        constexpr int SX = static_cast<int>(Coord::sx);
        constexpr int SY = static_cast<int>(Coord::sy);
        constexpr int SZ = static_cast<int>(Coord::sz);

        RawMoments m {};
        m.x_ms  = cov[X][X];   m.y_ms  = cov[Y][Y];   m.t_ms  = cov[T][T];
        m.px_ms = cov[PX][PX]; m.py_ms = cov[PY][PY]; m.pt_ms = cov[PT][PT];
        m.xpx = cov[X][PX]; m.ypy = cov[Y][PY]; m.tpt = cov[T][PT];
        m.xpt = cov[X][PT]; m.pxpt = cov[PX][PT]; m.ypt = cov[Y][PT]; m.pypt = cov[PY][PT];
        m.xy = cov[X][Y]; m.xpy = cov[X][PY]; m.xt = cov[X][T];
        m.pxy = cov[PX][Y]; m.pxpy = cov[PX][PY]; m.pxt = cov[PX][T];
        m.yt = cov[Y][T]; m.pyt = cov[PY][T];
        m.sx_ms = cov[SX][SX]; m.sy_ms = cov[SY][SY]; m.sz_ms = cov[SZ][SZ];
        m.mean_x  = shift[X]  + dmean[X];
        m.mean_y  = shift[Y]  + dmean[Y];
        m.mean_t  = shift[T]  + dmean[T];
        m.mean_px = shift[PX] + dmean[PX];
        m.mean_py = shift[PY] + dmean[PY];
        m.mean_pt = shift[PT] + dmean[PT];
        m.mean_sx = shift[SX] + dmean[SX];
        m.mean_sy = shift[SY] + dmean[SY];
        m.mean_sz = shift[SZ] + dmean[SZ];
        m.min_x = cmin[X]; m.min_y = cmin[Y]; m.min_t = cmin[T];
        m.min_px = cmin[PX]; m.min_py = cmin[PY]; m.min_pt = cmin[PT];
        m.max_x = cmax[X]; m.max_y = cmax[Y]; m.max_t = cmax[T];
        m.max_px = cmax[PX]; m.max_py = cmax[PY]; m.max_pt = cmax[PT];
        m.charge = q_C * w_sum;

        return m;
    }

    /** Dispatch the spin / min-max flags at a fixed compile-time profile. */
    template <MomentsProfile P>
    RawMoments dispatch_flags (
        ImpactXParticleContainer const & pc,
        std::array<amrex::ParticleReal, 9> const & shift,
        amrex::ParticleReal const q_C,
        MomentsSelection const & sel)
    {
        if (sel.spin && sel.minmax)  { return reduce_and_recover<P, true,  true >(pc, shift, q_C); }
        if (sel.spin && !sel.minmax) { return reduce_and_recover<P, true,  false>(pc, shift, q_C); }
        if (!sel.spin && sel.minmax) { return reduce_and_recover<P, false, true >(pc, shift, q_C); }
        return reduce_and_recover<P, false, false>(pc, shift, q_C);
    }

    /** Run the single-pass reduction for a selection, choosing the compile-time
     *  profile at runtime. All 4 x 2 x 2 instantiations are generated from the
     *  one reduce_and_recover kernel body.
     */
    RawMoments dispatch_reduce (
        ImpactXParticleContainer const & pc,
        std::array<amrex::ParticleReal, 9> const & shift,
        amrex::ParticleReal const q_C,
        MomentsSelection const & sel)
    {
        switch (sel.profile)
        {
            case MomentsProfile::Positions: return dispatch_flags<MomentsProfile::Positions>(pc, shift, q_C, sel);
            case MomentsProfile::Sizes:     return dispatch_flags<MomentsProfile::Sizes>(pc, shift, q_C, sel);
            case MomentsProfile::Twiss:     return dispatch_flags<MomentsProfile::Twiss>(pc, shift, q_C, sel);
            case MomentsProfile::Full:      return dispatch_flags<MomentsProfile::Full>(pc, shift, q_C, sel);
        }
        return dispatch_flags<MomentsProfile::Full>(pc, shift, q_C, sel);  // unreachable
    }
} // namespace

    MomentsSelection
    all_beam_moments_selection (bool const eigen)
    {
        MomentsSelection sel;
        sel.profile = MomentsProfile::Full;
        sel.spin = true;
        sel.minmax = true;
        sel.eigen = eigen;
        sel.keys.clear();  // empty => emit every covered key (i.e. all of them)
        return sel;
    }

    MomentsSelection
    default_beam_moments_selection (bool const spin_on, bool const eigen)
    {
        MomentsSelection sel;
        sel.profile = eigen ? MomentsProfile::Full : MomentsProfile::Twiss;
        sel.spin = spin_on;
        sel.minmax = false;
        sel.eigen = eigen;
        // empty => emit every covered key: all outputs except min/max, and except
        // the spin moments when spin tracking is off
        sel.keys.clear();
        return sel;
    }

    MomentsSelection
    resolve_beam_moments_selection (std::vector<std::string> const & names, bool const eigen_default)
    {
        // the "all" token requests the full set
        if (names.size() == 1 && names[0] == "all")
        {
            return all_beam_moments_selection(eigen_default);
        }

        MomentsSelection sel;
        sel.profile = MomentsProfile::Positions;
        sel.spin = false;
        sel.minmax = false;
        sel.eigen = false;
        sel.keys.reserve(names.size());
        for (auto const & name : names)
        {
            KeyReq const req = key_requirement(canonicalize(name));  // throws if unknown
            if (static_cast<int>(req.profile) > static_cast<int>(sel.profile))
            {
                sel.profile = req.profile;
            }
            sel.spin   = sel.spin   || req.spin;
            sel.minmax = sel.minmax || req.minmax;
            sel.eigen  = sel.eigen  || req.eigen;
            sel.keys.push_back(name);  // keep the requested spelling for the output
        }
        if (sel.eigen && static_cast<int>(sel.profile) < static_cast<int>(MomentsProfile::Full))
        {
            sel.profile = MomentsProfile::Full;
        }
        return sel;
    }

    std::unordered_map<std::string, amrex::ParticleReal>
    reduced_beam_characteristics (ImpactXParticleContainer const & pc)
    {
        BL_PROFILE("impactx::diagnostics::reduced_beam_characteristics(pc)");

        // full set, for backward-compatible callers (space charge, ASCII output,
        // deprecated Python binding); eigenemittances follow the diag flag
        amrex::ParmParse pp_diag("diag");
        bool eigen = false;
        pp_diag.queryAdd("eigenemittances", eigen);
        return reduced_beam_characteristics(pc, all_beam_moments_selection(eigen));
    }

    std::unordered_map<std::string, amrex::ParticleReal>
    reduced_beam_characteristics (
        ImpactXParticleContainer const & pc,
        MomentsSelection const & selection)
    {
        BL_PROFILE("impactx::diagnostics::reduced_beam_characteristics(pc, selection)");

        using namespace amrex::literals; // for _prt

        // preparing to access reference particle data: RefPart
        RefPart const ref_part = pc.GetRefParticle();
        // reference particle charge in C
        amrex::ParticleReal const q_C = ref_part.charge;
        // reference particle relativistic beta*gamma
        amrex::ParticleReal const bg = ref_part.beta_gamma();
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
         * Any global constant is exact under the parallel-axis theorem. Sampling the
         * first particle needs only offset/rms (i.e. to land within ~rms of the
         * centroid), not the true mean, and makes no assumption on the beam shape.
         * (In the offset >> rms limit the leading digits of the subtraction cancel
         * exactly per Sterbenz's lemma, but that is an asymptotic bonus, not the
         * mechanism.)
         */
        constexpr int comps[9] = {
            RealSoA::x, RealSoA::y, RealSoA::t,
            RealSoA::px, RealSoA::py, RealSoA::pt,
            RealSoA::sx, RealSoA::sy, RealSoA::sz
        };
        std::array shift = {0._prt, 0._prt, 0._prt, 0._prt, 0._prt, 0._prt, 0._prt, 0._prt, 0._prt};
        {
            // Sample the first particle held locally, then agree on a single global
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
                        for (int c = 0; c < 9; ++c) {
                            amrex::Gpu::dtoh_memcpy(
                                &shift[c], soa[comps[c]].dataPtr(), sizeof(amrex::ParticleReal));
                        }
                        found = true;
                        break;
                    }
                }
            }
            int src_rank = found ? amrex::ParallelDescriptor::MyProc()
                                 : amrex::ParallelDescriptor::NProcs();
            amrex::ParallelAllReduce::Min(src_rank, amrex::ParallelDescriptor::Communicator());
            if (src_rank < amrex::ParallelDescriptor::NProcs()) {
                amrex::ParallelDescriptor::Bcast(shift.data(), shift.size(), src_rank);
            }
        }
        // Single-pass reduction using the fastest profile that covers the
        // requested selection, then recover and assemble the requested outputs.
        RawMoments const raw = dispatch_reduce(pc, shift, q_C, selection);
        return derive_and_assemble(raw, bg, bg2, selection);
    }

    std::unordered_map<std::string, amrex::ParticleReal>
    reduced_beam_characteristics (Map6x6 const & cm, RefPart const & ref_part)
    {
        BL_PROFILE("impactx::diagnostics::reduced_beam_characteristics(cm)");

        using namespace amrex::literals; // for _prt

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
        auto const nan = std::numeric_limits<amrex::ParticleReal>::quiet_NaN();

        // A covariance matrix carries only the (central) second moments. The
        // means are zero by construction; per-coordinate extremes, spin moments
        // and beam charge are unavailable and reported as NaN (as before).
        RawMoments const raw {
            x_ms, y_ms, t_ms, px_ms, py_ms, pt_ms,
            xpx, ypy, tpt,
            xpt, pxpt, ypt, pypt,
            xy, xpy, xt, pxy, pxpy, pxt, yt, pyt,
            nan, nan, nan,
            0.0_prt, 0.0_prt, 0.0_prt, 0.0_prt, 0.0_prt, 0.0_prt,
            nan, nan, nan,
            nan, nan, nan, nan, nan, nan,
            nan, nan, nan, nan, nan, nan,
            nan  // charge_C (TODO: with space charge)
        };
        amrex::ParmParse pp_diag("diag");
        bool eigen = false;
        pp_diag.queryAdd("eigenemittances", eigen);
        return derive_and_assemble(raw, bg, bg2, all_beam_moments_selection(eigen));
    }

} // namespace impactx::diagnostics
