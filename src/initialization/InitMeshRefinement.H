/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#pragma once

#include "initialization/Algorithms.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>


namespace impactx::initialization
{
    inline amrex::Vector<amrex::Real>
    read_mr_prob_relative ()
    {
        amrex::ParmParse pp_algo("algo");
        amrex::ParmParse pp_amr("amr");
        amrex::ParmParse pp_geometry("geometry");

        auto space_charge = get_space_charge_algo();

        std::string poisson_solver = "fft";
        pp_algo.queryAdd("poisson_solver", poisson_solver);

        int max_level = 0;
        pp_amr.queryWithParser("max_level", max_level);

        if (max_level > 1 && space_charge != SpaceChargeAlgo::False)
            throw std::runtime_error(
                "Mesh-refinement (amr.max_level>=0) is only supported with "
                "space charge modeling (algo.space_charge=1).");

        // The box is expanded beyond the min and max of the particle beam.
        amrex::Vector<amrex::Real> prob_relative(max_level + 1, 1.0);
        prob_relative[0] = 3.0;  // top/bottom pad the beam on the lowest level by default by its width
        pp_geometry.queryarr("prob_relative", prob_relative);

        if (prob_relative[0] < 3.0 && space_charge != SpaceChargeAlgo::False && poisson_solver == "multigrid")
            ablastr::warn_manager::WMRecordWarning(
                    "ImpactX::read_mr_prob_relative",
                    "Dynamic resizing of the mesh uses a geometry.prob_relative "
                    "padding of less than 3 for level 0. This might result in boundary "
                    "artifacts for space charge calculation. "
                    "There is no minimum good value for this parameter, consider "
                    "doing a convergence test.",
                    ablastr::warn_manager::WarnPriority::high
            );

        if (prob_relative[0] < 1.0)
            throw std::runtime_error("geometry.prob_relative must be >= 1.0 (the beam size) on the coarsest level");

        // check that prob_relative[0] > prob_relative[1] > prob_relative[2] ...
        amrex::Real last_lev_rel = std::numeric_limits<amrex::Real>::max();
        for (int lev = 0; lev <= max_level; ++lev) {
            amrex::Real const prob_relative_lvl = prob_relative[lev];
            if (prob_relative_lvl <= 0.0)
                throw std::runtime_error("geometry.prob_relative must be strictly positive for all levels");
            if (prob_relative_lvl > last_lev_rel)
                throw std::runtime_error("geometry.prob_relative must be descending over refinement levels");

            last_lev_rel = prob_relative_lvl;
        }

        return prob_relative;
    }

    /** Settings that restrict the dynamically resized mesh to a discrete set of sizes.
     *
     * Fitting the box to the beam exactly means it moves with every fluctuation of the
     * particle extent, and the space-charge solver then has to rebuild its Green's
     * function each time. Allowing the box to be a little larger than requested, and
     * choosing it from a fixed set of lengths, makes a beam of nearly the same size land
     * on exactly the same mesh, which is what lets the Green's function be reused.
     */
    struct GridQuantization
    {
        bool enabled = false;              //!< restrict the box to the fit lengths
        int fit_lengths_per_doubling = 8;  //!< how many fit lengths per factor of two
    };

    /** Read the padding band, geometry.prob_relative and geometry.prob_relative_max */
    inline GridQuantization
    read_grid_quantization ()
    {
        amrex::ParmParse pp_algo("algo");
        amrex::ParmParse pp_geometry("geometry");

        std::string poisson_solver = "fft";
        pp_algo.queryAdd("poisson_solver", poisson_solver);

        amrex::Real const prob_relative_min = read_mr_prob_relative()[0];

        // By default allow a 10% band for the FFT solver, which reuses its Green's
        // function when the mesh holds still. The multigrid solver gains nothing from a
        // steady mesh, so there the box keeps fitting the beam exactly.
        amrex::Real prob_relative_max = (poisson_solver == "fft")
            ? prob_relative_min * amrex::Real(1.1) : prob_relative_min;
        pp_geometry.queryAddWithParser("prob_relative_max", prob_relative_max);

        if (prob_relative_max < prob_relative_min) {
            throw std::runtime_error(
                "geometry.prob_relative_max must be >= geometry.prob_relative[0]");
        }

        GridQuantization quant;

        amrex::Real const band = prob_relative_max / prob_relative_min;
        if (!(band > amrex::Real(1.))) {
            return quant;  // no band: fit the beam exactly, as if this were never set
        }

        // Pick the finest set of lengths whose spacing still fits inside the band, so
        // that rounding up can never exceed prob_relative_max. Keeping this a whole
        // number of lengths per doubling is what makes two boxes a factor of two apart
        // differ by exactly two, which the solver relies on to reuse a Green's function
        // built at another scale.
        quant.enabled = true;
        quant.fit_lengths_per_doubling =
            static_cast<int>(std::ceil(amrex::Real(1.) / std::log2(band)));

        return quant;
    }

    /** Smallest length not below @p extent that the mesh is allowed to take.
     *
     * The allowed lengths are 2^(k/m) meters for whole numbers k, with m fit lengths per
     * doubling. Anchoring them at 1 m, rather than at the first box seen, keeps them
     * reproducible across runs and across turns of a ring, so a beam that returns to a
     * size it had before lands on the very same mesh.
     *
     * @param[in] extent the length to round up, in meters
     * @param[in] fit_lengths_per_doubling number of allowed lengths m per factor of two
     * @return the smallest allowed length that is not below @p extent
     */
    inline amrex::Real
    smallest_fit_length (amrex::Real extent, int fit_lengths_per_doubling)
    {
        if (!(extent > amrex::Real(0.)) || !std::isfinite(extent)) { return extent; }

        auto const m = static_cast<amrex::Real>(fit_lengths_per_doubling);
        amrex::Real const k = std::ceil(std::log2(extent) * m);
        amrex::Real const length = std::exp2(k / m);

        // Guard against the rounding of log2/exp2 landing a hair below the request,
        // which would clip the beam.
        return (length < extent) ? std::exp2((k + amrex::Real(1.)) / m) : length;
    }
} // namespace impactx::initialization
