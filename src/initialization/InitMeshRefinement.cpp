/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Remi Lehe
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"
#include "initialization/Algorithms.H"
#include "initialization/InitAmrCore.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/distribution/Waterbag.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>


namespace impactx
{
namespace detail
{
    amrex::Vector<amrex::Real>
    read_mr_prob_relative ()
    {
        amrex::ParmParse pp_algo("algo");
        amrex::ParmParse pp_amr("amr");
        amrex::ParmParse pp_geometry("geometry");

        int max_level = 0;
        pp_amr.queryWithParser("max_level", max_level);

        std::string poisson_solver = "fft";
        pp_algo.queryAdd("poisson_solver", poisson_solver);

        // The box is expanded beyond the min and max of the particle beam.
        amrex::Vector<amrex::Real> prob_relative(max_level + 1, 1.0);
        prob_relative[0] = 3.0;  // top/bottom pad the beam on the lowest level by default by its width
        pp_geometry.queryarr("prob_relative", prob_relative);

        if (prob_relative[0] < 3.0 && poisson_solver == "multigrid")
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

    /** True if [lo, hi] has no extent, up to the rounding of the coordinates (8 ulp) */
    bool
    is_flat (amrex::ParticleReal lo, amrex::ParticleReal hi)
    {
        amrex::ParticleReal const eps = std::numeric_limits<amrex::ParticleReal>::epsilon();
        return hi - lo <= amrex::ParticleReal(8) * eps * std::max(std::abs(lo), std::abs(hi));
    }

    /** Throw if the beam has no extent in a direction that the space charge model needs.
     *
     * https://github.com/BLAST-ImpactX/impactx/issues/44
     * A direction counts as flat when its extent is zero up to the rounding of the particle
     * coordinates.  Such a direction gives the mesh no size, and the 2.5D and Gaussian
     * models divide by the beam extent.  The 2D model of an unbunched beam does not resolve
     * the longitudinal direction (the charge is projected to x-y), so it accepts a beam
     * without longitudinal extent, e.g. a DC beam with all particles at the same time of
     * flight.  An empty beam (min > max) is not checked, and an extent that is not
     * finite (an infinite position, or a NaN that reaches the min/max reduction) gets
     * its own error.
     */
    void
    check_beam_extent (
        amrex::ParticleReal x_min, amrex::ParticleReal y_min, amrex::ParticleReal z_min,
        amrex::ParticleReal x_max, amrex::ParticleReal y_max, amrex::ParticleReal z_max,
        SpaceChargeAlgo space_charge
    )
    {
        if (x_min > x_max) { return; }  // no particles

        for (amrex::ParticleReal const v : {x_min, y_min, z_min, x_max, y_max, z_max}) {
            if (!std::isfinite(v)) {
                throw std::runtime_error("ResizeMesh: the extent of the beam is not finite");
            }
        }
        std::string const see = " See https://github.com/BLAST-ImpactX/impactx/issues/44";
        if (is_flat(x_min, x_max) || is_flat(y_min, y_max)) {
            throw std::runtime_error("Flat beam detected: the beam has no transverse extent, "
                                     "which is not yet supported." + see);
        }
        if (is_flat(z_min, z_max) && space_charge != SpaceChargeAlgo::True_2D) {
            std::string msg = "Flat beam detected: the beam has no longitudinal extent, ";
            if (space_charge == SpaceChargeAlgo::False) {
                msg += "which sizing the mesh needs.";
            } else {
                msg += "which algo.space_charge = " + to_string(space_charge)
                       + " does not support yet. For an unbunched (e.g., DC) beam, use "
                         "algo.space_charge = 2D.";
            }
            throw std::runtime_error(msg + see);
        }
    }
}

    void ImpactX::ResizeMesh ()
    {
        BL_PROFILE("ImpactX::ResizeMesh");

        using namespace amrex::literals; // for _rt and _prt

        auto const space_charge = get_space_charge_algo();
        if (space_charge == SpaceChargeAlgo::False) {
            ablastr::warn_manager::WMRecordWarning(
                "ImpactX::ResizeMesh",
                "This is a simulation without space charge. "
                "ResizeMesh (and pc.Redistribute) should only be called "
                "in space charge simulations.",
                ablastr::warn_manager::WarnPriority::high
            );
        }

        // Extract the min and max of the particle positions
        auto [x_min, y_min, z_min, x_max, y_max, z_max] =
            amr_data->track_particles.m_particle_container->MinAndMaxPositions();

        // guard for flat beams (only the 2D model accepts a beam without longitudinal extent)
        detail::check_beam_extent(x_min, y_min, z_min, x_max, y_max, z_max, space_charge);

        amrex::ParmParse pp_geometry("geometry");
        bool dynamic_size = true;
        pp_geometry.query("dynamic_size", dynamic_size);

        amrex::Vector<amrex::RealBox> rb(amr_data->finestLevel() + 1);  // extent per level
        if (dynamic_size)
        {
            // The 2D model of an unbunched beam does not resolve the longitudinal
            // direction (the charge is projected to x-y), so a beam without longitudinal
            // extent only needs some finite box length: it gets the larger transverse
            // width of the beam.  A longitudinal extent below 1e-6 of that width (e.g.
            // jitter in the time of flight of a DC beam) is treated as none, which keeps
            // the cells from becoming degenerate; a short but resolved beam keeps its box.
            amrex::ParticleReal const width = std::max(x_max - x_min, y_max - y_min);
            bool const no_z_extent =
                detail::is_flat(z_min, z_max) || z_max - z_min < 1.0e-6_prt * width;
            if (space_charge == SpaceChargeAlgo::True_2D && no_z_extent) {
                amrex::ParticleReal const z_center = 0.5_prt * (z_min + z_max);
                z_min = z_center - 0.5_prt * width;
                z_max = z_center + 0.5_prt * width;
            }

            // The coarsest level is expanded (or reduced) relative the min and max of particles.
            auto const prob_relative = detail::read_mr_prob_relative();

            amrex::Real const frac = prob_relative[0];
            amrex::RealVect const beam_min(x_min, y_min, z_min);
            amrex::RealVect const beam_max(x_max, y_max, z_max);
            amrex::RealVect const beam_width(beam_max - beam_min);

            amrex::RealVect const beam_padding = beam_width * (frac - 1_rt) * 0.5_rt;
            //                           added to the beam extent --^         ^-- box half above/below the beam

            // In AMReX, all levels have the same problem domain, that of the
            // coarsest level, even if only partly covered.
            for (int lev = 0; lev <= amr_data->finestLevel(); ++lev)
            {
                rb[lev].setLo(beam_min - beam_padding);
                rb[lev].setHi(beam_max + beam_padding);
            }
        }
        else
        {
            // note: we read and set the size again because an interactive /
            //       Python user might have changed it between steps
            amrex::Vector<amrex::Real> prob_lo;
            amrex::Vector<amrex::Real> prob_hi;
            pp_geometry.getarr("prob_lo", prob_lo);
            pp_geometry.getarr("prob_hi", prob_hi);

            rb[0] = {prob_lo.data(), prob_hi.data()};

            if (amr_data->maxLevel() > 1)
                amrex::Abort("Did not implement ResizeMesh for static domains and >1 MR levels.");
        }

        // updating geometry.prob_lo/hi for consistency
        amrex::Vector<amrex::Real> const prob_lo = {rb[0].lo()[0], rb[0].lo()[1], rb[0].lo()[2]};
        amrex::Vector<amrex::Real> const prob_hi = {rb[0].hi()[0], rb[0].hi()[1], rb[0].hi()[2]};
        pp_geometry.addarr("prob_lo", prob_lo);
        pp_geometry.addarr("prob_hi", prob_hi);

        // Resize the domain size
        amrex::Geometry::ResetDefaultProbDomain(rb[0]);

        for (int lev = 0; lev <= amr_data->finestLevel(); ++lev)
        {
            amrex::Geometry g = amr_data->Geom(lev);
            g.ProbDomain(rb[lev]);
            amr_data->SetGeometry(lev, g);

            amr_data->track_particles.m_particle_container->SetParticleGeometry(lev, g);
        }
    }
} // namespace impactx
