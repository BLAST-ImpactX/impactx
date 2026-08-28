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
#include "initialization/InitMeshRefinement.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/distribution/Waterbag.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>


namespace impactx
{

    void ImpactX::ResizeMesh ()
    {
        BL_PROFILE("ImpactX::ResizeMesh");

        using namespace amrex::literals; // for _rt and _prt

        {
            auto space_charge = get_space_charge_algo();
            if (space_charge == SpaceChargeAlgo::False)
                ablastr::warn_manager::WMRecordWarning(
                    "ImpactX::ResizeMesh",
                    "This is a simulation without space charge. "
                    "ResizeMesh (and pc.Redistribute) should only be called "
                    "in space charge simulations.",
                    ablastr::warn_manager::WarnPriority::high
                );
        }

        // Extract the min and max of the particle positions
        auto const [x_min, y_min, z_min, x_max, y_max, z_max] = amr_data->track_particles.m_particle_container->MinAndMaxPositions();

        // guard for flat beams:
        //   https://github.com/BLAST-ImpactX/impactx/issues/44
        if (x_min == x_max || y_min == y_max || z_min == z_max)
            throw std::runtime_error("Flat beam detected. This is not yet supported: https://github.com/BLAST-ImpactX/impactx/issues/44");

        amrex::ParmParse pp_geometry("geometry");
        bool dynamic_size = true;
        pp_geometry.query("dynamic_size", dynamic_size);

        amrex::Vector<amrex::RealBox> rb(amr_data->finestLevel() + 1);  // extent per level
        if (dynamic_size)
        {
            // The coarsest level is expanded (or reduced) relative the min and max of particles.
            auto const prob_relative = initialization::read_mr_prob_relative();

            amrex::Real const frac = prob_relative[0];
            amrex::RealVect const beam_min(x_min, y_min, z_min);
            amrex::RealVect const beam_max(x_max, y_max, z_max);
            amrex::RealVect const beam_width(beam_max - beam_min);
            amrex::RealVect const beam_center = (beam_min + beam_max) * 0.5_rt;

            auto const quant = initialization::read_grid_quantization();

            amrex::RealVect box_lo;
            amrex::RealVect box_hi;

            if (!quant.enabled)
            {
                amrex::RealVect const beam_padding = beam_width * (frac - 1_rt) * 0.5_rt;
                //                       added to the beam extent --^         ^-- box half above/below the beam
                box_lo = beam_min - beam_padding;
                box_hi = beam_max + beam_padding;
            }
            else
            {
                amrex::RealVect box_width = beam_width * frac;

                // Round the box up to an allowed length, so that a beam of nearly the
                // same size lands on exactly the same grid and the space-charge solver
                // can reuse its Green's function. geometry.prob_relative stays the
                // minimum padding and geometry.prob_relative_max bounds the result.
                // The stretch the space-charge solver applies longitudinally. It reaches
                // the solver as a velocity and is turned back into a stretch there, which
                // for an ultrarelativistic beam loses about eps*gamma^2 of precision in
                // the round trip: 1 - 1/gamma^2 cancels away the leading digits. Deriving
                // it here the same way, rather than from gamma directly, makes that loss
                // identical on both sides so it cancels, and the mesh the solver sees is
                // the one we chose.
                amrex::ParticleReal const pt_ref =
                    amr_data->track_particles.m_particle_container->GetRefParticle().pt;
                amrex::ParticleReal const beta_s =
                    std::sqrt(1.0_prt - 1.0_prt / amrex::Math::powi<2>(pt_ref));
                auto const beta_z = static_cast<amrex::Real>(beta_s);
                amrex::Real const gamma_z = 1.0_rt / std::sqrt(1.0_rt - beta_z * beta_z);

                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                {
                    // The solver sees the longitudinal direction stretched, so that is the
                    // length which has to be an allowed one. A change of the reference
                    // energy then resizes the mesh by itself, with no separate tolerance.
                    amrex::Real const boost = (d == 2) ? gamma_z : 1.0_rt;
                    box_width[d] = initialization::smallest_fit_length(
                        box_width[d] * boost, quant.fit_lengths_per_doubling) / boost;
                }

                amrex::RealVect const box_half = box_width * 0.5_rt;
                box_lo = beam_center - box_half;
                box_hi = beam_center + box_half;
            }

            // In AMReX, all levels have the same problem domain, that of the
            // coarsest level, even if only partly covered.
            for (int lev = 0; lev <= amr_data->finestLevel(); ++lev)
            {
                rb[lev].setLo(box_lo);
                rb[lev].setHi(box_hi);
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
