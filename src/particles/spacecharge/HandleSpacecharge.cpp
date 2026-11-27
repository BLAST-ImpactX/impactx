/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "HandleSpacecharge.H"

#include "initialization/Algorithms.H"
#include "initialization/AmrCoreData.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/spacecharge/ForceFromSelfFields.H"
#include "particles/spacecharge/GatherAndPush.H"
#include "particles/spacecharge/Gauss3dPush.H"
#include "particles/spacecharge/Gauss2p5dPush.H"
#include "particles/spacecharge/PoissonSolve.H"
#include "particles/transformation/CoordinateTransformation.H"

#include <AMReX_REAL.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>


namespace impactx::particles::spacecharge
{
    void HandleSpacecharge (
        std::unique_ptr<initialization::AmrCoreData> & amr_data,
        std::function<void()> ResizeMesh,
        amrex::Real slice_ds
    )
    {
        if (slice_ds == 0.0) { return; }

        BL_PROFILE("impactx::particles::wakefields::HandleSpacecharge")

        auto space_charge = get_space_charge_algo();

        // turn off if disabled by user
        if (space_charge == SpaceChargeAlgo::False) { return; }

        // turn off if less than 2 particles
        if (amr_data->track_particles.m_particle_container->TotalNumberOfParticles(true, false) < 2) { return; }

        if (space_charge != SpaceChargeAlgo::True_2D && space_charge != SpaceChargeAlgo::True_2p5D && space_charge != SpaceChargeAlgo::Gauss2p5D)
        {
            // transform from x',y',t to x,y,z
            transformation::CoordinateTransformation(
                *amr_data->track_particles.m_particle_container,
                CoordSystem::t
            );
        }

        // the Gaussian models divide by the beam extent (the mesh models check it in ResizeMesh)
        if (space_charge == SpaceChargeAlgo::Gauss3D || space_charge == SpaceChargeAlgo::Gauss2p5D)
        {
            auto const [x_min, y_min, z_min, x_max, y_max, z_max] =
                amr_data->track_particles.m_particle_container->MinAndMaxPositions();
            CheckBeamExtent(x_min, y_min, z_min, x_max, y_max, z_max, space_charge);
        }

        if (space_charge == SpaceChargeAlgo::Gauss3D)
        {
            Gauss3dPush(*amr_data->track_particles.m_particle_container, slice_ds);
        }
        else if (space_charge == SpaceChargeAlgo::Gauss2p5D)
        {
            Gauss2p5dPush(*amr_data->track_particles.m_particle_container, slice_ds);
        }
        else if (space_charge == SpaceChargeAlgo::True_3D || space_charge == SpaceChargeAlgo::True_2D || space_charge == SpaceChargeAlgo::True_2p5D)
        {
            // Note: The following operations assume that
            // the particles are in x, y, z coordinates.

            // Resize the mesh, based on `amr_data->track_particles.m_particle_container` extent
            ResizeMesh();

            // Redistribute particles in the new mesh in:
            // 3D: x, y, z
            // 2D: x, y, t
            amr_data->track_particles.m_particle_container->Redistribute();

            // charge deposition
            amr_data->track_particles.m_particle_container->DepositCharge(
                amr_data->track_particles.m_rho,
                amr_data->refRatio()
            );

            // poisson solve in x,y,z
            spacecharge::PoissonSolve(
                *amr_data->track_particles.m_particle_container,
                amr_data->track_particles.m_rho,
                amr_data->track_particles.m_rho_2d,
                amr_data->track_particles.m_phi,
                amr_data->refRatio()
            );

            // calculate force in x,y,z
            spacecharge::ForceFromSelfFields(
                amr_data->track_particles.m_space_charge_field,
                amr_data->track_particles.m_phi,
                amr_data->Geom()
            );

            // gather and space-charge push in x,y,z , assuming the space-charge
            // field is the same before/after transformation
            // TODO: This is currently using linear order.
            spacecharge::GatherAndPush(
                *amr_data->track_particles.m_particle_container,
                amr_data->track_particles.m_space_charge_field,
                amr_data->track_particles.m_phi,
                amr_data->Geom(),
                slice_ds
            );
        }

        if (space_charge != SpaceChargeAlgo::True_2D && space_charge != SpaceChargeAlgo::True_2p5D && space_charge != SpaceChargeAlgo::Gauss2p5D)
        {
            // transform from x,y,z to x',y',t
            transformation::CoordinateTransformation(
                *amr_data->track_particles.m_particle_container,
                CoordSystem::s
            );
        }

        // for later: original Impact implementation as an option for 3D space charge to:
        // Redistribute particles in x',y',t
        //   TODO: only needed if we want to gather and push space charge
        //         in x',y',t
        //   TODO: change geometry beforehand according to transformation
        //m_amr_data->track_particles.m_particle_container->Redistribute();
        //
        // in original Impact, we gather and space-charge push in x',y',t ,
        // assuming that the distribution did not change
    }


    bool CheckBeamExtent (
        amrex::ParticleReal x_min, amrex::ParticleReal y_min, amrex::ParticleReal z_min,
        amrex::ParticleReal x_max, amrex::ParticleReal y_max, amrex::ParticleReal z_max,
        SpaceChargeAlgo space_charge
    )
    {
        // https://github.com/BLAST-ImpactX/impactx/issues/44
        auto const is_flat = [] (amrex::ParticleReal lo, amrex::ParticleReal hi) {
            amrex::ParticleReal const eps = std::numeric_limits<amrex::ParticleReal>::epsilon();
            return hi - lo <= amrex::ParticleReal(8) * eps * std::max(std::abs(lo), std::abs(hi));
        };
        std::string const unsupported =
            " This is not yet supported: https://github.com/BLAST-ImpactX/impactx/issues/44";
        if (is_flat(x_min, x_max) || is_flat(y_min, y_max)) {
            throw std::runtime_error(
                "Flat beam detected: the beam has no transverse extent." + unsupported);
        }
        bool const flat_longitudinal = is_flat(z_min, z_max);
        if (flat_longitudinal && space_charge != SpaceChargeAlgo::True_2D) {
            throw std::runtime_error(
                "Flat beam detected: the beam has no longitudinal extent. For an unbunched "
                "(e.g., DC) beam, use algo.space_charge = 2D." + unsupported);
        }
        return flat_longitudinal;
    }

} // namespace impactx::particles::spacecharge
