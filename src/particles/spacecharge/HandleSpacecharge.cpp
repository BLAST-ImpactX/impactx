/* Copyright 2022-2025 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "HandleSpacecharge.H"

#include "initialization/AmrCoreData.H"
#include "particles/ImpactXParticleContainer.H"
#include "particles/spacecharge/ForceFromSelfFields.H"
#include "particles/spacecharge/GatherAndPush.H"
#include "particles/spacecharge/PoissonSolve.H"
#include "particles/transformation/CoordinateTransformation.H"

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <memory>


namespace impactx::particles::spacecharge
{
    void HandleSpacecharge (
        std::unique_ptr<initialization::AmrCoreData> & amr_data,
        std::function<void()> ResizeMesh,
        amrex::Real slice_ds
    )
    {
        BL_PROFILE("impactx::particles::wakefields::HandleSpacecharge")

        amrex::ParmParse const pp_algo("algo");
        bool space_charge = false;
        pp_algo.query("space_charge", space_charge);

        // turn off if disabled by user
        if (!space_charge) { return; }

        // turn off if less than 2 particles
        if (amr_data->track_particles.m_particle_container->TotalNumberOfParticles(true, false) < 2) { return; }

        // transform from x',y',t to x,y,z
        transformation::CoordinateTransformation(
            *amr_data->track_particles.m_particle_container,
            CoordSystem::t
        );

        // Note: The following operation assume that
        // the particles are in x, y, z coordinates.

        // Resize the mesh, based on `amr_data->track_particles.m_particle_container` extent
        ResizeMesh();

        // Redistribute particles in the new mesh in x, y, z
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
            amr_data->Geom(),
            slice_ds
        );

        // transform from x,y,z to x',y',t
        transformation::CoordinateTransformation(
            *amr_data->track_particles.m_particle_container,
            CoordSystem::s
        );

        // for later: original Impact implementation as an option
        // Redistribute particles in x',y',t
        //   TODO: only needed if we want to gather and push space charge
        //         in x',y',t
        //   TODO: change geometry beforehand according to transformation
        //m_amr_data->track_particles.m_particle_container->Redistribute();
        //
        // in original Impact, we gather and space-charge push in x',y',t ,
        // assuming that the distribution did not change
    }

} // namespace impactx::particles::spacecharge
