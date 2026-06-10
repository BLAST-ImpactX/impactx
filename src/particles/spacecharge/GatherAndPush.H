/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Remi Lehe
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_GATHER_AND_PUSH_H
#define IMPACTX_GATHER_AND_PUSH_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_BaseFab.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <unordered_map>
#include <string>


namespace impactx::particles::spacecharge
{
    /** Gather force fields and push particles in x,y,z
     *
     * This gathers the space charge field with respect to particle position
     * and shape. The momentum of all particles is then pushed using a common
     * time step given by the reference particle speed and ds slice. The
     * position push is done in the lattice elements and not here.
     *
     * The mesh field arrays (space-charge force, potential) are stored at the
     * beam precision T_PC::ParticleRealType (amrex::MultiFab for double).
     *
     * @param[inout] pc container of the particles that deposited rho
     * @param[in] space_charge_field space charge force component in x,y,z per level
     * @param[in] space_charge_potential space charge potential in x,y,z per level
     * @param[in] geom geometry object
     * @param[in] slice_ds segment length in meters
     */
    template <class T_PC>
    void GatherAndPush (
        T_PC & pc,
        std::unordered_map<int, std::unordered_map<std::string, FieldMeshMF<typename T_PC::ParticleRealType>> > const & space_charge_field,
        std::unordered_map<int, FieldMeshMF<typename T_PC::ParticleRealType>> const & space_charge_potential,
        const amrex::Vector<amrex::Geometry>& geom,
        amrex::ParticleReal slice_ds
    );

}  // namespace impactx::particles::spacecharge

#endif // IMPACTX_GATHER_AND_PUSH_H
