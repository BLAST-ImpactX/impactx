/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Marco Garten, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_FORCEFROMSELFFIELDS_H
#define IMPACTX_FORCEFROMSELFFIELDS_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_BaseFab.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <string>
#include <unordered_map>


namespace impactx::particles::spacecharge
{
    /** Calculate the space charge force field from the electric potential
     *
     * This resets the values in scf_<component> to zero and then calculates the space
     * charge force field.
     *
     * The mesh field arrays (space-charge force, potential) are stored at the
     * beam precision as a canonical AMReX FabArray (amrex::MultiFab for double,
     * amrex::fMultiFab for float). The field-array type T_FieldMF is deduced from
     * the arguments, which is why it is templated on the FabArray type directly
     * (an alias template such as FieldMeshMF<T> is a non-deduced context).
     *
     * @tparam T_FieldMF the mesh field array type (amrex::MultiFab for double)
     * @param[inout] space_charge_field space charge force component in x,y,z per level
     * @param[in] phi scalar potential per level
     * @param[in] geom geometry object
     */
    template <class T_FieldMF = amrex::MultiFab>
    void ForceFromSelfFields (
        std::unordered_map<int, std::unordered_map<std::string, T_FieldMF> > & space_charge_field,
        std::unordered_map<int, T_FieldMF> const & phi,
        const amrex::Vector<amrex::Geometry>& geom
    );

}  // namespace impactx::particles::spacecharge

#endif // IMPACTX_FORCEFROMSELFFIELDS_H
