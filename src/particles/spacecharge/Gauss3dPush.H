/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACT_GAUSS3D_PUSH_H
#define IMPACT_GAUSS3D_PUSH_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_REAL.H>


namespace impactx::particles::spacecharge
{
    /** Apply a Space-Charge Kick according to a 3D Gaussian Distribution
     *
     * This calculates the space charge fields from a 3D Gaussian distribution
     * with respect to particle position from the following paper:
     * J. Qiang et al., "Two-and-a-half dimensional symplectic space-charge
     * solver", LBNL Report Number: LBNL-2001674 (2025).
     * The momentum of all particles is then pushed using a common
     * time step given by the reference particle speed and ds slice. The
     * position push is done in the lattice elements and not here.
     *
     * @param[inout] pc container of the particles that deposited rho
     * @param[in] slice_ds segment length in meters
     */
    void Gauss3dPush (
        ImpactXParticleContainer & pc,
        amrex::ParticleReal slice_ds
    );

}  // namespace impactx::particles::spacecharge

#endif // IMPACT_GAUSS3D_PUSH_H
