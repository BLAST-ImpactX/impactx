/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef PARTICLEBOUNDARY_H
#define PARTICLEBOUNDARY_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Enum.H>


namespace impactx::particles
{
    /** Boundary conditions for the particle bucket */
    AMREX_ENUM(ParticleBC,
        open,       //! particles might still be removed if they leave the simulation domain, e.g., w/ static grid geometry size for SC
        periodic,   //! wrap particles that leave the bucket around
        absorbing,  //! remove particles that leave the bucket
        reflecting  //! reflect particles that leave the bucket
    );

    /** Check boundary option and set default (open) if unset */
    ParticleBC
    get_particle_boundary_condition ();

    /** Function to apply longitudinal particle boundary (periodic, etc.)
     *
     * See @c ParticleBC.
     *
     * @param[in,out] pc the particle species container
     */
    template <class T_PC>
    void ParticleBoundary (
        T_PC & pc
    );

} // namespace impactx::particles

#endif // PARTICLEBOUNDARY_H
