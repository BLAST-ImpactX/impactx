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

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>


namespace impactx::particles
{
    /** Function to apply longitudinal particle boundary (periodic, etc.)
     *
     * @param[inout] pc the particle species container
     */
    void ParticleBoundary (
        ImpactXParticleContainer & pc
    );

} // namespace impactx::particles

#endif // PARTICLEBOUNDARY_H
