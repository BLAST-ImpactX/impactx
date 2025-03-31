/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef ISRPUSH_H
#define ISRPUSH_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>


namespace impactx::particles::wakefields
{

    /** Function to update particle momentum using a kick due to
     *  incoherent synchrotron radiation (ISR).
     *
     * @param[inout] pc the particle species container
     * @param[in] slice_ds slice spacing [m]
     * @param[in] rc bending radius [m]
     * @param[in] isr_order order of quantum ISR effects
     */
    void ISRPush (
        ImpactXParticleContainer & pc,
        amrex::ParticleReal slice_ds,
        amrex::ParticleReal rc,
        [[maybe_unused]] int isr_order
    );

} // namespace impactx::particles::wakefields

#endif // ISRPUSH_H
