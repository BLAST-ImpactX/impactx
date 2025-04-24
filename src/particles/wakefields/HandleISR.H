/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef HANDLE_ISR_H
#define HANDLE_ISR_H

#include "particles/ImpactXParticleContainer.H"
#include "ISRBendElement.H"
#include "ISRPush.H"

#include <cmath>
#include <vector>


namespace impactx::particles::wakefields
{
    /** Function to handle ISR in bending dipoles.
     *
     * @param[in] particle_container the particle species container
     * @param[in] element_variant variant type of the lattice element
     * @param[in] slice_ds slice spacing along s
     */
    template <typename T_Element>
    void HandleISR (
        impactx::ImpactXParticleContainer& particle_container,
        T_Element const& element_variant,
        amrex::Real slice_ds
    )
    {
        BL_PROFILE("impactx::particles::wakefields::HandleISR")

        amrex::ParmParse pp_algo("algo");
        bool isr = false;
        pp_algo.queryAdd("isr", isr);

        // Call the ISR bend function
        auto const [element_has_isr, R] = impactx::particles::wakefields::ISRBendElement(element_variant, particle_container.GetRefParticle());

        // Enter loop if lattice has bend element
        if (isr && element_has_isr)
        {

            int isr_order = 1;
            pp_algo.queryAddWithParser("isr_order", isr_order);

            // Call function to kick particles with wake
            impactx::particles::wakefields::ISRPush(particle_container, slice_ds, R, isr_order);
        }
    }

} // namespace impactx::particles::wakefields

#endif // HANDLE_ISR_H
