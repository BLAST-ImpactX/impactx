/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_OPENPMD_H
#define IMPACTX_OPENPMD_H

#ifdef ImpactX_USE_OPENPMD
#   include <openPMD/openPMD.hpp>
#endif

#include <string>
#include <utility>


namespace impactx::elements::diagnostics::detail
{
#ifdef ImpactX_USE_OPENPMD
    /** Unclutter a real_names to openPMD record
     *
     * @TODO move to ABLASTR
     *
     * @param fullName name as in real_names variable
     * @return pair of openPMD record and component name
     */
    std::pair< std::string, std::string >
    name2openPMD (const std::string& fullName);


    //! TODO@ move to ablastr
    openPMD::RecordComponent get_component_record (
        openPMD::ParticleSpecies & species,
        std::string comp_name
    );
#endif
} // namespace impactx::elements::diagnostics::detail

#endif // IMPACTX_OPENPMD_H
