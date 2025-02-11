/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "elements/diagnostics/openPMD.H"


namespace impactx::elements::diagnostics::detail
{
#ifdef ImpactX_USE_OPENPMD
    namespace io = openPMD;


    std::pair< std::string, std::string >
    name2openPMD (const std::string& fullName)
    {
        std::string record_name = fullName;
        std::string component_name = io::RecordComponent::SCALAR;

        // we use "_" as separator in names to group vector records
        std::size_t const startComp = fullName.find_last_of('_');
        if( startComp != std::string::npos ) {  // non-scalar
            record_name = fullName.substr(0, startComp);
            component_name = fullName.substr(startComp + 1u);
        }
        return make_pair(record_name, component_name);
    }


    io::RecordComponent get_component_record (
        io::ParticleSpecies & species,
        std::string comp_name
    )
    {
        // handle scalar and non-scalar records by name
        const auto [record_name, component_name] = name2openPMD(std::move(comp_name));
        return species[record_name][component_name];
    }
#endif
} // namespace impactx::elements::diagnostics::detail
