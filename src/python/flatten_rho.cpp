/* Copyright 2021-2023 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <ImpactX.H>
#include <particles/ChargeDeposition.H>

namespace py = pybind11;
using namespace impactx;


void init_flatten_rho(py::module& m) {
    m.def("flatten_charge_to_2D", [](
        ImpactX & sim
        ) {
            auto geom_3d = sim.amr_data->track_particles.m_particle_container->GetParGDB()->Geom();
            amrex::Box domain_3d = geom_3d[0].Domain();  // whole simulation index space (level 0)
            return flatten_charge_to_2D(
                sim.amr_data->track_particles.m_rho,
                domain_3d
            );
    });
}
