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
        ) -> py::object {
            // geometry is precision-agnostic
            amrex::Box const domain_3d = sim.amr_data->Geom(0).Domain();  // whole simulation index space (level 0)
            // dispatch on the runtime beam precision; the charge-density field
            // (and thus the flattened 2D result) is at the concrete precision
            return sim.amr_data_visit([&domain_3d](auto& data) -> py::object {
                return py::cast(
                    flatten_charge_to_2D(
                        data.track_particles.m_rho,
                        domain_3d
                    )
                );
            });
    });
}
