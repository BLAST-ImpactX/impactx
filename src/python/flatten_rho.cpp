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
            // geometry and the charge-density field are precision-agnostic
            amrex::Box const domain_3d = sim.amr_data->Geom(0).Domain();  // whole simulation index space (level 0)
            return flatten_charge_to_2D(
                sim.amr_data->rho(),
                domain_3d
            );
    });
}
