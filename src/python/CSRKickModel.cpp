/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/wakefields/CSRKickModel.H>

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>

#include <stdexcept>
#include <string>
#include <utility>

namespace py = pybind11;
using namespace impactx;
using namespace impactx::particles::wakefields;


impactx::particles::wakefields::CSRKickModel
make_csr_kick_model (py::object callable)
{
    return [callable] (CSRProfile const & profile, CSRElementContext const & context) -> CSRKick
    {
        py::gil_scoped_acquire gil;

        // borrowed views: only valid during the callable's execution
        py::object const py_profile = py::cast(profile, py::return_value_policy::reference);
        py::object const py_context = py::cast(context, py::return_value_policy::reference);

        py::object result = callable(py_profile, py_context);

        if (!py::isinstance<CSRKick>(result))
        {
            throw std::runtime_error(
                "The CSR kick model wrapper must return an impactx.CSRKick "
                "(install models via the csr_kick_model property) in element "
                "type '" + context.element_type + "' (name: '" +
                context.element_name + "')");
        }

        // take over the per-bin kick vectors from the temporary Python object
        return std::move(result.cast<CSRKick &>());
    };
}

void init_csr_kick_model (py::module & m)
{
    py::class_<CSRKick>(m, "CSRKick",
        "Per-bin CSR kick forces in Newtons, handed to ImpactX by a wrapped "
        "CSR kick model: the longitudinal component pt is required, the "
        "transverse components px and py are optional. Constructed from "
        "pyAMReX device vectors, which are moved in (the inputs are emptied)."
    )
        .def(py::init(
            [](
                amrex::Gpu::DeviceVector<amrex::Real> & pt,
                py::object px,
                py::object py_
            ) {
                CSRKick kick;
                kick.pt = std::move(pt);
                if (!px.is_none()) {
                    kick.px = std::move(px.cast<amrex::Gpu::DeviceVector<amrex::Real> &>());
                }
                if (!py_.is_none()) {
                    kick.py = std::move(py_.cast<amrex::Gpu::DeviceVector<amrex::Real> &>());
                }
                return kick;
            }),
            py::arg("pt"),
            py::arg("px") = py::none(),
            py::arg("py") = py::none()
        )
    ;

    py::class_<CSRProfile>(m, "CSRProfile",
        "Binned longitudinal beam profile passed to a user-provided CSR kick "
        "model. All arrays share the same binning: bin i spans "
        "[bin_min + i * bin_size, bin_min + (i+1) * bin_size) in the "
        "longitudinal coordinate t = ct (in meters) and have num_bins + 1 "
        "entries, where the last entry is a guard bin. The arrays live in "
        "device memory on GPU builds (use .to_xp() for zero-copy access) and "
        "are only valid during the model call: copy any data you retain."
    )
        .def_readonly("charge", &CSRProfile::charge,
            "line charge density lambda(t) [C/m] from nearest-grid-point "
            "deposition. sum(charge) * bin_size approximates the bunch charge [C]")
        .def_readonly("mean_x", &CSRProfile::mean_x,
            "charge-weighted mean of x per bin [m] (zero for empty bins)")
        .def_readonly("mean_y", &CSRProfile::mean_y,
            "charge-weighted mean of y per bin [m] (zero for empty bins)")
        .def_readonly("bin_min", &CSRProfile::bin_min,
            "lower edge of bin 0 in t = ct [m]")
        .def_readonly("bin_size", &CSRProfile::bin_size,
            "bin spacing [m]")
        .def_readonly("num_bins", &CSRProfile::num_bins,
            "number of kick bins (= csr_bins). The profile arrays have "
            "num_bins + 1 entries")
    ;

    py::class_<CSRElementContext>(m, "CSRElementContext",
        "Per-element, per-slice context passed to a user-provided CSR kick "
        "model. Only valid during the model call: copy any data you retain."
    )
        .def_readonly("element_name", &CSRElementContext::element_name,
            "user-given element name (empty if unnamed)")
        .def_readonly("element_type", &CSRElementContext::element_type,
            "element type, e.g. 'Sbend'")
        .def_readonly("rc", &CSRElementContext::rc,
            "radius-of-curvature magnitude |R| [m]")
        .def_readonly("signed_rc", &CSRElementContext::signed_rc,
            "signed radius of curvature R [m]")
        .def_readonly("ds", &CSRElementContext::ds,
            "element (arc) length [m]")
        .def_readonly("nslice", &CSRElementContext::nslice,
            "number of slices in the element")
        .def_readonly("slice", &CSRElementContext::slice,
            "current slice index, 0-based, in forward orientation")
        .def_readonly("s", &CSRElementContext::s,
            "distance from the element entrance at kick time [m] "
            "(= slice * slice_ds, also during back-tracking)")
        .def_readonly("slice_ds", &CSRElementContext::slice_ds,
            "slice length [m]")
        .def_readonly("ref", &CSRElementContext::ref,
            "reference particle (snapshot copy)")
    ;
}
