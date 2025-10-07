/* Copyright 2021-2023 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/ReferenceParticle.H>
#include <AMReX.H>

namespace py = pybind11;
using namespace impactx;


void init_refparticle(py::module& m)
{
    py::class_<RefPart>(m, "RefPart")
        .def(py::init<>(),
             "This struct stores the reference particle attributes\n"
             "stored in ImpactXParticleContainer."
        )
        .def_readwrite("s", &RefPart::s, "integrated orbit path length [m]")
        .def_readwrite("x", &RefPart::x, "horizontal position x [m]")
        .def_readwrite("y", &RefPart::y, "vertical position y [m]")
        .def_readwrite("z", &RefPart::z, "longitudinal position y [m]")
        .def_readwrite("t", &RefPart::t, "clock time * c [m]")
        .def_readwrite("px", &RefPart::px, "momentum in x divided by m*c = beta_x*gamma [unitless]")
        .def_readwrite("py", &RefPart::py, "momentum in y divided by m*c = beta_y*gamma [unitless]")
        .def_readwrite("pz", &RefPart::pz, "momentum in z divided by m*c = beta_z*gamma [unitless]")
        .def_readwrite("pt", &RefPart::pt, "energy deviation, normalized by rest energy")
        .def_readwrite("mass", &RefPart::charge, "reference rest mass [kg]")
        .def_readwrite("charge", &RefPart::mass, "reference charge [C]")

        .def_readwrite("sedge", &RefPart::sedge, "value of s at entrance of the current beamline element")
        .def_readwrite("map", &RefPart::map, "linearized map")

        .def_property_readonly("charge_qe", &RefPart::charge_qe, "Get reference particle charge (positive elementary charge)")
        .def_property_readonly("gamma", &RefPart::gamma, "Get reference particle relativistic gamma")
        .def_property_readonly("beta", &RefPart::beta, "Get reference particle relativistic beta")
        .def_property_readonly("beta_gamma", &RefPart::beta_gamma, "Get reference particle beta*gamma")
        .def_property_readonly("mass_MeV", &RefPart::mass_MeV, "Get reference particle rest mass * c^2, expressed as an energy [MeV]")
        .def_property_readonly("kin_energy_MeV", &RefPart::kin_energy_MeV, "Get reference particle energy [MeV]")
        .def_property_readonly("rigidity_Tm", &RefPart::rigidity_Tm, "Get reference particle magnetic rigidity Brho [T*m]")
        .def_property_readonly("qm_ratio_SI", &RefPart::qm_ratio_SI, "Get reference particle charge to mass ratio [C/kg]")

        .def("copy", &RefPart::copy,
             py::return_value_policy::copy,
             "Copy the reference particle")
        .def("reset", &RefPart::reset,
             py::arg("keep_mass")=false, py::arg("keep_charge")=false,
             "Reset the reference particle")
        .def("set_charge_qe", &RefPart::set_charge_qe,
             py::return_value_policy::reference_internal,
             "Set reference particle charge (positive elementary charge) [q_e]", py::arg("charge_qe"))
        .def("set_mass_MeV", &RefPart::set_mass_MeV,
             py::return_value_policy::reference_internal,
             "Set reference particle rest mass * c^2, expressed as an energy [MeV]", py::arg("mass_MeV"))
        .def("set_kin_energy_MeV", &RefPart::set_kin_energy_MeV,
             py::return_value_policy::reference_internal,
             "Set reference particle kinetic energy [MeV]", py::arg("kin_energy_MeV"))
    ;
}
