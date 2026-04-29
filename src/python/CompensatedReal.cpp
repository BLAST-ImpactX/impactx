/* Copyright 2022-2024 The ImpactX Community
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <CompensatedReal.H>
#include <AMReX_REAL.H>

#include <cmath>

namespace py = pybind11;
using namespace impactx;


void init_compensatedreal(py::module& m)
{
    // Expose the CompensatedParticleReal typedef
    py::class_<CompensatedParticleReal>(m, "CompensatedParticleReal")
        .def(py::init<>(),
             "Default (zero) constructor for compensated ParticleReal type"
        )
        .def(py::init<amrex::ParticleReal>(),
             py::arg("value"),
             "Constructor from a ParticleReal value"
        )
        .def(py::init<CompensatedParticleReal const&>(),
             py::arg("other"),
             "Copy constructor"
        )
        .def("assign", &CompensatedParticleReal::assign,
             py::arg("value"),
             py::return_value_policy::reference_internal,
             "Assign a new value and reset compensation"
        )

        // Arithmetic operators
        .def(py::self += amrex::ParticleReal())
        .def(py::self -= amrex::ParticleReal())
        .def(py::self + amrex::ParticleReal())
        .def(py::self - amrex::ParticleReal())
        .def(amrex::ParticleReal() + py::self)
        .def(amrex::ParticleReal() - py::self)
        .def(py::self += py::self)
        .def(py::self -= py::self)  // Clang <17 false-positive: -Wself-assign-overloaded
        .def(py::self + py::self)
        .def(py::self - py::self)  // Clang <17 false-positive: -Wself-assign-overloaded

        // Comparison operators
        .def(py::self == amrex::ParticleReal())
        .def(py::self != amrex::ParticleReal())
        .def(py::self < amrex::ParticleReal())
        .def(py::self <= amrex::ParticleReal())
        .def(py::self > amrex::ParticleReal())
        .def(py::self >= amrex::ParticleReal())
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self < py::self)
        .def(py::self <= py::self)
        .def(py::self > py::self)
        .def(py::self >= py::self)

        // Value access
        .def_property_readonly("value", &CompensatedParticleReal::value,
             "Get the final compensated value (main + first-order + second-order compensation)"
        )

        // conversion to float for compatibility
        .def("__float__", &CompensatedParticleReal::value,
             "Convert to Python float"
        )
        .def("__abs__", [](CompensatedParticleReal const& self) {
            return std::abs(self.value());
        }, "Convert to absolute value")

        // Utility methods
        .def("reset_compensation", &CompensatedParticleReal::reset_compensation,
             py::return_value_policy::reference_internal,
             "Reset compensation variables to zero"
        )

        // String representation
        .def("__repr__", [](CompensatedParticleReal const& self) {
            return "CompensatedParticleReal(" + std::to_string(self.value()) + ")";
        })
        .def("__str__", [](CompensatedParticleReal const& self) {
            return std::to_string(self.value());
        })

        .def_readonly_static("__doc__",
             "A compensated floating-point type specialized for ParticleReal values.\n"
             "\n"
             "This provides high-precision arithmetic for particle tracking calculations\n"
             "using the second-order iterative Kahan-Babuska algorithm.\n"
             "\n"
             "The algorithm follows Klein (2006) to provide high-precision arithmetic that avoids\n"
             "floating-point precision errors when dealing with large numbers of small values.\n"
             "\n"
             "References:\n"
             "- https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements\n"
             "- Klein (2006). \"A generalized Kahan–Babuška-Summation-Algorithm\". in\n"
             "  Computing. 76 (3–4). Springer-Verlag: 279–293. doi:10.1007/s00607-005-0139-x"
        )
    ;
}
