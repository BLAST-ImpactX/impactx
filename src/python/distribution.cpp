/* Copyright 2021-2023 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/distribution/All.H>
#include <initialization/InitDistribution.H>

#include <AMReX.H>
#include <AMReX_REAL.H>

#include <variant>

namespace py = pybind11;
using namespace impactx;


void init_distribution(py::module& m)
{
    py::module_ const md = m.def_submodule(
        "distribution",
        "Particle beam distributions in ImpactX"
    );

    py::class_<distribution::Gaussian>(md, "Gaussian")
        .def(py::init<
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal,
                amrex::ParticleReal, amrex::ParticleReal
             >(),
             py::arg("lambda_x"), py::arg("lambda_y"), py::arg("lambda_t"),
             py::arg("lambda_px"), py::arg("lambda_py"), py::arg("lambda_pt"),
             py::arg("mu_x_px")=0.0, py::arg("mu_y_py")=0.0, py::arg("mu_t_pt")=0.0,
             py::arg("mean_x")=0.0, py::arg("mean_y")=0.0, py::arg("mean_t")=0.0,
             py::arg("mean_px")=0.0, py::arg("mean_py")=0.0, py::arg("mean_pt")=0.0,
             py::arg("dispersion_x")=0.0, py::arg("dispersion_px")=0.0,
             py::arg("dispersion_y")=0.0, py::arg("dispersion_py")=0.0,
             "A 6D Gaussian distribution"
        );

    py::class_<distribution::Kurth4D>(md, "Kurth4D")
        .def(py::init<
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal
         >(),
             py::arg("lambda_x"), py::arg("lambda_y"), py::arg("lambda_t"),
             py::arg("lambda_px"), py::arg("lambda_py"), py::arg("lambda_pt"),
             py::arg("mu_x_px")=0.0, py::arg("mu_y_py")=0.0, py::arg("mu_t_pt")=0.0,
             py::arg("mean_x")=0.0, py::arg("mean_y")=0.0, py::arg("mean_t")=0.0,
             py::arg("mean_px")=0.0, py::arg("mean_py")=0.0, py::arg("mean_pt")=0.0,
             py::arg("dispersion_x")=0.0, py::arg("dispersion_px")=0.0,
             py::arg("dispersion_y")=0.0, py::arg("dispersion_py")=0.0,
             "A 4D Kurth distribution transversely + a uniform distribution\n"
             "in t + a Gaussian distribution in pt"
        );

    py::class_<distribution::Kurth6D>(md, "Kurth6D")
        .def(py::init<
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal
             >(),
             py::arg("lambda_x"), py::arg("lambda_y"), py::arg("lambda_t"),
             py::arg("lambda_px"), py::arg("lambda_py"), py::arg("lambda_pt"),
             py::arg("mu_x_px")=0.0, py::arg("mu_y_py")=0.0, py::arg("mu_t_pt")=0.0,
             py::arg("mean_x")=0.0, py::arg("mean_y")=0.0, py::arg("mean_t")=0.0,
             py::arg("mean_px")=0.0, py::arg("mean_py")=0.0, py::arg("mean_pt")=0.0,
             py::arg("dispersion_x")=0.0, py::arg("dispersion_px")=0.0,
             py::arg("dispersion_y")=0.0, py::arg("dispersion_py")=0.0,
             "A 6D Kurth distribution\n\n"
             "R. Kurth, Quarterly of Applied Mathematics vol. 32, pp. 325-329 (1978)\n"
             "C. Mitchell, K. Hwang and R. D. Ryne, IPAC2021, WEPAB248 (2021)"
        );

    py::class_<distribution::KVdist>(md, "KVdist")
        .def(py::init<
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal
             >(),
             py::arg("lambda_x"), py::arg("lambda_y"), py::arg("lambda_t"),
             py::arg("lambda_px"), py::arg("lambda_py"), py::arg("lambda_pt"),
             py::arg("mu_x_px")=0.0, py::arg("mu_y_py")=0.0, py::arg("mu_t_pt")=0.0,
             py::arg("mean_x")=0.0, py::arg("mean_y")=0.0, py::arg("mean_t")=0.0,
             py::arg("mean_px")=0.0, py::arg("mean_py")=0.0, py::arg("mean_pt")=0.0,
             py::arg("dispersion_x")=0.0, py::arg("dispersion_px")=0.0,
             py::arg("dispersion_y")=0.0, py::arg("dispersion_py")=0.0,
             "A K-V distribution transversely + a uniform distribution\n"
             "in t + a Gaussian distribution in pt"
        );

    py::class_<distribution::Empty>(md, "Empty")
        .def(py::init<>(),
             "Sets all values to zero."
        );

    py::class_<distribution::Semigaussian>(md, "Semigaussian")
        .def(py::init<
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal
             >(),
             py::arg("lambda_x"), py::arg("lambda_y"), py::arg("lambda_t"),
             py::arg("lambda_px"), py::arg("lambda_py"), py::arg("lambda_pt"),
             py::arg("mu_x_px")=0.0, py::arg("mu_y_py")=0.0, py::arg("mu_t_pt")=0.0,
             py::arg("mean_x")=0.0, py::arg("mean_y")=0.0, py::arg("mean_t")=0.0,
             py::arg("mean_px")=0.0, py::arg("mean_py")=0.0, py::arg("mean_pt")=0.0,
             py::arg("dispersion_x")=0.0, py::arg("dispersion_px")=0.0,
             py::arg("dispersion_y")=0.0, py::arg("dispersion_py")=0.0,
             "A 6D Semi-Gaussian distribution (uniform in position, Gaussian in momentum)."
        );

    py::class_<distribution::Thermal>(md, "Thermal")
        .def(py::init<
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal
             >(),
             py::arg("k"), py::arg("kT"), py::arg("kT_halo"),
             py::arg("normalize"), py::arg("normalize_halo"),
             py::arg("halo")=0.0,
             "A stationary thermal or bithermal distribution\n\n"
             "R. D. Ryne, J. Qiang, and A. Adelmann, in Proc. EPAC2004, pp. 1942-1944 (2004)"
        );

    py::class_<distribution::Triangle>(md, "Triangle")
        .def(py::init<
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal
             >(),
             py::arg("lambda_x"), py::arg("lambda_y"), py::arg("lambda_t"),
             py::arg("lambda_px"), py::arg("lambda_py"), py::arg("lambda_pt"),
             py::arg("mu_x_px")=0.0, py::arg("mu_y_py")=0.0, py::arg("mu_t_pt")=0.0,
             py::arg("mean_x")=0.0, py::arg("mean_y")=0.0, py::arg("mean_t")=0.0,
             py::arg("mean_px")=0.0, py::arg("mean_py")=0.0, py::arg("mean_pt")=0.0,
             py::arg("dispersion_x")=0.0, py::arg("dispersion_px")=0.0,
             py::arg("dispersion_y")=0.0, py::arg("dispersion_py")=0.0,
             "A triangle distribution for laser-plasma acceleration related applications.\n\n"
             "A ramped, triangular current profile with a Gaussian energy spread (possibly correlated).\n"
             "The transverse distribution is a 4D waterbag."
        );

    py::class_<distribution::Waterbag>(md, "Waterbag")
        .def(py::init<
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal,
                 amrex::ParticleReal, amrex::ParticleReal
             >(),
             py::arg("lambda_x"), py::arg("lambda_y"), py::arg("lambda_t"),
             py::arg("lambda_px"), py::arg("lambda_py"), py::arg("lambda_pt"),
             py::arg("mu_x_px")=0.0, py::arg("mu_y_py")=0.0, py::arg("mu_t_pt")=0.0,
             py::arg("mean_x")=0.0, py::arg("mean_y")=0.0, py::arg("mean_t")=0.0,
             py::arg("mean_px")=0.0, py::arg("mean_py")=0.0, py::arg("mean_pt")=0.0,
             py::arg("dispersion_x")=0.0, py::arg("dispersion_px")=0.0,
             py::arg("dispersion_y")=0.0, py::arg("dispersion_py")=0.0,
             "A 6D Waterbag distribution"
        );

    py::class_<Envelope>(m, "Envelope")
        .def(py::init<>())
        .def(py::init<CovarianceMatrix, amrex::ParticleReal>())
        .def_property("envelope", &Envelope::covariance_matrix, &Envelope::set_covariance_matrix)
        .def_property("beam_intensity", &Envelope::beam_intensity, &Envelope::set_beam_intensity)
    ;

    m.def("create_envelope", &initialization::create_envelope);
}
