#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "synergia/simulation/collective_operator_options.h"
#include "synergia/simulation/implemented_collective_options.h"

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE(collective, m)
{
  py::class_<Space_charge_2d_open_hockney_options>(
    m, "Space_charge_2d_open_hockney_options")
    .def(py::init<int, int, int>(),
         "Construct the space charge 2d open hockney solver.",
         "gridx"_a,
         "gridy"_a,
         "gridz"_a)
    .def_readwrite("comm_group_size",
                   &Space_charge_2d_open_hockney_options::comm_group_size,
                   "Communication group size (must be 1 on GPUs).");

  py::class_<Space_charge_3d_open_hockney_options>(
    m, "Space_charge_3d_open_hockney_options")
    .def(py::init<int, int, int>(),
         "Construct the space charge 3d open hockney solver.",
         "gridx"_a,
         "gridy"_a,
         "gridz"_a)
    .def_readwrite("comm_group_size",
                   &Space_charge_3d_open_hockney_options::comm_group_size,
                   "Communication group size (must be 1 on GPUs).");

#ifdef BUILD_FD_SPACE_CHARGE_SOLVER
  py::class_<Space_charge_3d_fd_options>(m, "Space_charge_3d_fd_options")
    .def(py::init<int, int, int>(),
         "Construct the space charge 3d finite difference solver solver.",
         "gridx"_a,
         "gridy"_a,
         "gridz"_a)
    .def_readwrite("comm_group_size",
                   &Space_charge_3d_fd_options::comm_group_size,
                   "Communication group size (must be 1 on GPUs).");
#endif

  py::class_<Space_charge_rectangular_options>(
    m, "Space_charge_rectangular_options")
    .def(py::init<std::array<int, 3> const&, std::array<double, 3> const&>(),
         "Construct the rectangular space charge solver.",
         "grid_shape"_a,
         "pipe_size"_a)
    .def_readwrite("comm_group_size",
                   &Space_charge_rectangular_options::comm_group_size,
                   "Communication group size (must be 1 on GPUs).");

  py::class_<Impedance_options>(m, "Impedance_options")
    .def(py::init<std::string const&, std::string const&>(),
         "Options for constructing the impedance operator.",
         "wake_file"_a,
         "wake_type"_a)

    .def(py::init<std::string const&, std::string const&, int>(),
         "Construct the impedance operator.",
         "wake_file"_a,
         "wake_type"_a,
         "z_grid"_a=1000)

    .def_readwrite("z_grid",
                   &Impedance_options::z_grid,
                   "Size of z-grid (int, default to 1000).")

    .def_readwrite("full_machine",
                   &Impedance_options::full_machine,
                   "Full machine (boolean, default to false).")

    .def_readwrite("nstored_turns",
                   &Impedance_options::nstored_turns,
                   "Number of stored turns (int, default to 15).")

    .def_readwrite("num_buckets",
                   &Impedance_options::num_buckets,
                   "Number of buckets (int, default to 1).")

    .def_readwrite("orbit_length",
                   &Impedance_options::orbit_length,
                   "Orbit length (double, default to 1.0).")

    .def_readwrite("bunch_spacing",
                   &Impedance_options::bunch_spacing,
                   "Bunch spacing (double, default to 1.0).")
                   
     .def_readwrite("mwf_xlead",
                    &Impedance_options::mwf_xlead,
                    "factor to scale x leading wake")

     .def_readwrite("mwf_xtrail",
                    &Impedance_options::mwf_xtrail,
                    "factor to scale x trailing wake")
                    
     .def_readwrite("mwf_ylead",
                    &Impedance_options::mwf_ylead,
                    "factor to scale y leading wake")
                    
     .def_readwrite("mwf_ytrail",
                    &Impedance_options::mwf_ytrail,
                    "factor to scale y trailing wake")
                    
     .def_readwrite("mwf_zwake",
                    &Impedance_options::mwf_zwake,
                    "factor to scale z (longitudinal) wake")
                    
                    ; // terminates Impedance_options class wrappings


  py::class_<Dummy_CO_options>(m, "Dummy_CO_options")
    .def(py::init<>(), "Construct a dummy collective operator.");
}
