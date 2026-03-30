
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "lsexpr.h"

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE(pylsexpr, m)
{
  py::class_<Lsexpr>(m, "Lsexpr")
    .def(py::init<>())
    .def(py::init<std::string const&>())
    .def(py::init<std::string const&, std::string const&>())
    .def(py::init<int>())
    .def(py::init<int, std::string const&>())
    .def(py::init<double>())
    .def(py::init<double, std::string const&>())
    .def(py::init<std::vector<std::string> const&>())
    .def(py::init<std::vector<std::string> const&, std::string const&>())
    .def(py::init<std::vector<int> const&>())
    .def(py::init<std::vector<int> const&, std::string const&>())
    .def(py::init<std::vector<double> const&>())
    .def(py::init<std::vector<double> const&, std::string const&>())
    .def("set_label", &Lsexpr::set_label)
    .def("is_labeled", &Lsexpr::is_labeled)
    .def("get_label", &Lsexpr::get_label)
    .def("is_atomic", &Lsexpr::is_atomic)
    .def("get_string", &Lsexpr::get_string)
    .def("get_int", &Lsexpr::get_int)
    .def("get_double", &Lsexpr::get_double)
    .def("get_string_vector", &Lsexpr::get_string_vector)
    .def("get_int_vector", &Lsexpr::get_int_vector)
    .def("get_double_vector", &Lsexpr::get_double_vector);

  m.def("read_lsexpr_file", read_lsexpr_file);
  m.def("write_lsexpr_file", write_lsexpr_file);
}
