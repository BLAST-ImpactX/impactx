#include "lsexpr.h"
#include "synergia/utils/container_conversions.h"
#include <boost/python.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(pylsexpr)
{
  class_<Lsexpr>("Lsexpr")
    .def(init<>())
    .def(init<std::string const&>())
    .def(init<std::string const&, std::string const&>())
    .def(init<int>())
    .def(init<int, std::string const&>())
    .def(init<double>())
    .def(init<double, std::string const&>())
    .def(init<std::vector<std::string> const&>())
    .def(init<std::vector<std::string> const&, std::string const&>())
    .def(init<std::vector<int> const&>())
    .def(init<std::vector<int> const&, std::string const&>())
    .def(init<std::vector<double> const&>())
    .def(init<std::vector<double> const&, std::string const&>())
    .def("set_label", &Lsexpr::set_label)
    .def("is_labeled", &Lsexpr::is_labeled)
    .def("get_label",
         &Lsexpr::get_label,
         return_value_policy<copy_const_reference>())
    .def("is_atomic", &Lsexpr::is_atomic)
    .def("get_string", &Lsexpr::get_string)
    .def("get_int", &Lsexpr::get_int)
    .def("get_double", &Lsexpr::get_double)
    .def("get_string_vector", &Lsexpr::get_string_vector)
    .def("get_int_vector", &Lsexpr::get_int_vector)
    .def("get_double_vector", &Lsexpr::get_double_vector);

  def("read_lsexpr_file", read_lsexpr_file);
  def("write_lsexpr_file", write_lsexpr_file);
}
