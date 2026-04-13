
#include <pybind11/pybind11.h>

#include "synergia/utils/commxx.h"
#include "synergia/utils/logger.h"
#include "synergia/utils/parallel_utils.h"
#include "synergia/utils/simple_timer.h"

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE(parallel_utils, m)
{
  py::enum_<LoggerV>(m, "LoggerV", py::arithmetic())
    .value("DEBUG", LoggerV::DEBUG)
    .value("DINFO", LoggerV::DINFO)

    .value("INFO", LoggerV::INFO)
    .value("INFO_OPN", LoggerV::INFO_OPN)
    .value("INFO_OPR", LoggerV::INFO_OPR)
    .value("INFO_STEP", LoggerV::INFO_STEP)
    .value("INFO_TURN", LoggerV::INFO_TURN)

    .value("WARNING", LoggerV::WARNING)
    .value("ERROR", LoggerV::ERROR);

  py::class_<Logger>(m, "Logger")
    .def(py::init<int, LoggerV, bool>(),
         "Construct a logger.",
         "rank"_a = 0,
         "verbosity"_a = LoggerV::DINFO,
         "log"_a = true)

    .def(py::init<int, std::string const&, LoggerV, bool>(),
         "Construct a logger.",
         "rank"_a = 0,
         "filename"_a = "logger",
         "verbosity"_a = LoggerV::DINFO,
         "log"_a = true)

    .def("write", &Logger::write, py::return_value_policy::reference_internal)
    .def("flush", &Logger::flush, py::return_value_policy::reference_internal);

  py::class_<Commxx, std::shared_ptr<Commxx>>(m, "Commxx")
    .def(py::init<>())
    .def("rank", &Commxx::rank)
    .def("size", &Commxx::size)
    .def("get_rank", &Commxx::rank)
    .def("get_size", &Commxx::size)
    .def("is_null", &Commxx::is_null)
    .def("is_root", &Commxx::is_root)

    .def_property_readonly_static("World",
                                  [](py::object) { return Commxx::World; })

    .def_property_readonly_static("Null",
                                  [](py::object) { return Commxx::Null; });

  m.def("simple_timer_print", &simple_timer_print, "logger"_a);
}
