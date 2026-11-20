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

#include <diagnostics/InSitu.H>
#include <particles/ImpactXParticleContainer.H>

#include <AMReX_Config.H>

#include <pybind11/numpy.h>

#ifdef AMREX_USE_CONDUIT
#   include <conduit/conduit.hpp>
#endif

#include <string>
#include <utility>

namespace py = pybind11;
using namespace impactx;


namespace
{
    /** JSON text of a Python object: strings are passed through */
    std::string
    to_json (py::object const & obj)
    {
        if (py::isinstance<py::str>(obj))
        {
            return obj.cast<std::string>();
        }
        return py::module_::import("json").attr("dumps")(obj).cast<std::string>();
    }

#ifdef AMREX_USE_CONDUIT
    /** Copy a Python object into a Conduit node, without JSON text in between
     *
     * Dicts become children, strings and numbers become values, and numpy arrays and
     * sequences of numbers, e.g., lists, become arrays.
     *
     * @param[in] obj Python object
     * @param[out] node Conduit node, empty on input
     * @throws py::type_error for other objects, e.g., booleans or arrays of other types
     */
    void
    to_node (py::handle obj, conduit::Node & node)
    {
        if (py::isinstance<py::dict>(obj))
        {
            for (auto const & item : py::reinterpret_borrow<py::dict>(obj))
            {
                if (!py::isinstance<py::str>(item.first))
                {
                    throw py::type_error("Ascent: the keys of a Conduit node are strings.");
                }
                to_node(item.second, node[item.first.cast<std::string>()]);
            }
            return;
        }
        if (py::isinstance<py::str>(obj))
        {
            node.set(obj.cast<std::string>());
            return;
        }
        if (py::isinstance<py::bool_>(obj))
        {
            throw py::type_error("Ascent: a Conduit node has no booleans.");
        }
        if (py::isinstance<py::int_>(obj))
        {
            node.set(obj.cast<conduit::int64>());
            return;
        }
        if (py::isinstance<py::float_>(obj))
        {
            node.set(obj.cast<conduit::float64>());
            return;
        }
        // numpy arrays and sequences of numbers, contiguous
        py::array const array = py::array::ensure(obj, py::array::c_style);
        if (!array)
        {
            throw py::type_error(
                "Ascent: cannot convert " + py::repr(obj).cast<std::string>().substr(0, 80)
                + " to a Conduit node.");
        }
        auto const num = static_cast<conduit::index_t>(array.size());
        py::dtype const dtype = array.dtype();
        if (dtype.kind() == 'f' && dtype.itemsize() == 8)
        {
            node.set_float64_ptr(static_cast<conduit::float64 const *>(array.data()), num);
        }
        else if (dtype.kind() == 'f' && dtype.itemsize() == 4)
        {
            node.set_float32_ptr(static_cast<conduit::float32 const *>(array.data()), num);
        }
        else if (dtype.kind() == 'i' && dtype.itemsize() == 8)
        {
            node.set_int64_ptr(static_cast<conduit::int64 const *>(array.data()), num);
        }
        else if (dtype.kind() == 'i' && dtype.itemsize() == 4)
        {
            node.set_int32_ptr(static_cast<conduit::int32 const *>(array.data()), num);
        }
        else
        {
            throw py::type_error(
                "Ascent: arrays of " + dtype.attr("name").cast<std::string>()
                + " are not supported, use float64, float32, int64 or int32.");
        }
    }
#endif
} // namespace

void init_insitu (py::module & m)
{
    py::module_ mi = m.def_submodule(
        "insitu",
        "In situ visualization of the particle beam"
    );

    using diagnostics::insitu::Ascent;

    // Ascent runs on the calling thread and keeps the GIL: Ascent Python
    // filters and extracts run in the same interpreter.
    py::class_<Ascent>(mi, "Ascent")
        .def(py::init(
                [](py::object const & options) { return Ascent(to_json(options)); }
             ),
             py::arg("options") = py::dict(),
             "Open Ascent.\n\n"
             "This is MPI-collective.\n\n"
             "Parameters\n"
             "----------\n"
             "options : dict or str, optional\n"
             "    Ascent options, as a dict or JSON text, e.g., {\"default_dir\": \"path\"}.\n"
             "    ImpactX sets the MPI communicator. Unless set here, errors raise exceptions,\n"
             "    only fields used by the actions are converted (\"field_filtering\": \"true\"),\n"
             "    the actions passed to execute run unchanged (\"actions_file\": \"\"), and\n"
             "    query results are not saved to a session file at exit (\"default_dir\" and\n"
             "    \"session_name\" are empty: images without a directory in their name are\n"
             "    written to the current directory).\n\n"
             "Raises\n"
             "------\n"
             "RuntimeError\n"
             "    If ImpactX is built without Ascent."
        )
        .def("execute",
             [](Ascent & self, py::object const & pc, py::object const & actions,
                int step, int period, py::object const & geometry)
             {
                 ImpactXParticleContainer const * beam =
                     pc.is_none() ? nullptr : &pc.cast<ImpactXParticleContainer const &>();
#ifdef AMREX_USE_CONDUIT
                 // copied into a Conduit node: large geometries skip JSON
                 conduit::Node node;
                 if (py::isinstance<py::str>(geometry))
                 {
                     node.parse(geometry.cast<std::string>(), "json");
                 }
                 else if (!geometry.is_none())
                 {
                     to_node(geometry, node);
                 }
                 std::string const info_json = self.execute(
                     beam, to_json(actions), step, period, geometry.is_none() ? nullptr : &node);
#else
                 std::string const geometry_json = geometry.is_none() ? "" : to_json(geometry);
                 std::string const info_json =
                     self.execute(beam, to_json(actions), step, period, geometry_json);
#endif
                 return py::module_::import("json").attr("loads")(info_json);
             },
             py::arg("pc"), py::arg("actions"), py::arg("step") = 0, py::arg("period") = 0,
             py::arg("geometry") = py::none(),
             "Publish the particle beam and execute Ascent actions.\n\n"
             "The beam is published without copies, in its coordinates relative to the\n"
             "reference particle at fixed s: the coordinates are position_x, position_y and\n"
             "position_t and every particle component is a field, e.g., momentum_t.\n"
             "Images are written by MPI rank 0.\n\n"
             "This is MPI-collective.\n\n"
             "Parameters\n"
             "----------\n"
             "pc : ImpactXParticleContainer or None\n"
             "    The particle beam, at fixed s. None publishes only the geometry.\n"
             "actions : list or str\n"
             "    Ascent actions, as a list of dicts or JSON text.\n"
             "step : int, optional\n"
             "    Global step, the Ascent cycle.\n"
             "period : int, optional\n"
             "    Current period (turn or cycle).\n"
             "geometry : dict or str, optional\n"
             "    A domain of a Conduit Mesh Blueprint, e.g., of the lattice, that MPI rank 0\n"
             "    publishes next to the beam: a dict of dicts, strings, numbers and arrays\n"
             "    (numpy arrays of float64, float32, int64 or int32, or lists of numbers), or\n"
             "    JSON text. Its coordinate sets, topologies and fields need names that\n"
             "    differ from the ones of the beam. Only the one of MPI rank 0 is used.\n\n"
             "Returns\n"
             "-------\n"
             "dict\n"
             "    The Ascent info, e.g., with the written images under \"images\".\n\n"
             "Raises\n"
             "------\n"
             "ValueError\n"
             "    Without a beam and without a geometry on MPI rank 0."
        )
        .def("close", &Ascent::close,
             "Close Ascent and release the published data. This is MPI-collective."
        )
        .def_property_readonly("is_open", &Ascent::is_open,
             "Is Ascent open?"
        )
        .def("__enter__", [](Ascent & self) -> Ascent & { return self; },
             py::return_value_policy::reference_internal)
        .def("__exit__",
             [](Ascent & self, py::object const &, py::object const &, py::object const &)
             {
                 self.close();
             })
    ;
}
