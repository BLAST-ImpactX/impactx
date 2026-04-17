/* Copyright 2021-2023 The ImpactX Community
 *
 * Authors: Axel Huebl, Eric G. Stern
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <elements/All.H>
#include <elements/mixin/lineartransport.H>
#include <elements/transformation/Insert.H>
#include <particles/CovarianceMatrix.H>
#include <particles/ReferenceParticle.H>

#include <string>
#include <type_traits>
#include <utility>
#include <variant>

namespace py = pybind11;
using namespace impactx;


void init_lattice(py::module& me)
{
    using elements::KnownElements;

    // all-element type list
    using KnownElementsList = std::list<KnownElements>;
    py::class_<KnownElementsList> kel(me, "KnownElementsList");
    kel
        .def(py::init<>())
        .def(py::init<KnownElements>())
        .def(py::init([](py::list const & l){
            KnownElementsList v;
            for (auto const & handle : l)
                v.push_back(handle.cast<KnownElements>());
            return v;  // return by value
        }))

        .def("append", [](KnownElementsList &v, KnownElements el) { v.emplace_back(std::move(el)); },
             "Add a single element to the list.")

        .def("extend",
             [](KnownElementsList &v, KnownElementsList const & l) {
                 for (auto const & el : l)
                     v.push_back(el);
                 return v;
             },
             "Add a list of elements to the list.")
        .def("extend",
             [](KnownElementsList &v, py::list const & l) {
                 for (auto const & handle : l)
                 {
                     auto el = handle.cast<KnownElements>();
                     v.push_back(el);
                 }
                 return v;
             },
             "Add a list of elements to the list."
        )

        .def("size", &KnownElementsList::size)
        .def("clear", &KnownElementsList::clear,
             "Clear the list to become empty.")
        .def("is_empty", &KnownElementsList::empty)
        .def("pop_back", &KnownElementsList::pop_back,
             "Return and remove the last element of the list.")
        .def("__len__", [](const KnownElementsList &v) { return v.size(); },
             "The length of the list.")
        .def("__iter__", [](KnownElementsList &v) {
            return py::make_iterator(v.begin(), v.end());
        }, py::keep_alive<0, 1>())  // Keep list alive while iterator is used
        .def("__getitem__", [](KnownElementsList &v, size_t index) -> elements::KnownElements& {
            if (index >= v.size()) {
                throw std::out_of_range("Index out of range");
            }
            auto it = std::next(v.begin(), index);
            return *it;  // return by reference
        }, py::return_value_policy::reference_internal)

        .def(
            "transfer_map",
            [](
                KnownElementsList &v,
                RefPart ref, // note: intentional copy
                std::string order,
                bool fallback_identity_map
            )
            {
                if (order != "linear") {
                    throw std::runtime_error("So far, only the calculation of linear transfer maps are supported in this function.");
                }
                Map6x6 linear_transfer_map = Map6x6::Identity();
                for (auto & el_v : v)
                {
                    // advance reference particle
                    std::visit([&ref](auto && el) {
                        el(ref);
                    }, el_v);

                    // extract element transport map, handle fallback
                    Map6x6 element_transport_map = Map6x6::Identity();
                    int element_nslice;
                    std::visit([&ref, &fallback_identity_map, &element_transport_map, &element_nslice](auto const & el) {
                        using Element = std::decay_t<decltype(el)>;
                        std::string not_impl_msg = "Undefined transfer map in lattice for element ";
                        if (el.has_name()) not_impl_msg += el.name() + " ";
                        not_impl_msg += std::string("of type ") + Element::type;

                        if constexpr (std::is_base_of_v<elements::mixin::LinearTransport<Element>, Element>) {
                            try {
                                element_transport_map = el.transport_map(ref);
                                element_nslice = el.nslice();
                            } catch (std::exception const &) {
                                if (!fallback_identity_map) {
                                    throw std::runtime_error(not_impl_msg);
                                }
                            }
                        } else {
                            if (!fallback_identity_map) {
                                throw std::runtime_error(not_impl_msg);
                            }
                        }
                    }, el_v);

                    // advance linear transfer map
                    for (int n = 0; n < element_nslice; n++) {
                        linear_transfer_map = element_transport_map * linear_transfer_map;
                    }
                }
                return linear_transfer_map;
            },
            py::arg("ref"),
            py::arg("order") = "linear",
            py::arg("fallback_identity_map") = false,
            "Calculate the transfer map of the elements in the list."
            )
    ;


    // lattice transformations
    py::module_ met = me.def_submodule(
        "transformation",
        "Transform and modify lattices"
    );

    met.def(
        "insert_element_every_ds",
        &impactx::elements::transformation::insert_element_every_ds,
        py::arg("list"),
        py::arg("ds"),
        py::arg("element"),
        "Insert an element every s into an element list"
    );
}
