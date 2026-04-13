// synmadx Python bindings
// Adapted from synergia/lattice/lattice_pywrap.cc and
// synergia/foundation/foundation_pywrap.cc

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "synergia/foundation/four_momentum.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/foundation/reference_particle.h"
#include "synergia/lattice/lattice.h"
#include "synergia/lattice/lattice_element.h"
#include "synergia/lattice/lattice_element_slice.h"
#include "synergia/lattice/madx_reader.h"

namespace py = pybind11;
using namespace py::literals;

// Custom type caster: convert std::list<Lattice_element> to Python tuple
// (from synergia/utils/container_conversions.h)
namespace pybind11 { namespace detail {
    template <typename Type, typename Value>
    struct py_tuple_caster {
        using value_conv = make_caster<Value>;

        bool load(handle src, bool convert) { return false; }

        template <typename T>
        static handle cast(T&& src, return_value_policy policy, handle parent)
        {
            if (!std::is_lvalue_reference<T>::value)
                policy = return_value_policy_override<Value>::policy(policy);

            tuple t(src.size());
            size_t index = 0;
            for (auto&& value : src) {
                auto value_ = reinterpret_steal<object>(
                    value_conv::cast(forward_like<T>(value), policy, parent));
                if (!value_) return handle();
                PyTuple_SET_ITEM(t.ptr(), (ssize_t)index++, value_.release().ptr());
            }
            return t.release();
        }

        PYBIND11_TYPE_CASTER(Type, _("Tuple[") + value_conv::name + _("]"));
    };
}} // namespace pybind11::detail

template <>
struct pybind11::detail::type_caster<std::list<Lattice_element>>
    : pybind11::detail::py_tuple_caster<std::list<Lattice_element>, Lattice_element> {};


PYBIND11_MODULE(synmadx_pybind, m)
{
    m.doc() = "synmadx: standalone MAD-X lattice parser (from Synergia)";

    // ---- physical constants submodule ----
    auto pc = m.def_submodule("pconstants", "PDG physical constants");
    pc.attr("pdg_year") = pconstants::pdg_year;
    pc.attr("mp") = pconstants::mp;
    pc.attr("proton_mass") = pconstants::proton_mass;
    pc.attr("me") = pconstants::me;
    pc.attr("electron_mass") = pconstants::electron_mass;
    pc.attr("mmu") = pconstants::mmu;
    pc.attr("muon_mass") = pconstants::muon_mass;
    pc.attr("e") = pconstants::e;
    pc.attr("c") = pconstants::c;
    pc.attr("mu0") = pconstants::mu0;
    pc.attr("epsilon0") = pconstants::epsilon0;
    pc.attr("re") = pconstants::re;
    pc.attr("rp") = pconstants::rp;
    pc.attr("rmu") = pconstants::rmu;
    pc.attr("proton_charge") = pconstants::proton_charge;
    pc.attr("antiproton_charge") = pconstants::antiproton_charge;
    pc.attr("electron_charge") = pconstants::electron_charge;
    pc.attr("positron_charge") = pconstants::positron_charge;
    pc.attr("muon_charge") = pconstants::muon_charge;
    pc.attr("antimuon_charge") = pconstants::antimuon_charge;
    pc.attr("kg_to_GeV") = pconstants::kg_to_GeV;

    // ---- enums ----
    py::enum_<element_format>(m, "element_format", py::arithmetic())
        .value("mad8", element_format::mad8)
        .value("madx", element_format::madx)
        ;

    py::enum_<element_type>(m, "element_type", py::arithmetic())
        .value("generic",     element_type::generic)
        .value("drift",       element_type::drift)
        .value("rbend",       element_type::rbend)
        .value("sbend",       element_type::sbend)
        .value("quadrupole",  element_type::quadrupole)
        .value("multipole",   element_type::multipole)
        .value("rfcavity",    element_type::rfcavity)
        .value("hkicker",     element_type::hkicker)
        .value("vkicker",     element_type::vkicker)
        .value("kicker",      element_type::kicker)
        .value("monitor",     element_type::monitor)
        .value("hmonitor",    element_type::hmonitor)
        .value("vmonitor",    element_type::vmonitor)
        .value("sextupole",   element_type::sextupole)
        .value("octupole",    element_type::octupole)
        .value("marker",      element_type::marker)
        .value("instrument",  element_type::instrument)
        .value("rcollimator", element_type::rcollimator)
        .value("nllens",      element_type::nllens)
        .value("solenoid",    element_type::solenoid)
        .value("elens",       element_type::elens)
        .value("foil",        element_type::foil)
        .value("dipedge",     element_type::dipedge)
        .value("matrix",      element_type::matrix)
        ;

    py::enum_<marker_type>(m, "marker_type", py::arithmetic())
        .value("h_tunes_corrector", marker_type::h_tunes_corrector)
        .value("v_tunes_corrector", marker_type::v_tunes_corrector)
        .value("h_chrom_corrector", marker_type::h_chrom_corrector)
        .value("v_chrom_corrector", marker_type::v_chrom_corrector)
        ;

    // ---- latt_func_t ----
    py::class_<latt_func_t::lf_val_t>(m, "lf_val_t")
        .def_readwrite("hor", &latt_func_t::lf_val_t::hor)
        .def_readwrite("ver", &latt_func_t::lf_val_t::ver)
        ;

    py::class_<latt_func_t>(m, "latt_func_t")
        .def_readwrite("arcLength",  &latt_func_t::arcLength)
        .def_readwrite("dispersion", &latt_func_t::dispersion)
        .def_readwrite("dPrime",     &latt_func_t::dPrime)
        .def_readwrite("beta",       &latt_func_t::beta)
        .def_readwrite("alpha",      &latt_func_t::alpha)
        .def_readwrite("psi",        &latt_func_t::psi)
        ;

    // ---- Four_momentum ----
    py::class_<Four_momentum>(m, "Four_momentum")
        .def(py::init<double>(),
             "Construct a Four_momentum in the rest frame."
             "\n\tparam: mass in GeV/c^2",
             "mass"_a)

        .def(py::init<double, double>(),
             "Construct a Four_momentum with the given total energy."
             "\n\tparam: mass in GeV/c^2"
             "\n\tparam: total_energy in GeV",
             "mass"_a,
             "total_energy"_a)

        .def("set_total_energy",
             &Four_momentum::set_total_energy,
             "Set the total energy (in GeV)",
             "total_energy"_a)

        .def("set_kinetic_energy",
             &Four_momentum::set_kinetic_energy,
             "Set the kinetic energy (in GeV)",
             "kinetic_energy"_a)

        .def("set_momentum",
             &Four_momentum::set_momentum,
             "Set the momentum (in GeV/c)",
             "momentum"_a)

        .def("set_gamma",
             &Four_momentum::set_gamma,
             "Set the relativistic gamma factor (unitless)",
             "gamma"_a)

        .def("set_beta",
             &Four_momentum::set_beta,
             "Set the relativistic beta factor (unitless)",
             "beta"_a)

        .def("get_mass", &Four_momentum::get_mass, "Get the mass in GeV/c^2")

        .def("get_total_energy",
             &Four_momentum::get_total_energy,
             "Get the total energy in GeV")

        .def("get_kinetic_energy",
             &Four_momentum::get_kinetic_energy,
             "Get the kinetic energy in GeV")

        .def("get_momentum",
             &Four_momentum::get_momentum,
             "Get the momentum in GeV/c")

        .def("get_gamma",
             &Four_momentum::get_gamma,
             "Get the relativistic gamma factor")

        .def("get_beta",
             &Four_momentum::get_beta,
             "Get the relativistic beta factor")

        .def("equal",
             &Four_momentum::equal,
             "Check equality to the given tolerance",
             "four_momentum"_a,
             "tolerance"_a)
        ;

    // ---- Reference_particle ----
    py::class_<Reference_particle>(m, "Reference_particle")
        .def(py::init<int, double, double>(),
             "Construct a Reference_particle with a given mass and total energy."
             "\n\tparam: charge in units of e"
             "\n\tparam: mass in GeV/c^2"
             "\n\tparam: total_energy in GeV in the lab frame",
             "charge"_a, "mass"_a, "total_energy"_a)

        .def(py::init<int, Four_momentum const&>(),
             "Construct a Reference_particle with a given four momentum",
             "charge"_a, "four_momentum"_a)

        .def(py::init<int, Four_momentum const&, std::array<double, 6> const&>(),
             "Construct a Reference_particle with a given four momentum and state",
             "charge"_a, "four_momentum"_a, "state"_a)

        .def("set_four_momentum",
             &Reference_particle::set_four_momentum,
             "Set the four momentum",
             "four_momentum"_a)

        .def("set_state",
             (void(Reference_particle::*)(std::array<double, 6> const&))
                 &Reference_particle::set_state,
             "Set the state vector in the reference frame",
             "state"_a)

        .def("set_state",
             (void(Reference_particle::*)(double, double, double, double, double, double))
                 &Reference_particle::set_state,
             "Set the state vector in the reference frame",
             "x"_a, "xp"_a, "y"_a, "yp"_a, "cdt"_a, "dpop"_a)

        .def("set_total_energy",
             &Reference_particle::set_total_energy,
             "Set the total energy in GeV in the lab frame",
             "total_energy"_a)

        .def("increment_trajectory",
             &Reference_particle::increment_trajectory,
             "Increment the trajectory length in m",
             "length"_a)

        .def("start_repetition",
             &Reference_particle::start_repetition,
             "Start a new repetition")

        .def("set_trajectory",
             &Reference_particle::set_trajectory,
             "Manually set the trajectory parameters",
             "repetition"_a, "repetition_length"_a, "s"_a)

        .def("set_bunch_abs_time",
             &Reference_particle::set_bunch_abs_time,
             "Set the bunch absolute time in [s]",
             "t"_a)

        .def("increment_bunch_abs_time",
             &Reference_particle::increment_bunch_abs_time,
             "Increment the bunch absolute time dt in [s]",
             "dt"_a)

        .def("set_bunch_abs_offset",
             &Reference_particle::set_bunch_abs_offset,
             "Set the bunch absolute offset time in [s]",
             "toffset"_a)

        .def("get_charge",
             &Reference_particle::get_charge,
             "Return the charge in units of e")

        .def("get_mass",
             &Reference_particle::get_mass,
             "Return the mass in units of GeV/c^2")

        .def("get_four_momentum",
             &Reference_particle::get_four_momentum,
             "Get the four momentum in the lab frame")

        .def("get_state",
             &Reference_particle::get_state,
             "Get the six-dimensional state vector in the reference frame")

        .def("get_beta",
             &Reference_particle::get_beta,
             "Get the relativistic beta in the lab frame")

        .def("get_gamma",
             &Reference_particle::get_gamma,
             "Get the relativistic gamma in the lab frame")

        .def("get_momentum",
             &Reference_particle::get_momentum,
             "Get the momentum in GeV/c in the lab frame")

        .def("get_total_energy",
             &Reference_particle::get_total_energy,
             "Get the total energy in GeV in the lab frame")

        .def("get_s",
             &Reference_particle::get_s,
             "Get the total path length in m")

        .def("get_s_n",
             &Reference_particle::get_s_n,
             "Get the distance traveled in m since the beginning of the current repetition")

        .def("get_repetition",
             &Reference_particle::get_repetition,
             "Get the number of repetitions")

        .def("get_repetition_length",
             &Reference_particle::get_repetition_length,
             "Get the repetition length in m")

        .def("get_bunch_abs_time",
             &Reference_particle::get_bunch_abs_time,
             "Get the bunch absolute time [s]")

        .def("get_bunch_abs_offset",
             &Reference_particle::get_bunch_abs_offset,
             "Get the bunch absolute offset [s]")

        .def("equal",
             &Reference_particle::equal,
             "Check equality to the given tolerance",
             "reference_particle"_a, "tolerance"_a)
        ;

    // ---- Lattice_element ----
    py::class_<Lattice_element>(m, "Lattice_element")
        .def(py::init<>(),
             "Construct a generic lattice element")

        .def(py::init<std::string const&, std::string const&, element_format>(),
             "Construct a lattice element with type, name, and format",
             "type"_a, "name"_a, "format"_a = element_format::madx)

        .def("get_type",
             &Lattice_element::get_type,
             "Returns lattice element type")

        .def("get_type_name",
             &Lattice_element::get_type_name,
             "Returns lattice element type string")

        .def("get_name",
             &Lattice_element::get_name,
             "Returns lattice element name")

        .def("get_format",
             &Lattice_element::get_format,
             "Returns lattice element format (madx or mad8)")

        .def("get_length",
             &Lattice_element::get_length)

        .def("get_bend_angle",
             &Lattice_element::get_bend_angle)

        .def("get_ancestors",
             &Lattice_element::get_ancestors)

        .def("add_ancestor",
             &Lattice_element::add_ancestor,
             "ancestor"_a)

        .def("copy_attributes_from",
             &Lattice_element::copy_attributes_from,
             "element"_a)

        // double attributes
        .def("has_double_attribute",
             &Lattice_element::has_double_attribute,
             "Check for existence of the named double attribute",
             "name"_a)

        .def("get_double_attribute",
             (double (Lattice_element::*)(std::string const&) const)
                 &Lattice_element::get_double_attribute,
             "Get the value of the named double attribute",
             "name"_a)

        .def("get_double_attribute",
             (double (Lattice_element::*)(std::string const&, double) const)
                 &Lattice_element::get_double_attribute,
             "Get the value of the named double attribute, or return the default",
             "name"_a, "default"_a)

        .def("set_double_attribute",
             (void (Lattice_element::*)(std::string const&, double, bool))
                 &Lattice_element::set_double_attribute,
             "Set the value of the named double attribute",
             "name"_a, "value"_a, "increment_revision"_a = true)

        .def("remove_double_attribute",
             (void (Lattice_element::*)(std::string const&))
                 &Lattice_element::remove_double_attribute,
             "Remove the named double attribute",
             "name"_a)

        .def("set_double_attribute",
             (void (Lattice_element::*)(std::string const&, std::string const&, bool))
                 &Lattice_element::set_double_attribute,
             "Set the value (as an expression) of the named double attribute",
             "name"_a, "value"_a, "increment_revision"_a = true)

        .def("get_double_attribute_names",
             &Lattice_element::get_double_attribute_names)

        // string attributes
        .def("has_string_attribute",
             &Lattice_element::has_string_attribute,
             "Check for existence of the named string attribute",
             "name"_a)

        .def("remove_string_attribute",
             &Lattice_element::remove_string_attribute,
             "Remove the named string attribute",
             "name"_a)

        .def("get_string_attribute",
             (std::string const& (Lattice_element::*)(std::string const&) const)
                 &Lattice_element::get_string_attribute,
             "Get the value of the named string attribute",
             "name"_a)

        .def("get_string_attribute",
             (std::string const& (Lattice_element::*)(std::string const&, std::string const&) const)
                 &Lattice_element::get_string_attribute,
             "Get the value of the named string attribute, or return the default",
             "name"_a, "default"_a)

        .def("set_string_attribute",
             &Lattice_element::set_string_attribute,
             "Set the value of the named string attribute",
             "name"_a, "value"_a, "increment_revision"_a = true)

        // vector attributes
        .def("has_vector_attribute",
             &Lattice_element::has_vector_attribute,
             "Check for existence of the named vector attribute",
             "name"_a)

        .def("remove_vector_attribute",
             &Lattice_element::remove_vector_attribute,
             "Remove the named vector attribute",
             "name"_a)

        .def("get_vector_attribute",
             (std::vector<double> (Lattice_element::*)(std::string const&) const)
                 &Lattice_element::get_vector_attribute,
             "Get the value of the named vector attribute",
             "name"_a)

        .def("get_vector_attribute",
             (std::vector<double> (Lattice_element::*)(std::string const&, std::vector<double> const&) const)
                 &Lattice_element::get_vector_attribute,
             "Get the value of the named vector attribute, or return the default",
             "name"_a, "default"_a)

        .def("set_vector_attribute",
             &Lattice_element::set_vector_attribute,
             "Set the value of the named vector attribute",
             "name"_a, "value"_a, "increment_revision"_a = true)

        .def("get_string_attributes",
             &Lattice_element::get_string_attributes)

        .def("get_double_attributes",
             &Lattice_element::get_double_attributes)

        // markers
        .def("set_marker",
             &Lattice_element::set_marker,
             "marker"_a)

        .def("reset_marker",
             &Lattice_element::reset_marker,
             "marker"_a)

        .def("reset_markers",
             &Lattice_element::reset_markers)

        .def("has_marker",
             &Lattice_element::has_marker,
             "marker"_a)

        // deposited charge
        .def("get_deposited_charge",
             &Lattice_element::get_deposited_charge,
             "bunch"_a = 0, "train"_a = 0)

        .def("set_deposited_charge",
             &Lattice_element::set_deposited_charge,
             "charge"_a, "bunch"_a = 0, "train"_a = 0)

        .def("deposit_charge",
             &Lattice_element::deposit_charge,
             "charge"_a, "bunch"_a = 0, "train"_a = 0)

        .def_static("get_all_type_names",
             &Lattice_element::get_all_type_names)

        .def_readwrite("lf",
             &Lattice_element::lf)

        .def("print_",
             &Lattice_element::print,
             "Print the lattice element")

        .def("as_madx",
             &Lattice_element::as_madx,
             "Get MAD-X representation",
             "sanitize"_a = false)

        .def("__repr__",
             &Lattice_element::as_string)
        ;

    // ---- Lattice_element_slice ----
    py::class_<Lattice_element_slice>(m, "Lattice_element_slice")
        .def(py::init<Lattice_element const&>(),
             "Construct a lattice element slice from an entire element",
             "element"_a)

        .def(py::init<Lattice_element const&, double, double>(),
             "Construct a lattice element slice with left and right from an element",
             "element"_a, "left"_a, "right"_a)

        .def("is_whole",
             &Lattice_element_slice::is_whole,
             "Is a whole element")

        .def("has_left_edge",
             &Lattice_element_slice::has_left_edge,
             "Does this slice include the left edge of the element")

        .def("has_right_edge",
             &Lattice_element_slice::has_right_edge,
             "Does this slice include the right edge of the element")

        .def("get_left",
             &Lattice_element_slice::get_left,
             "Get the start position of the slice")

        .def("get_right",
             &Lattice_element_slice::get_right,
             "Get the end position of the slice")

        .def("print_",
             &Lattice_element_slice::print,
             "Print the lattice element slice")

        .def("get_lattice_element",
             &Lattice_element_slice::get_lattice_element,
             "Get the lattice element corresponding to this slice")

        .def("__repr__",
             &Lattice_element_slice::as_string)
        ;

    // ---- Lattice_tree ----
    using lattice_tree_set_element_attribute_1 =
        void (Lattice_tree::*)(std::string const&, std::string const&, double);
    using lattice_tree_set_element_attribute_2 =
        void (Lattice_tree::*)(std::string const&, std::string const&, std::string const&);

    py::class_<Lattice_tree>(m, "Lattice_tree")
        .def("set_variable",
             (void (Lattice_tree::*)(std::string const&, double))
                 &Lattice_tree::set_variable,
             "name"_a, "val"_a)

        .def("set_variable",
             (void (Lattice_tree::*)(std::string const&, std::string const&))
                 &Lattice_tree::set_variable,
             "name"_a, "val"_a)

        .def("set_element_attribute",
             (lattice_tree_set_element_attribute_1)
                 &Lattice_tree::set_element_attribute,
             "label"_a, "attr"_a, "val"_a)

        .def("set_element_attribute",
             (lattice_tree_set_element_attribute_2)
                 &Lattice_tree::set_element_attribute,
             "label"_a, "attr"_a, "val"_a)

        .def("print",
             &Lattice_tree::print)
        ;

    // ---- Lattice ----
    py::class_<Lattice>(m, "Lattice")
        .def(py::init<>(),
             "Construct an unnamed empty lattice")

        .def(py::init<std::string const&>(),
             "Construct an empty lattice",
             "name"_a)

        .def(py::init<Lattice const&>(),
             "Construct a copy of a lattice")

        .def(py::init<std::string const&, Reference_particle const&>(),
             "Construct an empty lattice with given reference particle",
             "name"_a, "refpart"_a)

        .def(py::init<std::string const&, Lattice_tree const&>(),
             "Construct a dynamic lattice with given tree",
             "name"_a, "tree"_a)

        .def(py::init<Lsexpr const&>(),
             "Construct from the Lsexpr representation",
             "lsexpr"_a)

        .def("get_name",
             &Lattice::get_name,
             "Get the lattice name")

        .def("append",
             &Lattice::append,
             "Append a lattice element to the lattice",
             "element"_a)

        .def("get_reference_particle",
             (Reference_particle& (Lattice::*)()) &Lattice::get_reference_particle,
             py::return_value_policy::reference_internal,
             "Get the lattice reference particle")

        .def("set_reference_particle",
             &Lattice::set_reference_particle)

        .def("get_lattice_energy",
             &Lattice::get_lattice_energy,
             "Get the design energy for the lattice")

        .def("set_lattice_energy",
             &Lattice::set_lattice_energy,
             "Set the design energy of the lattice",
             "energy"_a)

        .def("get_elements",
             (std::list<Lattice_element>& (Lattice::*)())
                 &Lattice::get_elements,
             py::return_value_policy::reference_internal,
             "Get the list of all lattice elements")

        .def("get_lattice_tree",
             (Lattice_tree& (Lattice::*)())
                 &Lattice::get_lattice_tree,
             py::return_value_policy::reference_internal,
             "Get the Lattice_tree object for variable manipulation")

        .def("set_lattice_tree",
             &Lattice::set_lattice_tree,
             "Set the Lattice_tree object for variable manipulation",
             "lattice_tree"_a)

        .def("set_variable",
             (void (Lattice::*)(std::string const&, double))
                 &Lattice::set_variable,
             "Set a double variable in the variable table",
             "name"_a, "val"_a)

        .def("set_variable",
             (void (Lattice::*)(std::string const&, std::string const&))
                 &Lattice::set_variable,
             "Set a double variable expression in the variable table",
             "name"_a, "val"_a)

        .def("reset_all_markers",
             &Lattice::reset_all_markers,
             "Clear the h/v tunes and chromaticity corrector markers for all elements")

        .def("get_length",
             &Lattice::get_length,
             "Get the combined length of all elements in the lattice")

        .def("get_total_angle",
             &Lattice::get_total_angle,
             "Get the total angle in radians subtended by all elements in the lattice")

        .def("get_elements_const",
             (std::list<Lattice_element> const& (Lattice::*)() const) &Lattice::get_elements,
             py::return_value_policy::reference_internal,
             "Get the list of all lattice elements (const)")

        .def("set_all_double_attribute",
             &Lattice::set_all_double_attribute,
             "Set the value of the named double attribute on all elements",
             "name"_a, "value"_a, "increment_revision"_a = true)

        .def("set_all_string_attribute",
             &Lattice::set_all_string_attribute,
             "Set the value of the named string attribute on all elements",
             "name"_a, "value"_a, "increment_revision"_a = true)

        .def("export_madx_file",
             &Lattice::export_madx_file,
             "Export the lattice to a MadX file",
             "filename"_a, "sanitize"_a = false)

        .def_static("import_madx_file",
             &Lattice::import_madx_file,
             "Import a line or a sequence from MadX file",
             "filename"_a, "line"_a)

        .def("__repr__", &Lattice::as_string)
        ;

    // ---- MadX_reader ----
    using madx_reader_get_lattice_1 =
        Lattice (MadX_reader::*)(std::string const&);
    using madx_reader_get_lattice_2 =
        Lattice (MadX_reader::*)(std::string const&, std::string const&);

    py::class_<MadX_reader>(m, "MadX_reader")
        .def(py::init<>(), "Construct a MadX_reader")

        .def("get_lattice",
             (madx_reader_get_lattice_1) &MadX_reader::get_lattice,
             "Get the named lattice from an already parsed lattice",
             "line_name"_a)

        .def("get_lattice",
             (madx_reader_get_lattice_2) &MadX_reader::get_lattice,
             "Parse and get the named lattice",
             "line_name"_a, "filename"_a)

        .def("get_dynamic_lattice",
             (madx_reader_get_lattice_1) &MadX_reader::get_dynamic_lattice,
             "Get the named dynamic lattice from an already parsed lattice",
             "line_name"_a)

        .def("get_dynamic_lattice",
             (madx_reader_get_lattice_2) &MadX_reader::get_dynamic_lattice,
             "Parse and get the named dynamic lattice",
             "line_name"_a, "filename"_a)

        .def("parse",
             &MadX_reader::parse,
             "Parse a MAD-X lattice string",
             "lattice"_a)

        .def("parse_file",
             &MadX_reader::parse_file,
             "Parse a MAD-X lattice file",
             "filename"_a)

        .def("get_line_names",
             &MadX_reader::get_line_names,
             "Get all line names")

        .def("get_sequence_names",
             &MadX_reader::get_sequence_names,
             "Get all sequence names")

        .def("get_all_names",
             &MadX_reader::get_all_names,
             "Get all element names")

        .def("get_double_variable",
             &MadX_reader::get_double_variable,
             "Get variable value as double",
             "name"_a)

        .def("get_string_variable",
             &MadX_reader::get_string_variable,
             "Get variable value as string",
             "name"_a)
        ;
}
