
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "synergia/simulation/checkpoint.h"
#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/split_operator_stepper.h"
#include "synergia/simulation/split_operator_stepper_elements.h"
#include "synergia/simulation/stepper.h"

#include "synergia/bunch/bunch.h"
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/bunch/diagnostics_bulk_track.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/bunch/diagnostics_loss.h"
#include "synergia/bunch/diagnostics_particles.h"
#include "synergia/bunch/diagnostics_py.h"
#include "synergia/bunch/diagnostics_worker.h"
#include "synergia/bunch/populate.h"

#include "synergia/bunch/diagnostics_py.h"

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE(simulation, m)
{
    m.def("checkpoint_save",
          &syn::checkpoint_save,
          "propagator"_a,
          "bunch_simulator"_a);

    m.def("checkpoint_load", &syn::checkpoint_load);

    // Propagator::Lattice_element_slices
    py::class_<Propagator::Lattice_element_slices>(m, "Lattice_element_slices")
        .def(
            "__iter__",
            [](Propagator::Lattice_element_slices& s) {
                return py::make_iterator(s.begin(), s.end());
            },
            py::keep_alive<0, 1>());

    using halt_data_turn_t = std::function<bool(Bunch_simulator&, void*)>;
    using halt_turn_t = std::function<bool(Bunch_simulator&)>;

    // Propagator
    py::class_<Propagator>(m, "Propagator")
        .def(py::init<Lattice const&, Stepper const&>())

        .def("propagate", &Propagator::propagate)

        .def("print_steps", &Propagator::print_steps)

        .def("get_lattice_element_slices",
             &Propagator::get_lattice_element_slices,
             "Returns immutable copy of lattice element slices")

        .def("get_lattice_elements",
             &Propagator::get_lattice_elements,
             "Returns immutable copy of lattice elements")

        .def("get_lattice",
             py::overload_cast<>(&Propagator::get_lattice, py::const_),
             py::return_value_policy::reference_internal,
             "Returns immutable reference to the lattice, but the attributes "
             "of the elements contained with may be modified with "
             "set_<>_attribute member functions")

        .def("set_checkpoint_period",
             &Propagator::set_checkpoint_period,
             "period"_a)

        .def("get_checkpoint_period", &Propagator::get_checkpoint_period)

        .def("set_final_checkpoint", &Propagator::set_final_checkpoint, "val"_a)

        .def("get_final_checkpoint", &Propagator::get_final_checkpoint)

        .def("reg_halt_data_turn",
             py::overload_cast<halt_turn_t>(&Propagator::reg_halt_data_turn))

        ;

    // chormaticities_t
    py::class_<chromaticities_t>(m, "chromaticities_t")
        .def_readonly("momentum_compaction",
                      &chromaticities_t::momentum_compaction)

        .def_readonly("horizontal_chromaticity",
                      &chromaticities_t::horizontal_chromaticity)

        .def_readonly("horizontal_chromaticity_prime",
                      &chromaticities_t::horizontal_chromaticity_prime)

        .def_readonly("vertical_chromaticity",
                      &chromaticities_t::vertical_chromaticity)

        .def_readonly("vertical_chromaticity_prime",
                      &chromaticities_t::vertical_chromaticity_prime)

        .def_readonly("slip_factor", &chromaticities_t::slip_factor)

        .def_readonly("slip_factor_prime",
                      &chromaticities_t::slip_factor_prime);

    // Lattice simulator -- only a namespace
    m.def_submodule("Lattice_simulator")
        .def("set_closed_orbit_tolerance",
             &Lattice_simulator::set_closed_orbit_tolerance,
             "tolerance"_a)

        .def("get_closed_orbit_tolerance",
             &Lattice_simulator::get_closed_orbit_tolerance)

        .def("tune_linear_lattice",
             &Lattice_simulator::tune_linear_lattice,
             "Tune linear lattice.",
             "lattice"_a)

        .def("tune_circular_lattice",
             &Lattice_simulator::tune_circular_lattice,
             "Tune a circular lattice.",
             "lattice"_a)

        .def(
            "tune_rfcavities", &Lattice_simulator::tune_rfcavities, "lattice"_a)

        .def("calculate_closed_orbit",
             &Lattice_simulator::calculate_closed_orbit,
             "lattice"_a,
             "dpp"_a = 0.0)

        .def("calculate_tune_and_cdt",
             &Lattice_simulator::calculate_tune_and_cdt,
             "lattice"_a,
             "dpp"_a = 0.0)

        .def("get_chromaticities",
             &Lattice_simulator::get_chromaticities,
             "lattice"_a,
             "dpp"_a = 1.0e-5)

        .def("calculate_normal_form_o1",
             &Lattice_simulator::calculate_normal_form<1>,
             "lattice"_a)

        .def("calculate_normal_form_o2",
             &Lattice_simulator::calculate_normal_form<2>,
             "lattice"_a)

        .def("calculate_normal_form_o3",
             &Lattice_simulator::calculate_normal_form<3>,
             "lattice"_a)

        .def("calculate_normal_form_o4",
             &Lattice_simulator::calculate_normal_form<4>,
             "lattice"_a)

        .def("calculate_normal_form_o5",
             &Lattice_simulator::calculate_normal_form<5>,
             "lattice"_a)

        .def("calculate_normal_form_o6",
             &Lattice_simulator::calculate_normal_form<6>,
             "lattice"_a)

        .def("calculate_normal_form_o7",
             &Lattice_simulator::calculate_normal_form<7>,
             "lattice"_a)

        .def("adjust_tunes",
             &Lattice_simulator::adjust_tunes,
             "lattice"_a,
             "horizontal_tune"_a,
             "vertical_tune"_a,
             "tolerance"_a = 1e-5)

        .def("adjust_chromaticities",
             &Lattice_simulator::adjust_chromaticities,
             "lattice"_a,
             "horizontal_chromaticity"_a,
             "vertical_chromaticity"_a,
             "tolerance"_a = 1e-4,
             "max_steps"_a = 6)

        .def("CourantSnyderLatticeFunctions",
             (void (*)(Lattice&))(
                 &Lattice_simulator::CourantSnyderLatticeFunctions),
             "lattice"_a)

        .def("calc_dispersions",
             (void (*)(Lattice&)) & Lattice_simulator::calc_dispersions,
             "lattice"_a)

        .def("get_bucket_length",
             &Lattice_simulator::get_bucket_length,
             "lattice"_a)

        .def("get_rf_frequency",
             &Lattice_simulator::get_rf_frequency,
             "lattice"_a)

        .def(
            "get_linear_one_turn_map",
            [](Lattice const& lattice, int order) {
                using namespace Lattice_simulator;

                karray2d_row map;

                if (order == 1)
                    map = get_one_turn_map<1>(lattice).jacobian();
                else if (order == 2)
                    map = get_one_turn_map<2>(lattice).jacobian();
                else if (order == 3)
                    map = get_one_turn_map<3>(lattice).jacobian();
                else if (order == 4)
                    map = get_one_turn_map<4>(lattice).jacobian();
                else if (order == 5)
                    map = get_one_turn_map<5>(lattice).jacobian();
                else if (order == 6)
                    map = get_one_turn_map<6>(lattice).jacobian();
                else if (order == 7)
                    map = get_one_turn_map<7>(lattice).jacobian();
                else
                    throw std::runtime_error("invalid order");

                // get_linear_one_turn_map() is fixed to order 2
                // auto map = get_linear_one_turn_map(lattice);

                auto arr =
                    py::array_t<double>({map.extent(0), map.extent(1)},
                                        {map.stride(0) * sizeof(double),
                                         map.stride(1) * sizeof(double)});

                for (int i = 0; i < map.extent(0); ++i)
                    for (int j = 0; j < map.extent(1); ++j)
                        arr.mutable_at(i, j) = map(i, j);

                return arr;

#if 0
			    double *buf = new double[6*6];
			    for(int i=0; i<36; ++i) buf[i] = map.data()[i];

			    py::capsule free_when_done(buf, [](void* f) {
				double *b = reinterpret_cast<double*>(f);
				delete[] b;
			    });

			    return py::array_t<double>(
				{ 6, 6 },
				{ sizeof(double) * 6, sizeof(double) * 1 },
				buf,
				free_when_done
			    );
#endif
            },
            "lattice"_a,
            "order"_a = 2)

        .def("get_one_turn_map_o1",
             &Lattice_simulator::get_one_turn_map<1>,
             "lattice"_a,
             "dpp"_a = 0.0)

        .def("get_one_turn_map_o2",
             &Lattice_simulator::get_one_turn_map<2>,
             "lattice"_a,
             "dpp"_a = 0.0)

        .def("get_one_turn_map_o3",
             &Lattice_simulator::get_one_turn_map<3>,
             "lattice"_a,
             "dpp"_a = 0.0)

        .def("get_one_turn_map_o4",
             &Lattice_simulator::get_one_turn_map<4>,
             "lattice"_a,
             "dpp"_a = 0.0)

        .def("get_one_turn_map_o5",
             &Lattice_simulator::get_one_turn_map<5>,
             "lattice"_a,
             "dpp"_a = 0.0)

        .def("get_one_turn_map_o6",
             &Lattice_simulator::get_one_turn_map<6>,
             "lattice"_a,
             "dpp"_a = 0.0)

        .def("get_one_turn_map_o7",
             &Lattice_simulator::get_one_turn_map<7>,
             "lattice"_a,
             "dpp"_a = 0.0)

        .def(
            "map_to_twiss",
            [](py::array_t<double> map) {
                if (map.ndim() != 2 || map.shape(0) != 2 || map.shape(1) != 2)
                    throw std::runtime_error(
                        "Lattice_simulator::map_to_twiss(map): "
                        "map must be a numpy array of (2, 2)");

                auto accessor = map.unchecked<2>();
                karray2d_row ka_map("map_to_twiss_python_input", 2, 2);
                for (size_t idx1 = 0; idx1 < 2; idx1++) {
                    for (size_t idx2 = 0; idx2 < 2; idx2++) {
                        ka_map(idx1, idx2) = accessor(idx1, idx2);
                    }
                }

                return Lattice_simulator::map_to_twiss(ka_map);
            },
            "map"_a);

    // Collective operator options (base class)

    // Stepper base class
    py::class_<Stepper>(m, "Stepper")
        .def_readonly_static("fixed_step_tolerance",
                             &Stepper::fixed_step_tolerance);

    // Split_operator_stepper
    py::class_<Split_operator_stepper, Stepper>(m, "Split_operator_stepper")
        .def(py::init<CO_options const&, int>(),
             "collective_options"_a,
             "num_steps"_a)
        .def("append_collective_op",
             &Split_operator_stepper::append_collective_op,
             "Append a collective operator",
             "collective_options"_a);

    // Split_operator_stepper_elements
    py::class_<Split_operator_stepper_elements, Stepper>(
        m, "Split_operator_stepper_elements")
        .def(py::init<CO_options const&, int>(),
             "collective_options"_a,
             "num_steps"_a);

    // Independent_stepper_elements
    py::class_<Independent_stepper_elements, Stepper>(
        m, "Independent_stepper_elements")
        .def(py::init<int>(), "steps_per_element"_a = 1);

    // Bunch_simulator
    using action_step_t =
        std::function<void(Bunch_simulator&, Lattice&, int, int)>;
    using action_turn_t = std::function<void(Bunch_simulator&, Lattice&, int)>;

    py::class_<Bunch_simulator>(m, "Bunch_simulator")
        .def_static("create_single_bunch_simulator",
                    &Bunch_simulator::create_single_bunch_simulator,
                    py::return_value_policy::move,
                    "Create a Bunch_simulator with a single bunch.",
                    "reference_particle"_a,
                    "num_particles"_a,
                    "num_real_particles"_a,
                    "commxx"_a = Commxx(),
                    "num_spectators"_a = 0)

        .def_static("create_bunch_train_simulator",
                    &Bunch_simulator::create_bunch_train_simulator,
                    py::return_value_policy::move,
                    "Create a Bunch_simulator with a single bunch train.",
                    "reference_particle"_a,
                    "num_particles"_a,
                    "num_real_particles"_a,
                    "num_bunches"_a,
                    "spacing"_a,
                    "commxx"_a = Commxx(),
                    "num_spectators"_a = 0)

        .def_static("create_two_trains_simulator",
                    &Bunch_simulator::create_two_trains_simulator,
                    py::return_value_policy::move,
                    "Create a Bunch_simulator with two bunch trains.",
                    "reference_particle_pri"_a,
                    "reference_particle_sec"_a,
                    "num_particles"_a,
                    "num_real_particles"_a,
                    "num_bunches_pri"_a = 1,
                    "num_bunch_sec"_a = 1,
                    "spacing_pri"_a = 0,
                    "spacing_sec"_a = 0,
                    "commxx"_a = Commxx(),
                    "num_spectators"_a = 0)

#if 0
		.def( "set_turns",
			&Bunch_simulator::set_turns,
			"Set the simulation start turn and total number of turns.",
			"first_turn"_a, "num_turns"_a )
#endif

        .def("get_bunch",
             (Bunch const& (Bunch_simulator::*)(size_t, size_t) const) &
                 Bunch_simulator::get_bunch,
             // py::overload_cast<size_t, size_t>(&Bunch_simulator::get_bunch,
             // py::const_),
             py::return_value_policy::reference_internal,
             "Get the bunch reference from the bunch_simulator.",
             "train"_a = 0,
             "bunch"_a = 0)

        .def("has_local_bunch",
             &Bunch_simulator::has_local_bunch,
             "Check if the bunch is present on current rank.",
             "train"_a,
             "bunch"_a)

        .def(
            "reg_prop_action_step_end",
            (void(Bunch_simulator::*)(action_step_t)) &
                Bunch_simulator::reg_prop_action_step_end,
            // py::overload_cast<action_step_t>(&Bunch_simulator::reg_prop_action_step_end),
            "Register the step end propagate action (callback)."
            "\nthe callback should have a signature of (Bunch_simulator, "
            "Lattice, "
            "turn_num, step_num)",
            "action"_a)

        .def(
            "reg_prop_action_turn_end",
            (void(Bunch_simulator::*)(action_turn_t)) &
                Bunch_simulator::reg_prop_action_turn_end,
            // py::overload_cast<action_turn_t>(&Bunch_simulator::reg_prop_action_turn_end),
            "Register the turn end propagate action (callback)."
            "\nthe callback should have a signature of (Bunch_simulator, "
            "Lattice, "
            "turn_num)",
            "action"_a)

        .def(
            "reg_diag_per_turn",
            [](Bunch_simulator& self,
               std::shared_ptr<Diagnostics> const& diag,
               int period,
               int bunch_idx,
               int train_idx) {
                PyDiagnostics* p = dynamic_cast<PyDiagnostics*>(diag.get());

                if (p) p->reg_self();

                return self.reg_diag_per_turn<std::shared_ptr<Diagnostics>>(
                    diag, period, bunch_idx, train_idx);
            },
            "Register a per turn diagnostics.",
            "diag"_a,
            "period"_a = 1,
            "bunch_idx"_a = 0,
            "train_idx"_a = 0)

        .def("reg_diag_per_step",
             &Bunch_simulator::reg_diag_per_step<std::shared_ptr<Diagnostics>>,
             "Register a per step diagnostics.",
             "diag"_a,
             "period"_a = 1,
             "bunch_idx"_a = 0,
             "train_idx"_a = 0)

        .def("reg_diag_turn_listed",
             &Bunch_simulator::reg_diag_turn_listed<
                 std::shared_ptr<Diagnostics>>,
             "Register a per step diagnostics.",
             "diag"_a,
             "turns"_a = std::vector<int>(),
             "bunch_idx"_a = 0,
             "train_idx"_a = 0)

        .def(
            "reg_diag_at_element",
            &Bunch_simulator::reg_diag_at_element<std::shared_ptr<Diagnostics>>,
            "Register a per element diagnostics.",
            "diag"_a,
            "element"_a,
            "turn_period"_a = 1,
            "bunch_idx"_a = 0,
            "train_idx"_a = 0)

        .def("reg_diag_loss_aperture",
             &Bunch_simulator::reg_diag_loss_aperture,
             "Register an aperture loss diagnostics to the bunch.",
             "filename"_a,
             "bunch_idx"_a = 0,
             "train_idx"_a = 0)

        .def("reg_diag_loss_zcut",
             &Bunch_simulator::reg_diag_loss_zcut,
             "Register a zcut loss diagnostics to the bunch.",
             "filename"_a,
             "bunch_idx"_a = 0,
             "train_idx"_a = 0)

        .def("current_turn",
             &Bunch_simulator::current_turn,
             "Get the current simulation turn.")

        .def("max_turns",
             &Bunch_simulator::max_turns,
             "Get the max simulation turns.")

        .def("set_max_turns",
             &Bunch_simulator::set_max_turns,
             "Set the max simulation turns.",
             "max_turns"_a)

        .def("set_longitudinal_boundary",
             &Bunch_simulator::set_longitudinal_boundary,
             "Set the longitudinal boundary for each bunch in the simulator",
             "boundary"_a,
             "param"_a = 0.0)

        .def("dump", &Bunch_simulator::dump, "Dump.")

        .def("populate_6d",
             [](Bunch_simulator& sim,
                uint64_t seed,
                py::buffer p_means,
                py::buffer p_covars) {
                 using ka1d_unmanaged =
                     Kokkos::View<double*,
                                  Kokkos::HostSpace,
                                  Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

                 using ka2d_unmanaged =
                     Kokkos::View<double**,
                                  Kokkos::LayoutRight,
                                  Kokkos::HostSpace,
                                  Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

                 py::buffer_info pm_info = p_means.request();
                 py::buffer_info pc_info = p_covars.request();

                 ka1d_unmanaged means((double*)pm_info.ptr, pm_info.shape[0]);

                 ka2d_unmanaged covars(
                     (double*)pc_info.ptr, pc_info.shape[0], pc_info.shape[1]);

                 sim.populate_6d(seed, means, covars);
             });
}
