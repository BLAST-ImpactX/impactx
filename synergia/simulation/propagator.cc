#include "synergia/simulation/propagator.h"
#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/checkpoint.h"

#include "synergia/utils/digits.h"

void
Propagator::do_before_start(Bunch_simulator& simulator, Logger& logger)
{
    if (simulator.current_turn() == 0) {
        Kokkos::Profiling::pushRegion("diagnostic-actions");
        simulator.diag_action_step_and_turn(PRE_TURN, FINAL_STEP);
        Kokkos::Profiling::popRegion();

        simulator.prop_action_first(lattice);
        simulator.set_lattice_reference_particle(
            lattice.get_reference_particle());
    }

    auto updates = lattice.update();
    // if (updates.structure) steps = stepper.apply(lattice);
}

void
Propagator::do_start_repetition(Bunch_simulator& simulator)
{

    for (auto& bunch : simulator[0].get_bunches())
        bunch.get_reference_particle().start_repetition();

    for (auto& bunch : simulator[1].get_bunches())
        bunch.get_reference_particle().start_repetition();
}

void
Propagator::do_step(Bunch_simulator& simulator,
                    Step& step,
                    int step_count,
                    int turn_count,
                    Logger& logger)
{
    double t_step0 = MPI_Wtime();

    // make sure the lattice is up-to-date
    // e.g., update the chef_lattice after any of the
    // lattice elements has been updated
    lattice.update();

    // ensure that particles are on the device in case any of the
    // following might have transferred them onto the host
    // custom diagnostics routines, turn_end_action, etc
    Kokkos::Profiling::pushRegion("memory-location check before step-apply");
    for (auto& train : simulator.get_trains()) {
        for (auto& bunch : train.get_bunches()) {
            auto bparts = bunch.get_bunch_particles();
            if (bparts.get_memory_location() == MemoryLocation::Host) {
                std::runtime_error(
                    "Bunch particles are active on Host memory space!");
            }
        }
    }
    Kokkos::Profiling::popRegion();

    // propagate through the step
    step.apply(simulator, logger);

    // t = simple_timer_show(t, "propagate-step_apply");

    // operations associated with bunches
    // e.g., bunch longitudinal operations (periodic, zcut, etc)
    // these operations are not from lattices so are not included
    // in the steps
    // now I believe this can be included as part of the
    // prop_action_step_end() method
    // simulator.bunch_operation_step_end();
    // t = simple_timer_show(t, "propagate-bunch_operations_step");

    // general diagnostics
    Kokkos::Profiling::pushRegion("diagnostic-actions");
    simulator.diag_action_step_and_turn(turn_count, step_count);
    Kokkos::Profiling::popRegion();

    // t = simple_timer_show(t, "propagate-diagnostics_actions_step");

    // propagate action
    simulator.prop_action_step_end(lattice, turn_count, step_count);
    // t = simple_timer_show(t, "propagate-general_actions-step");

    double t_step1 = MPI_Wtime();

    logger(LoggerV::INFO_STEP)
        << "Propagator: step " << std::setw(digits(steps.size())) << step_count
        << "/" << steps.size()

        << ", s_n = " << std::fixed << std::setprecision(4)
        << simulator[0][0].get_reference_particle().get_s_n()

        << ", time = " << std::fixed << std::setprecision(3)
        << t_step1 - t_step0 << "s, macroparticles = ";

    for (auto const& train : simulator.get_trains()) {
        logger << "(";

        for (auto const& bunch : train.get_bunches()) {
            logger << bunch.get_total_num();
            if (bunch.get_array_index() != train.get_bunch_array_size() - 1)
                logger << ", ";
        }

        logger << ")";
        if (train.get_index() == 0) logger << " / ";
    }

    logger << "\n";

    logger(LoggerV::INFO_OPR) << "\n\n";
}

bool
Propagator::check_out_of_particles(Bunch_simulator const& simulator,
                                   Logger& logger)
{
    return false;
}

bool
Propagator::do_turn_end(Bunch_simulator& simulator,
                        int turn_count,
                        Logger& logger)
{
    // t = simple_timer_current();

    // diagnostic actions
    Kokkos::Profiling::pushRegion("diagnostic-actions");
    simulator.diag_action_step_and_turn(turn_count, Propagator::FINAL_STEP);
    Kokkos::Profiling::popRegion();

    // propagate actions
    simulator.prop_action_turn_end(lattice, turn_count);

    // update lattice in case it has been chaged in the propagate action
    auto updates = lattice.update();
    // if (updates.structure()) steps = stepper.apply(lattice);

    // increment the turn number
    simulator.inc_turn();

    if (halt_func_turn) {
        bool res = halt_func_turn(simulator);
        return res;
    } else {
        return false;
    }
}

void
Propagator::propagate(Bunch_simulator& sim, Logger& logger, int max_turns)
{
    Kokkos::Profiling::pushRegion("propagate-execute");
    const int total_turns = sim.max_turns();

    // parameter check
    if (max_turns == -1 && total_turns == -1) {
        throw std::runtime_error(
            "Max number of turns must be set either in the Bunch_simulator "
            "(Bunch_simulator::set_num_turns(int)), or in the Propagator "
            "(Propagator::propagate(..., int max_turns))");
    }

    try {
        do_before_start(sim, logger);

        // first turn is always the current turn from simulator
        int turn = sim.current_turn();

        // decide the last turn
        int last_turn = (max_turns == -1) ? total_turns : turn + max_turns;
        if ((last_turn > total_turns) && (total_turns != -1))
            last_turn = total_turns;

        int turns_since_checkpoint = 0;
        bool out_of_particles = false;
        double t_prop0 = MPI_Wtime();

        logger(LoggerV::INFO_TURN) << "Propagator: starting turn " << turn + 1
                                   << ", final turn " << last_turn << "\n\n";

        for (; turn < last_turn; ++turn) {
            double t_turn0 = MPI_Wtime();

            do_start_repetition(sim);

            int step_count = 0;
            for (auto& step : steps) {
                ++step_count;
                do_step(sim, step, step_count, turn, logger);

                out_of_particles = check_out_of_particles(sim, logger);
                if (out_of_particles) break;
            }

            double t_turn1 = MPI_Wtime();

            // turn end log
            logger(LoggerV::INFO_TURN)
                << "Propagator: turn " << std::setw(digits(total_turns))
                << turn + 1 << "/";

            if (total_turns == -1)
                logger << "inf.";
            else
                logger << total_turns;

            logger << ", time = " << std::fixed << std::setprecision(3)
                   << t_turn1 - t_turn0 << "s"

                   << ", macroparticles = ";

            for (auto const& train : sim.get_trains()) {
                logger << "(";

                for (auto const& bunch : train.get_bunches()) {
                    logger << bunch.get_total_num();
                    if (bunch.get_array_index() !=
                        train.get_bunch_array_size() - 1)
                        logger << ", ";
                }

                logger << ")";
                if (train.get_index() == 0) logger << " / ";
            }

            logger << "\n";
            logger(LoggerV::INFO_STEP) << "\n";

            // out of particles
            if (out_of_particles) break;

            ++turns_since_checkpoint;
            bool halt = do_turn_end(sim, turn, logger);

            if (halt) {
                logger(LoggerV::INFO_TURN) << "Propagator: halting condition "
                                              "satisfied, ending propagation\n";
                break;
            }

            // checkpoint save
            // syn::checkpoint_save(*this, sim);

            if ((turns_since_checkpoint == checkpoint_period) ||
                ((turn == (sim.max_turns() - 1)) && final_checkpoint)) {
                // t = simple_timer_current();
                syn::checkpoint_save(*this, sim);
                // t = simple_timer_show(t, "propagate-checkpoint_period");
                turns_since_checkpoint = 0;
            }
        }

        if (last_turn != total_turns) {
            logger(LoggerV::INFO_TURN)
                << "Propagator: maximum number of turns reached\n";

            // TODO: checkpoint
        }

        if (out_of_particles) {
            logger(LoggerV::WARNING) << "Propagator: no particles left\n";
        }

        double t_prop1 = MPI_Wtime();

        logger(LoggerV::INFO_TURN)
            << "Propagator: total time = " << std::fixed << std::setprecision(3)
            << t_prop1 - t_prop0 << "s\n";
    }
    catch (std::exception const& e) {
        std::cerr << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 888);
    }

    Kokkos::Profiling::popRegion();
}
