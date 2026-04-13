#ifndef SYNERGIA_SIMULATION_CHECKPOINT_H
#define SYNERGIA_SIMULATION_CHECKPOINT_H

#include <utility>

class Propagator;
class Bunch_simulator;

namespace syn {
    // save the current state in the checkpoint
    void checkpoint_save(Propagator const& prop, Bunch_simulator const& sim);

    // load from the most recent checkpoint
    std::pair<Propagator, Bunch_simulator> checkpoint_load();

    // load and resume simulation from most recent checkpoint
    void resume();
}

#endif
