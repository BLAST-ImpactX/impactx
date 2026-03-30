#ifndef PROPAGATE_ACTIONS_H
#define PROPAGATE_ACTIONS_H

#include <cereal/archives/json.hpp>
#include <cereal/types/polymorphic.hpp>
// #include <cereal/archives/binary.hpp>
// #include <cereal/archives/xml.hpp>

class Bunch_simulator;
class Lattice;

class Propagate_actions {

  public:
    Propagate_actions() {}
    virtual ~Propagate_actions() = default;

    virtual void
    first(Bunch_simulator& sim, Lattice& l)
    {}
    virtual void
    turn_end(Bunch_simulator& sim, Lattice& l, int turn)
    {}
    virtual void
    step_end(Bunch_simulator& sim, Lattice& l, int turn, int step)
    {}
    virtual void
    before_resume(Bunch_simulator& sim, Lattice& l)
    {}

    template <class AR>
    void
    serialize(AR& ar)
    {}
};

#endif
