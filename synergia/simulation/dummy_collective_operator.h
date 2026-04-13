#ifndef DUMMY_COLLECTIVE_OPERATOR_H
#define DUMMY_COLLECTIVE_OPERATOR_H

#include "synergia/simulation/collective_operator.h"

class Dummy_collective_operator : public Collective_operator {
private:
  virtual void
  apply_impl(Bunch_simulator& simulator,
             double time_step,
             Logger& logger) override
  {}

public:
  Dummy_collective_operator() : Collective_operator("dummy collective", 1.0) {}
};

#endif
