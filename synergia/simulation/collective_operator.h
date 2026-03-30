#ifndef COLLECTIVE_OPERATOR_H_
#define COLLECTIVE_OPERATOR_H_

#include "synergia/simulation/operator.h"

// collective
class Collective_operator : public Operator {
private:
  void
  create_operations_impl(Lattice const& lattice) final
  {}
  void
  print_impl(Logger& logger) const override
  {}

public:
  Collective_operator(std::string const& name, double time)
    : Operator(name, "collective", time)
  {}
};

#endif
