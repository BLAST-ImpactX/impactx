#ifndef OPERATOR_H_
#define OPERATOR_H_

#include <string>

#include "synergia/lattice/lattice_element_slice.h"

#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/independent_operation.h"

#include "synergia/utils/logger.h"

class Propagator;

class Operator {

private:
  std::string name;
  std::string type;

  double time_fraction;

private:
  virtual void create_operations_impl(Lattice const& lattice) = 0;

  virtual void apply_impl(Bunch_simulator& simulator,
                          double time_step,
                          Logger& logger) = 0;

  virtual void print_impl(Logger& logger) const = 0;

public:
  Operator(std::string const& name, std::string const& type, double time)
    : name(name), type(type), time_fraction(time)
  {}

  virtual ~Operator() = default;

  std::string const&
  get_name() const
  {
    return name;
  }
  std::string const&
  get_type() const
  {
    return type;
  }
  double
  get_time_fraction() const
  {
    return time_fraction;
  }

  void
  create_operations(Lattice const& lattice)
  {
    create_operations_impl(lattice);
  }

  void
  apply(Bunch_simulator& simulator, double time_step, Logger& logger)
  {
    apply_impl(simulator, time_step * time_fraction, logger);
  }

  void
  print(Logger& logger) const
  {
    logger(LoggerV::DEBUG) << "operator name = " << name << ", type = " << type
                           << ", time_fraction = " << time_fraction << "\n";

    print_impl(logger);
  }
};

#endif /* OPERATOR_H_ */
