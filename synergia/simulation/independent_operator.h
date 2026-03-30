#ifndef INDEPENDENT_OPERATOR_H_
#define INDEPENDENT_OPERATOR_H_

#include "synergia/simulation/operator.h"

// Independent
class Independent_operator : public Operator {
private:
  void apply_impl(Bunch_simulator& simulator,
                  double time_step,
                  Logger& logger) override;

  void create_operations_impl(Lattice const& lattice) override;

  void print_impl(Logger& logger) const override;

  bool need_update(Reference_particle const& ref, Logger& logger);

  void update_operations(Reference_particle const& ref);

private:
  std::vector<Lattice_element_slice> slices;
  std::vector<std::unique_ptr<Independent_operation>> operations;

public:
  Independent_operator(std::string const& name, double time);

  template <class... Args>
  Independent_operator&
  append_slice(Args&&... args)
  {
    slices.emplace_back(std::forward<Args>(args)...);
    return *this;
  }

  std::vector<Lattice_element_slice> const&
  get_slices() const
  {
    return slices;
  }

  friend class Propagator;
};

#endif
