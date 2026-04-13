#ifndef STEP_H_
#define STEP_H_

#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/collective_operator_options.h"
#include "synergia/simulation/independent_operator.h"
#include "synergia/simulation/operator.h"
#include "synergia/utils/logger.h"

class Propagator;

class Step {

private:
  std::vector<std::shared_ptr<Operator>> operators;

  double length;
  std::vector<double> step_betas;

public:
  explicit Step(double length);

  template <class... Args>
  Independent_operator&
  append_independent(Args&&... args)
  {
    operators.emplace_back(
      std::make_shared<Independent_operator>(std::forward<Args>(args)...));
    return *(dynamic_cast<Independent_operator*>(operators.back().get()));
  }

  void
  append_collective(std::unique_ptr<CO_options> const& co_ops)
  {
    operators.emplace_back(std::visit(create_collective_operator{}, *co_ops));
  }

  void
  append(std::shared_ptr<Operator> col_op)
  {
    operators.push_back(col_op);
  }

  void apply(Bunch_simulator& simulator, Logger& logger) const;
  void create_operations(Lattice const& lattice);

  double
  get_length() const
  {
    return length;
  }

  void
  print(Logger& logger) const
  {
    logger(LoggerV::DEBUG) << "step length = " << length << "\n";
    for (auto const& opr : operators) opr->print(logger);
    logger(LoggerV::DEBUG) << "\n";
  }

#if 0
    void set_betas(double betax, double betay);
    std::vector<double > get_betas();

    void print(int index) const;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif

  friend class Propagator;
};

#endif /* STEP_H_ */
