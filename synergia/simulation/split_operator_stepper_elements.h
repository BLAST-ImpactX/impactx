#ifndef SPLIT_OPERATOR_STEPPER_ELEMENTS_H_
#define SPLIT_OPERATOR_STEPPER_ELEMENTS_H_

#include "synergia/simulation/collective_operator_options.h"
#include "synergia/simulation/stepper.h"

/// The Split_operator_stepper_elements class generates a constant number of
/// split-operator steps per thick element. Thin elements are assigned
/// a single step each. One or more collective effects are included per
/// step.
class Split_operator_stepper_elements : public Stepper {
private:
  int steps_per_element;
  std::shared_ptr<const CO_options> co_ops;

  std::vector<Step> apply_impl(Lattice const& lattice) const override;

public:
  Split_operator_stepper_elements(
    CO_options const& coo = CO_options{std::in_place_type<Dummy_CO_options>,
                                       Dummy_CO_options()},
    int steps_per_element = 1)
    : steps_per_element(steps_per_element)
    , co_ops(std::make_shared<const CO_options>(coo))
  {}

  std::unique_ptr<Stepper>
  clone() const override
  {
    return std::make_unique<Split_operator_stepper_elements>(*this);
  }

private:
  friend class cereal::access;

  template <class Archive>
  void
  serialize(Archive& ar)
  {
    ar(cereal::base_class<Stepper>(this));
    ar(steps_per_element);
    ar(co_ops);
  }
};

CEREAL_REGISTER_TYPE(Split_operator_stepper_elements)

#endif /* SPLIT_OPERATOR_STEPPER_ELEMENTS_H_ */
