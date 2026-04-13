#ifndef SPLIT_OPERATOR_STEPPER_H_
#define SPLIT_OPERATOR_STEPPER_H_

#include "synergia/utils/cereal.h"

#include "synergia/simulation/collective_operator_options.h"
#include "synergia/simulation/stepper.h"

/// The Split_operator_stepper class generates evenly-spaced split-operator
/// steps through a Lattice. One or more collective effects are included per
/// step.
class Split_operator_stepper : public Stepper {
  private:
    int num_steps;
    std::vector<std::shared_ptr<const CO_options>> co_ops;

    std::vector<Step> apply_impl(Lattice const& lattice) const override;

  public:
    Split_operator_stepper(
        CO_options const& coo = CO_options(Dummy_CO_options()),
        int num_steps = 1)
        : num_steps(num_steps)
    {
        co_ops.emplace_back(std::make_shared<const CO_options>(coo));
    }

    template <class... ARGS>
    Split_operator_stepper(int num_steps, ARGS const&... args)
        : num_steps(num_steps)
        , co_ops({(std::make_shared<const CO_options>(args))...})
    {}

    void
    append_collective_op(CO_options const& coo)
    {
        co_ops.push_back(std::make_shared<const CO_options>(coo));
    }

    std::unique_ptr<Stepper>
    clone() const override
    {
        return std::make_unique<Split_operator_stepper>(*this);
    }

  private:
    friend class cereal::access;

    template <class Archive>
    void
    serialize(Archive& ar)
    {
        ar(cereal::virtual_base_class<Stepper>(this));
        ar(num_steps);
        ar(co_ops);
    }
};

CEREAL_REGISTER_TYPE(Split_operator_stepper)

#endif /* SPLIT_OPERATOR_STEPPER_H_ */
