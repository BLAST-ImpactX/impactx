#include "independent_stepper_elements.h"

std::vector<Step>
Independent_stepper_elements::apply_impl(Lattice const& lattice) const
{
  if (steps_per_element < 1) {
    throw std::runtime_error(
      "Split_operator_stepper_elements: steps_per_element must be >= 1");
  }

  std::vector<Step> steps;

  for (auto const& ele : lattice.get_elements()) {
    double length = ele.get_length();

    // zero-length element
    if (length == 0.0) {
      steps.emplace_back(0.0);
      steps.back().append_independent("step", 1.0).append_slice(ele);
    } else {
      double step_length = length / steps_per_element;

      for (int i = 0; i < steps_per_element; ++i) {
        double left = i * step_length;
        double right = (i + 1) * step_length;

        steps.emplace_back(step_length);

        steps.back()
          .append_independent("step", 1.0)
          .append_slice(ele, left, right);
      }
    }
  }

  return steps;
}
