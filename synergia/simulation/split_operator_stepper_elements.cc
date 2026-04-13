#include "split_operator_stepper_elements.h"

std::vector<Step>
Split_operator_stepper_elements::apply_impl(Lattice const& lattice) const
{
  if (steps_per_element < 1) {
    throw std::runtime_error(
      "Split_operator_stepper_elements: steps_per_element must be >= 1");
  }

  std::vector<Step> steps;

  auto col_op_ptr = std::visit(create_collective_operator{}, *co_ops);

  for (auto const& ele : lattice.get_elements()) {
    double length = ele.get_length();

    // zero-length element
    if (length == 0.0) {
      steps.emplace_back(0.0);
      steps.back().append_independent("step", 1.0).append_slice(ele);
    } else {
      double step_length = length / steps_per_element;
      double half_step_length = 0.5 * step_length;

      for (int i = 0; i < steps_per_element; ++i) {
        double left = i * step_length;
        double middle = left + half_step_length;
        double right = (i + 1) * step_length;

        steps.emplace_back(step_length);

        // 1st Half
        steps.back()
          .append_independent("first_half", 0.5)
          .append_slice(ele, left, middle);

        // Collective Effects
        // steps.back().append_collective(co_ops);
        steps.back().append(col_op_ptr);

        // 2nd Half
        steps.back()
          .append_independent("second_half", 0.5)
          .append_slice(ele, middle, right);
      }
    }
  }

  return steps;
}
