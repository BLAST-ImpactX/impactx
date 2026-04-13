#include "split_operator_stepper.h"
#include "synergia/utils/string_utils.h"

namespace {
    const std::string force_diag_attr = "force_diagnostics";
    const double fixed_step_tolerance = 1.0e-8;

    // Return an Independent_operator for a half step, starting at the
    // lattice_element given by lattice_it at position left. Both lattice_it
    // and left are updated by the function.

    template <typename ITER>
    std::vector<Lattice_element_slice>
    get_fixed_step_slices(ITER& it,
                          ITER const& end,
                          double& left,
                          double step_length,
                          double& offset_fudge,
                          bool end_on_force_diagnostics)
    {
        std::vector<Lattice_element_slice> slices;

        double length = offset_fudge;
        bool complete = false;

        while (!complete) {
            double right = it->get_length();

            if (length + (right - left) - fixed_step_tolerance <= step_length) {
                // The rest of the element fits in the half step
                Lattice_element_slice slice(*it, left, right);
                slices.push_back(slice);
                length += (right - left);

                if (end_on_force_diagnostics && slice.has_right_edge() &&
                    it->has_string_attribute(force_diag_attr) &&
                    !false_string(it->get_string_attribute(force_diag_attr))) {
                    complete = true;
                }

                ++it;
                left = 0.0;

                if (std::abs(length - step_length) < fixed_step_tolerance) {
                    if ((it == end) || (it->get_length() != 0.0))
                        complete = true;
                } else {
                    if (it == end)
                        throw std::runtime_error(
                            "get_step stepped beyond end of lattice");
                }
            } else {
                // Need to take a portion of the element...
                bool end_within_error = false;
                double old_right = right;
                right = step_length - length + left;

                if (std::abs(old_right - right) < fixed_step_tolerance) {
                    // ... unless we are within an accumulated tolerance of the
                    // end
                    right = old_right;
                    end_within_error = true;
                }

                Lattice_element_slice slice(*it, left, right);
                slices.push_back(slice);
                length += (right - left);

                if (end_within_error) {
                    ++it;
                    left = 0.0;
                } else {
                    left = right;
                }

                complete = true;
            }
        }

        offset_fudge = length - step_length;
        return slices;
    }

    template <typename ITER>
    void
    create_substep(std::vector<Step>& steps,
                   std::vector<std::shared_ptr<Operator>> const& col_op_ptrs,
                   ITER& it,
                   ITER const& end,
                   double left,
                   double length,
                   double& offset_fudge)
    {
        // create the substep
        steps.emplace_back(length);

        if (length == 0.0) {
            auto slices =
                get_fixed_step_slices(it, end, left, 0.0, offset_fudge, true);

            auto& op = steps.back().append_independent("zero_length", 1.0);
            for (auto const& s : slices)
                op.append_slice(s);
        } else {
            double half_length = 0.5 * length;

            // slices for first half substep
            auto fst_half_slices = get_fixed_step_slices(
                it, end, left, half_length, offset_fudge, false);

            // slices for second half substep
            auto snd_half_slices = get_fixed_step_slices(
                it, end, left, half_length, offset_fudge, false);

            auto& op1 = steps.back().append_independent("first_half", 0.5);
            for (auto const& s : fst_half_slices)
                op1.append_slice(s);

            // collective
            for (auto col_op_ptr : col_op_ptrs)
                steps.back().append(col_op_ptr);

            // 2nd half
            auto& op2 = steps.back().append_independent("second_half", 0.5);
            for (auto const& s : snd_half_slices)
                op2.append_slice(s);
        }
    }
}

std::vector<Step>
Split_operator_stepper::apply_impl(Lattice const& lattice) const
{
    std::vector<Step> steps;
    std::vector<std::shared_ptr<Operator>> col_op_ptrs;

    for (auto co_op : co_ops) {
        CO_options _co_op = *co_op.get();
        std::shared_ptr<Operator> _op =
            std::visit(create_collective_operator{}, _co_op);
        col_op_ptrs.push_back(_op);
    }

    double step_length = lattice.get_length() / num_steps;
    double half_step_length = 0.5 * step_length;

    auto lattice_it = lattice.get_elements().begin();
    auto lattice_end = lattice.get_elements().end();

    // empty lattice
    if (lattice_it == lattice_end) return steps;

    double left = 0.0;
    double offset_fudge = 0.0;

    for (int i = 0; i < num_steps; ++i) {
        auto substep_lattice_it = lattice_it;

        double substep_left = left;
        double substep_offset_fudge = offset_fudge;

        // slices for first half step
        auto fst_half_slices = get_fixed_step_slices(lattice_it,
                                                     lattice_end,
                                                     left,
                                                     half_step_length,
                                                     offset_fudge,
                                                     false);

        // slices for second half step
        auto snd_half_slices = get_fixed_step_slices(lattice_it,
                                                     lattice_end,
                                                     left,
                                                     half_step_length,
                                                     offset_fudge,
                                                     false);

        // merge first half slices and second half slices
        std::vector<Lattice_element_slice> all_slices;
        all_slices.reserve(fst_half_slices.size() + snd_half_slices.size());
        all_slices.insert(
            all_slices.end(), fst_half_slices.begin(), fst_half_slices.end());
        all_slices.insert(
            all_slices.end(), snd_half_slices.begin(), snd_half_slices.end());

        double substep_length = 0.0;
        double all_substeps_length = 0.0;
        bool found_force = false;

        // skip the last slice
        for (int it = 0; it < all_slices.size() - 1; ++it) {
            auto const& s = all_slices[it];
            auto const& ele = s.get_lattice_element();

            substep_length += s.get_right() - s.get_left();

            // force diagnostics?
            if (s.has_right_edge() &&
                ele.has_string_attribute(force_diag_attr) &&
                !false_string(ele.get_string_attribute(force_diag_attr))) {
                found_force = true;
                all_substeps_length += substep_length;

                // create the substep
                create_substep(steps,
                               col_op_ptrs,
                               substep_lattice_it,
                               lattice_end,
                               substep_left,
                               substep_length,
                               substep_offset_fudge);

                // reset length
                substep_length = 0.0;
            }
        }

        if (found_force) {
            double remain_length = step_length - all_substeps_length;
            if (remain_length <= fixed_step_tolerance) remain_length = 0.0;

            // create the substep
            create_substep(steps,
                           col_op_ptrs,
                           substep_lattice_it,
                           lattice_end,
                           substep_left,
                           remain_length,
                           substep_offset_fudge);

            if (substep_lattice_it != lattice_it)
                throw std::runtime_error(
                    "internal error: Split_operator_stepper "
                    "created an inconsistent force_diagnostics step");
        } else {
            // new step
            steps.emplace_back(step_length);

            // 1st half
            auto& op1 = steps.back().append_independent("first_half", 0.5);
            for (auto const& s : fst_half_slices)
                op1.append_slice(s);

            // collective
            for (auto col_op_ptr : col_op_ptrs)
                steps.back().append(col_op_ptr);

            // 2nd half
            auto& op2 = steps.back().append_independent("second_half", 0.5);
            for (auto const& s : snd_half_slices)
                op2.append_slice(s);
        }
    }

    if (lattice_it != lattice_end) {
        throw std::runtime_error("internal error: split_operator_stepper "
                                 "did not make it to the end of the lattice\n");
    }

    return steps;
}
