#include <iostream>

#include "synergia/foundation/reference_particle.h"

#include "synergia/simulation/aperture_operation.h"
#include "synergia/simulation/independent_operator.h"
#include "synergia/simulation/operation_extractor.h"

Independent_operator::Independent_operator(std::string const& name, double time)
  : Operator(name, "independent", time), slices(), operations()
{}

bool
Independent_operator::need_update(Reference_particle const& ref, Logger& logger)
{
  return false;
}

void
Independent_operator::update_operations(
  Reference_particle const& reference_particle)
{}

void
Independent_operator::create_operations_impl(Lattice const& lattice)
{
  operations.clear();
  // operations_revisions.clear();

  std::string aperture_type("");
  bool need_left_aperture, need_right_aperture;

  std::string extractor_type(""), last_extractor_type("");

  // Group slices of equal extractor_type and pass to operation_extractor
  // to get operations.
  std::vector<Lattice_element_slice> group;
  for (auto const& slice : slices) {
    auto const& element = slice.get_lattice_element();

    if (element.has_string_attribute("aperture_type")) {
      aperture_type = element.get_string_attribute("aperture_type");
      need_left_aperture = slice.has_left_edge();
      need_right_aperture = slice.has_right_edge();
    } else {
      need_left_aperture = false;
      need_right_aperture = false;
    }

    extractor_type = element.get_string_attribute("extractor_type", "default");

    if (((extractor_type != last_extractor_type) || need_left_aperture) &&
        (!group.empty())) {
      extract_independent_operations(
        extractor_type, lattice, group, operations);
      group.clear();
    }

    if (need_left_aperture) {
      operations.emplace_back(extract_aperture_operation(aperture_type, slice));
    }

    group.push_back(slice);
    last_extractor_type = extractor_type;

    if (need_right_aperture) {
      extract_independent_operations(
        extractor_type, lattice, group, operations);
      operations.emplace_back(extract_aperture_operation(aperture_type, slice));
      group.clear();
    }

    // operations_revisions.push_back(element.get_revision());
  }

  if (!group.empty()) {
    extract_independent_operations(extractor_type, lattice, group, operations);
  }

  // always attach a finite aperture and a circular aperture by default
  operations.emplace_back(
    extract_aperture_operation(Finite_aperture::type, slices.back()));

  operations.emplace_back(
    extract_aperture_operation(Circular_aperture::type, slices.back()));

#if 0
    have_operations = true;
    operations_reference_particle = reference_particle;
#endif
}

void
Independent_operator::apply_impl(Bunch_simulator& simulator,
                                 double time_step,
                                 Logger& logger)
{
  using LV = LoggerV;

  // double t_total = simple_timer_current();

  logger(LV::DINFO) << "    Independent_operator: slice(s) = ";

  for (auto const& slice : slices)
    logger(LV::DINFO) << slice.as_string() << ", ";

  logger(LV::DINFO) << "\n";

  for (auto& train : simulator.get_trains()) {
    for (auto& bunch : train.get_bunches()) {
      for (auto const& opn : operations) {
        logger(LV::INFO_OPN)
          << "    Independent_operator: operation type = " << opn->get_type()
          << "\n";

        opn->apply(bunch, logger);

        simulator.diag_action_operation(*opn);
      }

      // per element diagnostics
      // the former "forced diagnostics"
      for (auto const& slice : slices) {
        if (slice.has_right_edge())
          simulator.diag_action_element(slice.get_lattice_element());
      }

      // update per-bunch per-independent-operator
      bunch.update_total_num();
    }
  }
}

void
Independent_operator::print_impl(Logger& logger) const
{
  logger(LoggerV::DEBUG) << "\tslices: "
                         << "\n"
                         << "\toperations: ";

  for (auto const& opn : operations) opn->print(logger);

  logger(LoggerV::DEBUG) << "\n";
}
