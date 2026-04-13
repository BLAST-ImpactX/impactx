
#include "independent_operation.h"
#include "synergia/libFF/ff_element.h"

void
LibFF_operation::apply_impl(Bunch& bunch, Logger& logger) const
{
  for (auto const& slice : slices) FF_element::apply(slice, bunch);
}
