#ifndef IMPLEMENTED_COLLECTIVE_OPERATORS_H
#define IMPLEMENTED_COLLECTIVE_OPERATORS_H

#include "synergia/collective/impedance.h"
#include "synergia/collective/space_charge_2d_kv.h"
#include "synergia/collective/space_charge_2d_open_hockney.h"

#ifdef BUILD_FD_SPACE_CHARGE_SOLVER
#include "synergia/collective/space_charge_3d_fd.h"
#endif

#include "synergia/collective/space_charge_3d_open_hockney.h"
#include "synergia/collective/space_charge_rectangular.h"

#endif
