#ifndef SYNERGIA_SIMULATION_LATTICE_SIMULATOR_HOST_H
#define SYNERGIA_SIMULATION_LATTICE_SIMULATOR_HOST_H

#include <iomanip>
#include <iostream>
#include <sstream>

#include "Eigen/Eigen"

namespace Lattice_simulator {
    std::array<double, 2> filter_transverse_tunes(double const* jac_arr);
}

#endif
