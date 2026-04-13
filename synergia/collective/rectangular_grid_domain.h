#ifndef RECTANGULAR_GRID_DOMAIN_H_
#define RECTANGULAR_GRID_DOMAIN_H_

#include "synergia/collective/rectangular_grid.h"
#include "synergia/utils/fast_int_floor.h"
#include <array>

class Rectangular_grid_domain {

private:
  std::array<int, 3> grid_shape;
  std::array<double, 3> physical_size;
  std::array<double, 3> physical_offset;
  std::array<double, 3> left;
  std::array<double, 3> cell_size;

  bool periodic_z;

public:
  Rectangular_grid_domain(std::array<int, 3> const& grid_shape,
                          std::array<double, 3> const& physical_size,
                          std::array<double, 3> const& physical_offset = {0,
                                                                          0,
                                                                          0},
                          bool periodic_z = false)
    : grid_shape(grid_shape)
    , physical_size(physical_size)
    , physical_offset(physical_offset)
    , left{physical_offset[0] - physical_size[0] / 2.0,
           physical_offset[1] - physical_size[1] / 2.0,
           physical_offset[2] - physical_size[2] / 2.0}
    , cell_size{physical_size[0] / (1.0 * grid_shape[0]),
                physical_size[1] / (1.0 * grid_shape[1]),
                physical_size[2] / (1.0 * grid_shape[2])}
    , periodic_z(periodic_z)
  {}

  bool
  is_periodic() const
  {
    return periodic_z;
  }

  std::array<double, 3> const&
  get_physical_size() const
  {
    return physical_size;
  }

  std::array<double, 3> const&
  get_physical_offset() const
  {
    return physical_offset;
  }

  std::array<int, 3> const&
  get_grid_shape() const
  {
    return grid_shape;
  }

  std::array<double, 3> const&
  get_cell_size() const
  {
    return cell_size;
  }

  std::array<double, 3> const&
  get_left() const
  {
    return left;
  }

  Rectangular_grid_2dc
  make_2dc_grid_xy(bool zero = true) const
  {
    return Rectangular_grid_2dc(grid_shape[0], grid_shape[1], zero);
  }

  Rectangular_grid_1d
  make_1d_grid_z(bool zero = true) const
  {
    return Rectangular_grid_1d(grid_shape[2], zero);
  }

  // returns cell location and fractional offset
  inline bool
  get_leftmost_indices_offsets(double x,
                               double y,
                               double z,
                               int& ix,
                               int& iy,
                               int& iz,
                               double& offx,
                               double& offy,
                               double& offz) const
  {
    bool retval;
    double scaled_location;

    scaled_location = (x - left[0]) / cell_size[0] - 0.5;
    ix = fast_int_floor(scaled_location);
    offx = scaled_location - ix;

    scaled_location = (y - left[1]) / cell_size[1] - 0.5;
    iy = fast_int_floor(scaled_location);
    offy = scaled_location - iy;

    scaled_location = (z - left[2]) / cell_size[2] - 0.5;
    iz = fast_int_floor(scaled_location);
    offz = scaled_location - iz;

    if (grid_shape[2] == 1) {
      // csp: For grid_shape = 1, iz and offz are not used in deposit
      //      and interpolation. These are just for the reference.
      //      iz is 0 or 1, so that all particles are in domain, i.e.,
      //      no cutting edge.
      iz += 1;
      if (iz == 0) offz = -0.5 + offz;
      if (iz == 1) offz = 0.5 - offz;

      retval = ((ix >= 0) && (ix < grid_shape[0] - 1) && (iy >= 0) &&
                (iy < grid_shape[1] - 1)) &&
               (periodic_z || ((iz >= 0) && (iz <= grid_shape[2])));
    } else {
      retval = ((ix >= 0) && (ix < grid_shape[0] - 1) && (iy >= 0) &&
                (iy < grid_shape[1] - 1)) &&
               (periodic_z || ((iz >= 0) && (iz < grid_shape[2] - 1)));
    }

    return retval;
  }

  void
  get_cell_coordinates(int ix, int iy, int iz, double& x, double& y, double& z)
    const
  {
    x = left[0] + cell_size[0] * (0.5 + ix);
    y = left[1] + cell_size[1] * (0.5 + iy);
    z = left[2] + cell_size[2] * (0.5 + iz);
  }
};

#endif /* RECTANGULAR_GRID_DOMAIN_H_ */
