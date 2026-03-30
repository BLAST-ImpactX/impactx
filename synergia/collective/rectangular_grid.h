#ifndef RECTANGULAR_GRID_H_
#define RECTANGULAR_GRID_H_

#include "synergia/utils/kokkos_views.h"

class Rectangular_grid_1d {

private:
  using kokkos_noinit = Kokkos::ViewAllocateWithoutInitializing;

  karray1d_row_dev grid_;
  karray1d_row_hst hgrid_;

  double normalization;

public:
  Rectangular_grid_1d(size_t shape_x, bool zero = true)
    : grid_(zero ? karray1d_row_dev("g1d", shape_x) :
                   karray1d_row_dev(kokkos_noinit("g1d"), shape_x))
    , hgrid_(Kokkos::create_mirror_view(grid_))
    , normalization(1.0)
  {}

  size_t
  shape(size_t dim) const
  {
    return grid_.extent(dim);
  }

  size_t
  span() const
  {
    return grid_.span();
  }

  double*
  data(size_t x = 0) const
  {
    return &grid_(x);
  }

  // void set_zero()
  //{ points_.setZero(); }

  void
  set_normalization(double val)
  {
    normalization = val;
  }

  double
  get_normalization() const
  {
    return normalization;
  }
};

class Rectangular_grid_2dc {

private:
  using kokkos_noinit = Kokkos::ViewAllocateWithoutInitializing;

  karray2dc_row_dev grid_;
  karray2dc_row_hst hgrid_;

  double normalization;

public:
  Rectangular_grid_2dc(size_t shape_x, size_t shape_y, bool zero = true)
    : grid_(zero ? karray2dc_row_dev("g2dc", shape_x, shape_y) :
                   karray2dc_row_dev(kokkos_noinit("g2dc"), shape_x, shape_y))
    , hgrid_(Kokkos::create_mirror_view(grid_))
    , normalization(1.0)
  {}

  size_t
  shape(size_t dim) const
  {
    return grid_.extent(dim);
  }

  size_t
  span() const
  {
    return grid_.span();
  }

  Kokkos::complex<double>*
  data(size_t x = 0, size_t y = 0) const
  {
    return &grid_(x, y);
  }

  // void set_zero()
  //{ points_.setZero(); }

  void
  set_normalization(double val)
  {
    normalization = val;
  }

  double
  get_normalization() const
  {
    return normalization;
  }
};

#endif /* RECTANGULAR_GRID_EIGEN_H_ */
