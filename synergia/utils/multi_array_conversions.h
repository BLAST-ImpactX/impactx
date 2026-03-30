#ifndef SYNERGIA_UTILS_MULTI_ARRAY_CONVERSIONS_H
#define SYNERGIA_UTILS_MULTI_ARRAY_CONVERSIONS_H

#include "Eigen/Eigen"
#include "synergia/utils/kokkos_views.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
  MatrixD;

inline MatrixD
karray_to_matrix(karray2d_row arr)
{
  auto d0 = arr.extent(0);
  auto d1 = arr.extent(1);
  MatrixD m(d0, d1);

  for (auto i = 0; i < d0; ++i)
    for (auto j = 0; j < d1; ++j) m(i, j) = arr(i, j);

  return m;
}

#endif
