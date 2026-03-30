#ifndef DISTRIBUTED_FFT3D_RECT_H
#define DISTRIBUTED_FFT3D_RECT_H

// discrete FFT on 3D grid.
// 2D DST-II on x-y plane and DFT along the z-axis

#include <array>

#include "synergia/utils/commxx.h"
#include "synergia/utils/kokkos_views.h"

class Distributed_fft3d_rect_base {
protected:
  std::array<int, 3> shape;
  Commxx comm;

  int lower;
  int nx;

public:
  static int
  get_padded_shape_real(int s)
  {
    return 2 * (s / 2 + 1);
  }

  static int
  get_padded_shape_cplx(int s)
  {
    return s / 2 + 1;
  }

  Distributed_fft3d_rect_base() : shape(), comm(Commxx::Null), lower(0), nx(0)
  {}

  virtual ~Distributed_fft3d_rect_base() = default;

  int
  get_lower() const
  {
    return lower;
  }
  int
  get_upper() const
  {
    return lower + nx;
  }

  int
  padded_nz_real() const
  {
    return get_padded_shape_real(shape[2]);
  }
  int
  padded_nz_cplx() const
  {
    return get_padded_shape_cplx(shape[2]);
  }

  std::array<int, 3> const&
  get_shape() const
  {
    return shape;
  }

  Commxx const&
  get_comm() const
  {
    return comm;
  }

  double
  get_roundtrip_normalization() const
  {
    return 1.0 / (shape[0] * shape[1] * shape[2]);
  }

  virtual void construct(std::array<int, 3> const& shape,
                         Commxx const& comm) = 0;

  virtual void transform(karray1d_dev& in, karray1d_dev& out) = 0;

  virtual void inv_transform(karray1d_dev& in, karray1d_dev& out) = 0;
};

#ifdef SYNERGIA_ENABLE_CUDA
#include "synergia/utils/distributed_fft3d_rect_cuda.h"
#else
#include "synergia/utils/distributed_fft3d_rect_fftw.h"
#endif

#endif
