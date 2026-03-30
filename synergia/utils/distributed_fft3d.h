#ifndef DISTRIBUTED_FFT3D_H
#define DISTRIBUTED_FFT3D_H

#include <array>

#include "synergia/utils/commxx.h"
#include "synergia/utils/kokkos_views.h"

class Distributed_fft3d_base {
  protected:
    std::array<int, 3> shape;
    Commxx comm;

    int lower;
    int nz;

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

    Distributed_fft3d_base() : shape(), comm(Commxx::Null), lower(0), nz(0) {}

    virtual ~Distributed_fft3d_base() = default;

    int
    get_lower() const
    {
        return lower;
    }
    int
    get_upper() const
    {
        return lower + nz;
    }

    int
    padded_nx_real() const
    {
        return get_padded_shape_real(shape[0]);
    }
    int
    padded_nx_cplx() const
    {
        return get_padded_shape_cplx(shape[0]);
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
};

#ifdef SYNERGIA_ENABLE_CUDA
#include "synergia/utils/distributed_fft3d_cuda.h"
#else
#include "synergia/utils/distributed_fft3d_fftw.h"
#endif

#endif
