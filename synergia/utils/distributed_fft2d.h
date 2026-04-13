#ifndef DISTRIBUTED_FFT2D_H
#define DISTRIBUTED_FFT2D_H

#include <array>

#include "synergia/utils/commxx.h"
#include "synergia/utils/kokkos_views.h"

class Distributed_fft2d_base {

  protected:
    std::array<int, 2> shape;
    Commxx comm;

    int lower;
    int nx;

  public:
    Distributed_fft2d_base() : shape(), comm(Commxx::Null), lower(0), nx(0) {}

    virtual ~Distributed_fft2d_base() {}

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

    std::array<int, 2> const&
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
        return 1.0 / (shape[0] * shape[1]);
    }
};

#ifdef SYNERGIA_ENABLE_CUDA
#include "synergia/utils/distributed_fft2d_cuda.h"
#else
#include "synergia/utils/distributed_fft2d_fftw.h"
#endif

#endif
