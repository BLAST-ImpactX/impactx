#include "distributed_fft2d.h"
#include <stdexcept>

Distributed_fft2d::Distributed_fft2d() : Distributed_fft2d_base(), plan() {}

void
Distributed_fft2d::construct(std::array<int, 2> const& new_shape,
                             Commxx const& new_comm)
{
  cufftDestroy(plan);

  if (new_comm.is_null()) return;

  if (new_comm.size() != 1) {
    throw std::runtime_error("Distributed_fft2d: number of processor must be 1 "
                             "for CUDA implementation");
  }

  shape = new_shape;
  comm = new_comm;

  lower = 0;
  nx = shape[0];

  cufftPlan2d(&plan, shape[0], shape[1], CUFFT_Z2Z);
}

void
Distributed_fft2d::transform(karray1d_dev& in, karray1d_dev& out)
{
  cufftExecZ2Z(plan,
               (cufftDoubleComplex*)in.data(),
               (cufftDoubleComplex*)out.data(),
               CUFFT_FORWARD);
}

void
Distributed_fft2d::inv_transform(karray1d_dev& in, karray1d_dev& out)
{
  cufftExecZ2Z(plan,
               (cufftDoubleComplex*)in.data(),
               (cufftDoubleComplex*)out.data(),
               CUFFT_INVERSE);
}

Distributed_fft2d::~Distributed_fft2d()
{
  cufftDestroy(plan);
}
