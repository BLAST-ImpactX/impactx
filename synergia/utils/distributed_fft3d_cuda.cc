#include "distributed_fft3d.h"
#include <cstring>
#include <stdexcept>

Distributed_fft3d::Distributed_fft3d()
    : Distributed_fft3d_base(), plan(), invplan()
{}

void
Distributed_fft3d::construct(std::array<int, 3> const& new_shape,
                             Commxx const& new_comm)
{
    cufftDestroy(plan);

    if (new_comm.is_null()) return;

    if (new_comm.size() != 1) {
        throw std::runtime_error(
            "Distributed_fft3d: number of processor must be 1 "
            "for CUDA implementation");
    }

    shape = new_shape;
    comm = new_comm;

    lower = 0;
    nz = shape[2];

    auto res = cufftPlan3d(&plan, shape[2], shape[1], shape[0], CUFFT_D2Z);
    res = cufftPlan3d(&invplan, shape[2], shape[1], shape[0], CUFFT_Z2D);
}

void
Distributed_fft3d::transform(karray1d_dev& in, karray1d_dev& out)
{
    auto res = cufftExecD2Z(
        plan, (cufftDoubleReal*)in.data(), (cufftDoubleComplex*)out.data());
}

void
Distributed_fft3d::inv_transform(karray1d_dev& in, karray1d_dev& out)
{
    cufftExecZ2D(
        invplan, (cufftDoubleComplex*)in.data(), (cufftDoubleReal*)out.data());
}

Distributed_fft3d::~Distributed_fft3d()
{
    cufftDestroy(plan);
    cufftDestroy(invplan);
}
