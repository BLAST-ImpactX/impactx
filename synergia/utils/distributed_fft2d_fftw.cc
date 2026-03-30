#include <cstring>
#include <stdexcept>

#include "distributed_fft2d.h"

Distributed_fft2d::Distributed_fft2d()
    : Distributed_fft2d_base()
    , plan(nullptr)
    , inv_plan(nullptr)
    , data(nullptr)
    , workspace(nullptr)
{
    fftw_mpi_init();
}

void
Distributed_fft2d::construct(std::array<int, 2> const& new_shape,
                             Commxx const& new_comm)
{
    if (data || workspace) {
        fftw_destroy_plan(plan);
        fftw_destroy_plan(inv_plan);
        fftw_free(data);
        fftw_free(workspace);
    }

    plan = nullptr;
    inv_plan = nullptr;
    data = nullptr;
    workspace = nullptr;

    if (new_comm.is_null()) return;

    if (new_comm.size() / 2 >= new_shape[0] / 2) {
        throw std::runtime_error("Distributed_fft2d: (number of processors)/2 "
                                 "must be <= shape[0]/2");
    }

    shape = new_shape;
    comm = new_comm;

    ptrdiff_t local_nx, local_x_start;
    ptrdiff_t fftw_local_size = fftw_mpi_local_size_2d(
        shape[0], shape[1], comm, &local_nx, &local_x_start);

    int local_size_real = local_nx * shape[1];

    // MEDIUM HACK. fftw often (always?) returns a local size that is
    // impossibly small. Adjust to the smallest possible size.
    if (fftw_local_size < local_size_real) fftw_local_size = local_size_real;

    data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftw_local_size);
    workspace =
        (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftw_local_size);

    plan = fftw_mpi_plan_dft_2d(
        shape[0], shape[1], data, workspace, comm, FFTW_FORWARD, FFTW_ESTIMATE);

    inv_plan = fftw_mpi_plan_dft_2d(shape[0],
                                    shape[1],
                                    workspace,
                                    data,
                                    comm,
                                    FFTW_BACKWARD,
                                    FFTW_ESTIMATE);

    lower = local_x_start;
    nx = local_nx;
}

void
Distributed_fft2d::transform(karray1d_dev& in, karray1d_dev& out)
{
    if (!data || !workspace)
        throw std::runtime_error(
            "Distributed_fft2d::transform() uninitialized");

    memcpy((void*)data,
           (void*)&in(lower * shape[1] * 2),
           nx * shape[1] * sizeof(double) * 2);

    fftw_execute(plan);

    memcpy((void*)&out(lower * shape[1] * 2),
           (void*)(workspace),
           nx * shape[1] * sizeof(double) * 2);
}

void
Distributed_fft2d::inv_transform(karray1d_dev& in, karray1d_dev& out)
{
    if (!data || !workspace)
        throw std::runtime_error(
            "Distributed_fft2d::transform() uninitialized");

    memcpy((void*)workspace,
           (void*)&in(lower * shape[1] * 2),
           nx * shape[1] * sizeof(double) * 2);

    fftw_execute(inv_plan);

    memcpy((void*)&out(lower * shape[1] * 2),
           (void*)data,
           nx * shape[1] * sizeof(double) * 2);
}

Distributed_fft2d::~Distributed_fft2d()
{
    if (data || workspace) {
        fftw_destroy_plan(plan);
        fftw_destroy_plan(inv_plan);
        fftw_free(data);
        fftw_free(workspace);
    }

    // fftw_mpi_cleanup();
}
