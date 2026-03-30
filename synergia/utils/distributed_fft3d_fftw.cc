#include <cstring>
#include <stdexcept>
#include <thread>

#include "distributed_fft3d.h"

Distributed_fft3d::Distributed_fft3d()
    : Distributed_fft3d_base()
    , plan(nullptr)
    , inv_plan(nullptr)
    , data(nullptr)
    , workspace(nullptr)
{
    fftw_init_threads();
    fftw_mpi_init();
    fftw_plan_with_nthreads(std::thread::hardware_concurrency());
}

void
Distributed_fft3d::construct(std::array<int, 3> const& new_shape,
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

    if (new_comm.size() > new_shape[2]) {
        throw std::runtime_error(
            "Distributed_fft3d: (number of processors) must be "
            "<= shape[2]");
    }

    shape = new_shape;
    comm = new_comm;

    int padded_cplx_s0 = get_padded_shape_cplx(shape[0]);
    int padded_real_s0 = get_padded_shape_real(shape[0]);

    ptrdiff_t local_n, local_start;
    ptrdiff_t fftw_local_size = fftw_mpi_local_size_3d(
        shape[2], shape[1], shape[0], comm, &local_n, &local_start);

    int local_size_real = local_n * shape[1] * padded_real_s0;
    int local_size_cplx = local_n * shape[1] * padded_cplx_s0;

#if 0
    // MEDIUM HACK. fftw often (always?) returns a local size that is
    // impossibly small. Adjust to the smallest possible size.
    if (fftw_local_size < local_size_real)
        fftw_local_size = local_size_real;
#endif

    data = (double*)fftw_malloc(sizeof(double) * local_size_real);
    workspace =
        (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * local_size_cplx);

    plan = fftw_mpi_plan_dft_r2c_3d(
        shape[2], shape[1], shape[0], data, workspace, comm, FFTW_ESTIMATE);

    inv_plan = fftw_mpi_plan_dft_c2r_3d(
        shape[2], shape[1], shape[0], workspace, data, comm, FFTW_ESTIMATE);

    lower = local_start;
    nz = local_n;
}

void
Distributed_fft3d::transform(karray1d_dev& in, karray1d_dev& out)
{
    if (!data || !workspace)
        throw std::runtime_error(
            "Distributed_fft3d::transform() uninitialized");

    // plane size = padded_nx * ny
    int plane_real = padded_nx_real() * shape[1]; // padded_nx * ny
    int plane_cplx = padded_nx_cplx() * shape[1]; // padded_nx * ny

    memcpy((void*)data,
           (void*)&in(lower * plane_real),
           nz * plane_real * sizeof(double));

    fftw_execute(plan);

    memcpy((void*)&out(lower * plane_cplx * 2),
           (void*)(workspace),
           nz * plane_cplx * sizeof(double) * 2);
}

void
Distributed_fft3d::inv_transform(karray1d_dev& in, karray1d_dev& out)
{
    if (!data || !workspace)
        throw std::runtime_error(
            "Distributed_fft3d::transform() uninitialized");

    // plane size = padded_nx * ny
    int plane_real = padded_nx_real() * shape[1];
    int plane_cplx = padded_nx_cplx() * shape[1];

    memcpy((void*)workspace,
           (void*)&in(lower * plane_cplx * 2),
           nz * plane_cplx * sizeof(double) * 2);

    fftw_execute(inv_plan);

    memcpy((void*)&out(lower * plane_real),
           (void*)data,
           nz * plane_real * sizeof(double));
}

Distributed_fft3d::~Distributed_fft3d()
{
    if (data || workspace) {
        fftw_destroy_plan(plan);
        fftw_destroy_plan(inv_plan);
        fftw_free(data);
        fftw_free(workspace);
    }

    // fftw_mpi_cleanup();
}
