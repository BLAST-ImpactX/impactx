#include "distributed_fft3d_rect.h"
#include <cstring>
#include <stdexcept>
#include <thread>

Distributed_fft3d_rect::Distributed_fft3d_rect()
  : Distributed_fft3d_rect_base()
  , plan_xy(nullptr)
  , plan_z(nullptr)
  , inv_plan_xy(nullptr)
  , inv_plan_z(nullptr)
  , data(nullptr)
  , workspace(nullptr)
{
  fftw_init_threads();
  fftw_mpi_init();
  fftw_plan_with_nthreads(std::thread::hardware_concurrency());
}

void
Distributed_fft3d_rect::construct(std::array<int, 3> const& new_shape,
                                  Commxx const& new_comm)
{
  if (data || workspace) {
    fftw_destroy_plan(plan_xy);
    fftw_destroy_plan(plan_z);
    fftw_destroy_plan(inv_plan_xy);
    fftw_destroy_plan(inv_plan_z);
    fftw_free(data);
    fftw_free(workspace);
  }

  plan_xy = nullptr;
  plan_z = nullptr;
  inv_plan_xy = nullptr;
  inv_plan_z = nullptr;

  data = nullptr;
  workspace = nullptr;

  if (new_comm.is_null()) return;

  if (new_comm.size() > new_shape[0]) {
    throw std::runtime_error(
      "Distributed_fft3d_rect: (number of processors) must be "
      "<= shape[0]");
  }

  shape = new_shape;
  comm = new_comm;

  const ptrdiff_t ndim_xy[] = {shape[0], shape[1]};
  const int ndim_z[] = {shape[2]};
  ptrdiff_t local_n, local_start;

  ptrdiff_t fftw_local_size = fftw_mpi_local_size_many(
    2, ndim_xy, shape[2], FFTW_MPI_DEFAULT_BLOCK, comm, &local_n, &local_start);

  lower = local_start;
  nx = local_n;

  // plans for x-y DST
  data = (double*)fftw_malloc(sizeof(double) * fftw_local_size);

  fftw_r2r_kind kind_direct[] = {FFTW_RODFT10, FFTW_RODFT10};
  plan_xy = fftw_mpi_plan_many_r2r(2,
                                   ndim_xy,
                                   shape[2],
                                   FFTW_MPI_DEFAULT_BLOCK,
                                   FFTW_MPI_DEFAULT_BLOCK,
                                   data,
                                   data,
                                   comm,
                                   kind_direct,
                                   // FFTW_EXHAUSTIVE);
                                   FFTW_ESTIMATE);

  fftw_r2r_kind kind_inv[] = {FFTW_RODFT01, FFTW_RODFT01};
  inv_plan_xy = fftw_mpi_plan_many_r2r(2,
                                       ndim_xy,
                                       shape[2],
                                       FFTW_MPI_DEFAULT_BLOCK,
                                       FFTW_MPI_DEFAULT_BLOCK,
                                       data,
                                       data,
                                       comm,
                                       kind_inv,
                                       // FFTW_EXHAUSTIVE);
                                       FFTW_ESTIMATE);

  // plans for z DFT (r to c)
  int padded_cplx_s2 = get_padded_shape_cplx(shape[2]);
  int local_size_cplx = local_n * shape[1] * padded_cplx_s2;

  workspace =
    (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * local_size_cplx);

  plan_z = fftw_plan_many_dft_r2c(1,
                                  ndim_z,
                                  nx * shape[1], // rank, ndim, howmany
                                  data,
                                  NULL,
                                  1,
                                  shape[2], // in, inembed, stride, dist
                                  workspace,
                                  NULL,
                                  1,
                                  padded_cplx_s2,
                                  FFTW_ESTIMATE);

  inv_plan_z = fftw_plan_many_dft_c2r(1,
                                      ndim_z,
                                      nx * shape[1], // rank, ndim, howmany
                                      workspace,
                                      NULL,
                                      1,
                                      padded_cplx_s2,
                                      data,
                                      NULL,
                                      1,
                                      shape[2], // in, inembed, stride, dist
                                      FFTW_ESTIMATE);
}

void
Distributed_fft3d_rect::transform(karray1d_dev& in, karray1d_dev& out)
{
  if (!data || !workspace)
    throw std::runtime_error(
      "Distributed_fft3d_rect::transform() uninitialized");

  int plane_yz_real = shape[1] * shape[2];
  int plane_yz_cplx = shape[1] * (shape[2] / 2 + 1);

  memcpy((void*)data,
         (void*)&in(lower * plane_yz_real),
         nx * plane_yz_real * sizeof(double));

  fftw_execute(plan_xy);
  fftw_execute(plan_z);

  memcpy((void*)&out(lower * plane_yz_cplx * 2),
         (void*)workspace,
         nx * plane_yz_cplx * sizeof(double) * 2);
}

void
Distributed_fft3d_rect::inv_transform(karray1d_dev& in, karray1d_dev& out)
{
  if (!data || !workspace)
    throw std::runtime_error(
      "Distributed_fft3d_rect::transform() uninitialized");

  int plane_yz_real = shape[1] * shape[2];
  int plane_yz_cplx = shape[1] * (shape[2] / 2 + 1);

  memcpy((void*)workspace,
         (void*)&in(lower * plane_yz_cplx * 2),
         nx * plane_yz_cplx * sizeof(double) * 2);

  fftw_execute(inv_plan_z);
  fftw_execute(inv_plan_xy);

  memcpy((void*)&out(lower * plane_yz_real),
         (void*)data,
         nx * plane_yz_real * sizeof(double));
}

Distributed_fft3d_rect::~Distributed_fft3d_rect()
{
  if (data || workspace) {
    fftw_destroy_plan(plan_xy);
    fftw_destroy_plan(plan_z);
    fftw_destroy_plan(inv_plan_xy);
    fftw_destroy_plan(inv_plan_z);

    fftw_free(data);
    fftw_free(workspace);
  }

  // fftw_mpi_cleanup();
}
