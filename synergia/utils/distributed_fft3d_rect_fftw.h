#ifndef DISTRIBUTED_FFT3D_RECT_FFTW_H_
#define DISTRIBUTED_FFT3D_RECT_FFTW_H_

#include <fftw3-mpi.h>
#include <fftw3.h>

class Distributed_fft3d_rect : public Distributed_fft3d_rect_base {

private:
  fftw_plan plan_xy, plan_z;
  fftw_plan inv_plan_xy, inv_plan_z;

  double* data;
  fftw_complex* workspace;

public:
  Distributed_fft3d_rect();
  virtual ~Distributed_fft3d_rect();

  void construct(std::array<int, 3> const& shape, Commxx const& comm) override;

  void transform(karray1d_dev& in, karray1d_dev& out) override;

  void inv_transform(karray1d_dev& in, karray1d_dev& out) override;
};

#endif /* DISTRIBUTED_FFT3D_RECT_FFTW_H_ */
