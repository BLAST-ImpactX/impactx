#ifndef DISTRIBUTED_FFT3D_RECT_CUDA_H_
#define DISTRIBUTED_FFT3D_RECT_CUDA_H_

#include <cufft.h>

class Distributed_fft3d_rect : public Distributed_fft3d_rect_base {

private:
  karray1d_dev data1;
  karray1d_dev data2;

  cufftHandle plan_x;
  cufftHandle plan_y;
  cufftHandle plan_z;

  cufftHandle inv_plan_x;
  cufftHandle inv_plan_y;
  cufftHandle inv_plan_z;

public:
  Distributed_fft3d_rect();
  virtual ~Distributed_fft3d_rect();

  void construct(std::array<int, 3> const& shape, Commxx const& comm) override;

  void transform(karray1d_dev& in, karray1d_dev& out) override;

  void inv_transform(karray1d_dev& in, karray1d_dev& out) override;
};

#endif /* DISTRIBUTED_FFT3D_RECT_CUDA_H_ */
