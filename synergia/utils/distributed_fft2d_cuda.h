#ifndef DISTRIBUTED_FFT2D_CUDA_H_
#define DISTRIBUTED_FFT2D_CUDA_H_

#include <cufft.h>

class Distributed_fft2d : public Distributed_fft2d_base {

private:
  cufftHandle plan;

public:
  Distributed_fft2d();
  virtual ~Distributed_fft2d();

  void construct(std::array<int, 2> const& shape, Commxx const& comm);

  void transform(karray1d_dev& in, karray1d_dev& out);

  void inv_transform(karray1d_dev& in, karray1d_dev& out);
};

#endif /* DISTRIBUTED_FFT2D_H_ */
