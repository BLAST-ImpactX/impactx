#ifndef DISTRIBUTED_FFT3D_CUDA_H_
#define DISTRIBUTED_FFT3D_CUDA_H_

#include <cufft.h>

class Distributed_fft3d : public Distributed_fft3d_base {

  private:
    cufftHandle plan;
    cufftHandle invplan;

  public:
    Distributed_fft3d();
    virtual ~Distributed_fft3d();

    void construct(std::array<int, 3> const& shape, Commxx const& comm);

    void transform(karray1d_dev& in, karray1d_dev& out);

    void inv_transform(karray1d_dev& in, karray1d_dev& out);
};

#endif /* DISTRIBUTED_FFT2D_H_ */
