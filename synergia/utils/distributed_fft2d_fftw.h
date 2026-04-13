#ifndef DISTRIBUTED_FFT2D_FFTW_H_
#define DISTRIBUTED_FFT2D_FFTW_H_

#include <fftw3-mpi.h>
#include <fftw3.h>

class Distributed_fft2d : public Distributed_fft2d_base {

  private:
    fftw_plan plan;
    fftw_plan inv_plan;
    fftw_complex* data;
    fftw_complex* workspace;

  public:
    Distributed_fft2d();
    virtual ~Distributed_fft2d();

    void construct(std::array<int, 2> const& shape, Commxx const& comm);

    void transform(karray1d_dev& in, karray1d_dev& out);

    void inv_transform(karray1d_dev& in, karray1d_dev& out);
};

#endif /* DISTRIBUTED_FFT2D_H_ */
