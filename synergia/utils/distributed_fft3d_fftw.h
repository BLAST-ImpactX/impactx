#ifndef DISTRIBUTED_FFT3D_FFTW_H_
#define DISTRIBUTED_FFT3D_FFTW_H_

#include <fftw3-mpi.h>
#include <fftw3.h>

class Distributed_fft3d : public Distributed_fft3d_base {

  private:
    fftw_plan plan;
    fftw_plan inv_plan;
    double* data;
    fftw_complex* workspace;

  public:
    Distributed_fft3d();
    virtual ~Distributed_fft3d();

    void construct(std::array<int, 3> const& shape, Commxx const& comm);

    void transform(karray1d_dev& in, karray1d_dev& out);

    void inv_transform(karray1d_dev& in, karray1d_dev& out);
};

#endif /* DISTRIBUTED_FFT3D_H_ */
