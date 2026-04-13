#include "distributed_fft3d_rect.h"
#include <Kokkos_MathematicalConstants.hpp>
#include <cstring>
#include <stdexcept>

namespace {
    void
    print(karray1d_dev const& arr)
    {
        karray1d_hst harr = Kokkos::create_mirror_view(arr);
        Kokkos::deep_copy(harr, arr);

        std::cout << harr.label() << "\n  ";

        for (int i = 0; i < harr.extent(0); ++i)
            std::cout << harr(i) << ", ";

        std::cout << "\n";
    }

    struct alg_zeroer {
        karray1d_dev arr;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            arr(i) = 0.0;
        }
    };

    void
    zero(karray1d_dev const& arr)
    {
        alg_zeroer alg{arr};
        Kokkos::parallel_for(arr.extent(0), alg);
    }

    struct alg_pad_x {
        // in: (x, y, z) of (gx, gy, gz) double
        // out: (z, y, x) of (gz, gy, gx*2) complex
        karray1d_dev in;
        karray1d_dev out;

        int dgx, gy, gz;
        double igygz;
        double igz;

        alg_pad_x(karray1d_dev const& in,
                  karray1d_dev const& out,
                  std::array<int, 3> const& g)
            : in(in)
            , out(out)
            , dgx(2 * g[0])
            , gy(g[1])
            , gz(g[2])
            , igygz(1.0 / (gy * gz))
            , igz(1.0 / gz)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            int ix = i * igygz;
            int iy = (i - ix * gy * gz) * igz;
            int iz = i - ix * gy * gz - iy * gz;

            // { 1, 2, 3 } ->
            // { 1, 0, 2, 0, 3, 0, -3, 0, -2, 0, -1, 0 }
            out((iz * gy * dgx + iy * dgx + ix) * 2) = in(i);
            out((iz * gy * dgx + iy * dgx + dgx - 1 - ix) * 2) = -in(i);
        }
    };

    struct alg_pad_y {
        // in: (z, y, x) of (gz, gy, gx*2) complex
        // out: (z, x, y) of (gz, gx, gy*2) complex
        karray1d_dev in;
        karray1d_dev out;

        int gx, gy;
        int dgx, dgy;
        double igxgy;
        double igx;

        alg_pad_y(karray1d_dev const& in,
                  karray1d_dev const& out,
                  std::array<int, 3> const& g)
            : in(in)
            , out(out)
            , gx(g[0])
            , gy(g[1])
            , dgx(2 * g[0])
            , dgy(2 * g[1])
            , igxgy(1.0 / (gy * gx))
            , igx(1.0 / gx)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            int iz = i * igxgy;
            int iy = (i - iz * gy * gx) * igx;
            int ix = i - iz * gy * gx - iy * gx;

            // { 1, 2, 3 } ->
            // { 1, 0, 2, 0, 3, 0, -3, 0, -2, 0, -1, 0 }
            out((iz * gx * dgy + ix * dgy + iy) * 2) =
                -in((iz * gy * dgx + iy * dgx + ix + 1) * 2 + 1);

            out((iz * gx * dgy + ix * dgy + dgy - 1 - iy) * 2) =
                in((iz * gy * dgx + iy * dgx + ix + 1) * 2 + 1);
        }
    };

    struct alg_pad_z {
        // in: (z, x, y) of (gz, gx, gy*2) complex
        // out: (x, y, z) of (gx, gy, gz) double
        karray1d_dev in;
        karray1d_dev out;

        int gx, gy, gz;
        int dgy;
        double igxgy;
        double igx;

        alg_pad_z(karray1d_dev const& in,
                  karray1d_dev const& out,
                  std::array<int, 3> const& g)
            : in(in)
            , out(out)
            , gx(g[0])
            , gy(g[1])
            , gz(g[2])
            , dgy(g[1] * 2)
            , igxgy(1.0 / (gx * gy))
            , igx(1.0 / gx)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            int iz = i * igxgy;
            int iy = (i - iz * gy * gx) * igx;
            int ix = i - iz * gy * gx - iy * gx;

            out(ix * gy * gz + iy * gz + iz) =
                -in((iz * gx * dgy + ix * dgy + iy + 1) * 2 + 1);
        }
    };

    struct alg_inv_pad_y {
        // in: (x, y, z) of (gx, gy, gz) double
        // out: (z, x, y) of (gz, gx, gy*2) complex
        karray1d_dev in;
        karray1d_dev out;

        int gx, gy, gz;
        int dgy;
        double igygz;
        double igz;

        alg_inv_pad_y(karray1d_dev const& in,
                      karray1d_dev const& out,
                      std::array<int, 3> const& g)
            : in(in)
            , out(out)
            , gx(g[0])
            , gy(g[1])
            , gz(g[2])
            , dgy(2 * g[1])
            , igygz(1.0 / (gy * gz))
            , igz(1.0 / gz)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            int ix = i * igygz;
            int iy = (i - ix * gy * gz) * igz;
            int iz = i - ix * gy * gz - iy * gz;

            // { 1, 2, 3 } ->
            // { 0, 0, 0, 1, 0, 2, 0, 3, 0, 2, 0, 1 }
            int out_idx = iz * gx * dgy + ix * dgy + iy + 1;

            out(out_idx * 2 + 0) = 0.0;
            out(out_idx * 2 + 1) = in(ix * gy * gz + iy * gz + iz);

            // mirror
            if (iy == gy - 1) {
                out_idx = iz * gx * dgy + ix * dgy + 0;
                out(out_idx * 2 + 0) = 0.0;
                out(out_idx * 2 + 1) = 0.0;
            } else {
                out_idx = iz * gx * dgy + ix * dgy + dgy - 1 - iy;

                out(out_idx * 2 + 0) = 0.0;
                out(out_idx * 2 + 1) = in(ix * gy * gz + iy * gz + iz);
            }
        }
    };

    struct alg_inv_pad_x {
        // in: (z, x, y) of (gz, gx, gy*2) complex
        // out: (z, y, x) of (gz, gy, gx*2) complex
        karray1d_dev in;
        karray1d_dev out;

        int gx, gy;
        int dgx, dgy;
        double igxgy;
        double igy;

        alg_inv_pad_x(karray1d_dev const& in,
                      karray1d_dev const& out,
                      std::array<int, 3> const& g)
            : in(in)
            , out(out)
            , gx(g[0])
            , gy(g[1])
            , dgx(2 * g[0])
            , dgy(2 * g[1])
            , igxgy(1.0 / (gx * gy))
            , igy(1.0 / gy)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            int iz = i * igxgy;
            int ix = (i - iz * gx * gy) * igy;
            int iy = i - iz * gx * gy - ix * gy;

            // { 1, 2, 3 } ->
            // { 0, 0, 0, 1, 0, 2, 0, 3, 0, 2, 0, 1 }
            int out_idx = iz * gy * dgx + iy * dgx + ix + 1;
            int in_idx = iz * gx * dgy + ix * dgy + iy;

            out(out_idx * 2 + 0) = 0.0;
            out(out_idx * 2 + 1) = -in(in_idx * 2 + 0);

            // mirror
            out_idx = (ix == gx - 1) ?
                          (iz * gy * dgx + iy * dgx + 0) :
                          (iz * gy * dgx + iy * dgx + dgx - 1 - ix);

            out(out_idx * 2 + 0) = 0.0;
            out(out_idx * 2 + 1) = -in(in_idx * 2 + 0);
        }
    };

    struct alg_inv_pad_z {
        // in: (z, y, x) of (gz, gy, gx*2) complex
        // out: (x, y, z) of (gx, gy, gz) double
        karray1d_dev in;
        karray1d_dev out;

        int gx, gy, gz;
        int dgx;
        double igxgy;
        double igx;

        alg_inv_pad_z(karray1d_dev const& in,
                      karray1d_dev const& out,
                      std::array<int, 3> const& g)
            : in(in)
            , out(out)
            , gx(g[0])
            , gy(g[1])
            , gz(g[2])
            , dgx(g[0] * 2)
            , igxgy(1.0 / (gx * gy))
            , igx(1.0 / gx)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            int iz = i * igxgy;
            int iy = (i - iz * gy * gx) * igx;
            int ix = i - iz * gy * gx - iy * gx;

            out(ix * gy * gz + iy * gz + iz) =
                -in((iz * gy * dgx + iy * dgx + ix) * 2 + 0);
        }
    };

    struct alg_shift_forward {
        karray1d_dev in;
        int gx;
        double igx;

        alg_shift_forward(karray1d_dev const& in, int gx)
            : in(in), gx(gx), igx(1.0 / gx)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            const int idx = i * igx;
            const int ix = i - idx * gx;

            double dk = (ix + 1) * Kokkos::numbers::pi_v<double> * igx * 0.5;

            double imag = -in(idx * gx * 4 + (ix + 1) * 2 + 0) * sin(dk) +
                          in(idx * gx * 4 + (ix + 1) * 2 + 1) * cos(dk);

            in(idx * gx * 4 + (ix + 1) * 2 + 1) = imag;
        }
    };

    struct alg_shift_backward {
        karray1d_dev in;
        int gx;
        double igx;

        alg_shift_backward(karray1d_dev const& in, int gx)
            : in(in), gx(gx), igx(1.0 / gx)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            const int idx = i * igx;
            const int ix = i - idx * gx;

            double dk = ix * Kokkos::numbers::pi_v<double> * igx;

            double real = in((idx * gx + ix) * 2 + 0) * cos(dk) -
                          in((idx * gx + ix) * 2 + 1) * sin(dk);

            double imag = in((idx * gx + ix) * 2 + 0) * sin(dk) +
                          in((idx * gx + ix) * 2 + 1) * cos(dk);

            in((idx * gx + ix) * 2 + 0) = real;
            in((idx * gx + ix) * 2 + 1) = imag;
        }
    };
}

Distributed_fft3d_rect::Distributed_fft3d_rect()
    : Distributed_fft3d_rect_base()
    , data1()
    , data2()
    , plan_x()
    , plan_y()
    , plan_z()
    , inv_plan_x()
    , inv_plan_y()
    , inv_plan_z()
{}

void
Distributed_fft3d_rect::construct(std::array<int, 3> const& new_shape,
                                  Commxx const& new_comm)
{
    cufftDestroy(plan_x);
    cufftDestroy(plan_y);
    cufftDestroy(plan_z);

    cufftDestroy(inv_plan_x);
    cufftDestroy(inv_plan_y);
    cufftDestroy(inv_plan_z);

    if (new_comm.is_null()) return;

    if (new_comm.size() > new_shape[0]) {
        throw std::runtime_error(
            "Distributed_fft3d_rect: (number of processors) must be "
            "<= shape[0]");
    }

    shape = new_shape;
    comm = new_comm;

    lower = 0;
    nx = shape[0];

    // (x*2, y, z) of a complex grid for batched DST along x
    // (x, y*2, z) of a complex grid for batched DST along y
    data1 = karray1d_dev("data1", shape[0] * 4 * shape[1] * shape[2]);
    data2 = karray1d_dev("data2", shape[0] * shape[1] * 4 * shape[2]);

    int rank = 1;
    int istride = 1;
    int ostride = 1;

    // plan x
    int n_x[] = {shape[0] * 2};
    int inembed_x[] = {shape[0] * 2};
    int onembed_x[] = {shape[0] * 2};
    int idist_x = shape[0] * 2;
    int odist_x = shape[0] * 2;

    cufftPlanMany(&plan_x,
                  rank,
                  n_x,
                  inembed_x,
                  istride,
                  idist_x,
                  onembed_x,
                  ostride,
                  odist_x,
                  CUFFT_Z2Z,
                  shape[1] * shape[2]);

    cufftPlanMany(&inv_plan_x,
                  rank,
                  n_x,
                  onembed_x,
                  ostride,
                  odist_x,
                  inembed_x,
                  istride,
                  idist_x,
                  CUFFT_Z2Z,
                  shape[1] * shape[2]);

    // plan y
    int n_y[] = {shape[1] * 2};
    int inembed_y[] = {shape[1] * 2};
    int onembed_y[] = {shape[1] * 2};
    int idist_y = shape[1] * 2;
    int odist_y = shape[1] * 2;

    cufftPlanMany(&plan_y,
                  rank,
                  n_y,
                  inembed_y,
                  istride,
                  idist_y,
                  onembed_y,
                  ostride,
                  odist_y,
                  CUFFT_Z2Z,
                  shape[0] * shape[2]);

    cufftPlanMany(&inv_plan_y,
                  rank,
                  n_y,
                  onembed_y,
                  ostride,
                  odist_y,
                  inembed_y,
                  istride,
                  idist_y,
                  CUFFT_Z2Z,
                  shape[0] * shape[2]);

    // plan z
    int n_z[] = {shape[2]};
    int inembed_z[] = {shape[2]};
    int onembed_z[] = {shape[2]};
    int idist_z = shape[2];
    int odist_z = shape[2] / 2 + 1;

    cufftPlanMany(&plan_z,
                  rank,
                  n_z,
                  inembed_z,
                  istride,
                  idist_z,
                  onembed_z,
                  ostride,
                  odist_z,
                  CUFFT_D2Z,
                  shape[0] * shape[1]);

    cufftPlanMany(&inv_plan_z,
                  rank,
                  n_z,
                  onembed_z,
                  ostride,
                  odist_z,
                  inembed_z,
                  istride,
                  idist_z,
                  CUFFT_Z2D,
                  shape[0] * shape[1]);
}

void
Distributed_fft3d_rect::transform(karray1d_dev& in, karray1d_dev& out)
{
    const int gx = shape[0];
    const int gy = shape[1];
    const int gz = shape[2];

    zero(data1);
    zero(data2);

    // pad x
    alg_pad_x padx(in, data1, shape);
    Kokkos::parallel_for(gx * gy * gz, padx);

    // dft x
    auto res = cufftExecZ2Z(plan_x,
                            (cufftDoubleComplex*)data1.data(),
                            (cufftDoubleComplex*)data1.data(),
                            CUFFT_FORWARD);

    // shift x
    alg_shift_forward shiftx(data1, gx);
    Kokkos::parallel_for(gx * gy * gz, shiftx);

    // pad y
    alg_pad_y pady(data1, data2, shape);
    Kokkos::parallel_for(gx * gy * gz, pady);

    // dft y
    res = cufftExecZ2Z(plan_y,
                       (cufftDoubleComplex*)data2.data(),
                       (cufftDoubleComplex*)data2.data(),
                       CUFFT_FORWARD);

    // shift y
    alg_shift_forward shifty(data2, gy);
    Kokkos::parallel_for(gx * gy * gz, shifty);

    // pad z
    alg_pad_z padz(data2, data1, shape);
    Kokkos::parallel_for(gx * gy * gz, padz);

    // dft z
    res = cufftExecD2Z(plan_z,
                       (cufftDoubleReal*)data1.data(),
                       (cufftDoubleComplex*)out.data());
}

void
Distributed_fft3d_rect::inv_transform(karray1d_dev& in, karray1d_dev& out)
{
    const int gx = shape[0];
    const int gy = shape[1];
    const int gz = shape[2];

    // dft z
    auto res = cufftExecZ2D(inv_plan_z,
                            (cufftDoubleComplex*)in.data(),
                            (cufftDoubleReal*)data1.data());

    // pady
    alg_inv_pad_y pady(data1, data2, shape);
    Kokkos::parallel_for(gx * gy * gz, pady);

    // shifty
    alg_shift_backward shifty(data2, gy * 2);
    Kokkos::parallel_for(gx * gy * gz * 2, shifty);

    // z2z along y
    res = cufftExecZ2Z(inv_plan_y,
                       (cufftDoubleComplex*)data2.data(),
                       (cufftDoubleComplex*)data2.data(),
                       CUFFT_INVERSE);

    // padx
    alg_inv_pad_x padx(data2, data1, shape);
    Kokkos::parallel_for(gx * gy * gz, padx);

    // shiftx
    alg_shift_backward shiftx(data1, gx * 2);
    Kokkos::parallel_for(gx * gy * gz * 2, shiftx);

    // z2z along x
    res = cufftExecZ2Z(inv_plan_x,
                       (cufftDoubleComplex*)data1.data(),
                       (cufftDoubleComplex*)data1.data(),
                       CUFFT_INVERSE);

    // re-arrange for final array
    alg_inv_pad_z padz(data1, out, shape);
    Kokkos::parallel_for(gx * gy * gz, padz);
}

Distributed_fft3d_rect::~Distributed_fft3d_rect()
{
    cufftDestroy(plan_x);
    cufftDestroy(plan_y);
    cufftDestroy(plan_z);

    cufftDestroy(inv_plan_x);
    cufftDestroy(inv_plan_y);
    cufftDestroy(inv_plan_z);
}
