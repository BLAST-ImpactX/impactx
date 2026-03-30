#include <Kokkos_ScatterView.hpp>

#include "deposit.h"
#include "utils.h"

#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/distributed_fft3d.h"

namespace deposit_impl {
    using scatter_t =
        Kokkos::Experimental::ScatterView<double*, Kokkos::LayoutLeft>;

    KOKKOS_INLINE_FUNCTION
    bool
    ingrid(int i, int g)
    {
        return i >= 0 && i < g;
    }

    KOKKOS_INLINE_FUNCTION
    bool
    ingrid(int ix, int iy, int iz, int gx, int gy, int gz)
    {
        // return ix>=0 && ix<gx && iy>=0 && iy<gy && iz>=0 && iz<gz;

        // exclude edges
        return ix > 0 && ix < gx - 1 && iy > 0 && iy < gy - 1 && iz > 0 &&
               iz < gz - 1;
    }

    struct rho_reducer {
        typedef double value_type[];

        const int value_count;
        ConstParticles p;
        karray2d_dev bin;
        int gx, gy, gz;
        double ihx, ihy, ihz;
        double lx, ly, lz;
        double w0;

        rho_reducer(ConstParticles const& p,
                    karray2d_dev const& bin,
                    std::array<int, 3> const& g,
                    std::array<double, 3> const& h,
                    std::array<double, 3> const& l,
                    double w0)
            : value_count(g[0] * g[1] + g[2])
            , p(p)
            , bin(bin)
            , gx(g[0])
            , gy(g[1])
            , gz(g[2])
            , ihx(1.0 / h[0])
            , ihy(1.0 / h[1])
            , ihz(1.0 / h[2])
            , lx(l[0])
            , ly(l[1])
            , lz(l[2])
            , w0(w0)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i, value_type sum) const
        {
            int ix, iy, iz;
            double offx, offy, offz;

            get_leftmost_indices_offset(p(i, 0), lx, ihx, ix, offx);
            get_leftmost_indices_offset(p(i, 2), ly, ihy, iy, offy);
            get_leftmost_indices_offset(p(i, 4), lz, ihz, iz, offz);

            bin(i, 0) = ix;
            bin(i, 1) = offx;
            bin(i, 2) = iy;
            bin(i, 3) = offy;
            bin(i, 4) = iz;
            bin(i, 5) = offz;

            int cellz1 = iz;
            int cellz2 = cellz1 + 1;

            if (ingrid(cellz1, gz)) sum[gx * gy + cellz1] += (1.0 - offz) * ihz;
            if (ingrid(cellz2, gz)) sum[gx * gy + cellz2] += offz * ihz;

            if (ix < 0 || ix > gx - 1 || iy < 0 || iy > gy - 1) return;

            int cellx1, cellx2, celly1, celly2;
            cellx1 = ix;
            cellx2 = ix + 1;
            celly1 = iy;
            celly2 = iy + 1;

            double aoffx, aoffy;
            aoffx = 1. - offx;
            aoffy = 1. - offy;

            sum[cellx1 * gy + celly1] += w0 * aoffx * aoffy;
            sum[cellx1 * gx + celly2] += w0 * aoffx * offy;
            sum[cellx2 * gy + celly1] += w0 * offx * aoffy;
            sum[cellx2 * gx + celly2] += w0 * offx * offy;
        }
    };

    // move the rho data from double array to complex array
    struct rho_mover {
        karray1d_dev r0;
        karray1d_dev r1;
        int gx, gy;

        rho_mover(karray1d_dev const& rho_double,
                  karray1d_dev const& rho_complex,
                  std::array<int, 3> const& g)
            : r0(rho_double), r1(rho_complex), gx(g[0]), gy(g[1])
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            if (i < gx * gy) {
                r1(i * 2) = r0(i);
            } else {
                int z = i - gx * gy;
                r1(gx * gy * 2 + z) = r0(gx * gy + z);
            }
        }
    };

    // use atomic add
    struct atomic_rho_reducer {
        ConstParticles p;
        karray1d_atomic_dev rho;
        karray2d_dev bin;
        int gx, gy, gz;
        double ihx, ihy, ihz;
        double lx, ly, lz;
        double w0;

        atomic_rho_reducer(ConstParticles const& p,
                           karray1d_atomic_dev const& rho,
                           karray2d_dev const& bin,
                           std::array<int, 3> const& g,
                           std::array<double, 3> const& h,
                           std::array<double, 3> const& l,
                           double w0)
            : p(p)
            , rho(rho)
            , bin(bin)
            , gx(g[0])
            , gy(g[1])
            , gz(g[2])
            , ihx(1.0 / h[0])
            , ihy(1.0 / h[1])
            , ihz(1.0 / h[2])
            , lx(l[0])
            , ly(l[1])
            , lz(l[2])
            , w0(w0)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            int ix, iy, iz;
            double offx, offy, offz;

            get_leftmost_indices_offset(p(i, 0), lx, ihx, ix, offx);
            get_leftmost_indices_offset(p(i, 2), ly, ihy, iy, offy);
            get_leftmost_indices_offset(p(i, 4), lz, ihz, iz, offz);

            bin(i, 0) = ix;
            bin(i, 1) = offx;
            bin(i, 2) = iy;
            bin(i, 3) = offy;
            bin(i, 4) = iz;
            bin(i, 5) = offz;

            int cellz1 = iz;
            int cellz2 = cellz1 + 1;

            if (cellz1 >= 0 && cellz1 < gz)
                rho(gx * gy * 2 + cellz1) += (1.0 - offz) * ihz;
            if (cellz2 >= 0 && cellz2 < gz)
                rho(gx * gy * 2 + cellz2) += offz * ihz;

            if (ix < 0 || ix > gx - 1 || iy < 0 || iy > gy - 1) return;

            int cellx1, cellx2, celly1, celly2;
            cellx1 = ix;
            cellx2 = ix + 1;
            celly1 = iy;
            celly2 = iy + 1;

            double aoffx, aoffy;
            aoffx = 1. - offx;
            aoffy = 1. - offy;

            rho((cellx1 * gy + celly1) * 2) += w0 * aoffx * aoffy;
            rho((cellx1 * gx + celly2) * 2) += w0 * aoffx * offy;
            rho((cellx2 * gy + celly1) * 2) += w0 * offx * aoffy;
            rho((cellx2 * gx + celly2) * 2) += w0 * offx * offy;
        }
    };

    struct rho_zeroer {
        karray1d_dev rho;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            rho(i) = 0.0;
        }
    };

    // use scatter view
    struct sv_rho_reducer {
        ConstParticles p;
        ConstParticleMasks masks;

        Kokkos::Experimental::ScatterView<double*, Kokkos::LayoutLeft> scatter;
        karray2d_dev bin;
        int gx, gy, gz;
        double ihx, ihy, ihz;
        double lx, ly, lz;
        double w0;

        sv_rho_reducer(
            ConstParticles const& p,
            ConstParticleMasks const& masks,
            Kokkos::Experimental::ScatterView<double*,
                                              Kokkos::LayoutLeft> const&
                scatter,
            karray2d_dev const& bin,
            std::array<int, 3> const& g,
            std::array<double, 3> const& h,
            std::array<double, 3> const& l,
            double w0)
            : p(p)
            , masks(masks)
            , scatter(scatter)
            , bin(bin)
            , gx(g[0])
            , gy(g[1])
            , gz(g[2])
            , ihx(1.0 / h[0])
            , ihy(1.0 / h[1])
            , ihz(1.0 / h[2])
            , lx(l[0])
            , ly(l[1])
            , lz(l[2])
            , w0(w0)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            if (masks(i)) {
                auto access = scatter.access();

                int ix, iy, iz;
                double offx, offy, offz;

                get_leftmost_indices_offset(p(i, 0), lx, ihx, ix, offx);
                get_leftmost_indices_offset(p(i, 2), ly, ihy, iy, offy);
                get_leftmost_indices_offset(p(i, 4), lz, ihz, iz, offz);

                bin(i, 0) = ix;
                bin(i, 1) = offx;
                bin(i, 2) = iy;
                bin(i, 3) = offy;
                bin(i, 4) = iz;
                bin(i, 5) = offz;

                int cellz1 = iz;
                int cellz2 = cellz1 + 1;

                if (cellz1 >= 0 && cellz1 < gz)
                    access(gx * gy * 2 + cellz1) += (1.0 - offz) * ihz;
                if (cellz2 >= 0 && cellz2 < gz)
                    access(gx * gy * 2 + cellz2) += offz * ihz;

                if (ix < 0 || ix > gx - 1 || iy < 0 || iy > gy - 1) return;

                int cellx1, cellx2, celly1, celly2;
                cellx1 = ix;
                cellx2 = ix + 1;
                celly1 = iy;
                celly2 = iy + 1;

                double aoffx, aoffy;
                aoffx = 1. - offx;
                aoffy = 1. - offy;

                access((cellx1 * gy + celly1) * 2) += w0 * aoffx * aoffy;
                access((cellx1 * gx + celly2) * 2) += w0 * aoffx * offy;
                access((cellx2 * gy + celly1) * 2) += w0 * offx * aoffy;
                access((cellx2 * gx + celly2) * 2) += w0 * offx * offy;
            }
        }
    };

    // use scatter view
    struct sv_zyx_rho_reducer_non_periodic {
        ConstParticles p;
        ConstParticleMasks masks;

        scatter_t scatter;

        int gx, gy, gz; // original grid size
        int dx, dy, dz; // dimensions of the grid
        double ihx, ihy, ihz;
        double lx, ly, lz;
        double w0;

        sv_zyx_rho_reducer_non_periodic(ConstParticles const& p,
                                        ConstParticleMasks const& masks,
                                        scatter_t const& scatter,
                                        std::array<int, 3> const& g,
                                        std::array<int, 3> const& d,
                                        std::array<double, 3> const& h,
                                        std::array<double, 3> const& l,
                                        double w0)
            : p(p)
            , masks(masks)
            , scatter(scatter)
            , gx(g[0])
            , gy(g[1])
            , gz(g[2])
            , dx(d[0])
            , dy(d[1])
            , dz(d[2])
            , ihx(1.0 / h[0])
            , ihy(1.0 / h[1])
            , ihz(1.0 / h[2])
            , lx(l[0])
            , ly(l[1])
            , lz(l[2])
            , w0(w0)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            if (masks(i)) {
                auto access = scatter.access();

                int ix, iy, iz;
                double ox, oy, oz;

                get_leftmost_indices_offset(p(i, 0), lx, ihx, ix, ox);
                get_leftmost_indices_offset(p(i, 2), ly, ihy, iy, oy);
                get_leftmost_indices_offset(p(i, 4), lz, ihz, iz, oz);

                double aox = 1.0 - ox;
                double aoy = 1.0 - oy;
                double aoz = 1.0 - oz;

                int base = iz * dx * dy;

                if (ingrid(ix, iy, iz, gx, gy, gz))
                    access(base + iy * dx + ix) += w0 * aox * aoy * aoz;

                if (ingrid(ix + 1, iy, iz, gx, gy, gz))
                    access(base + iy * dx + ix + 1) += w0 * ox * aoy * aoz;

                if (ingrid(ix, iy + 1, iz, gx, gy, gz))
                    access(base + (iy + 1) * dx + ix) += w0 * aox * oy * aoz;

                if (ingrid(ix + 1, iy + 1, iz, gx, gy, gz))
                    access(base + (iy + 1) * dx + ix + 1) += w0 * ox * oy * aoz;

                base = (iz + 1) * dx * dy;

                if (ingrid(ix, iy, iz + 1, gx, gy, gz))
                    access(base + iy * dx + ix) += w0 * aox * aoy * oz;

                if (ingrid(ix + 1, iy, iz + 1, gx, gy, gz))
                    access(base + iy * dx + ix + 1) += w0 * ox * aoy * oz;

                if (ingrid(ix, iy + 1, iz + 1, gx, gy, gz))
                    access(base + (iy + 1) * dx + ix) += w0 * aox * oy * oz;

                if (ingrid(ix + 1, iy + 1, iz + 1, gx, gy, gz))
                    access(base + (iy + 1) * dx + ix + 1) += w0 * ox * oy * oz;
            }
        }
    };

    // use scatter view
    struct sv_xyz_rho_reducer_non_periodic {
        ConstParticles p;
        ConstParticleMasks masks;

        scatter_t scatter;

        int gx, gy, gz; // original grid size
        int dx, dy, dz; // dimensions of the grid
        double ihx, ihy, ihz;
        double lx, ly, lz;
        double w0;

        sv_xyz_rho_reducer_non_periodic(ConstParticles const& p,
                                        ConstParticleMasks const& masks,
                                        scatter_t const& scatter,
                                        std::array<int, 3> const& g,
                                        std::array<int, 3> const& d,
                                        std::array<double, 3> const& h,
                                        std::array<double, 3> const& l,
                                        double w0)
            : p(p)
            , masks(masks)
            , scatter(scatter)
            , gx(g[0])
            , gy(g[1])
            , gz(g[2])
            , dx(d[0])
            , dy(d[1])
            , dz(d[2])
            , ihx(1.0 / h[0])
            , ihy(1.0 / h[1])
            , ihz(1.0 / h[2])
            , lx(l[0])
            , ly(l[1])
            , lz(l[2])
            , w0(w0)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            if (masks(i)) {
                auto access = scatter.access();

                int ix, iy, iz;
                double ox, oy, oz;

                get_leftmost_indices_offset(p(i, 0), lx, ihx, ix, ox);
                get_leftmost_indices_offset(p(i, 2), ly, ihy, iy, oy);
                get_leftmost_indices_offset(p(i, 4), lz, ihz, iz, oz);

                double aox = 1.0 - ox;
                double aoy = 1.0 - oy;
                double aoz = 1.0 - oz;

                int base = ix * dy * dz;

                if (ingrid(ix, iy, iz, gx, gy, gz))
                    access(base + iy * dz + iz) += w0 * aox * aoy * aoz;

                if (ingrid(ix, iy, iz + 1, gx, gy, gz))
                    access(base + iy * dz + iz + 1) += w0 * aox * aoy * oz;

                if (ingrid(ix, iy + 1, iz, gx, gy, gz))
                    access(base + (iy + 1) * dz + iz) += w0 * aox * oy * aoz;

                if (ingrid(ix, iy + 1, iz + 1, gx, gy, gz))
                    access(base + (iy + 1) * dz + iz + 1) += w0 * aox * oy * oz;

                base = (ix + 1) * dy * dz;

                if (ingrid(ix + 1, iy, iz, gx, gy, gz))
                    access(base + iy * dz + iz) += w0 * ox * aoy * aoz;

                if (ingrid(ix + 1, iy, iz + 1, gx, gy, gz))
                    access(base + iy * dz + iz + 1) += w0 * ox * aoy * oz;

                if (ingrid(ix + 1, iy + 1, iz, gx, gy, gz))
                    access(base + (iy + 1) * dz + iz) += w0 * ox * oy * aoz;

                if (ingrid(ix + 1, iy + 1, iz + 1, gx, gy, gz))
                    access(base + (iy + 1) * dz + iz + 1) += w0 * ox * oy * oz;
            }
        }
    };

}

karray1d_dev
deposit_charge_rectangular_2d_kokkos(Rectangular_grid_domain& domain,
                                     karray2d_dev& particle_bin,
                                     Bunch const& bunch)
{
    using deposit_impl::rho_mover;
    using deposit_impl::rho_reducer;

    auto g = domain.get_grid_shape();
    auto h = domain.get_cell_size();
    auto l = domain.get_left();

    auto parts = bunch.get_local_particles();
    int nparts = bunch.get_local_num();

    double weight0 = (bunch.get_real_num() / bunch.get_total_num()) *
                     bunch.get_particle_charge() * pconstants::e /
                     (h[0] * h[1]); // * h[2]);

    karray1d_dev rho_dbl("rho", g[0] * g[1] + g[2]);
    rho_reducer rr(parts, particle_bin, g, h, l, weight0);
    Kokkos::parallel_reduce(nparts, rr, rho_dbl);

    karray1d_dev rho_cplx("rho_cplx", g[0] * g[1] * 2 + g[2]);
    rho_mover rm(rho_dbl, rho_cplx, g);
    Kokkos::parallel_for(g[0] * g[1] + g[2], rm);
    Kokkos::fence();

    return rho_cplx;
}

karray1d_dev
deposit_charge_rectangular_2d_kokkos_atomic(Rectangular_grid_domain& domain,
                                            karray2d_dev& particle_bin,
                                            Bunch const& bunch)
{
    using deposit_impl::atomic_rho_reducer;

    auto g = domain.get_grid_shape();
    auto h = domain.get_cell_size();
    auto l = domain.get_left();

    auto parts = bunch.get_local_particles();
    int nparts = bunch.get_local_num();

    double weight0 = (bunch.get_real_num() / bunch.get_total_num()) *
                     bunch.get_particle_charge() * pconstants::e /
                     (h[0] * h[1]); // * h[2]);

    // double[x][y][2] + double[z]
    karray1d_dev rho_dev("rho", g[0] * g[1] * 2 + g[2]);
    karray1d_atomic_dev rho_atomic = rho_dev;

    atomic_rho_reducer rr(parts, rho_atomic, particle_bin, g, h, l, weight0);
    Kokkos::parallel_for(nparts, rr);
    Kokkos::fence();

    return rho_dev;
}

void
deposit_charge_rectangular_2d_kokkos_scatter_view(
    karray1d_dev& rho_dev,
    Rectangular_grid_domain& domain,
    karray2d_dev& particle_bin,
    Bunch const& bunch)
{
    using deposit_impl::rho_zeroer;
    using deposit_impl::sv_rho_reducer;

    auto g = domain.get_grid_shape();
    auto h = domain.get_cell_size();
    auto l = domain.get_left();

    auto parts = bunch.get_local_particles();
    auto masks = bunch.get_local_particle_masks();
    int nparts = bunch.size();

    double weight0 = (bunch.get_real_num() / bunch.get_total_num()) *
                     bunch.get_particle_charge() * pconstants::e /
                     (h[0] * h[1]); // * h[2]);

    // double[x][y][2] + double[z]
    if (rho_dev.extent(0) != g[0] * g[1] * 2 + g[2])
        throw std::runtime_error("wrong size for rho in deposit charge");

    // zero first
    rho_zeroer rz{rho_dev};
    Kokkos::parallel_for(rho_dev.extent(0), rz);
    Kokkos::fence();

    // deposit
    Kokkos::Experimental::ScatterView<double*, Kokkos::LayoutLeft> scatter(
        rho_dev);

    sv_rho_reducer rr(parts, masks, scatter, particle_bin, g, h, l, weight0);
    Kokkos::parallel_for(nparts, rr);
    Kokkos::Experimental::contribute(rho_dev, scatter);

    Kokkos::fence();
}

void
deposit_charge_rectangular_3d_kokkos_scatter_view(
    karray1d_dev& rho_dev,
    Rectangular_grid_domain& domain,
    std::array<int, 3> const& dims,
    Bunch const& bunch)
{
    using namespace deposit_impl;

    auto g = domain.get_grid_shape();
    auto h = domain.get_cell_size();
    auto l = domain.get_left();

    // g[0] (nx) needs to be padded to (2*(g[0]/2+1))
    // int padded_gx = Distributed_fft3d::get_padded_shape_real(g[0]);

    auto parts = bunch.get_local_particles();
    auto masks = bunch.get_local_particle_masks();
    int nparts = bunch.size();

    double weight0 = (bunch.get_real_num() / bunch.get_total_num()) *
                     bunch.get_particle_charge() * pconstants::e /
                     (h[0] * h[1] * h[2]);

    if (rho_dev.extent(0) < g[0] * g[1] * g[2])
        throw std::runtime_error("insufficient size for rho in deposit charge");

    // zero first
    rho_zeroer rz{rho_dev};
    Kokkos::parallel_for(rho_dev.extent(0), rz);
    Kokkos::fence();

    // deposit
    scatter_t scatter(rho_dev);
    sv_zyx_rho_reducer_non_periodic rr(
        parts, masks, scatter, g, dims, h, l, weight0);

    Kokkos::parallel_for(nparts, rr);
    Kokkos::Experimental::contribute(rho_dev, scatter);

    Kokkos::fence();
}

void
deposit_charge_rectangular_3d_kokkos_scatter_view_xyz(
    karray1d_dev& rho_dev,
    Rectangular_grid_domain& domain,
    std::array<int, 3> const& dims,
    Bunch const& bunch)
{
    using namespace deposit_impl;

    auto g = domain.get_grid_shape();
    auto h = domain.get_cell_size();
    auto l = domain.get_left();

    // g[0] (nx) needs to be padded to (2*(g[0]/2+1))
    // int padded_gx = Distributed_fft3d::get_padded_shape_real(g[0]);

    auto parts = bunch.get_local_particles();
    auto masks = bunch.get_local_particle_masks();
    int nparts = bunch.size();

    double weight0 = (bunch.get_real_num() / bunch.get_total_num()) *
                     bunch.get_particle_charge() * pconstants::e /
                     (h[0] * h[1] * h[2]);

    if (rho_dev.extent(0) < g[0] * g[1] * g[2])
        throw std::runtime_error("insufficient size for rho in deposit charge");

    // zero first
    rho_zeroer rz{rho_dev};
    Kokkos::parallel_for(rho_dev.extent(0), rz);
    Kokkos::fence();

    // deposit
    scatter_t scatter(rho_dev);
    sv_xyz_rho_reducer_non_periodic rr(
        parts, masks, scatter, g, dims, h, l, weight0);

    Kokkos::parallel_for(nparts, rr);
    Kokkos::Experimental::contribute(rho_dev, scatter);

    Kokkos::fence();
}

#ifdef SYNERGIA_ENABLE_OPENMP
void
deposit_charge_rectangular_2d_omp_reduce(karray1d_dev& rho_dev,
                                         Rectangular_grid_domain& domain,
                                         karray2d_dev& bin,
                                         Bunch const& bunch)
{
    using namespace deposit_impl;

    auto g = domain.get_grid_shape();
    auto h = domain.get_cell_size();
    auto l = domain.get_left();

    auto parts = bunch.get_local_particles();
    auto masks = bunch.get_local_particle_masks();
    int npart = bunch.size();

    double w0 = (bunch.get_real_num() / bunch.get_total_num()) *
                bunch.get_particle_charge() * pconstants::e / (h[0] * h[1]);

    if (rho_dev.extent(0) != g[0] * g[1] * 2 + g[2])
        throw std::runtime_error("insufficient size for rho in deposit charge");

    // zero first
    rho_zeroer rz{rho_dev};
    Kokkos::parallel_for(rho_dev.extent(0), rz);
    Kokkos::fence();

    int gx = g[0];
    int gy = g[1];
    int gz = g[2];

    int nc = gx * gy + gz; // num of cells

    static int nt = 0;
    static int ncc = 0;
    static double* rl = 0;

    if (nt == 0) {
#pragma omp parallel
        {
            nt = omp_get_num_threads();
        }

        ncc = nc;

        rl = new double[nt * nc]; // nt copies of +1 cells
    }

    if (nc != ncc) {
        // bunch geometry has been changed
        delete[] rl;

        ncc = nc;

        rl = new double[nt * nc]; // nt copies of +1 cells
    }

    double lx = l[0];
    double ly = l[1];
    double lz = l[2];

    double ihx = 1.0 / h[0];
    double ihy = 1.0 / h[1];
    double ihz = 1.0 / h[2];

#pragma omp parallel shared(npart,                                             \
                            parts,                                             \
                            masks,                                             \
                            bin,                                               \
                            lx,                                                \
                            ly,                                                \
                            lz,                                                \
                            ihx,                                               \
                            ihy,                                               \
                            ihz,                                               \
                            w0,                                                \
                            gx,                                                \
                            gy,                                                \
                            gz,                                                \
                            rl,                                                \
                            nc,                                                \
                            rho_dev)
    {
        int nt = omp_get_num_threads();
        int it = omp_get_thread_num();

        int np = npart;
        int le = npart / nt;
        int ps = it * le;
        int pe = (it == nt - 1) ? np : (it + 1) * le;

        // zero the worksheet
        std::memset(rl + it * nc, 0, sizeof(double) * nc);

        double ox, oy, oz, w;
        int ix, iy, iz;

        for (int n = ps; n < pe; ++n) {
            if (!masks(n)) continue;

            get_leftmost_indices_offset(parts(n, 0), lx, ihx, ix, ox);
            get_leftmost_indices_offset(parts(n, 2), ly, ihy, iy, oy);
            get_leftmost_indices_offset(parts(n, 4), lz, ihz, iz, oz);

            bin(n, 0) = ix;
            bin(n, 1) = ox;
            bin(n, 2) = iy;
            bin(n, 3) = oy;
            bin(n, 4) = iz;
            bin(n, 5) = oz;

            int base = it * nc;

            int cellz1 = iz;
            int cellz2 = cellz1 + 1;

            if (cellz1 >= 0 && cellz1 < gz)
                rl[base + gx * gy + cellz1] += (1.0 - oz) * ihz;

            if (cellz2 >= 0 && cellz2 < gz)
                rl[base + gx * gy + cellz2] += oz * ihz;

            if (ix < 0 || ix > gx - 1 || iy < 0 || iy > gy - 1) continue;

            int cellx1 = ix;
            int cellx2 = ix + 1;
            int celly1 = iy;
            int celly2 = iy + 1;

            double aox, aoy;
            aox = 1. - ox;
            aoy = 1. - oy;

            rl[base + cellx1 * gy + celly1] += w0 * aox * aoy;
            rl[base + cellx1 * gy + celly2] += w0 * aox * oy;
            rl[base + cellx2 * gy + celly1] += w0 * ox * aoy;
            rl[base + cellx2 * gy + celly2] += w0 * ox * oy;
        }

#pragma omp barrier

        // reduction
        le = gy / nt;
        ps = it * le;
        pe = (it == nt - 1) ? gy : (ps + le);

        for (int y = ps; y < pe; ++y) {
            for (int x = 0; x < gx; ++x) {
                w = 0.0;

                for (int n = 0; n < nt; ++n)
                    w += rl[n * nc + (x * gy + y)];

                rho_dev((x * gy + y) * 2) = w;
            }
        }

#pragma omp barrier

        le = gz / nt;
        ps = it * le;
        pe = (it == nt - 1) ? gz : (ps + le);

        for (int z = ps; z < pe; ++z) {
            w = 0.0;

            for (int n = 0; n < nt; ++n)
                w += rl[n * nc + gx * gy + z];

            rho_dev(gx * gy * 2 + z) = w;
        }

#pragma omp barrier

    } //  end of #pragma parallel
}

void
deposit_charge_rectangular_3d_omp_reduce(karray1d_dev& rho_dev,
                                         Rectangular_grid_domain& domain,
                                         std::array<int, 3> const& dims,
                                         Bunch const& bunch)
{
    using namespace deposit_impl;

    auto g = domain.get_grid_shape();
    auto h = domain.get_cell_size();
    auto l = domain.get_left();

    auto parts = bunch.get_local_particles();
    auto masks = bunch.get_local_particle_masks();
    int npart = bunch.size();

    double w0 = (bunch.get_real_num() / bunch.get_total_num()) *
                bunch.get_particle_charge() * pconstants::e /
                (h[0] * h[1] * h[2]);

    if (rho_dev.extent(0) < g[0] * g[1] * g[2])
        throw std::runtime_error("insufficient size for rho in deposit charge");

    // zero first
    rho_zeroer rz{rho_dev};
    Kokkos::parallel_for(rho_dev.extent(0), rz);
    Kokkos::fence();

    int gx = g[0];
    int gy = g[1];
    int gz = g[2];

    int dx = dims[0];
    int dy = dims[1];
    int dz = dims[2];

    int nc = gx * gy * gz; // num of cells

    static int nt = 0;
    static int ncc = 0;
    static double* rl = 0;

    if (nt == 0) {
#pragma omp parallel
        {
            nt = omp_get_num_threads();
        }

        ncc = nc;

        rl = new double[nt * nc]; // nt copies of +1 cells
    }

    if (nc != ncc) {
        // bunch geometry has been changed
        delete[] rl;

        ncc = nc;

        rl = new double[nt * nc]; // nt copies of +1 cells
    }

    double lx = l[0];
    double ly = l[1];
    double lz = l[2];

    double ihx = 1.0 / h[0];
    double ihy = 1.0 / h[1];
    double ihz = 1.0 / h[2];

#pragma omp parallel shared(npart,                                             \
                            parts,                                             \
                            masks,                                             \
                            lx,                                                \
                            ly,                                                \
                            lz,                                                \
                            ihx,                                               \
                            ihy,                                               \
                            ihz,                                               \
                            w0,                                                \
                            gx,                                                \
                            gy,                                                \
                            gz,                                                \
                            dx,                                                \
                            dy,                                                \
                            dz,                                                \
                            rl,                                                \
                            nc,                                                \
                            rho_dev)
    {
        int nt = omp_get_num_threads();
        int it = omp_get_thread_num();

        int np = npart;
        int le = npart / nt;
        int ps = it * le;
        int pe = (it == nt - 1) ? np : (it + 1) * le;

        // zero the worksheet
        std::memset(rl + it * nc, 0, sizeof(double) * nc);

        double ox, oy, oz, w;
        int ix, iy, iz;

        for (int n = ps; n < pe; ++n) {
            if (!masks(n)) continue;

            get_leftmost_indices_offset(parts(n, 0), lx, ihx, ix, ox);
            get_leftmost_indices_offset(parts(n, 2), ly, ihy, iy, oy);
            get_leftmost_indices_offset(parts(n, 4), lz, ihz, iz, oz);

            int base = it * nc + iz * gx * gy;

            if (ingrid(ix, iy, iz, gx, gy, gz))
                // rl[(base + (iy  )*gx + ix  )*nt+it] +=
                // w0*(1-ox)*(1-oy)*(1-oz);
                rl[(base + (iy)*gx + ix)] +=
                    w0 * (1 - ox) * (1 - oy) * (1 - oz);

            if (ingrid(ix + 1, iy, iz, gx, gy, gz))
                // rl[(base + (iy  )*gx + ix+1)*nt+it] += w0*(
                // ox)*(1-oy)*(1-oz);
                rl[(base + (iy)*gx + ix + 1)] +=
                    w0 * (ox) * (1 - oy) * (1 - oz);

            if (ingrid(ix, iy + 1, iz, gx, gy, gz))
                // rl[(base + (iy+1)*gx + ix  )*nt+it] += w0*(1-ox)*(
                // oy)*(1-oz);
                rl[(base + (iy + 1) * gx + ix)] +=
                    w0 * (1 - ox) * (oy) * (1 - oz);

            if (ingrid(ix + 1, iy + 1, iz, gx, gy, gz))
                // rl[(base + (iy+1)*gx + ix+1)*nt+it] += w0*(  ox)*(
                // oy)*(1-oz);
                rl[(base + (iy + 1) * gx + ix + 1)] +=
                    w0 * (ox) * (oy) * (1 - oz);

            base = it * nc + (iz + 1) * gx * gy;

            if (ingrid(ix, iy, iz + 1, gx, gy, gz))
                // rl[(base + (iy  )*gx + ix  )*nt+it] += w0*(1-ox)*(1-oy)*(oz);
                rl[(base + (iy)*gx + ix)] += w0 * (1 - ox) * (1 - oy) * (oz);

            if (ingrid(ix + 1, iy, iz + 1, gx, gy, gz))
                // rl[(base + (iy  )*gx + ix+1)*nt+it] += w0*(  ox)*(1-oy)*(oz);
                rl[(base + (iy)*gx + ix + 1)] += w0 * (ox) * (1 - oy) * (oz);

            if (ingrid(ix, iy + 1, iz + 1, gx, gy, gz))
                // rl[(base + (iy+1)*gx + ix  )*nt+it] += w0*(1-ox)*(  oy)*(oz);
                rl[(base + (iy + 1) * gx + ix)] += w0 * (1 - ox) * (oy) * (oz);

            if (ingrid(ix + 1, iy + 1, iz + 1, gx, gy, gz))
                // rl[(base + (iy+1)*gx + ix+1)*nt+it] += w0*(  ox)*(  oy)*(oz);
                rl[(base + (iy + 1) * gx + ix + 1)] += w0 * (ox) * (oy) * (oz);
        }

#pragma omp barrier

        // reduction
        le = gz / nt;
        ps = it * le;
        pe = (it == nt - 1) ? gz : (ps + le);

        for (int z = ps; z < pe; ++z) {
            for (int y = 0; y < gy; ++y) {
                for (int x = 0; x < gx; ++x) {
                    w = 0.0;
                    for (int n = 0; n < nt; ++n)
                        // w += rl[(z*gx*gy + y*gx + x)*nt+n];
                        w += rl[n * nc + (z * gx * gy + y * gx + x)];

                    rho_dev(z * dx * dy + y * dx + x) = w;
                }
            }
        }

#pragma omp barrier

    } //  end of #pragma parallel
}

void
deposit_charge_rectangular_3d_omp_reduce_xyz(karray1d_dev& rho_dev,
                                             Rectangular_grid_domain& domain,
                                             std::array<int, 3> const& dims,
                                             Bunch const& bunch)
{
    using namespace deposit_impl;

    auto g = domain.get_grid_shape();
    auto h = domain.get_cell_size();
    auto l = domain.get_left();

    auto parts = bunch.get_local_particles();
    auto masks = bunch.get_local_particle_masks();
    int npart = bunch.size();

    double w0 = (bunch.get_real_num() / bunch.get_total_num()) *
                bunch.get_particle_charge() * pconstants::e /
                (h[0] * h[1] * h[2]);

    if (rho_dev.extent(0) < g[0] * g[1] * g[2])
        throw std::runtime_error("insufficient size for rho in deposit charge");

    // zero first
    rho_zeroer rz{rho_dev};
    Kokkos::parallel_for(rho_dev.extent(0), rz);
    Kokkos::fence();

    int gx = g[0];
    int gy = g[1];
    int gz = g[2];

    int dx = dims[0];
    int dy = dims[1];
    int dz = dims[2];

    int nc = gx * gy * gz; // num of cells

    static int nt = 0;
    static int ncc = 0;
    static double* rl = 0;

    if (nt == 0) {
#pragma omp parallel
        {
            nt = omp_get_num_threads();
        }

        ncc = nc;

        rl = new double[nt * nc]; // nt copies of +1 cells
    }

    if (nc != ncc) {
        // bunch geometry has been changed
        delete[] rl;

        ncc = nc;

        rl = new double[nt * nc]; // nt copies of +1 cells
    }

    double lx = l[0];
    double ly = l[1];
    double lz = l[2];

    double ihx = 1.0 / h[0];
    double ihy = 1.0 / h[1];
    double ihz = 1.0 / h[2];

#pragma omp parallel shared(npart,                                             \
                            parts,                                             \
                            masks,                                             \
                            lx,                                                \
                            ly,                                                \
                            lz,                                                \
                            ihx,                                               \
                            ihy,                                               \
                            ihz,                                               \
                            w0,                                                \
                            gx,                                                \
                            gy,                                                \
                            gz,                                                \
                            dx,                                                \
                            dy,                                                \
                            dz,                                                \
                            rl,                                                \
                            nc,                                                \
                            rho_dev)
    {
        int nt = omp_get_num_threads();
        int it = omp_get_thread_num();

        int np = npart;
        int le = npart / nt;
        int ps = it * le;
        int pe = (it == nt - 1) ? np : (it + 1) * le;

        // zero the worksheet
        std::memset(rl + it * nc, 0, sizeof(double) * nc);

        double ox, oy, oz, w;
        int ix, iy, iz;

        for (int n = ps; n < pe; ++n) {
            if (!masks(n)) continue;

            get_leftmost_indices_offset(parts(n, 0), lx, ihx, ix, ox);
            get_leftmost_indices_offset(parts(n, 2), ly, ihy, iy, oy);
            get_leftmost_indices_offset(parts(n, 4), lz, ihz, iz, oz);

            int base = it * nc + ix * gy * gz;

            if (ingrid(ix, iy, iz, gx, gy, gz))
                rl[(base + (iy)*gz + iz)] +=
                    w0 * (1 - ox) * (1 - oy) * (1 - oz);

            if (ingrid(ix, iy, iz + 1, gx, gy, gz))
                rl[(base + (iy)*gz + iz + 1)] +=
                    w0 * (1 - ox) * (1 - oy) * (oz);

            if (ingrid(ix, iy + 1, iz, gx, gy, gz))
                rl[(base + (iy + 1) * gz + iz)] +=
                    w0 * (1 - ox) * (oy) * (1 - oz);

            if (ingrid(ix, iy + 1, iz + 1, gx, gy, gz))
                rl[(base + (iy + 1) * gz + iz + 1)] +=
                    w0 * (1 - ox) * (oy) * (oz);

            base = it * nc + (ix + 1) * gy * gz;

            if (ingrid(ix + 1, iy, iz, gx, gy, gz))
                rl[(base + (iy)*gz + iz)] += w0 * (ox) * (1 - oy) * (1 - oz);

            if (ingrid(ix + 1, iy, iz + 1, gx, gy, gz))
                rl[(base + (iy)*gz + iz + 1)] += w0 * (ox) * (1 - oy) * (oz);

            if (ingrid(ix + 1, iy + 1, iz, gx, gy, gz))
                rl[(base + (iy + 1) * gz + iz)] += w0 * (ox) * (oy) * (1 - oz);

            if (ingrid(ix + 1, iy + 1, iz + 1, gx, gy, gz))
                rl[(base + (iy + 1) * gz + iz + 1)] += w0 * (ox) * (oy) * (oz);
        }

#pragma omp barrier

        // reduction
        le = gz / nt;
        ps = it * le;
        pe = (it == nt - 1) ? gz : (ps + le);

        for (int z = ps; z < pe; ++z) {
            for (int y = 0; y < gy; ++y) {
                for (int x = 0; x < gx; ++x) {
                    w = 0.0;
                    for (int n = 0; n < nt; ++n)
                        // w += rl[(z*gx*gy + y*gx + x)*nt+n];
                        w += rl[n * nc + (x * gy * gz + y * gz + z)];

                    rho_dev(x * dy * dz + y * dz + z) = w;
                }
            }
        }

#pragma omp barrier

    } //  end of #pragma parallel
}
#endif // SYNERGIA_ENABLE_OPENMP
