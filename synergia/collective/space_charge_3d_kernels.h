#ifndef SPACE_CHARGE_3D_KERNELS_H_
#define SPACE_CHARGE_3D_KERNELS_H_

#include "utils.h"

#include "synergia/bunch/bunch.h"
#include "synergia/utils/distributed_fft3d.h"
#include "synergia/utils/kokkos_views.h"

/*! common kernels for 3D space charge solvers */
namespace sc3d_kernels {

    /*! Kernels where the data has zyx ordering */
    namespace zyx {

        struct alg_force_extractor {
            karray1d_dev phi2;
            karray1d_dev enx;
            karray1d_dev eny;
            karray1d_dev enz;

            int gx, gy, gz;
            double ihx, ihy, ihz;
            double igxgy;
            double igx;

            alg_force_extractor(karray1d_dev const& phi2,
                                karray1d_dev const& enx,
                                karray1d_dev const& eny,
                                karray1d_dev const& enz,
                                std::array<int, 3> const& g,
                                std::array<double, 3> const& h)
                : phi2(phi2)
                , enx(enx)
                , eny(eny)
                , enz(enz)
                , gx(g[0])
                , gy(g[1])
                , gz(g[2])
                , ihx(0.5 / h[0])
                , ihy(0.5 / h[1])
                , ihz(0.5 / h[2])
                , igxgy(1.0 / (gx * gy))
                , igx(1.0 / gx)
            {}

            KOKKOS_INLINE_FUNCTION
            void
            operator()(const int i) const
            {
                int iz = i * igxgy;
                int iy = (i - iz * gx * gy) * igx;
                int ix = i - iz * gx * gy - iy * gx;

                int ixl, ixr, iyl, iyr, izl, izr;

                double idx, idy, idz;

                // all boundaries will be skipped (set to 0)
                if (ix == 0 || ix == gx - 1 || iy == 0 || iy == gy - 1 ||
                    iz == 0 || iz == gz - 1)
                    return;

                ixl = ix - 1;
                ixr = ix + 1;
                iyl = iy - 1;
                iyr = iy + 1;
                izl = iz - 1;
                izr = iz + 1;

                idx = ihx;
                idy = ihy;
                idz = ihz;

                int idx_r = iz * gx * gy + iy * gx + ixr;
                int idx_l = iz * gx * gy + iy * gx + ixl;
                enx(i) = -(phi2(idx_r) - phi2(idx_l)) * idx;

                idx_r = iz * gx * gy + iyr * gx + ix;
                idx_l = iz * gx * gy + iyl * gx + ix;
                eny(i) = -(phi2(idx_r) - phi2(idx_l)) * idy;

                idx_r = izr * gx * gy + iy * gx + ix;
                idx_l = izl * gx * gy + iy * gx + ix;
                auto val_right = phi2(idx_r) * 1.0;
                auto val_left = phi2(idx_l) * 1.0;
                enz(i) = -(phi2(idx_r) - phi2(idx_l)) * idz;
            }
        };

        struct alg_kicker {
            Particles parts;
            ConstParticleMasks masks;

            karray1d_dev enx;
            karray1d_dev eny;
            karray1d_dev enz;

            int gx, gy, gz;
            double ihx, ihy, ihz;
            double lx, ly, lz;
            double factor, pref, m;

            alg_kicker(Particles parts,
                       ConstParticleMasks masks,
                       karray1d_dev const& enx,
                       karray1d_dev const& eny,
                       karray1d_dev const& enz,
                       std::array<int, 3> const& g,
                       std::array<double, 3> const& h,
                       std::array<double, 3> const& l,
                       double factor,
                       double pref,
                       double m)
                : parts(parts)
                , masks(masks)
                , enx(enx)
                , eny(eny)
                , enz(enz)
                , gx(g[0])
                , gy(g[1])
                , gz(g[2])
                , ihx(1.0 / h[0])
                , ihy(1.0 / h[1])
                , ihz(1.0 / h[2])
                , lx(l[0])
                , ly(l[1])
                , lz(l[2])
                , factor(factor)
                , pref(pref)
                , m(m)
            {}

            KOKKOS_INLINE_FUNCTION
            void
            operator()(const int i) const
            {
                if (masks(i)) {
                    int ix, iy, iz;
                    double ox, oy, oz;

                    get_leftmost_indices_offset(parts(i, 0), lx, ihx, ix, ox);
                    get_leftmost_indices_offset(parts(i, 2), ly, ihy, iy, oy);
                    get_leftmost_indices_offset(parts(i, 4), lz, ihz, iz, oz);

                    double aox = 1.0 - ox;
                    double aoy = 1.0 - oy;
                    double aoz = 1.0 - oz;

                    if ((ix >= 0 && ix < gx - 1) && (iy >= 0 && iy < gy - 1) &&
                        ((iz >= 0 && iz < gz - 1) /* || periodic_z */)) {
                        double val = 0;
                        int base = 0;

                        // enz
                        base = iz * gx * gy + iy * gx + ix;
                        val = aoz * aoy * aox * enz(base);      // z, y, x
                        val += aoz * aoy * ox * enz(base + 1);  // z, y, x+1
                        val += aoz * oy * aox * enz(base + gx); // z, y+1, x
                        val +=
                            aoz * oy * ox * enz(base + gx + 1); // z, y+1, x+1

                        base = (iz + 1) * gx * gy + iy * gx + ix;
                        val += oz * aoy * aox * enz(base);     // z+1, y, x
                        val += oz * aoy * ox * enz(base + 1);  // z+1, y, x+1
                        val += oz * oy * aox * enz(base + gx); // z+1, y+1, x
                        val +=
                            oz * oy * ox * enz(base + gx + 1); // z+1, y+1, x+1

                        double p = pref + parts(i, 5) * pref;
                        double Eoc_i = std::sqrt(p * p + m * m);
                        double Eoc_f = Eoc_i + factor * (-pref) * val;
                        double dpop = (std::sqrt(Eoc_f * Eoc_f - m * m) -
                                       std::sqrt(Eoc_i * Eoc_i - m * m)) /
                                      pref;

                        parts(i, 5) += dpop;

                        // eny
                        base = iz * gx * gy + iy * gx + ix;
                        val = aoz * aoy * aox * eny(base);      // z, y, x
                        val += aoz * aoy * ox * eny(base + 1);  // z, y, x+1
                        val += aoz * oy * aox * eny(base + gx); // z, y+1, x
                        val +=
                            aoz * oy * ox * eny(base + gx + 1); // z, y+1, x+1

                        base = (iz + 1) * gx * gy + iy * gx + ix;
                        val += oz * aoy * aox * eny(base);     // z+1, y, x
                        val += oz * aoy * ox * eny(base + 1);  // z+1, y, x+1
                        val += oz * oy * aox * eny(base + gx); // z+1, y+1, x
                        val +=
                            oz * oy * ox * eny(base + gx + 1); // z+1, y+1, x+1

                        parts(i, 3) += factor * val;

                        // enx
                        base = iz * gx * gy + iy * gx + ix;
                        val = aoz * aoy * aox * enx(base);      // z, y, x
                        val += aoz * aoy * ox * enx(base + 1);  // z, y, x+1
                        val += aoz * oy * aox * enx(base + gx); // z, y+1, x
                        val +=
                            aoz * oy * ox * enx(base + gx + 1); // z, y+1, x+1

                        base = (iz + 1) * gx * gy + iy * gx + ix;
                        val += oz * aoy * aox * enx(base);     // z+1, y, x
                        val += oz * aoy * ox * enx(base + 1);  // z+1, y, x+1
                        val += oz * oy * aox * enx(base + gx); // z+1, y+1, x
                        val +=
                            oz * oy * ox * enx(base + gx + 1); // z+1, y+1, x+1

                        parts(i, 1) += factor * val;
                    }
                }
            }
        };

        struct alg_force_extractor_doubled_domain {
            karray1d_dev phi2;
            karray1d_dev enx;
            karray1d_dev eny;
            karray1d_dev enz;

            int gx, gy, gz;
            int dgx, dgy;
            double ihx, ihy, ihz;
            double igxgy;
            double igx;

            alg_force_extractor_doubled_domain(karray1d_dev const& phi2,
                                               karray1d_dev const& enx,
                                               karray1d_dev const& eny,
                                               karray1d_dev const& enz,
                                               std::array<int, 3> const& g,
                                               std::array<int, 3> const& dg,
                                               std::array<double, 3> const& h)
                : phi2(phi2)
                , enx(enx)
                , eny(eny)
                , enz(enz)
                , gx(g[0])
                , gy(g[1])
                , gz(g[2])
                , dgx(Distributed_fft3d::get_padded_shape_real(dg[0]))
                , dgy(dg[1])
                , ihx(0.5 / h[0])
                , ihy(0.5 / h[1])
                , ihz(0.5 / h[2])
                , igxgy(1.0 / (gx * gy))
                , igx(1.0 / gx)
            {}

            KOKKOS_INLINE_FUNCTION
            void
            operator()(const int i) const
            {
                int iz = i * igxgy;
                int iy = (i - iz * gx * gy) * igx;
                int ix = i - iz * gx * gy - iy * gx;

                int ixl, ixr, iyl, iyr, izl, izr;

                double idx, idy, idz;

                // all boundaries will be skipped (set to 0)
                if (ix == 0 || ix == gx - 1 || iy == 0 || iy == gy - 1 ||
                    iz == 0 || iz == gz - 1)
                    return;

                ixl = ix - 1;
                ixr = ix + 1;
                iyl = iy - 1;
                iyr = iy + 1;
                izl = iz - 1;
                izr = iz + 1;

                idx = ihx;
                idy = ihy;
                idz = ihz;

#if 0
            // or, calculate the full domain
            if (ix==0) { ixl = 0; ixr = 1; idx = ihx*2; }
            else if (ix==gx-1) { ixl = gx-2; ixr = gx-1; idx = ihx; }
            else { ixl = ix-1; ixr = ix+1; idx = ihx; }

            if (iy==0) { iyl = 0; iyr = 1; idy = ihy*2; }
            else if (iy==gy-1) { iyl = gy-2; iyr = gy-1; idy = ihy; }
            else { iyl = iy-1; iyr = iy+1; idy = ihy; }

            if (iz==0) { izl = 0; izr = 1; idz = ihz*2; }
            else if (iz==gz-1) { izl = gz-2; izr = gz-1; idz = ihz; }
            else { izl = iz-1; izr = iz+1; idz = ihz; }
#endif

                int idx_r = iz * dgx * dgy + iy * dgx + ixr;
                int idx_l = iz * dgx * dgy + iy * dgx + ixl;
                enx(i) = -(phi2(idx_r) - phi2(idx_l)) * idx;

                idx_r = iz * dgx * dgy + iyr * dgx + ix;
                idx_l = iz * dgx * dgy + iyl * dgx + ix;
                eny(i) = -(phi2(idx_r) - phi2(idx_l)) * idy;

                idx_r = izr * dgx * dgy + iy * dgx + ix;
                idx_l = izl * dgx * dgy + iy * dgx + ix;
                enz(i) = -(phi2(idx_r) - phi2(idx_l)) * idz;
            }
        };

    }

    /*! Kernels where the data has xyz ordering */
    namespace xyz {
        struct alg_force_extractor {
            karray1d_dev phi;
            karray1d_dev enx;
            karray1d_dev eny;
            karray1d_dev enz;

            int gx, gy, gz;
            int dgx, dgy;
            double ihx, ihy, ihz;
            double igygz;
            double igz;

            alg_force_extractor(karray1d_dev const& phi,
                                karray1d_dev const& enx,
                                karray1d_dev const& eny,
                                karray1d_dev const& enz,
                                std::array<int, 3> const& g,
                                std::array<double, 3> const& h)
                : phi(phi)
                , enx(enx)
                , eny(eny)
                , enz(enz)
                , gx(g[0])
                , gy(g[1])
                , gz(g[2])
                , ihx(0.5 / h[0])
                , ihy(0.5 / h[1])
                , ihz(0.5 / h[2])
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

                int ixl, ixr, iyl, iyr, izl, izr;
                double idx, idy, idz;

                // all x-y boundaries will be skipped (set to 0)
                if (ix == 0 || ix == gx - 1 || iy == 0 || iy == gy - 1) return;

                ixl = ix - 1;
                ixr = ix + 1;
                iyl = iy - 1;
                iyr = iy + 1;
                izl = iz - 1;
                izr = iz + 1;

                idx = ihx;
                idy = ihy;
                idz = ihz;

                // periodic z
                if (iz == 0) {
                    izl = gz - 1;
                } else if (iz == gz - 1) {
                    izr = 0;
                }

                int idx_r = ixr * gy * gz + iy * gz + iz;
                int idx_l = ixl * gy * gz + iy * gz + iz;
                enx(i) = -(phi(idx_r) - phi(idx_l)) * idx;

                idx_r = ix * gy * gz + iyr * gz + iz;
                idx_l = ix * gy * gz + iyl * gz + iz;
                eny(i) = -(phi(idx_r) - phi(idx_l)) * idy;

                idx_r = ix * gy * gz + iy * gz + izr;
                idx_l = ix * gy * gz + iy * gz + izl;
                enz(i) = -(phi(idx_r) - phi(idx_l)) * idz;
            }
        };

        struct alg_kicker {
            Particles parts;
            ConstParticleMasks masks;

            karray1d_dev enx;
            karray1d_dev eny;
            karray1d_dev enz;

            int gx, gy, gz;
            double ihx, ihy, ihz;
            double lx, ly, lz;
            double factor, pref, m;

            alg_kicker(Particles parts,
                       ConstParticleMasks masks,
                       karray1d_dev const& enx,
                       karray1d_dev const& eny,
                       karray1d_dev const& enz,
                       std::array<int, 3> const& g,
                       std::array<double, 3> const& h,
                       std::array<double, 3> const& l,
                       double factor,
                       double pref,
                       double m)
                : parts(parts)
                , masks(masks)
                , enx(enx)
                , eny(eny)
                , enz(enz)
                , gx(g[0])
                , gy(g[1])
                , gz(g[2])
                , ihx(1.0 / h[0])
                , ihy(1.0 / h[1])
                , ihz(1.0 / h[2])
                , lx(l[0])
                , ly(l[1])
                , lz(l[2])
                , factor(factor)
                , pref(pref)
                , m(m)
            {}

            KOKKOS_INLINE_FUNCTION
            void
            operator()(const int i) const
            {
                if (masks(i)) {
                    int ix, iy, iz;
                    double ox, oy, oz;

                    get_leftmost_indices_offset(parts(i, 0), lx, ihx, ix, ox);
                    get_leftmost_indices_offset(parts(i, 2), ly, ihy, iy, oy);
                    get_leftmost_indices_offset(parts(i, 4), lz, ihz, iz, oz);

                    double aox = 1.0 - ox;
                    double aoy = 1.0 - oy;
                    double aoz = 1.0 - oz;

                    if (ix >= 0 && ix < gx - 1 && iy >= 0 && iy < gy - 1) {
                        while (iz > gz - 1)
                            iz -= gz;
                        while (iz < 0)
                            iz += gz;

                        int izp1 = (iz == (gz - 1)) ? 0 : iz + 1;

                        double val = 0;
                        int base = 0;

                        // enz
                        base = ix * gy * gz + iy * gz;
                        val = aox * aoy * aoz * enz(base + iz);   // x, y, z
                        val += aox * aoy * oz * enz(base + izp1); // x, y, z+1
                        val +=
                            aox * oy * aoz * enz(base + gz + iz); // x, y+1, z
                        val += aox * oy * oz *
                               enz(base + gz + izp1); // x, y+1, z+1

                        base = (ix + 1) * gy * gz + iy * gz;
                        val += ox * aoy * aoz * enz(base + iz);  // x+1, y, z
                        val += ox * aoy * oz * enz(base + izp1); // x+1, y, z+1
                        val +=
                            ox * oy * aoz * enz(base + gz + iz); // x+1, y+1, z
                        val += ox * oy * oz *
                               enz(base + gz + izp1); // x+1, y+1, z+1

                        double p = pref + parts(i, 5) * pref;
                        double Eoc_i = std::sqrt(p * p + m * m);
                        double Eoc_f = Eoc_i - factor * pref * val;
                        double dpop = (std::sqrt(Eoc_f * Eoc_f - m * m) -
                                       std::sqrt(Eoc_i * Eoc_i - m * m)) /
                                      pref;

                        parts(i, 5) += dpop;

                        // eny
                        base = ix * gy * gz + iy * gz;
                        val = aox * aoy * aoz * eny(base + iz);   // x, y, z
                        val += aox * aoy * oz * eny(base + izp1); // x, y, z+1
                        val +=
                            aox * oy * aoz * eny(base + gz + iz); // x, y+1, z
                        val += aox * oy * oz *
                               eny(base + gz + izp1); // x, y+1, z+1

                        base = (ix + 1) * gy * gz + iy * gz;
                        val += ox * aoy * aoz * eny(base + iz);  // x+1, y, z
                        val += ox * aoy * oz * eny(base + izp1); // x+1, y, z+1
                        val +=
                            ox * oy * aoz * eny(base + gz + iz); // x+1, y+1, z
                        val += ox * oy * oz *
                               eny(base + gz + izp1); // x+1, y+1, z+1

                        parts(i, 3) += factor * val;

                        // enx
                        base = ix * gy * gz + iy * gz;
                        val = aox * aoy * aoz * enx(base + iz);   // x, y, z
                        val += aox * aoy * oz * enx(base + izp1); // x, y, z+1
                        val +=
                            aox * oy * aoz * enx(base + gz + iz); // x, y+1, z
                        val += aox * oy * oz *
                               enx(base + gz + izp1); // x, y+1, z+1

                        base = (ix + 1) * gy * gz + iy * gz;
                        val += ox * aoy * aoz * enx(base + iz);  // x+1, y, z
                        val += ox * aoy * oz * enx(base + izp1); // x+1, y, z+1
                        val +=
                            ox * oy * aoz * enx(base + gz + iz); // x+1, y+1, z
                        val += ox * oy * oz *
                               enx(base + gz + izp1); // x+1, y+1, z+1

                        parts(i, 1) += factor * val;
                    }
                }
            }
        };

    }

}
#endif
