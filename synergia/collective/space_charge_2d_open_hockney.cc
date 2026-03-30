
#include "space_charge_2d_open_hockney.h"
#include "deposit.h"

#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/kokkos_utils.h"
#include "synergia/utils/simple_timer.h"

constexpr auto pi = Kokkos::numbers::pi_v<double>;

namespace {
    void
    print_grid(Logger& logger,
               karray1d_dev const& grid,
               int x0,
               int x1,
               int y0,
               int y1,
               int gx,
               int gy,
               int off = 0)
    {
        karray1d_hst hgrid = Kokkos::create_mirror_view(grid);
        Kokkos::deep_copy(hgrid, grid);

        double sum = 0;

        int dim = grid.extent(0);
        for (int i = 0; i < dim; ++i)
            sum += hgrid(i);

#if 0
        for(int x=0; x<gx; ++x)
            for(int y=0; y<gy; ++y)
                sum += hgrid((x*gy + y)*2 + off);
#endif

        logger
            //<< std::resetiosflags(std::ios::fixed)
            //<< std::setiosflags(std::ios::showpos | std::ios::scientific)
            << std::setprecision(12);

        logger << "      " << grid.label() << " = " << sum << "\n";

        for (int x = x0; x < x1; ++x) {
            logger << "        " << x << " | ";

            for (int y = y0; y < y1; ++y) {
                logger << hgrid((x * gy + y) * 2 + off) << ", ";
            }

            logger << "\n";
        }
    }

    double
    get_smallest_non_tiny(double val, double other1, double other2, double tiny)
    {
        double retval;
        if (val > tiny) {
            retval = val;
        } else {
            if ((other1 > tiny) && (other2 > tiny)) {
                retval = std::min(other1, other2);
            } else {
                retval = std::max(other1, other2);
            }
        }

        return retval;
    }

    KOKKOS_INLINE_FUNCTION
    int
    fast_int_floor_kokkos(const double x)
    {
        int ix = static_cast<int>(x);
        return x > 0.0 ? ix : ((x - ix == 0) ? ix : ix - 1);
    }

#if 0
    inline std::complex<double >
    interpolate_rectangular_2d(double * bin, bool periodic_z,
            MArray2dc_ref const& a, MArray1d_ref const& b)
    {
        // bi-linear interpolation
        int ix, iy, iz;
        double offx, offy, offz;
        ix = fast_int_floor(bin[0]);
        iy = fast_int_floor(bin[2]);
        iz = fast_int_floor(bin[4]);

        offx = bin[1];
        offy = bin[3];
        offz = bin[5];

        std::vector<int > grid_shape(3);
        grid_shape[0] = a.shape()[0];
        grid_shape[1] = a.shape()[1];
        grid_shape[2] = b.shape()[0];

        std::complex<double > val;
        if ( (grid_shape[2] == 1)
                && ( (ix >= 0)
                    && (ix < grid_shape[0] - 1)
                    && (iy >= 0)
                    && (iy < grid_shape[1] - 1) )
                && ( periodic_z
                    || ( (iz >= 0) && (iz <= grid_shape[2]) ) ) )
        {
            double line_density = b[0];
            val = line_density * ((1.0 - offx) * (1.0 - offy) * a[ix][iy]
                    + offx * (1.0 - offy) * a[ix + 1][iy]
                    + (1.0 - offx) * offy * a[ix][iy + 1]
                    + offx * offy * a[ix + 1][iy + 1]);
        }
        else if ( (grid_shape[2] > 1)
                && ( (ix >= 0)
                    && (ix < grid_shape[0] - 1)
                    && (iy >= 0)
                    && (iy < grid_shape[1] - 1) )
                && ( periodic_z
                    || ( (iz >= 0) && (iz < grid_shape[2] - 1) ) ) )
        {
            double line_density = (1.0 - offz) * b[iz] + offz * b[iz + 1];
            val = line_density * ((1.0 - offx) * (1.0 - offy) * a[ix][iy]
                    + offx * (1.0 - offy) * a[ix + 1][iy]
                    + (1.0 - offx) * offy * a[ix][iy + 1]
                    + offx * offy * a[ix + 1][iy + 1]);
        }

        return val;
    }
#endif

    struct alg_g2_pointlike {
        const double epsilon = 0.01;

        karray1d_dev g2;
        int gx, gy;
        int dgx, dgy;
        double hx, hy;

        alg_g2_pointlike(karray1d_dev const& g2,
                         std::array<int, 3> const& g,
                         std::array<int, 3> const& dg,
                         std::array<double, 3> const& h)
            : g2(g2)
            , gx(g[0])
            , gy(g[1])
            , dgx(dg[0])
            , dgy(dg[1])
            , hx(h[0])
            , hy(h[1])
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            int ix = i / dgy;
            int iy = i - ix * dgy;

            double dx = (ix > gx) ? (ix - dgx) * hx : ix * hx;
            double dy = (iy > gy) ? (iy - dgy) * hy : iy * hy;

            double Gx, Gy;

            if (ix == gx || iy == gy) {
                Gx = 0.0;
                Gy = 0.0;
            } else {
                double inv =
                    1.0 / (dx * dx + dy * dy + hx * hy * epsilon * epsilon);
                Gx = dx * inv;
                Gy = dy * inv;
            }

            g2(i * 2) = Gx;
            g2(i * 2 + 1) = Gy;
        }
    };

    struct alg_cplx_multiplier {
        karray1d_dev prod;
        karray1d_dev m1, m2;
        int off;

        alg_cplx_multiplier(karray1d_dev const& prod,
                            karray1d_dev const& m1,
                            karray1d_dev const& m2,
                            int offset)
            : prod(prod), m1(m1), m2(m2), off(offset)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            const int real = (off + i) * 2;
            const int imag = (off + i) * 2 + 1;

            prod[real] = m1[real] * m2[real] - m1[imag] * m2[imag];
            prod[imag] = m1[real] * m2[imag] + m1[imag] * m2[real];
        }
    };

    struct alg_kicker {
        Particles p;
        ConstParticleMasks masks;

        karray1d_dev fn;
        karray1d_dev rho;
        karray2d_dev bin;

        int gx, gy, gz;
        double factor;

        alg_kicker(Particles p,
                   ConstParticleMasks masks,
                   karray1d_dev const& fn,
                   karray1d_dev const& rho,
                   karray2d_dev const& bin,
                   std::array<int, 3> const& g,
                   double factor)
            : p(p)
            , masks(masks)
            , fn(fn)
            , rho(rho)
            , bin(bin)
            , gx(g[0])
            , gy(g[1])
            , gz(g[2])
            , factor(factor)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            if (masks(i)) {
                int ix = fast_int_floor_kokkos(bin(i, 0));
                int iy = fast_int_floor_kokkos(bin(i, 2));
                int iz = fast_int_floor_kokkos(bin(i, 4));

                double ox = bin(i, 1);
                double oy = bin(i, 3);
                double oz = bin(i, 5);

                double aox = 1.0 - ox;
                double aoy = 1.0 - oy;
                double aoz = 1.0 - oz;

                int zbase = gx * gy * 2;

                if ((ix >= 0 && ix < gx - 1) && (iy >= 0 && iy < gy - 1) &&
                    (/*periodic_z ||*/ (iz >= 0 && iz < gz - 1))) {
                    double line_density =
                        aoz * rho(zbase + iz) + oz * rho(zbase + iz + 1);

                    double vx = line_density *
                                (aox * aoy * fn(((ix)*gy + (iy)) * 2) +
                                 ox * aoy * fn(((ix + 1) * gy + (iy)) * 2) +
                                 aox * oy * fn(((ix)*gy + (iy + 1)) * 2) +
                                 ox * oy * fn(((ix + 1) * gy + (iy + 1)) * 2));

                    double vy =
                        line_density *
                        (aox * aoy * fn(((ix)*gy + (iy)) * 2 + 1) +
                         ox * aoy * fn(((ix + 1) * gy + (iy)) * 2 + 1) +
                         aox * oy * fn(((ix)*gy + (iy + 1)) * 2 + 1) +
                         ox * oy * fn(((ix + 1) * gy + (iy + 1)) * 2 + 1));

                    p(i, 1) += factor * vx;
                    p(i, 3) += factor * vy;
                }
            }
        }
    };
}

Space_charge_2d_open_hockney::Space_charge_2d_open_hockney(
    Space_charge_2d_open_hockney_options const& ops)
    : Collective_operator("sc_2d_open_hockney", 1.0)
    , options(ops)
    , bunch_sim_id()
    , domain(ops.shape, {1.0, 1.0, 1.0})
    , doubled_domain(ops.doubled_shape, {1.0, 1.0, 1.0})
    , particle_bin()
    , ffts()
{}

void
Space_charge_2d_open_hockney::apply_impl(Bunch_simulator& sim,
                                         double time_step,
                                         Logger& logger)
{
    logger << "    Space charge 2d open hockney\n";

    scoped_simple_timer timer("sc2d_total");

    // construct the workspace for a new bunch simulator
    if (bunch_sim_id != sim.id()) {
        construct_workspaces(sim);
        bunch_sim_id = sim.id();
    }

    // apply to bunches
    for (size_t t = 0; t < 2; ++t) {
        for (size_t b = 0; b < sim[t].get_bunch_array_size(); ++b) {
            apply_bunch(sim[t][b], ffts[t][b], time_step, logger);
        }
    }
}

void
Space_charge_2d_open_hockney::apply_bunch(Bunch& bunch,
                                          Distributed_fft2d& fft,
                                          double time_step,
                                          Logger& logger)
{
    update_domain(bunch);

    get_local_charge_density(bunch); // [C/m^3]
    get_global_charge_density(bunch);

    get_green_fn2_pointlike();

    get_local_force2(fft);
    get_global_force2(fft.get_comm());

    auto fn_norm = get_normalization_force(bunch, fft);

    apply_kick(bunch, fn_norm, time_step);
}

void
Space_charge_2d_open_hockney::construct_workspaces(Bunch_simulator const& sim)
{
    scoped_simple_timer timer("sc2d_workspaces");

    auto const& s = options.doubled_shape;

    rho2 = karray1d_dev("rho2", s[0] * s[1] * 2 + s[2]);
    g2 = karray1d_dev("g2", s[0] * s[1] * 2);
    phi2 = karray1d_dev("phi2", s[0] * s[1] * 2);

    h_rho2 = Kokkos::create_mirror_view(rho2);
    h_phi2 = Kokkos::create_mirror_view(phi2);

    for (size_t t = 0; t < 2; ++t) {
        int num_local_bunches = sim[t].get_bunch_array_size();
        ffts[t] = std::vector<Distributed_fft2d>(num_local_bunches);

        for (size_t b = 0; b < num_local_bunches; ++b) {
            auto comm = sim[t][b].get_comm().divide(options.comm_group_size);

            ffts[t][b].construct({s[0], s[1]}, comm);
        }
    }
}

void
Space_charge_2d_open_hockney::update_domain(Bunch const& bunch)
{
    scoped_simple_timer timer("sc2d_domain");

    auto mean = Core_diagnostics::calculate_mean(bunch);
    auto std = Core_diagnostics::calculate_std(bunch, mean);

    const double tiny = 1.0e-10;

    const auto ix = Bunch::x;
    const auto iy = Bunch::y;
    const auto iz = Bunch::z;

    if ((std[ix] < tiny) && (std[iy] < tiny) && (std[iz] < tiny)) {
        throw std::runtime_error(
            "Space_charge_3d_open_hockney_eigen::update_domain: "
            "all three spatial dimensions have neglible extent");
    }

    std::array<double, 3> offset{mean[0], mean[2], mean[4]};

    std::array<double, 3> size{
        options.n_sigma * get_smallest_non_tiny(std[0], std[2], std[4], tiny),
        options.n_sigma * get_smallest_non_tiny(std[2], std[0], std[4], tiny),
        options.n_sigma * get_smallest_non_tiny(std[4], std[0], std[2], tiny)};

    std::array<double, 3> doubled_size{size[0] * 2.0, size[1] * 2.0, size[2]};

    domain = Rectangular_grid_domain(options.shape, size, offset, false);

    doubled_domain = Rectangular_grid_domain(
        options.doubled_shape, doubled_size, offset, false);
}

void
Space_charge_2d_open_hockney::get_local_charge_density(Bunch const& bunch)
{
    scoped_simple_timer timer("sc2d_local_rho");

    if (bunch.size() > particle_bin.extent(0))
        Kokkos::resize(particle_bin, bunch.size(), 6);

#ifdef SYNERGIA_ENABLE_CUDA
    deposit_charge_rectangular_2d_kokkos_scatter_view(
        rho2, doubled_domain, particle_bin, bunch);
#else
    deposit_charge_rectangular_2d_omp_reduce(
        rho2, doubled_domain, particle_bin, bunch);
#endif
}

void
Space_charge_2d_open_hockney::get_global_charge_density(Bunch const& bunch)
{
    // do nothing if the solver only has a single rank
    if (bunch.get_comm().size() == 1) return;

    scoped_simple_timer timer("sc2d_global_rho");

    auto dg = doubled_domain.get_grid_shape();

    simple_timer_start("sc2d_global_rho_copy");
    Kokkos::deep_copy(h_rho2, rho2);
    simple_timer_stop("sc2d_global_rho_copy");

    simple_timer_start("sc2d_global_rho_reduce");
    int err = MPI_Allreduce(MPI_IN_PLACE,
                            (void*)h_rho2.data(),
                            dg[0] * dg[1] * 2 + dg[2],
                            MPI_DOUBLE,
                            MPI_SUM,
                            bunch.get_comm());
    simple_timer_stop("sc2d_global_rho_reduce");

    if (err != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI error in Space_charge_2d_open_hockney"
            "(MPI_Allreduce in get_global_charge_density)");
    }

    simple_timer_start("sc2d_global_rho_copy");
    Kokkos::deep_copy(rho2, h_rho2);
    simple_timer_stop("sc2d_global_rho_copy");
}

void
Space_charge_2d_open_hockney::get_green_fn2_pointlike()
{
    scoped_simple_timer timer("sc2d_green_fn2");

    auto g = domain.get_grid_shape();
    auto h = doubled_domain.get_cell_size();
    auto dg = doubled_domain.get_grid_shape();

    alg_g2_pointlike alg(g2, g, dg, h);

    Kokkos::parallel_for(dg[0] * dg[1], alg);
    Kokkos::fence();
}

void
Space_charge_2d_open_hockney::get_local_force2(Distributed_fft2d& fft)
{
    scoped_simple_timer timer("sc2d_local_f");

    auto dg = doubled_domain.get_grid_shape();

    // FFT
    fft.transform(rho2, rho2);
    Kokkos::fence();

    fft.transform(g2, g2);
    Kokkos::fence();

    // zero phi2 when using multiple ranks
    if (fft.get_comm().size() > 1) {
        ku::alg_zeroer az{phi2};
        Kokkos::parallel_for(dg[0] * dg[1] * 2, az);
    }

    int lower = fft.get_lower();
    int upper = fft.get_upper();
    int nx = upper - lower;
    int offset = lower * dg[1];

    alg_cplx_multiplier alg(phi2, rho2, g2, offset);
    Kokkos::parallel_for(nx * dg[1], alg);
    Kokkos::fence();

    // inv fft
    fft.inv_transform(phi2, phi2);
    Kokkos::fence();
}

void
Space_charge_2d_open_hockney::get_global_force2(Commxx const& comm)
{
    // do nothing if the solver only has a single rank
    if (comm.size() == 1) return;

    scoped_simple_timer timer("sc2d_global_f");

    auto dg = doubled_domain.get_grid_shape();

    Kokkos::deep_copy(h_phi2, phi2);

    int err = MPI_Allreduce(MPI_IN_PLACE,
                            (void*)h_phi2.data(),
                            dg[0] * dg[1] * 2,
                            MPI_DOUBLE,
                            MPI_SUM,
                            comm);

    if (err != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI error in Space_charge_2d_open_hockney"
            "(MPI_Allreduce in get_global_electric_force2_allreduce)");
    }

    Kokkos::deep_copy(phi2, h_phi2);
}

double
Space_charge_2d_open_hockney::get_normalization_force(
    Bunch const& bunch,
    Distributed_fft2d const& fft)
{
    auto h = domain.get_cell_size();

    double hx = h[0];
    double hy = h[1];
    double hz = h[2];

    // volume element in integral
    double normalization = hx * hy;

    // dummy factor from weight0 of deposit.cc
    normalization *= 1.0 / (4.0 * pi * pconstants::epsilon0);
    normalization *= hz;

    // average line density for real particles
    normalization *= 1.0 / doubled_domain.get_physical_size()[2];

    // average line density for macro particles
    normalization *= 2.0 * doubled_domain.get_physical_size()[2] /
                     (hz * bunch.get_total_num());

    normalization *= bunch.get_particle_charge() * pconstants::e;

    normalization *= 1.0; // charge_density2.get_normalization();
    normalization *= 1.0; // green_fn2.get_normalization();

    normalization *= fft.get_roundtrip_normalization();

    return normalization;
}

void
Space_charge_2d_open_hockney::apply_kick(Bunch& bunch,
                                         double fn_norm,
                                         double time_step)
{
    scoped_simple_timer timer("sc2d_kick");

    // EGS ported AM changes for kicks in lab frame from 3D solver
    // AM: kicks are done in the z_lab frame
    // $\delta \vec{p} = \vec{F} \delta t = q \vec{E} \delta t$
    // delta_t_beam: [s] in beam frame
    //  See chapter 11, jackson electrodynamics, for field transformation
    //  from bunch frame (BF) to the lab frame (LF). Keep in mind that
    //  \vec{B}_BF=0.
    //  Ex_LF=gamma*Ex_BF, Ey_LF=gamma*Ey_BF, Ez_LF=Ez_BF
    //  Bx_LF=gamma*beta*Ey_BF, By_LF=-gamma*beta*Ex_BF, Bz_LF=Bz_BF=0
    //  Transverse Lorentz force in the lab frame:
    //        Fx_LF = q*(Ex_L-beta_z*By_LF)
    //              = q*gamma*(1-beta*beta_z)*Ex_BF
    //  Longitudinal Lorentz force in the lab frame:
    //        Fz = q*(Ez_LF+beta_x*By_LF-beta_y*Bx_LF)
    //           = q*(Ez_BF-gamma*beta*(beta_x*Ex_BF+beta_y*Ey_BF ))
    // In order to get a conservative approximation!:
    // The following approximations are done: beta_z=beta, beta_x=beta_y=0,
    // thus suppresing the particles' movement relative to the reference
    // particle. The same approximation was employed when the field in
    // the bunch frame was calculated.
    // Thus: Fx_LF=q*Ex_BF/gamma, Fz=q*Ez_BF

    double delta_t_beam =
        time_step / bunch.get_reference_particle().get_gamma();

    // unit_conversion: [N] = [kg m/s^2] to [Gev/c]
    double unit_conversion = pconstants::c / (1.0e9 * pconstants::e);

    // scaled p = p/p_ref
    double gamma = bunch.get_reference_particle().get_gamma();
    double beta = bunch.get_reference_particle().get_beta();
    double p_scale = 1.0 / bunch.get_reference_particle().get_momentum();

    // gamma*beta factor introduced here when we are no longer going to
    // the t_bunch frame. That factor was introduced in the fixed_z_lab
    // to fixed_t_bunch conversion.  Physically, gamma factor comes from
    // the lorentz expansion longitudinally in the bunch frame and beta
    // comes because the stored coordinate is c*dt whereas the actual
    // domain is beta*c*dt.
    double factor =
        unit_conversion * delta_t_beam * fn_norm * p_scale / (gamma * beta);

    auto parts = bunch.get_local_particles();
    auto masks = bunch.get_local_particle_masks();

    alg_kicker kicker(parts,
                      masks,
                      phi2,
                      rho2,
                      particle_bin,
                      doubled_domain.get_grid_shape(),
                      factor);

    Kokkos::parallel_for(bunch.size(), kicker);
    Kokkos::fence();
}
