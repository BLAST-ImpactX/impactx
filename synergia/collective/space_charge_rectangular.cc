#include "synergia/collective/space_charge_rectangular.h"

#include "deposit.h"
#include "space_charge_3d_kernels.h"

#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/kokkos_utils.h"
#include "synergia/utils/simple_timer.h"

namespace {
    void
    print_grid(Logger& logger,
               karray1d_dev const& grid,
               int x0,
               int x1,
               int y0,
               int y1,
               int z0,
               int z1,
               int gx,
               int gy,
               int gz,
               int off = 0)
    {
        karray1d_hst hgrid = Kokkos::create_mirror_view(grid);
        Kokkos::deep_copy(hgrid, grid);

        double sum = 0;

        int dim = grid.extent(0);
        for (int i = 0; i < dim; ++i) {
            sum += fabs(hgrid(i));
        }

#if 0
        for(int x=0; x<gx; ++x)
            for(int y=0; y<gy; ++y)
                sum += hgrid((x*gy + y)*2 + off);
#endif

        logger << std::setprecision(12) << std::scientific;
        logger << "      " << grid.label() << ".sum = " << sum << "\n";

        for (int z = z0; z < z1; ++z) {
            logger << "        " << z << ", ";

            for (int y = y0; y < y1; ++y) {
                logger << y << ", " << x0 << " | ";

                for (int x = x0; x < x1; ++x) {
                    logger << std::setprecision(12)
                           //<< hgrid(z*gx*gy + y*gx + x) << ", ";
                           << hgrid(x * gy * gz + y * gz + z) << ", ";
                }

                logger << "\n";
            }
        }
    }

    void
    print_statistics(Bunch& bunch, Logger& logger)
    {

        logger << "Bunch statistics: "
               << "num_valid = " << bunch.get_local_num()
               << ", size = " << bunch.size()
               << ", capacity = " << bunch.capacity()
               << ", total_num = " << bunch.get_total_num()
               << "\nMean and std: ";

        // print particles after propagate
        auto mean = Core_diagnostics::calculate_mean(bunch);
        auto std = Core_diagnostics::calculate_std(bunch, mean);

        logger << std::setprecision(16) << std::showpos << std::scientific
               << "\n"
            //<< "\nmean\tstd\n"
            ;

        for (int i = 0; i < 6; ++i)
            logger << mean[i] << ", " << std[i] << "\n";

        logger << "\n";

        for (int p = 0; p < 4; ++p)
            bunch.print_particle(p, logger);

        logger << "\n";
    }

    struct alg_phi {
        karray1d_dev phi;
        int lower, gy, gz;
        double igygz, igz;
        double px, py, pz;
        double gamma;

        alg_phi(karray1d_dev const& phi,
                int lower,
                int gy,
                int gz,
                double px,
                double py,
                double pz,
                double gamma)
            : phi(phi)
            , lower(lower)
            , gy(gy)
            , gz(gz)
            , igygz(1.0 / (gy * gz))
            , igz(1.0 / gz)
            , px(px)
            , py(py)
            , pz(pz)
            , gamma(gamma)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            const int ix = i * igygz;
            const int iy = (i - ix * gy * gz) * igz;
            const int iz = i - ix * gy * gz - iy * gz;

            int xt = ix + lower + 1;
            int yt = iy + 1;

            double denominator = Kokkos::numbers::pi_v<double> *
                                 Kokkos::numbers::pi_v<double> *
                                 (xt * xt / (px * px) + yt * yt / (py * py) +
                                  4.0 * iz * iz / (pz * pz * gamma * gamma));

            int base = ix * gy * gz + iy * gz + iz;

            phi(base * 2 + 0) /= denominator;
            phi(base * 2 + 1) /= denominator;
        }
    };

}

Space_charge_rectangular::Space_charge_rectangular(
    Space_charge_rectangular_options const& ops)
    : Collective_operator("sc_rectangular", 1.0)
    , options(ops)
    , domain(ops.shape, {1.0, 1.0, 1.0})
    , ffts()
{}

void
Space_charge_rectangular::apply_impl(Bunch_simulator& sim,
                                     double time_step,
                                     Logger& logger)
{
    logger << "    Space charge 3d open hockney\n";

    scoped_simple_timer timer("sc_rect_total");

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
Space_charge_rectangular::apply_bunch(Bunch& bunch,
                                      Distributed_fft3d_rect& fft,
                                      double time_step,
                                      Logger& logger)
{
    update_domain(bunch);

    get_local_charge_density(bunch);
    get_global_charge_density(bunch);

    double gamma = bunch.get_reference_particle().get_gamma();

    get_local_phi(fft, gamma);
    get_global_phi(fft);

    auto fn_norm = get_normalization_force();

    extract_force();
    apply_kick(bunch, fn_norm, time_step);
}

void
Space_charge_rectangular::construct_workspaces(Bunch_simulator const& sim)
{
    // shape
    auto const& s = options.shape;

    // fft objects
    for (size_t t = 0; t < 2; ++t) {
        int num_local_bunches = sim[t].get_bunch_array_size();
        ffts[t] = std::vector<Distributed_fft3d_rect>(num_local_bunches);

        for (size_t b = 0; b < num_local_bunches; ++b) {
            auto comm = sim[t][b].get_comm().divide(options.comm_group_size);

            ffts[t][b].construct(s, comm);
        }
    }

    // local workspaces
    int nz_cplx = Distributed_fft3d_rect::get_padded_shape_cplx(s[2]);

    rho = karray1d_dev("rho", s[0] * s[1] * s[2]);
    phi = karray1d_dev("phi", s[0] * s[1] * s[2]);

    phihat = karray1d_dev("phihat", s[0] * s[1] * nz_cplx * 2);

    h_rho = Kokkos::create_mirror_view(rho);
    h_phi = Kokkos::create_mirror_view(phi);

    // En is in the original domain
    enx = karray1d_dev("enx", s[0] * s[1] * s[2]);
    eny = karray1d_dev("eny", s[0] * s[1] * s[2]);
    enz = karray1d_dev("enz", s[0] * s[1] * s[2]);
}

void
Space_charge_rectangular::update_domain(Bunch const& bunch)
{
    double beta = bunch.get_reference_particle().get_beta();
    auto dsize = options.pipe_size;
    dsize[2] /= beta; // size in z_lab frame, longitudinal cdt coordinate

    // A.M physical_offsets of the domain should be rescaled too,
    // but in this case they are zero
    domain =
        Rectangular_grid_domain(options.shape, dsize, {0.0, 0.0, 0.0}, true);
}

void
Space_charge_rectangular::get_local_charge_density(Bunch const& bunch)
{
    scoped_simple_timer timer("sc_rect_local_rho");

    auto g = domain.get_grid_shape();
    // g[2] = Distributed_fft3d::get_padded_shape_real(g[2]);
    // g[2] = (g[2]/2+1)*2;

#ifdef SYNERGIA_ENABLE_OPENMP
    deposit_charge_rectangular_3d_omp_reduce_xyz(rho, domain, g, bunch);
#else
    deposit_charge_rectangular_3d_kokkos_scatter_view_xyz(
        rho, domain, g, bunch);
#endif
}

void
Space_charge_rectangular::get_global_charge_density(Bunch const& bunch)
{
    // do nothing if the bunch occupis a single rank
    if (bunch.get_comm().size() == 1) return;

    scoped_simple_timer timer("sc_rect_global_rho");

    auto g = domain.get_grid_shape();

    simple_timer_start("sc_rect_global_rho_copy");
    Kokkos::deep_copy(h_rho, rho);
    simple_timer_stop("sc_rect_global_rho_copy");

    simple_timer_start("sc_rect_global_rho_reduce");
    int err = MPI_Allreduce(MPI_IN_PLACE,
                            (void*)h_rho.data(),
                            h_rho.extent(0),
                            MPI_DOUBLE,
                            MPI_SUM,
                            bunch.get_comm());
    simple_timer_stop("sc_rect_global_rho_reduce");

    if (err != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI error in Space_charge_rectangular"
            "(MPI_Allreduce in get_global_charge_density)");
    }

    simple_timer_start("sc_rect_global_rho_copy");
    Kokkos::deep_copy(rho, h_rho);
    simple_timer_stop("sc_rect_global_rho_copy");
}

void
Space_charge_rectangular::get_local_phi(Distributed_fft3d_rect& fft,
                                        double gamma)
{
    scoped_simple_timer timer("sc_rect_local_phi");

    ku::alg_zeroer az{phihat};
    Kokkos::parallel_for(phihat.extent(0), az);

    fft.transform(rho, phihat);

    int lower = fft.get_lower();
    int upper = fft.get_upper();
    auto ps = options.pipe_size;
    auto g = domain.get_grid_shape();
    int gy = g[1];
    int gz_padded_cplx = g[2] / 2 + 1;

    alg_phi aphi(phihat, lower, gy, gz_padded_cplx, ps[0], ps[1], ps[2], gamma);

    Kokkos::parallel_for((upper - lower) * gy * gz_padded_cplx, aphi);

    fft.inv_transform(phihat, phi);
}

void
Space_charge_rectangular::get_global_phi(Distributed_fft3d_rect const& fft)
{
    // do nothing if the solver only has a single rank
    if (fft.get_comm().size() == 1) return;

    scoped_simple_timer timer("sc_rect_global_phi");

    Kokkos::deep_copy(h_phi, phi);

    int err = MPI_Allreduce(MPI_IN_PLACE,
                            (void*)h_phi.data(),
                            h_phi.extent(0),
                            MPI_DOUBLE,
                            MPI_SUM,
                            fft.get_comm());

    if (err != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI error in Space_charge_3d_open_hockney"
            "(MPI_Allreduce in get_global_electric_force2_allreduce)");
    }

    Kokkos::deep_copy(phi, h_phi);
}

double
Space_charge_rectangular::get_normalization_force()
{
    auto g = domain.get_grid_shape();
    return 1.0 / (4.0 * g[0] * g[1] * g[2] * pconstants::epsilon0);
}

void
Space_charge_rectangular::extract_force()
{
    scoped_simple_timer timer("sc_rect_get_en");

    auto g = domain.get_grid_shape();
    auto h = domain.get_cell_size();

    // phi is in (gx, gy, gz)
    // en{x|y|z} is in (gx, gy, gz)
    sc3d_kernels::xyz::alg_force_extractor alg(phi, enx, eny, enz, g, h);
    Kokkos::parallel_for(g[0] * g[1] * g[2], alg);
    Kokkos::fence();
}

void
Space_charge_rectangular::apply_kick(Bunch& bunch,
                                     double fn_norm,
                                     double time_step)
{
    scoped_simple_timer timer("sc_rect_kick");

    auto ref = bunch.get_reference_particle();

    double q = bunch.get_particle_charge() * pconstants::e;
    double m = bunch.get_mass();

    double gamma = ref.get_gamma();
    double beta = ref.get_beta();
    double pref = ref.get_momentum();

    double unit_conversion = pconstants::c / (1e9 * pconstants::e);
    double factor = unit_conversion * q * time_step * fn_norm /
                    (pref * gamma * gamma * beta);

    auto parts = bunch.get_local_particles();
    auto masks = bunch.get_local_particle_masks();

    auto g = domain.get_grid_shape();
    auto h = domain.get_cell_size();
    auto l = domain.get_left();

    sc3d_kernels::xyz::alg_kicker kicker(
        parts, masks, enx, eny, enz, g, h, l, factor, pref, m);

    Kokkos::parallel_for(bunch.size(), kicker);
    Kokkos::fence();
}
