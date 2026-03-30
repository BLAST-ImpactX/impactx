#include "space_charge_2d_kv.h"

#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/physical_constants.h"

#include "synergia/utils/simple_timer.h"

namespace {
    // returns the "normalized" electric field in the rest frame of the bunch,
    // in inverse meters.  To get the field [V/m], this must be multiplied
    // by Q/(2 pi epsilon_o), where Q is the line density of charge [C/m]
    // (in rest frame).

    // returns Ex, Ey field for a unit kv charge
    //
    // from Miguel A. Furman, Compact Complex Expressions for the Electric
    // Field of 2-D Elliptical Charge Distributions: LBL-34682, CBP Note 014,
    // PEP-II/AP Note 34-93
    //
    // If strictly_linear is true, use the linear field over all space.
    // Otherwise use the Furman expression outside of the ellipse boundaries.

    KOKKOS_INLINE_FUNCTION
    Kokkos::complex<double>
    unit_efield(double x, double y, bool linear, double sigma_x, double sigma_y)
    {
        Kokkos::complex<double> E;

        // projected std on 1 axis of 2d uniform distribution of
        // radius R = R/2.
        double a = 2.0 * sigma_x;
        double b = 2.0 * sigma_y;

        bool inside = (x / a) * (x / a) + (y / b) * (y / b) < 1.0;

        int xneg = x < 0.0 ? -1 : 1;
        int yneg = y < 0.0 ? -1 : 1;

        if (linear || inside) {
            E.real() = 4 * x / (a * (a + b));
            E.imag() = 4 * y / (b * (a + b));
        } else {
            // for the exterior, we have to reflect to the
            // first quadrant
            x *= xneg;
            y *= yneg;

            Kokkos::complex<double> zbar(x, -y);
            E = 4.0 / (zbar + sqrt(zbar * zbar - a * a + b * b));

            E.real() *= xneg;
            E.imag() *= yneg;
        }

        return E;
    }

    struct sc_kv_uniform {
        Particles p;
        ConstParticleMasks masks;

        double offx;
        double offy;

        double sigma_x;
        double sigma_y;

        double factor;
        double line_charge_density;

        bool linear;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            if (masks(i)) {
                double x = p(i, 0) - offx;
                double y = p(i, 2) - offy;

                auto E = unit_efield(x, y, linear, sigma_x, sigma_y);

                p(i, 1) += E.real() * factor * line_charge_density;
                p(i, 3) += E.imag() * factor * line_charge_density;
            }
        }
    };

    struct sc_kv_gaussian {
        Particles p;
        ConstParticleMasks masks;

        double offx;
        double offy;
        double offz;

        double sigma_x;
        double sigma_y;
        double sigma_cdt;

        double factor;
        double coeff;

        bool linear;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            if (masks(i)) {
                double x = p(i, 0) - offx;
                double y = p(i, 2) - offy;
                double z = p(i, 4) - offz;

                // coeff = total_q / (sqrt(2.0*pi)*sigma_cdt*beta);
                double line_charge_density =
                    exp(-z * z / (2.0 * sigma_cdt * sigma_cdt)) * coeff;

                auto E = unit_efield(x, y, linear, sigma_x, sigma_y);

                p(i, 1) += E.real() * factor * line_charge_density;
                p(i, 3) += E.imag() * factor * line_charge_density;
            }
        }
    };
}

Space_charge_2d_kv::Space_charge_2d_kv(Space_charge_2d_kv_options const& opts)
    : Collective_operator("sc_2d_kv", 1.0), opts(opts)
{}

void
Space_charge_2d_kv::apply_impl(Bunch_simulator& sim,
                               double time_step,
                               Logger& logger)
{
    logger << "    Space charge 2d kv\n";

    scoped_simple_timer timer("sc2d_kv_total");

    // apply to bunches
    for (size_t t = 0; t < 2; ++t) {
        for (size_t b = 0; b < sim[t].get_bunch_array_size(); ++b) {
            apply_bunch(sim[t][b], time_step, logger);
        }
    }
}

void
Space_charge_2d_kv::apply_bunch(Bunch& bunch, double time_step, Logger& logger)
{
    auto mean = Core_diagnostics::calculate_mean(bunch);
    auto std = Core_diagnostics::calculate_std(bunch, mean);

    // dp/p kick =
    //
    //   N r_p * (1/gamma**2) *        delta-t   *         (1/(beta*gamma)) *
    //   unit_E_field
    //              E-B cancellation                               ^
    //                                                             |
    //                                                       1/p to get dp/p
    //

    double beta = bunch.get_reference_particle().get_beta();
    double gamma = bunch.get_reference_particle().get_gamma();

    double factor = pconstants::rp * pconstants::c * time_step /
                    (gamma * gamma * gamma * beta);

    // set longitudinal density depending on the longitudinal flag
    double total_q = bunch.get_real_num() * bunch.get_particle_charge();

    double offx = 0;
    double offy = 0;
    double offz = 0;

    if (!opts.strictly_centered) {
        offx = mean(0);
        offy = mean(2);
        offz = mean(4);
    }

    auto parts = bunch.get_local_particles();
    auto masks = bunch.get_local_particle_masks();

    if (opts.longitudinal_distribution ==
        Space_charge_2d_kv_options::LD::uniform) {
        // if longitudinally uniform, set the length based on
        // sqrt(12)*sigma_z because you might not necessarily
        // have a periodic bunch with a z_period_length and your
        // actual bunch may not fill the entire bucket.
        double bunch_length = std::sqrt(12) * std(4) * beta;
        double line_charge_density = total_q / bunch_length;

        sc_kv_uniform alg{parts,
                          masks,
                          offx,
                          offy,
                          std(0),
                          std(2),
                          factor,
                          line_charge_density,
                          opts.strictly_linear};

        Kokkos::parallel_for(bunch.size(), alg);
    } else {
        double coeff = total_q / (sqrt(2.0 * Kokkos::numbers::pi_v<double>) *
                                  std(4) * beta);

        sc_kv_gaussian alg{parts,
                           masks,
                           offx,
                           offy,
                           offz,
                           std(0),
                           std(2),
                           std(4),
                           factor,
                           coeff,
                           opts.strictly_linear};

        Kokkos::parallel_for(bunch.size(), alg);
    }

    Kokkos::fence();
}
