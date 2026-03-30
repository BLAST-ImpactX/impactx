#ifndef FF_SEXTUPOLE_H
#define FF_SEXTUPOLE_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/libFF/ff_patterned_propagator.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/foundation/physical_constants.h"


namespace FF_sextupole
{
    template<class T>
    KOKKOS_INLINE_FUNCTION
    void kick(T const&x, T& xp, T const& y, T& yp, T const&, double const* kL)
    { FF_algorithm::thin_sextupole_unit(x, xp, y, yp, kL); }

    template<class BunchT>
    void apply(Lattice_element_slice const& slice, BunchT& bunch)
    {
        auto const& element = slice.get_lattice_element();
        double length = slice.get_right() - slice.get_left();

        // strength
        double k[2] = { 
            element.get_double_attribute("k2", 0.0),
            element.get_double_attribute("k2s", 0.0)
        };

        // tilting
        double tilt = element.get_double_attribute("tilt", 0.0);
        if (tilt != 0.0)
        {
            Kokkos::complex<double> ck2(k[0], k[1]);
            ck2 = ck2 * exp(Kokkos::complex<double>(0.0, -3.0*tilt));
            k[0] = ck2.real();
            k[1] = ck2.imag();
        }

        // scaling
        Reference_particle       & ref_l = bunch.get_design_reference_particle();
        Reference_particle const & ref_b = bunch.get_reference_particle();

        double brho_l = ref_l.get_momentum() / ref_l.get_charge();  // GV/c
        double brho_b = ref_b.get_momentum() 
                        * (1.0 + ref_b.get_state()[Bunch::dpop]) 
                        / ref_l.get_charge();  // GV/c

        double scale = brho_l / brho_b;

        k[0] *= scale;
        k[1] *= scale;

        using gsv_t = typename BunchT::gsv_t;
        using pp = FF_patterned_propagator<BunchT, gsv_t,
              kick<gsv_t>, kick<double>>;

        if ( close_to_zero(length) )
        {
            pp::get_reference_cdt_zero(ref_l, k);

            pp::apply_thin_kick(bunch, ParticleGroup::regular, k);
            pp::apply_thin_kick(bunch, ParticleGroup::spectator, k);
        }
        else
        {
            double pref = bunch.get_reference_particle().get_momentum();
            double mass = bunch.get_mass();

            // yoshida steps
            int steps = (int)element.get_double_attribute(
                    "yoshida_steps", 4.0);

            // reference cdt
            double ref_cdt = pp::get_reference_cdt_yoshida(
                    ref_l, length, k, steps);

            // propagate
            pp::apply_yoshida_kick(bunch, ParticleGroup::regular,
                    pref, mass, ref_cdt, length, k, steps);

            pp::apply_yoshida_kick(bunch, ParticleGroup::spectator,
                    pref, mass, ref_cdt, length, k, steps);

            // trajectory
            bunch.get_reference_particle().increment_trajectory(length);
            // absolute time
            double const velocity = bunch.get_reference_particle().get_beta()*pconstants::c;
            bunch.get_reference_particle().increment_bunch_abs_time(length/velocity);

        }

    }
}

#endif // FF_SEXTUPOLE_H
