#ifndef FF_HKICKER_H
#define FF_HKICKER_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/libFF/ff_patterned_propagator.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/foundation/physical_constants.h"

namespace FF_kicker
{
    template<class T>
    KOKKOS_INLINE_FUNCTION
    void kick(T const&x, T& xp, T const& y, T& yp, T const&, double const* kL)
    { FF_algorithm::thin_kicker_unit(xp, yp, kL); }


    template<class BunchT>
    void apply(Lattice_element_slice const& slice, BunchT& bunch)
    {
        scoped_simple_timer timer("libFF_kicker");

        auto const& elem = slice.get_lattice_element();
        const double length = slice.get_right() - slice.get_left();

        // hk and vk are the hk/vk under lattice reference momentum
        double    l = elem.get_double_attribute("l");
        double  hk0 = elem.get_double_attribute("hkick");
        double  vk0 = elem.get_double_attribute("vkick");
        double tilt = elem.get_double_attribute("tilt");

        // tilt
        double hk = cos(tilt)*hk0 - sin(tilt)*vk0;
        double vk = sin(tilt)*hk0 + cos(tilt)*vk0;

        //double k = elem.get_double_attribute(k_attr, 0.0);

        auto& ref_lattice = bunch.get_design_reference_particle();
        auto const& ref_bunch = bunch.get_reference_particle();

        double plattice = ref_lattice.get_momentum();
        double pbunch = ref_bunch.get_momentum();

        // scale is to scale the kick strength defined relative to the lattice momentum to
        // the scale of the bunch particles defined relative to the bunch momentum
        double scale = plattice/pbunch;

        // kick strength is defined as momentum change/reference momentum
        double b_hk = hk * ( ref_bunch.get_charge() / ref_lattice.get_charge() ) * scale;
        double b_vk = vk * ( ref_bunch.get_charge() / ref_lattice.get_charge() ) * scale;

        double k[2] = {hk, vk};
        double sk[2] = {b_hk, b_vk};

        using gsv_t = typename BunchT::gsv_t;
        using pp = FF_patterned_propagator<BunchT, gsv_t,
              kick<gsv_t>, kick<double>>;

        if ( close_to_zero(length) )
        {
            // the reference time is calculated with the design reference 
            // particle which is relative to the p-lattice.
            // also update the reference particle
            pp::get_reference_cdt_zero(ref_lattice, k);

            pp::apply_thin_kick(bunch, ParticleGroup::regular, sk);
            pp::apply_thin_kick(bunch, ParticleGroup::spectator, sk);
        }
        else
        {
            // simple drift-kick-drift scheme
            double simple_d = slice
                .get_lattice_element()
                .get_double_attribute("simple", 0.0);

            bool simple = fabs(simple_d) > 1e-16;

            // strength per unit length
            k[0] = k[0]/l;
            k[1] = k[1]/l;

            sk[0] = sk[0]/l;
            sk[1] = sk[1]/l;

            double pref = bunch.get_reference_particle().get_momentum();
            double mass = bunch.get_mass();

            if (simple)
            {
                // use un-scaled k_pul
                double ref_cdt = pp::get_reference_cdt_simple(
                        ref_lattice, length, k);

                pp::apply_simple_kick(bunch, ParticleGroup::regular, 
                        pref, mass, ref_cdt, length, sk);

                pp::apply_simple_kick(bunch, ParticleGroup::spectator, 
                        pref, mass, ref_cdt, length, sk);
            }
            else
            {
                // yoshida steps
                int steps = (int)elem.get_double_attribute(
                        "yoshida_steps", 4.0);

                // use drift for reference particle:
                // MadX uses drift for reference particle if a closed
                // orbit is not found. libFF also gives an option of
                // using the drift, but mostly for testing purposes
                int use_drift = (int)elem.get_double_attribute(
                        "cdt_use_drift", 0.0);

                // use un-scaled k_pul
                double ref_cdt = pp::get_reference_cdt_yoshida(
                        ref_lattice, length, k, steps, use_drift);

                pp::apply_yoshida_kick(bunch, ParticleGroup::regular,
                        pref, mass, ref_cdt, length, sk, steps);

                pp::apply_yoshida_kick(bunch, ParticleGroup::spectator,
                        pref, mass, ref_cdt, length, sk, steps);
            }

            bunch.get_reference_particle().increment_trajectory(length);
            // absolute time
            double const velocity = bunch.get_reference_particle().get_beta()*pconstants::c;
            bunch.get_reference_particle().increment_bunch_abs_time(length/velocity);

        }
    }
}

#endif // FF_HKICKER_H
