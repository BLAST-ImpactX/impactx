#ifndef FF_RFCAVITY_H
#define FF_RFCAVITY_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/libFF/ff_patterned_propagator.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/foundation/physical_constants.h"

namespace rfcavity_impl {
    struct RFCavityParams {
        int nh;
        double mhp[12];

        double length;
        double str;
        double phi_s;
        double w_rf;

        double pref_b;
        double m_b;
        double new_pref_b;

        double ref_cdt_1;
        double ref_cdt_2;
    };

    template <class BunchT>
    struct PropThinRFCavity {
        using gsv_t = typename BunchT::gsv_t;

        typename BunchT::bp_t::parts_t p;
        typename BunchT::bp_t::const_masks_t m;
        const RFCavityParams rp;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int idx) const
        {
            int i = idx * gsv_t::size();
            int mask = 0;
            for (int x = i; x < i + gsv_t::size(); ++x)
                mask |= m(x);

            if (mask) {
                gsv_t p1(&p(i, 1));
                gsv_t p3(&p(i, 3));
                gsv_t p4(&p(i, 4));
                gsv_t p5(&p(i, 5));

                FF_algorithm::thin_rfcavity_unit(p1,
                                                 p3,
                                                 p4,
                                                 p5,
                                                 rp.w_rf,
                                                 rp.str,
                                                 rp.phi_s,
                                                 rp.m_b,
                                                 rp.pref_b,
                                                 rp.new_pref_b,
                                                 rp.mhp,
                                                 rp.nh);

                p1.store(&p(i, 1));
                p3.store(&p(i, 3));
                p4.store(&p(i, 4));
                p5.store(&p(i, 5));
            }
        }
    };

    template <class BunchT>
    struct PropRFCavity {
        using gsv_t = typename BunchT::gsv_t;

        typename BunchT::bp_t::parts_t p;
        typename BunchT::bp_t::const_masks_t m;
        const RFCavityParams rp;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int idx) const
        {
            int i = idx * gsv_t::size();
            int mask = 0;
            for (int x = i; x < i + gsv_t::size(); ++x)
                mask |= m(x);

            if (mask) {
                gsv_t p0(&p(i, 0));
                gsv_t p1(&p(i, 1));
                gsv_t p2(&p(i, 2));
                gsv_t p3(&p(i, 3));
                gsv_t p4(&p(i, 4));
                gsv_t p5(&p(i, 5));

                FF_algorithm::drift_unit(p0,
                                         p1,
                                         p2,
                                         p3,
                                         p4,
                                         p5,
                                         0.5 * rp.length,
                                         rp.pref_b,
                                         rp.m_b,
                                         rp.ref_cdt_1);

                FF_algorithm::thin_rfcavity_unit(p1,
                                                 p3,
                                                 p4,
                                                 p5,
                                                 rp.w_rf,
                                                 rp.str,
                                                 rp.phi_s,
                                                 rp.m_b,
                                                 rp.pref_b,
                                                 rp.new_pref_b,
                                                 rp.mhp,
                                                 rp.nh);

                FF_algorithm::drift_unit(p0,
                                         p1,
                                         p2,
                                         p3,
                                         p4,
                                         p5,
                                         0.5 * rp.length,
                                         rp.pref_b,
                                         rp.m_b,
                                         rp.ref_cdt_2);

                p0.store(&p(i, 0));
                p1.store(&p(i, 1));
                p2.store(&p(i, 2));
                p3.store(&p(i, 3));
                p4.store(&p(i, 4));
                p5.store(&p(i, 5));
            }
        }
    };

}

namespace FF_rfcavity {
    template <class BunchT>
    void
    apply(Lattice_element_slice const& slice, BunchT& bunch)
    {
        using namespace rfcavity_impl;

        scoped_simple_timer timer("libFF_rfcavity");

        RFCavityParams rp;

        rp.length = slice.get_right() - slice.get_left();
        Lattice_element const& elm = slice.get_lattice_element();

        int harmonic_number = elm.get_double_attribute("harmon", -1.0);
        double volt = elm.get_double_attribute("volt", 0.0);
        double lag = elm.get_double_attribute("lag", 0.0);
        double shunt = elm.get_double_attribute("shunt", 0.0);
        // freq is the synchronous frequency of the cavity in MHz
        double freq = elm.get_double_attribute("freq", -1.0);
        // delta_freq is the frequency offset from synchronous in MHz like freq
        double delta_freq = elm.get_double_attribute("delta_freq", 0.0);

        // harmonics,
        // mhp[h*3]   = harmonic multiple
        // mhp[h*3+1] = relative strength
        // mhp[h*3+2] = phase shift
        // nh = number of harmonics
        rp.nh = 1;
        rp.mhp[0] = 1.0;
        rp.mhp[1] = 1.0;
        rp.mhp[2] = 0.0;

        // harmonic multiple = 2
        if (elm.has_double_attribute("h2_str")) {
            rp.mhp[rp.nh * 3 + 0] = 2;
            rp.mhp[rp.nh * 3 + 1] = elm.get_double_attribute("h2_str", 0.0);
            rp.mhp[rp.nh * 3 + 2] = elm.get_double_attribute("h2_phase", 0.0);
            ++rp.nh;
        }

        // harmonic multiple = 3
        if (elm.has_double_attribute("h3_str")) {
            rp.mhp[rp.nh * 3 + 0] = 3;
            rp.mhp[rp.nh * 3 + 1] = elm.get_double_attribute("h3_str", 0.0);
            rp.mhp[rp.nh * 3 + 2] = elm.get_double_attribute("h3_phase", 0.0);
            ++rp.nh;
        }

        // harmonic multiple = 4
        if (elm.has_double_attribute("h4_str")) {
            rp.mhp[rp.nh * 3 + 0] = 4;
            rp.mhp[rp.nh * 3 + 1] = elm.get_double_attribute("h4_str", 0.0);
            rp.mhp[rp.nh * 3 + 2] = elm.get_double_attribute("h4_phase", 0.0);
            ++rp.nh;
        }

        rp.str = volt * 1.0e-3;

        // keep lag within the range of [0, 1).
        while (lag < 0.0) {
            lag += 1.0;
        }
        while (lag >= 1.0) {
            lag -= 1.0;
        }

        // elm.set_double_attribute("lag", lag);
        rp.phi_s = 2.0 * Kokkos::numbers::pi_v<double> * lag;
        rp.w_rf =
            2.0 * Kokkos::numbers::pi_v<double> * (freq + delta_freq) * 1.0e6;

        Reference_particle& ref_l = bunch.get_design_reference_particle();
        Reference_particle& ref_b = bunch.get_reference_particle();

        // The bunch particles momentum is with respect to the bunch reference
        // particle
        rp.pref_b = ref_b.get_momentum();
        rp.m_b = bunch.get_mass();

        rp.new_pref_b = FF_algorithm::thin_rfcavity_pnew(
            rp.pref_b, rp.m_b, rp.str, rp.phi_s);

        // reference_cdt uses the lattice reference particle
        // double reference_cdt = get_reference_cdt(length, ref_l);
        double ref_l_x = ref_l.get_state()[Bunch::x];
        double ref_l_xp = ref_l.get_state()[Bunch::xp];
        double ref_l_y = ref_l.get_state()[Bunch::y];
        double ref_l_yp = ref_l.get_state()[Bunch::yp];
        double ref_l_cdt = 0.0;
        // EGS: Don't use reference dpop, for calculation of cdt
        double ref_l_dpop = 0.0;

        double ref_l_p = ref_l.get_momentum();
        double ref_l_m = ref_l.get_mass();

        // double new_ref_l_p =
        //     FF_algorithm::thin_rfcavity_pnew(ref_l_p, ref_l_m, str, phi_s);

        // first half drift
        FF_algorithm::drift_unit(ref_l_x,
                                 ref_l_xp,
                                 ref_l_y,
                                 ref_l_yp,
                                 ref_l_cdt,
                                 ref_l_dpop,
                                 0.5 * rp.length,
                                 ref_l_p,
                                 ref_l_m,
                                 0.0);

        double total_ref_cdt = ref_l_cdt;
        rp.ref_cdt_1 = ref_l_cdt;
        ref_l_cdt = 0.0;

        // do not give it the new_ref_l_p because the xp and yp dont get scaled,
        // the momentum of the lattice reference particle remains unchanged,
        // only the dpop of the state has been changed
        FF_algorithm::thin_rfcavity_unit(ref_l_xp,
                                         ref_l_yp,
                                         ref_l_cdt,
                                         ref_l_dpop,
                                         rp.w_rf,
                                         rp.str,
                                         rp.phi_s,
                                         ref_l_m,
                                         ref_l_p,
                                         ref_l_p,
                                         rp.mhp,
                                         rp.nh);

        // second half drift
        FF_algorithm::drift_unit(ref_l_x,
                                 ref_l_xp,
                                 ref_l_y,
                                 ref_l_yp,
                                 ref_l_cdt,
                                 ref_l_dpop,
                                 0.5 * rp.length,
                                 ref_l_p,
                                 ref_l_m,
                                 0.0);

        total_ref_cdt += ref_l_cdt;
        rp.ref_cdt_2 = ref_l_cdt;

        // save the new state. The cdt has to be saved so that the transit time
	// in the RF cavity is accounted for. We don't save the new dp/p here but all
	// elements should not use it as their starting value in reference time
	// calculation.
        ref_l.set_state(
            ref_l_x, ref_l_xp, ref_l_y, ref_l_yp, total_ref_cdt, 0.0);

        // bunch particles
        auto apply = [&](ParticleGroup pg) {
            if (!bunch.get_local_num(pg)) return;

            auto parts = bunch.get_local_particles(pg);
            auto masks = bunch.get_local_particle_masks(pg);

            using exec = typename BunchT::exec_space;
            auto range = Kokkos::RangePolicy<exec>(0, bunch.size_in_gsv(pg));

            if (close_to_zero(rp.length)) {
                PropThinRFCavity<BunchT> rfcavity{parts, masks, rp};
                Kokkos::parallel_for(range, rfcavity);
            } else {
                PropRFCavity<BunchT> rfcavity{parts, masks, rp};
                Kokkos::parallel_for(range, rfcavity);
            }
        };

        apply(ParticleGroup::regular);
        apply(ParticleGroup::spectator);

        // updated four momentum
        Four_momentum fm = ref_b.get_four_momentum();
        fm.set_momentum(rp.new_pref_b);

        // update the bunch reference particle with the updated ref_p
        ref_b.set_four_momentum(fm);
        ref_b.increment_trajectory(rp.length);
        // absolute time
        double const velocity = ref_b.get_beta()*pconstants::c;
        ref_b.increment_bunch_abs_time(rp.length/velocity);

    }
}

#endif // FF_RFCAVITY_H
