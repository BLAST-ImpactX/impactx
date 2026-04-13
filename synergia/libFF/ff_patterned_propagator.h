#ifndef SYNERGIA_LIBFF_PATTERNED_PROPAGATOR_H
#define SYNERGIA_LIBFF_PATTERNED_PROPAGATOR_H

namespace pp_impl
{
    template<class T>
    using kf_t = void(*)(T const&, T&, T const&, T&, T const&, double const*);
}

template< class BUNCH,
          class gsv_t,
          pp_impl::kf_t<gsv_t> kf_gsv,
          pp_impl::kf_t<double> kf_d,
          int COMP = 1 >
struct FF_patterned_propagator
{
    using bp_t = typename BUNCH::bp_t;
    using parts_t = typename bp_t::parts_t;
    using const_masks_t = typename bp_t::const_masks_t;

    struct thin_kicker
    {
        parts_t p;
        const_masks_t m;

        double k[COMP*2];

        thin_kicker(bp_t & bp, double const* str)
        : p(bp.parts), m(bp.masks)
        { for(int i=0; i<COMP*2; ++i) k[i] = str[i]; }

        KOKKOS_INLINE_FUNCTION
        void operator() (const int idx) const
        { 
            int i = idx * gsv_t::size();
            int mask = 0;
            for(int x=i; x<i+gsv_t::size(); ++x) mask |= m(x);

            if(mask) 
            {
                gsv_t p0(&p(i, 0));
                gsv_t p1(&p(i, 1));
                gsv_t p2(&p(i, 2));
                gsv_t p3(&p(i, 3));
                gsv_t p5(&p(i, 5));

                kf_gsv(p0, p1, p2, p3, p5, k);

                p1.store(&p(i, 1));
                p3.store(&p(i, 3));
            }
        }
    };

    struct simple_kicker
    {
        parts_t p;
        const_masks_t m;

        double pref, mass, ref_cdt, len;
        double k[COMP*2];

        simple_kicker(bp_t & bp,
                double pref, double mass, double ref_cdt, double len,
                double const* str)
        : p(bp.parts), m(bp.masks), pref(pref), mass(mass)
        , ref_cdt(ref_cdt), len(len)
        { for(int i=0; i<COMP*2; ++i) k[i] = str[i]; }

        KOKKOS_INLINE_FUNCTION
        void operator() (const int idx) const
        { 
            int i = idx * gsv_t::size();
            int mask = 0;
            for(int x=i; x<i+gsv_t::size(); ++x) mask |= m(x);

            if(mask) 
            {
                gsv_t p0(&p(i, 0));
                gsv_t p1(&p(i, 1));
                gsv_t p2(&p(i, 2));
                gsv_t p3(&p(i, 3));
                gsv_t p4(&p(i, 4));
                gsv_t p5(&p(i, 5));

                FF_algorithm::drift_unit(
                        p0, p1, p2, p3, p4, p5,
                        len*0.5, pref, mass, ref_cdt*0.5 
                );

                kf_gsv(p0, p1, p2, p3, p5, k);

                FF_algorithm::drift_unit(
                        p0, p1, p2, p3, p4, p5,
                        len*0.5, pref, mass, ref_cdt*0.5 
                );

                // drift updates x, y, and cdt
                // kick function updates xp and yp
                p0.store(&p(i, 0));
                p1.store(&p(i, 1));
                p2.store(&p(i, 2));
                p3.store(&p(i, 3));
                p4.store(&p(i, 4));
            }
        }
    };


    struct yoshida_kicker
    {
        parts_t p;
        const_masks_t m;

        double pref, mass, step_ref_cdt, step_len;
        double step_k[COMP*2];
        int steps;

        yoshida_kicker(bp_t & bp,
                double pref, double mass, double step_ref_cdt, 
                double step_len, double const* step_str, int steps)
        : p(bp.parts), m(bp.masks), pref(pref), mass(mass)
        , step_ref_cdt(step_ref_cdt), step_len(step_len)
        , steps(steps)
        { for(int i=0; i<COMP*2; ++i) step_k[i] = step_str[i]; }

        KOKKOS_INLINE_FUNCTION
        void operator() (const int idx) const
        { 
            int i = idx * gsv_t::size();
            int mask = 0;
            for(int x=i; x<i+gsv_t::size(); ++x) mask |= m(x);

            if(mask) 
            {
                gsv_t p0(&p(i, 0));
                gsv_t p1(&p(i, 1));
                gsv_t p2(&p(i, 2));
                gsv_t p3(&p(i, 3));
                gsv_t p4(&p(i, 4));
                gsv_t p5(&p(i, 5));

                FF_algorithm::yoshida6<gsv_t, kf_gsv, COMP>( 
                        p0, p1, p2, p3, p4, p5,
                        pref, mass, step_ref_cdt,
                        step_len, step_k, steps );

                p0.store(&p(i, 0));
                p1.store(&p(i, 1));
                p2.store(&p(i, 2));
                p3.store(&p(i, 3));
                p4.store(&p(i, 4));
            }
        }
    };

    static void apply_thin_kick(BUNCH& bunch, ParticleGroup pg, 
            double const* k)
    {
        if(!bunch.get_local_num(pg)) return;

        using exec = typename BUNCH::exec_space;
        auto range = Kokkos::RangePolicy<exec>(0, bunch.size_in_gsv(pg));

        thin_kicker tk(bunch.get_bunch_particles(pg), k);
        Kokkos::parallel_for(range, tk);
    }

    static void apply_simple_kick(BUNCH& bunch, ParticleGroup pg,
            double pref, double mass, double ref_cdt, 
            double len, double const* k)
    {
        if(!bunch.get_local_num(pg)) return;

        using exec = typename BUNCH::exec_space;
        auto range = Kokkos::RangePolicy<exec>(0, bunch.size_in_gsv(pg));

        simple_kicker sk(
            bunch.get_bunch_particles(pg),
            pref, mass, ref_cdt, len, k
        );

        Kokkos::parallel_for(range, sk);
    }


    static void apply_yoshida_kick(BUNCH& bunch, ParticleGroup pg,
            double pref, double mass, double ref_cdt,
            double len, double const* str, int steps)
    {
        if(!bunch.get_local_num(pg)) return;

        double step_len = len / steps;
        double step_ref_cdt = ref_cdt / steps;
        double step_str[2*COMP];
        for(int i=0; i<2*COMP; ++i) step_str[i] = str[i]*step_len;

        using exec = typename BUNCH::exec_space;
        auto range = Kokkos::RangePolicy<exec>(0, bunch.size_in_gsv(pg));

        yoshida_kicker yk( 
            bunch.get_bunch_particles(pg),
            pref, mass, step_ref_cdt, 
            step_len, step_str, steps 
        );

        Kokkos::parallel_for(range, yk);
    }

    static void get_reference_cdt_zero(
            Reference_particle& ref, 
            double const* k)
    {
        // propagate the bunch design reference particle
        double x  = ref.get_state()[0];
        double xp = ref.get_state()[1];
        double y  = ref.get_state()[2];
        double yp = ref.get_state()[3];
        // Don't use dp/p for getting reference time
        // double dpop = ref.get_state()[5];
        double dpop = 0.0;

        kf_d(x, xp, y, yp, dpop, k);

        ref.set_state_xp(xp);
        ref.set_state_yp(yp);
        ref.set_state_cdt(0.0);
    }

    static double get_reference_cdt_simple(
            Reference_particle& ref, 
            double len, double const* k)
    {
        double pref = ref.get_momentum();
        double mass = ref.get_mass();
        double cdt = 0.0;

        auto  st = ref.get_state();
        // EGS: use dp/p=0 for reference time calculation bacuse it might
        //   have been changed by acceleration
        st[5] = 0.0;

        // for >0 length, hk,vk is the strength/length of the kick
        double str[2*COMP];
        for(int i=0; i<2*COMP; ++i) str[i] = k[i]*len;

        FF_algorithm::drift_unit(
                st[0], st[1], st[2], 
                st[3], cdt, st[5], 
                len* 0.5, pref, mass, 0.0);

        kf_d(st[0], st[1], st[2], st[3], st[5], str);

        FF_algorithm::drift_unit(
                st[0], st[1], st[2], 
                st[3], cdt, st[5], 
                len* 0.5, pref, mass, 0.0);

        st[4] = cdt;
        ref.set_state(st);
        return cdt;
    }


    // non zero length
    static double get_reference_cdt_yoshida(
            Reference_particle& ref, 
            double len, double const* k, int steps, 
            bool use_drift = false)
    {
        double pref = ref.get_momentum();
        double mass = ref.get_mass();
        double cdt = 0.0;

        auto  st = ref.get_state();
        // EGS: use dp/p=0 for reference time calculation bacuse it might
        //   have been changed by acceleration
        st[5] = 0.0;

        // steps comes from base class, set in apply method
        double step_len = len / steps;

        // for >0 length, hk,vk is the strength/length of the kick
        double step_str[2*COMP];
        for(int i=0; i<2*COMP; ++i) step_str[i] = k[i]*step_len;

        // propagate
        if (use_drift)
        {
            FF_algorithm::drift_unit(
                    st[0], st[1], st[2], 
                    st[3], cdt, st[5], 
                    len, pref, mass, 0.0);

        }
        else
        {
            // hardcoded to use 6th order
            FF_algorithm::yoshida<double, kf_d, 6, COMP>( 
                    st[0], st[1], st[2], 
                    st[3], cdt, st[5],
                    pref, mass, 0.0,
                    step_len, step_str, steps );
        }

        st[4] = cdt;
        ref.set_state(st);
        return cdt;
    }

};

#endif
