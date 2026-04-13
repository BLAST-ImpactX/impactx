#ifndef FF_QUADRUPOLE_H
#define FF_QUADRUPOLE_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/foundation/physical_constants.h"


namespace quad_impl
{
    // max multipole moments
    constexpr const int max_mp_order = 8;

    // pure thin quadrupole kick
    template<class T>
    KOKKOS_INLINE_FUNCTION
    void kick(T const&x, T& xp, T const& y, T& yp, T const&, double const* kL)
    { FF_algorithm::thin_quadrupole_unit(x, xp, y, yp, kL); }

    // thin quadrupole kick with high order corrections
    template<class T>
    KOKKOS_INLINE_FUNCTION
    void cf_kick(T const&x, T& xp, T const& y, T& yp, T const&, double const* kL)
    { 
        FF_algorithm::thin_quadrupole_unit(x, xp, y, yp, kL); 

        if (kL[2] || kL[3])   FF_algorithm::thin_quadrupole_unit(x, xp, y, yp, kL+2);
        if (kL[4] || kL[5])   FF_algorithm::thin_sextupole_unit(x, xp, y, yp, kL+4);
        if (kL[6] || kL[7])   FF_algorithm::thin_octupole_unit(x, xp, y, yp, kL+6);
        if (kL[8] || kL[9])   FF_algorithm::thin_magnet_unit(x, xp, y, yp, kL+8, 5);
        if (kL[10] || kL[11]) FF_algorithm::thin_magnet_unit(x, xp, y, yp, kL+10, 6);
        if (kL[12] || kL[13]) FF_algorithm::thin_magnet_unit(x, xp, y, yp, kL+12, 7);
        if (kL[14] || kL[15]) FF_algorithm::thin_magnet_unit(x, xp, y, yp, kL+14, 8);
    }

    template<class BP>
    struct PropQuadThin
    {
        typename BP::parts_t p;
        typename BP::const_masks_t masks;
        double k[2];
        double xoff, yoff;

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            if (masks(i))
            {
                auto x = p(i, 0) - xoff;
                auto y = p(i, 2) - yoff;

                FF_algorithm::thin_quadrupole_unit(
                        x, p(i, 1), y, p(i, 3), k);
            }
        }
    };

    template<class BP>
    struct PropQuadThinSimd
    {
        using gsv_t = typename BP::gsv_t;

        typename BP::parts_t p;
        typename BP::const_masks_t masks;
        double k[2];
        gsv_t xoff, yoff;

        KOKKOS_INLINE_FUNCTION
        void operator()(const int idx) const
        {
            int i = idx * gsv_t::size();

            int m = 0;
            for(int x=i; x<i+gsv_t::size(); ++x) m |= masks(x);

            if (m)
            {
                gsv_t p0(&p(i, 0));
                gsv_t p1(&p(i, 1));
                gsv_t p2(&p(i, 2));
                gsv_t p3(&p(i, 3));

                p0 = p0 - xoff;
                p2 = p2 - yoff;

                FF_algorithm::thin_quadrupole_unit(
                        p0, p1, p2, p3, k);

                p1.store(&p(i, 1));
                p3.store(&p(i, 3));
            }
        }
    };

    template<class BP>
    struct PropCFQuadThinSimd
    {
        using gsv_t = typename BP::gsv_t;

        typename BP::parts_t p;
        typename BP::const_masks_t masks;
        kt::arr_t<double, 2*max_mp_order> k;
        gsv_t xoff, yoff;

        KOKKOS_INLINE_FUNCTION
        void operator()(const int idx) const
        {
            int i = idx * gsv_t::size();

            int m = 0;
            for(int x=i; x<i+gsv_t::size(); ++x) m |= masks(x);

            if (m)
            {
                gsv_t p0(&p(i, 0));
                gsv_t p1(&p(i, 1));
                gsv_t p2(&p(i, 2));
                gsv_t p3(&p(i, 3));

                p0 = p0 - xoff;
                p2 = p2 - yoff;

                cf_kick(p0, p1, p2, p3, gsv_t(0.0), k.data);

                p1.store(&p(i, 1));
                p3.store(&p(i, 3));
            }
        }
    };


    template<class BP>
    struct PropQuad
    {
        typename BP::parts_t p;
        typename BP::const_masks_t masks;
        int steps;
        double xoff, yoff;
        double ref_p, ref_m, step_ref_t, step_l;
        double step_k[2];

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            using part_t = typename BP::part_t;

            if (masks(i))
            {
                p(i, 0) -= xoff;
                p(i, 2) -= yoff;

                FF_algorithm::yoshida6<part_t, kick<part_t>, 1> (
                            p(i,0), p(i,1), p(i,2), 
                            p(i,3), p(i,4), p(i,5), 
                            ref_p, ref_m, step_ref_t, 
                            step_l, step_k, steps );

                p(i, 0) += xoff;
                p(i, 2) += yoff;
            }
        }
    };

    template<class BP>
    struct PropQuadSimd
    {
        using gsv_t = typename BP::gsv_t;

        typename BP::parts_t p;
        typename BP::const_masks_t masks;
        int steps;
        gsv_t xoff, yoff;
        double ref_p, ref_m, step_ref_t, step_l;
        double step_k[2];

        KOKKOS_INLINE_FUNCTION
        void operator()(const int idx) const
        {
            int i = idx * gsv_t::size();

            int m = 0;
            for(int x=i; x<i+gsv_t::size(); ++x) m |= masks(x);

            if (m)
            {
                gsv_t p0(&p(i, 0));
                gsv_t p1(&p(i, 1));
                gsv_t p2(&p(i, 2));
                gsv_t p3(&p(i, 3));
                gsv_t p4(&p(i, 4));
                gsv_t p5(&p(i, 5));

                p0 = p0 - xoff;
                p2 = p2 - yoff;

                FF_algorithm::yoshida6<gsv_t, kick<gsv_t>, 1>( 
                        p0, p1, p2, p3, p4, p5,
                        ref_p, ref_m, step_ref_t, 
                        step_l, step_k, steps );

                p0 = p0 + xoff;
                p2 = p2 + yoff;

                p0.store(&p(i, 0));
                p1.store(&p(i, 1));
                p2.store(&p(i, 2));
                p3.store(&p(i, 3));
                p4.store(&p(i, 4));
            }
        }
    };

    template<class BP>
    struct PropCFQuadSimd
    {
        using gsv_t = typename BP::gsv_t;

        typename BP::parts_t p;
        typename BP::const_masks_t masks;
        int steps;
        gsv_t xoff, yoff;
        double ref_p, ref_m, step_ref_t, step_l;
        kt::arr_t<double, 2*max_mp_order> step_k;

        KOKKOS_INLINE_FUNCTION
        void operator()(const int idx) const
        {
            int i = idx * gsv_t::size();

            int m = 0;
            for(int x=i; x<i+gsv_t::size(); ++x) m |= masks(x);

            if (m)
            {
                gsv_t p0(&p(i, 0));
                gsv_t p1(&p(i, 1));
                gsv_t p2(&p(i, 2));
                gsv_t p3(&p(i, 3));
                gsv_t p4(&p(i, 4));
                gsv_t p5(&p(i, 5));

                p0 = p0 - xoff;
                p2 = p2 - yoff;

                FF_algorithm::yoshida6<gsv_t, cf_kick<gsv_t>, max_mp_order>( 
                        p0, p1, p2, p3, p4, p5,
                        ref_p, ref_m, step_ref_t, 
                        step_l, step_k.data, steps );

                p0 = p0 + xoff;
                p2 = p2 + yoff;

                p0.store(&p(i, 0));
                p1.store(&p(i, 1));
                p2.store(&p(i, 2));
                p3.store(&p(i, 3));
                p4.store(&p(i, 4));
            }
        }
    };


    inline double get_reference_cdt(
            double length, int steps, 
            double const* kn,
            double xoff, double yoff,
            Reference_particle& ref)
    {
        double x  = ref.get_state()[Bunch::x];
        double xp = ref.get_state()[Bunch::xp];
        double y  = ref.get_state()[Bunch::y];
        double yp = ref.get_state()[Bunch::yp];
        double cdt  = 0.0;
        // EGS: don't use incomoing dp/p for reference time calculation
        //double dpop = ref.get_state()[Bunch::dpop];
        double dpop = 0.0;

        x -= xoff;
        y -= yoff;

        if (close_to_zero(length))
        {
            FF_algorithm::thin_quadrupole_unit(x, xp,  y, yp, kn);
        }
        else
        {
            double ref_p = ref.get_momentum();
            double ref_m = ref.get_mass();
            double ref_cdt = 0.0;

            FF_algorithm::yoshida6<double, cf_kick<double>, max_mp_order>
                    ( x, xp, y, yp, cdt, dpop,
                      ref_p, ref_m, ref_cdt,
                      length/steps, kn, steps );
        }

        x += xoff;
        y += yoff;

        ref.set_state(x, xp, y, yp, cdt, dpop);

        return cdt;
    }
}


namespace FF_quadrupole
{
    template<class BunchT>
    inline void apply(Lattice_element_slice const& slice, BunchT & bunch)
    {
        using namespace quad_impl;

        scoped_simple_timer timer("libFF_quad");

        // element
        auto const& ele = slice.get_lattice_element();

        // length
        double length = slice.get_right() - slice.get_left();

        // offsets
        const double xoff = ele.get_double_attribute("hoffset", 0.0);
        const double yoff = ele.get_double_attribute("voffset", 0.0);

        // tilt
        double tilt = ele.get_double_attribute("tilt", 0.0);

        // quadrupole strength
        // kn[0], kn[1] are the quadrupole strength
        kt::arr_t<double, 2*max_mp_order> kn;

        kn[0] = ele.get_double_attribute("k1", 0.0);
        kn[1] = ele.get_double_attribute("k1s", 0.0);

        // quad strength is strictly k1 as it appears in the element
        double quad_str = kn[0];

        // rot_angle is phase of ck1/(order+1) (order = 1 for quadrupole)
        // for pure skew quad, phase is -pi/2, order=1 so rotation = -pi/4
        // so rot_angle is the rotation needed to turn a normal quad
        // into one with k1 and k1s components.
        std::complex<double> ck1(kn[0], -kn[1]);
        double rot = std::arg(ck1)/2.0 + tilt;

        if (tilt != 0.0)
        {
            std::complex<double> ck2(kn[0], kn[1]);
            ck2 = ck2 * exp(std::complex<double>(0.0, -2.0*tilt));

            kn[0] = ck2.real();
            kn[1] = ck2.imag();
        }

        // multipole moments 
        // kn[2] ... kn[2*max_mp_order-1]
        bool has_mp = false;

        std::string a_attr = "a0";
        std::string b_attr = "b0";

        for(int i=1; i<max_mp_order; ++i)
        {
            a_attr[1] = '0' + i;
            b_attr[1] = '0' + i;

            kn[i*2+0] = ele.get_double_attribute(b_attr, 0.0);
            kn[i*2+1] = ele.get_double_attribute(a_attr, 0.0);

            if (kn[i*2+0] || kn[i*2+1]) 
            {
                has_mp = true;

                double coeff = FF_algorithm::factorial(i);

                kn[i*2+0] *= quad_str * coeff;
                kn[i*2+1] *= quad_str * coeff;

                if (rot || tilt)
                {
                    // the first rotate counts for the tilt in this mp
                    // second rotate is to compensate the rotate in quad
                    auto ck = std::complex<double>(kn[i*2+0], kn[i*2+1]);
                    ck *= exp(std::complex<double>(0.0, -(i+1)*tilt));
                    ck *= exp(std::complex<double>(0.0, -(i+1)*rot));

                    kn[i*2+0] = ck.real();
                    kn[i*2+1] = ck.imag();
                }
            }
        }

        // scaling
        Reference_particle       & ref_l = bunch.get_design_reference_particle();
        Reference_particle const & ref_b = bunch.get_reference_particle();

        double brho_l = ref_l.get_momentum() / ref_l.get_charge();  // GV/c
        double brho_b = ref_b.get_momentum()
                        * (1.0 + ref_b.get_state()[Bunch::dpop])
                        / ref_l.get_charge();  // GV/c

        double scale = brho_l / brho_b;
        for(auto & k : kn) k *= scale;

        if (close_to_zero(length))
        {
            // propagate the design reference particle
            get_reference_cdt(length, 0, kn.data, xoff, yoff, ref_l);

            // propagate the bunch particles
            auto apply = [&](ParticleGroup pg) {
                auto bp = bunch.get_bunch_particles(pg);
                if (!bp.num_valid()) return;

                using exec = typename BunchT::exec_space;

#if LIBFF_USE_GSV
                // thin quadrupole cant have mp components
                PropQuadThinSimd<typename BunchT::bp_t> pqt{ 
                    bp.parts, bp.masks, {kn[0], kn[1]}, xoff, yoff };

                auto range = Kokkos::RangePolicy<exec>(0, bp.size_in_gsv());
                Kokkos::parallel_for(range, pqt);
#else
                auto range = Kokkos::RangePolicy<exec>(0, bp.size());
                PropQuadThin<typename BunchT::bp_t> pqt{ 
                    bp.parts, bp.masks, {kn[0], kn[1]}, xoff, yoff };
                Kokkos::parallel_for(range, pqt);
#endif
            };

            // apply
            apply(ParticleGroup::regular);
            apply(ParticleGroup::spectator);
        }
        else
        {
            // yoshida steps
            int steps = (int)slice
                .get_lattice_element()
                .get_double_attribute("yoshida_steps", 4.0);

            //auto k2 = kn;

            // perstep k
            for(auto & k : kn) k *= length/steps;

            // params
            double ref_p = ref_b.get_momentum();
            double ref_m = ref_b.get_mass();
            double ref_t = get_reference_cdt(
                    length, steps, kn.data, xoff, yoff, ref_l);

            // bunch particles
            auto apply = [&](ParticleGroup pg) {
                auto bp = bunch.get_bunch_particles(pg);
                if (!bp.num_valid()) return;

                using exec = typename BunchT::exec_space;

#if LIBFF_USE_GSV
                if (has_mp)
                {
#if 1
                    PropCFQuadSimd<typename BunchT::bp_t> pq{
                        bp.parts, bp.masks, steps,
                        xoff, yoff, ref_p, ref_m,
                        ref_t/steps, length/steps, kn
                    };

                    auto range = Kokkos::RangePolicy<exec>(0, bp.size_in_gsv());
                    Kokkos::parallel_for(range, pq);
#endif

#if 0
                    // drift - kick - drift for testing
                    // also needs to make a copy of kn before scaling
                    // and comment out quad_unit in cf_kick
                    ref_t *= 0.5;
                    length *= 0.5;
                    k2[0] *= length/steps;
                    k2[1] *= length/steps;

                    auto range = Kokkos::RangePolicy<exec>(
                            0, bp.size_in_gsv());

                    PropQuadSimd<typename BunchT::bp_t> pq1{
                        bp.parts, bp.masks, steps,
                        xoff, yoff, ref_p, ref_m,
                        ref_t/steps, length/steps, 
                        {k2[0], k2[1]}
                    };

                    Kokkos::parallel_for(range, pq1);

                    PropCFQuadThinSimd<typename BunchT::bp_t> pqt{ 
                        bp.parts, bp.masks, k2, xoff, yoff };

                    Kokkos::parallel_for(range, pqt);

                    PropQuadSimd<typename BunchT::bp_t> pq2{
                        bp.parts, bp.masks, steps,
                        xoff, yoff, ref_p, ref_m,
                        ref_t/steps, length/steps, 
                        {k2[0], k2[1]}
                    };

                    Kokkos::parallel_for(range, pq2);
#endif
                }
                else
                {
                    PropQuadSimd<typename BunchT::bp_t> pq{
                        bp.parts, bp.masks, steps,
                        xoff, yoff, ref_p, ref_m,
                        ref_t/steps, length/steps, 
                        {kn[0], kn[1]}
                    };

                    auto range = Kokkos::RangePolicy<exec>(0, bp.size_in_gsv());
                    Kokkos::parallel_for(range, pq);
                }
#else
                if (has_mp)
                {
                    PropQuad<typename BunchT::bp_t> pq{
                        bp.parts, bp.masks, steps,
                        xoff, yoff, ref_p, ref_m,
                        ref_t/steps, length/steps,
                        {kn[0], kn[1]}
                    };

                    auto range = Kokkos::RangePolicy<exec>(0, bp.size());
                    Kokkos::parallel_for(range, pq);
                }
                else
                {
                    PropQuad<typename BunchT::bp_t> pq{
                        bp.parts, bp.masks, steps,
                        xoff, yoff, ref_p, ref_m,
                        ref_t/steps, length/steps,
                        {kn[0], kn[1]}
                    };

                    auto range = Kokkos::RangePolicy<exec>(0, bp.size());
                    Kokkos::parallel_for(range, pq);
                }
#endif
            };

            // apply
            apply(ParticleGroup::regular);
            apply(ParticleGroup::spectator);

            // advance the ref_part
            bunch.get_reference_particle().increment_trajectory(length);

            // absolute time
            double const velocity = bunch.get_reference_particle().get_beta()*pconstants::c;
            bunch.get_reference_particle().increment_bunch_abs_time(length/velocity);

        }

        Kokkos::fence();
    }
}

#endif // FF_QUADRUPOLE_H
