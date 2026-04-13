#ifndef FF_SOLENOID_H
#define FF_SOLENOID_H

#include "synergia/libFF/ff_drift.h"
#include "synergia/libFF/ff_algorithm.h"
#include "synergia/libFF/ff_patterned_propagator.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/foundation/physical_constants.h"

namespace solenoid_impl
{
    template<class T>
    using kf_t = void(*)(T const&, T&, T const&, T&, double);

    template<class BunchT, kf_t<typename BunchT::gsv_t> KF>
    struct solenoid_edge_kicker
    {
        using gsv_t = typename BunchT::gsv_t;

        typename BunchT::parts_t p;
        typename BunchT::const_masks_t m;

        double kse;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int idx) const
        { 
            int i = idx * gsv_t::size();
            int mask = 0;
            for(int x=i; x<i+gsv_t::size(); ++x) mask |= m(x);

            if (mask) 
            {
                gsv_t p0(&p(i, 0));
                gsv_t p1(&p(i, 1));
                gsv_t p2(&p(i, 2));
                gsv_t p3(&p(i, 3));

                KF(p0, p1, p2, p3, kse); 

                p1.store(&p(i, 1));
                p3.store(&p(i, 3));
            }
        }
    };

    template<class BunchT, kf_t<typename BunchT::gsv_t> KF>
    void apply_edge_kick(BunchT& bunch, ParticleGroup pg, double kse)
    {
        if(!bunch.get_local_num(pg)) return;

        auto parts = bunch.get_local_particles(pg);
        auto masks = bunch.get_local_particle_masks(pg);

        using exec = typename BunchT::exec_space;
        auto range = Kokkos::RangePolicy<exec>(0, bunch.size_in_gsv(pg));

        solenoid_edge_kicker<BunchT, KF> sk{parts, masks, kse};
        Kokkos::parallel_for(range, sk);
    }

    template<class BunchT>
    struct solenoid_unit_kicker
    {
        using gsv_t = typename BunchT::gsv_t;

        typename BunchT::parts_t p;
        typename BunchT::const_masks_t m;

        double ksl;
        double ks;
        double length;
        double ref_p;
        double mass;
        double ref_cdt;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int idx) const
        { 
            int i = idx * gsv_t::size();
            int mask = 0;
            for(int x=i; x<i+gsv_t::size(); ++x) mask |= m(x);

            if (mask) 
            {
                gsv_t p0(&p(i, 0));
                gsv_t p1(&p(i, 1));
                gsv_t p2(&p(i, 2));
                gsv_t p3(&p(i, 3));
                gsv_t p4(&p(i, 4));
                gsv_t p5(&p(i, 5));

                FF_algorithm::solenoid_unit(
                        p0, p1, p2, p3, p4, p5,
                        ksl, ks, length, ref_p, mass, ref_cdt); 

                p0.store(&p(i, 0));
                p1.store(&p(i, 1));
                p2.store(&p(i, 2));
                p3.store(&p(i, 3));
                p4.store(&p(i, 4));
            }
        }
    };

    template<class BunchT>
    void apply_solenoid_unit(BunchT& bunch, ParticleGroup pg,
            double ksl, double ks, double length, double ref_p,
            double mass, double ref_cdt)
    {
        if(!bunch.get_local_num(pg)) return;

        auto parts = bunch.get_local_particles(pg);
        auto masks = bunch.get_local_particle_masks(pg);

        using exec = typename BunchT::exec_space;
        auto range = Kokkos::RangePolicy<exec>(0, bunch.size_in_gsv(pg));

        solenoid_unit_kicker<BunchT> sk{
            parts, masks, ksl, ks, length, 
            ref_p, mass, ref_cdt
        };

        Kokkos::parallel_for(range, sk);
    }


    inline double get_reference_cdt_solenoid(
            double length, Reference_particle & reference_particle,
            bool in_edge, bool out_edge, double ks, double kse, double ksl)
    {
        double x(reference_particle.get_state()[Bunch::x]);
        double xp(reference_particle.get_state()[Bunch::xp]);
        double y(reference_particle.get_state()[Bunch::y]);
        double yp(reference_particle.get_state()[Bunch::yp]);
        double cdt(0.0);
        // For calculating cdt, the dp/p of the reference particle should be 0
        // double dpop(reference_particle.get_state()[Bunch::dpop]);
        double dpop=0.0;

        double ref_p = reference_particle.get_momentum()*(1+dpop);
        double mass  = reference_particle.get_mass();

        // in edge
        if (in_edge)
        {
            FF_algorithm::solenoid_in_edge_kick(x, xp, y, yp, kse);
        }

        // body
        FF_algorithm::solenoid_unit( x, xp, y, yp, cdt, dpop,
                ksl, ks, length, ref_p, mass, 0.0);

        // out edge
        if (out_edge)
        {
            FF_algorithm::solenoid_out_edge_kick(x, xp, y, yp, kse);
        }

        // propagate and update the bunch design reference particle state
        reference_particle.set_state(x, xp, y, yp, cdt, dpop);

        return cdt;
    }
}


namespace FF_solenoid
{
    template<class BunchT>
    void apply(Lattice_element_slice const& slice, BunchT& bunch)
    {
        using namespace solenoid_impl;

        const double  length = slice.get_right() - slice.get_left();
        const double    mass = bunch.get_mass();

        // in edge, out edge
        const bool has_in_edge  = slice.has_left_edge();
        const bool has_out_edge = slice.has_right_edge();

        // ks here is the field
        double ks = slice.get_lattice_element().get_double_attribute("ks");

        Reference_particle       & ref_l = bunch.get_design_reference_particle();
        Reference_particle const & ref_b = bunch.get_reference_particle();
        const double ref_p = ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]);

        // if ks=0 it is effectively a drift
        if (fabs(ks) < 1e-12)
        {
            FF_drift::apply(slice, bunch);
            return;
        }

        // scaling
        double brho_l = ref_l.get_momentum() 
            * (1.0 + ref_l.get_state()[Bunch::dpop]) / ref_l.get_charge();  // GV/c

        double brho_b = ref_b.get_momentum() 
            * (1.0 + ref_b.get_state()[Bunch::dpop]) / ref_l.get_charge();  // GV/c

        double scale = brho_l / brho_b;

                ks = ks * scale;
        double ksl = ks * length;
        double kse = ks * 0.5;

        // reference cdt
        const double ref_cdt = get_reference_cdt_solenoid(length, ref_l, 
                has_in_edge, has_out_edge, ks, kse, ksl);

        using gsv_t = typename BunchT::gsv_t;

        // in-edge
        if (has_in_edge)
        {
            apply_edge_kick<BunchT, 
                FF_algorithm::solenoid_in_edge_kick<gsv_t>>(
                    bunch, ParticleGroup::regular, kse);

            apply_edge_kick<BunchT, 
                FF_algorithm::solenoid_in_edge_kick<gsv_t>>(
                    bunch, ParticleGroup::spectator, kse);
        }

        // body
        apply_solenoid_unit(bunch, ParticleGroup::regular,
                ksl, ks, length, ref_p, mass, ref_cdt);

        apply_solenoid_unit(bunch, ParticleGroup::spectator,
                ksl, ks, length, ref_p, mass, ref_cdt);

        // out-edge
        if (has_out_edge)
        {
            apply_edge_kick<BunchT,
                FF_algorithm::solenoid_out_edge_kick<gsv_t>>(
                    bunch, ParticleGroup::regular, kse);

            apply_edge_kick<BunchT,
                FF_algorithm::solenoid_out_edge_kick<gsv_t>>(
                    bunch, ParticleGroup::spectator, kse);
        }

        // trajectory
        bunch.get_reference_particle().increment_trajectory(length);
        // absolute time
        double const velocity = bunch.get_reference_particle().get_beta()*pconstants::c;
        bunch.get_reference_particle().increment_bunch_abs_time(length/velocity);

    }
}


#endif // FF_SOLENOID_H
