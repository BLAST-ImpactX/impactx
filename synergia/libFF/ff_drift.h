#ifndef FF_DRIFT_H
#define FF_DRIFT_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/foundation/physical_constants.h"


namespace drift_impl
{
    template<class BP>
    struct PropDrift
    {
        typename BP::parts_t p;
        typename BP::const_masks_t masks;

        double len, ref_p, mass, ref_t;

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            if (masks(i))
                FF_algorithm::drift_unit(
                        p(i, 0), p(i, 1), p(i, 2),
                        p(i, 3), p(i, 4), p(i, 5),
                        len, ref_p, mass, ref_t);
        }
    };

    template<class BP>
    struct PropDriftSimd
    {
        using gsv_t = typename BP::gsv_t;

        typename BP::parts_t p;
        typename BP::const_masks_t masks;

        double len, ref_p, mass, ref_t;

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

                FF_algorithm::drift_unit(
                        p0, p1, p2, p3, p4, p5,
                        len, ref_p, mass, ref_t);

                p0.store(&p(i, 0));
                p2.store(&p(i, 2));
                p4.store(&p(i, 4));
            }
        }
    };


    inline double get_reference_cdt(double length, Reference_particle & ref)
    {
        double x(ref.get_state()[Bunch::x]);
        double xp(ref.get_state()[Bunch::xp]);
        double y(ref.get_state()[Bunch::y]);
        double yp(ref.get_state()[Bunch::yp]);
        double cdt(0.0);
        // EGS: don't use refpart state dpop for calculation of cdt
        // double dpop(ref.get_state()[Bunch::dpop]);
        double dpop(0.0);
        double ref_p = ref.get_momentum();
        double m = ref.get_mass();

        FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, ref_p, m, 0.0);

        // propagate and update the bunch design reference particle state
        ref.set_state(x, xp, y, yp, cdt, dpop);

        return cdt;
    }
}


namespace FF_drift
{
    template<class BunchT>
    inline void apply(Lattice_element_slice const& slice, BunchT & bunch)
    {
        using namespace drift_impl;

        scoped_simple_timer timer("libFF_drift");

        const double  length = slice.get_right() - slice.get_left();
        const double    mass = bunch.get_mass();

        // zero-length drift does nothing
        if (close_to_zero(length)) 
        {
            bunch.get_design_reference_particle().set_state_cdt(0.0);
            return;
        }

        Reference_particle       & ref_l = bunch.get_design_reference_particle();
        Reference_particle const & ref_b = bunch.get_reference_particle();

        const double ref_p   = ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]);
        const double ref_cdt = get_reference_cdt(length, ref_l);

        // apply method
        auto apply_impl = [&](ParticleGroup pg) {
            auto bp = bunch.get_bunch_particles(pg);
            if (!bp.num_valid()) return;

            using namespace Kokkos;
            using bp_t = typename BunchT::bp_t;
            using exec = typename BunchT::exec_space;

#if LIBFF_USE_GSV
            auto range = RangePolicy<exec>(0, bp.size_in_gsv());
            PropDriftSimd<bp_t> drift{
                bp.parts, bp.masks, length, ref_p, mass, ref_cdt};
            parallel_for(range, drift);
#else
            auto range = RangePolicy<exec>(0, bp.size());
            PropDrift<bp_t> drift{
                bp.parts, bp.masks, length, ref_p, mass, ref_cdt};
            parallel_for(range, drift);
#endif
        };

        // apply on bunch
        apply_impl(ParticleGroup::regular);
        apply_impl(ParticleGroup::spectator);

        // trajectory

        bunch.get_reference_particle().increment_trajectory(length);

        // absolute time
        double const velocity = bunch.get_reference_particle().get_beta()*pconstants::c;
        bunch.get_reference_particle().increment_bunch_abs_time(length/velocity);

        Kokkos::fence();
    }
}

#endif // FF_DRIFT_H
