#ifndef FF_REF_COORDS_H
#define FF_REF_COORDS_H

#include "synergia/foundation/physical_constants.h"
#include "synergia/libFF/ff_algorithm.h"
#include "synergia/libFF/ff_patterned_propagator.h"
#include "synergia/utils/simple_timer.h"

namespace ref_energy_impl {
    struct CoordParams {
        double mass;   // mass
        double old_E; // old energy
        double new_E; // new energy
        double old_p; // old momentum
        double new_p; // new momentum
    };

    template <class BunchT>
    struct AdjustCoords {
        using gsv_t = typename BunchT::gsv_t;

        typename BunchT::bp_t::parts_t p;
        typename BunchT::bp_t::const_masks_t m;
        const CoordParams cp;

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

                FF_algorithm::adjust_ref_unit(p1,
                                                 p3,
                                                 p4,
                                                 p5,
                                                 cp.mass,
                                                 cp.old_p,
                                                 cp.new_p,
                                                 cp.old_E,
                                                 cp.new_E);

                p1.store(&p(i, 1));
                p3.store(&p(i, 3));
                p4.store(&p(i, 4));
                p5.store(&p(i, 5));
            }
        }
    };

    template <class BunchT>
    struct PropAdjust {
        using gsv_t = typename BunchT::gsv_t;

        typename BunchT::bp_t::parts_t p;
        typename BunchT::bp_t::const_masks_t m;
        const CoordParams cp;

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

                FF_algorithm::adjust_ref_unit(p1,
                                                 p3,
                                                 p4,
                                                 p5,
                                                 cp.mass,
                                                 cp.old_p,
                                                 cp.new_p,
                                                 cp.old_E,
                                                 cp.new_E);

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

namespace FF_adjust_bunch_ref_coords {
    template <class BunchT>
    void
    apply(double new_E, BunchT& bunch)
    {
        using namespace ref_energy_impl;

        scoped_simple_timer timer("libFF_ref_energy");

        CoordParams cp;

        // actually we don't care about the slice at all. Maybe
        // we can eliminate the slice argument
 
        Reference_particle& ref_l = bunch.get_design_reference_particle();
        Reference_particle& ref_b = bunch.get_reference_particle();

        // The bunch particles momentum is with respect to the bunch reference
        // particle
        cp.old_p = ref_l.get_momentum();
        cp.old_E = ref_l.get_total_energy();
        cp.mass = bunch.get_mass();

        cp.new_E = new_E;
        cp.new_p = sqrt(new_E*new_E - cp.mass*cp.mass);

        // bunch particles
        auto apply = [&](ParticleGroup pg) {
            if (!bunch.get_local_num(pg)) return;

            auto parts = bunch.get_local_particles(pg);
            auto masks = bunch.get_local_particle_masks(pg);

            using exec = typename BunchT::exec_space;
            auto range = Kokkos::RangePolicy<exec>(0, bunch.size_in_gsv(pg));

            AdjustCoords<BunchT> adjuster{parts, masks, cp};
            Kokkos::parallel_for(range, adjuster);
        };
        apply(ParticleGroup::regular);
        apply(ParticleGroup::spectator);

        // now set the bunch energy to match
        ref_b.set_total_energy(cp.new_E);
    }
}

#endif // FF_REF_COORDS_H
