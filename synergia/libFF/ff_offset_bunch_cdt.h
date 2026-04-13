#ifndef FF_OFFSET_CDT_H
#define FF_OFFSET_CDT_H

#include "synergia/foundation/physical_constants.h"
#include "synergia/libFF/ff_algorithm.h"
#include "synergia/libFF/ff_patterned_propagator.h"
#include "synergia/utils/simple_timer.h"

namespace offset_cdt_impl {
    struct CoordParams {
        double offset; // offset for cdt
    };

    template <class BunchT>
    struct OffsetCoords {
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
                gsv_t p4(&p(i, 4));

                FF_algorithm::adjust_cdt_unit(p4, cp.offset);

                p4.store(&p(i, 4));
            }
        }
    };

    template <class BunchT>
    struct PropOffset {
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
                gsv_t p4(&p(i, 4));

                FF_algorithm::adjust_cdt_unit(p4, cp.offset);

                p4.store(&p(i, 4));
            }
        }
    };

}

namespace FF_offset_cdt {
    template <class BunchT>
    void
    apply(double offset, BunchT& bunch)
    {
        using namespace offset_cdt_impl;

        scoped_simple_timer timer("libFF_offset_cdt");

        CoordParams cp;

	cp.offset = offset;
        // actually we don't care about the slice at all. Maybe
        // we can eliminate the slice argument
 
        Reference_particle& ref_l = bunch.get_design_reference_particle();
        Reference_particle& ref_b = bunch.get_reference_particle();

        // This operation does not affect any reference particle coordinates

        // bunch particles
        auto apply = [&](ParticleGroup pg) {
            if (!bunch.get_local_num(pg)) return;

            auto parts = bunch.get_local_particles(pg);
            auto masks = bunch.get_local_particle_masks(pg);

            using exec = typename BunchT::exec_space;
            auto range = Kokkos::RangePolicy<exec>(0, bunch.size_in_gsv(pg));

            PropOffset<BunchT> adjuster{parts, masks, cp};
            Kokkos::parallel_for(range, adjuster);
        };
        apply(ParticleGroup::regular);
        apply(ParticleGroup::spectator);

    }
}

#endif // FF_OFFSET_CDT_H
