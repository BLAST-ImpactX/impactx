#ifndef DIAGNOSTICS_FULL2_H_
#define DIAGNOSTICS_FULL2_H_

#include "synergia/bunch/diagnostics.h"
#include "synergia/foundation/reference_particle.h"
#include "synergia/utils/kokkos_views.h"

/// Diagnostics_full2 provides the full set of statistical
/// quantities to be calculated for a Bunch up to the second moments.
class Diagnostics_full2 : public Diagnostics {
  private:
    Reference_particle ref;

    int num_particles;
    double real_num_particles;

    karray1d mean;
    karray1d std;
    karray1d min;
    karray1d max;
    karray2d_row mom2;
    karray2d_row corr;

    double emitx, emity, emitz, emitxy, emitxyz;

  public:
    Diagnostics_full2(std::string const& filename = "diag_full2.h5")
        : Diagnostics("diagnostics_full2", filename, true)
        , ref()
        , num_particles(0)
        , real_num_particles(0.0)
        , mean("mean", 6)
        , std("std", 6)
        , min("min", 3)
        , max("max", 3)
        , mom2("mom2", 6, 6)
        , corr("corr", 6, 6)
    {}

  private:
    void do_update(Bunch const& bunch) override;
    void
    do_reduce(Commxx const& comm, int root) override
    {}
    void do_first_write(io_device& file) override;
    void do_write(io_device& file, size_t iteration) override;

    friend class cereal::access;

    template <class AR>
    void
    serialize(AR& ar)
    {
        ar(cereal::base_class<Diagnostics>(this));
    }
};

CEREAL_REGISTER_TYPE(Diagnostics_full2)

#endif /* DIAGNOSTICS_FULL2_H_ */
