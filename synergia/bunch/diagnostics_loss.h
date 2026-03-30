#ifndef DIAGNOSTICS_LOSS_H_
#define DIAGNOSTICS_LOSS_H_

#include "synergia/bunch/diagnostics.h"
#include "synergia/utils/kokkos_views.h"

class Diagnostics_loss : public Diagnostics {
  private:
    int bucket_index;
    int repetition;
    double s_ref_particle;
    double sn_ref_particle;

    karray2d_row coords;

#ifdef SYNERGIA_HAVE_OPENPMD
    size_t discards_local_num;
    size_t discards_offset_num;
    size_t discards_total_num;
#endif

  public:
    Diagnostics_loss(std::string const& filename = "diag_loss.h5")
        : Diagnostics("diagnostics_loss", filename, true)
        , bucket_index(-1)
        , repetition(0)
        , s_ref_particle(0.0)
        , sn_ref_particle(0.0)
        , coords()
    {}

  private:
    void do_update(Bunch const& bunch) override;
    void
    do_reduce(Commxx const& comm, int root) override
    {}
    void do_first_write(io_device& file) override;
    void do_write(io_device& file, size_t iteration) override;

    friend class cereal::access;

    template <class Archive>
    void
    serialize(Archive& ar)
    {
        ar(cereal::base_class<Diagnostics>(this));
    }
};

CEREAL_REGISTER_TYPE(Diagnostics_loss)

#endif /* LOSS_DIAGNOSTICS_H_ */
