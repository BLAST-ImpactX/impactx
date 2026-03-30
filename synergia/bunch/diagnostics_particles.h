#ifndef DIAGNOSTICS_PARTICLES_H_
#define DIAGNOSTICS_PARTICLES_H_

#include "synergia/bunch/bunch_particles.h"
#include "synergia/bunch/diagnostics.h"
#include <functional>
#include <optional>

/// Diagnostics_particles dumps the state of particles in a bunch
class Diagnostics_particles : public Diagnostics {

  private:
    // Bunch const& bunch_ptr;
    std::optional<std::reference_wrapper<const Bunch>> bunch_ref;

    int num_part, offset;
    int num_spec_part, spec_offset;
    size_t local_num, local_offset, file_offset;
    size_t spec_local_num, spec_local_offset, spec_file_offset;

#ifdef SYNERGIA_HAVE_OPENPMD
    bunch_particles_t<double>::host_parts_t parts_subset;
    bunch_particles_t<double>::host_masks_t masks_subset;

    bunch_particles_t<double>::host_parts_t spec_parts_subset;
    bunch_particles_t<double>::host_masks_t spec_masks_subset;
#else
#endif

  public:
    Diagnostics_particles(std::string const& filename = "diag_particles.h5",
                          int num_part = -1,
                          int offset = 0,
                          int num_spec_part = 0,
                          int spec_offset = 0);

  private:
    void
    do_reduce(Commxx const& comm, int root) override
    {}
    void do_first_write(io_device& file) override;
    void
    do_update(Bunch const& bunch) override
    {
        bunch_ref.emplace(std::cref<Bunch>(bunch));
    }
    void do_write(io_device& file, const size_t iteration) override;

    friend class cereal::access;

    template <class AR>
    void
    serialize(AR& ar)
    {
        ar(cereal::base_class<Diagnostics>(this));
        ar(num_part);
        ar(offset);
        ar(num_spec_part);
        ar(spec_offset);
    }
};

CEREAL_REGISTER_TYPE(Diagnostics_particles)

#endif /* DIAGNOSTICS_PARTICLES_H_ */
