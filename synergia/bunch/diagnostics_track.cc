
#include "diagnostics_track.h"
#include "synergia/bunch/bunch.h"

Diagnostics_track::Diagnostics_track(int particle_id,
                                     std::string const& filename,
                                     std::string const& local_dir)
    : Diagnostics(diag_type, diag_write_serial, filename, local_dir)
    , found(false)
    , first_search(true)
    , first_write(true)
    , index(Bunch::particle_index_null)
    , particle_id(particle_id)
    , s(0.0)
    , s_n(0.0)
    , repetition(0)
    , coords()
{}

void
Diagnostics_track::do_update(Bunch const& bunch)
{
    auto const& ref = bunch.get_reference_particle();

    repetition = ref.get_repetition();
    s = ref.get_s();

    if (first_search) {
        index = bunch.search_particle(particle_id);
        found = (index != Bunch::particle_index_null);
        first_search = false;
    } else if (found) {
        index = bunch.search_particle(particle_id, index);
        found = (index != Bunch::particle_index_null);
    }

    if (found) {
        s = ref.get_s();
        s_n = ref.get_s_n();
        repetition = ref.get_repetition();

        coords = bunch.get_particle(index);
    }
}

void
Diagnostics_track::do_write(Bunch const& bunch)
{
    if (found && get_write_helper(bunch).write_locally()) {
        auto& file = get_write_helper(bunch).get_hdf5_file();

        if (first_write) {
            auto const& ref = bunch.get_reference_particle();

            file.write("charge", ref.get_charge());
            file.write("mass", ref.get_four_momentum().get_mass());
            file.write("pz", ref.get_four_momentum().get_momentum());

            first_write = false;
        }

        file.write_serial("coords", coords);
        file.write_serial("s_n", s_n);
        file.write_serial("repetition", repetition);
        file.write_serial("s", s);

        get_write_helper(bunch).finish_write();
    }
}
