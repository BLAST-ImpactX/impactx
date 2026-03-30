#ifndef DIAGNOSTICS_TRACK_H_
#define DIAGNOSTICS_TRACK_H_

#include "synergia/bunch/diagnostics.h"

/// Diagnostics_track records the phase space coordinates of a single particle.
/// Particles will only be tracked if they stay on the same processor.
/// Lost particles that are somehow restored or particles not available when
/// the first update is called will also not be tracked.
class Diagnostics_track : public Diagnostics {
  public:
    /// Multiple serial diagnostics can be written to a single file.
    /// The Diagnostics_track class is serial.
    constexpr static const char* diag_type = "diagnostics_track";
    constexpr static const bool diag_write_serial = true;

  private:
    bool found;
    bool first_search;
    bool first_write;

    int index;
    int particle_id;

    double s;
    double s_n;
    int repetition;

    karray1d_row coords;

  private:
    /// Update the diagnostics
    void do_update(Bunch const& bunch) override;
    void do_write(Bunch const& bunch) override;

    std::unique_ptr<Diagnostics>
    do_pilfer() override
    {
        return std::make_unique<Diagnostics_track>(std::move(*this));
    }

  public:
    /// Create an empty Diagnostics_track object
    /// @param bunch_sptr the Bunch
    /// @param filename the base name for file to write to (base names will have
    ///        a numerical index inserted
    /// @param particle_id the particle ID to track
    /// @param local_dir local directory to use for temporary scratch
    Diagnostics_track(int particle_id,
                      std::string const& filename,
                      std::string const& local_dir = "");

    template <class Archive>
    void serialize(Archive& ar, const unsigned int version);
};

#endif /* DIAGNOSTICS_TRACK_H_ */
