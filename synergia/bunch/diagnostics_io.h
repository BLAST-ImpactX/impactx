#ifndef DIAGNOSTICS_IO_H
#define DIAGNOSTICS_IO_H

#include <cereal/types/memory.hpp>
#include <memory>
#include <optional>

#include "synergia/utils/commxx.h"
#include "synergia/utils/kokkos_views.h"
#include "synergia/utils/synergia_config.h"

class Diagnostics_io {
  private:
    std::optional<io_device> file;

    const std::shared_ptr<Commxx> comm;

    bool single_file;
    size_t iteration_count;

    std::string filename;

    std::string filename_base;
    std::string filename_suffix;

    static const int flush_period = 100;

  public:
    Diagnostics_io(std::string filename,
                   bool single_file,
                   std::shared_ptr<Commxx> const& comm_input);

    ~Diagnostics_io();

    Diagnostics_io(Diagnostics_io const&) = delete;
    Diagnostics_io(Diagnostics_io&&) = default;

    io_device& get_io_device();

    int get_root_rank();

    // count for non-single-file writer
    int
    get_count() const
    {
        return iteration_count;
    }
    void
    set_count(int count_input)
    {
        iteration_count = count_input;
    }
    void
    increment_count()
    {
        ++iteration_count;
    }

    void finish_write();

    void open_file();
    void flush_file();
    void close_file();

  private:
    std::string get_filename();

    void move_file_overwrite_if_exists(std::string const& src,
                                       std::string const& dst);

    friend class Diagnostics_worker;
    friend class cereal::access;

    Diagnostics_io();

    template <class AR>
    void
    serialize(AR& ar)
    {
#ifdef SYNERGIA_HAVE_OPENPMD
#else
        try {
            ar(CEREAL_NVP(file.value()));
        }
        catch (const std::bad_optional_access& e) {
            // consistent with old behavior! previously we let empty unique_ptrs
            // be serialized
        }

#endif
        ar(CEREAL_NVP(single_file));
        ar(CEREAL_NVP(iteration_count));
        ar(CEREAL_NVP(filename));
        ar(CEREAL_NVP(filename_base));
        ar(CEREAL_NVP(filename_suffix));
    }
};

#endif
