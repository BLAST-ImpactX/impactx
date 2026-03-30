#include "synergia/bunch/diagnostics_io.h"

#include <iomanip>
#include <string>
#include <synergia_version.h>

Diagnostics_io::Diagnostics_io(std::string filename,
                               bool single_file,
                               std::shared_ptr<Commxx> const& comm_input)
    : file()
    , comm(comm_input)
    , single_file(single_file)
    , iteration_count(0)
    , filename(filename)
    , filename_base()
    , filename_suffix()
{

#ifdef SYNERGIA_HAVE_OPENPMD
    if (!single_file) {
        auto idx = filename.rfind('.');
        if (idx == std::string::npos) {
            filename_suffix = "";
        } else {
            filename_base = filename.substr(0, idx);
            filename_suffix = filename.substr(idx);
        }
        filename = "";
        filename.append(filename_base).append("_04%T").append(filename_suffix);
    }
#else
    auto idx = filename.rfind('.');
    if (idx == std::string::npos) {
        filename_base = filename;
        filename_suffix = "";
    } else {
        filename_base = filename.substr(0, idx);
        filename_suffix = filename.substr(idx);
    }
#endif

#ifdef SYNERGIA_HAVE_OPENPMD
    file.emplace(filename, openPMD::Access::CREATE, MPI_Comm(*comm));

    std::string version =
        std::to_string(synergia_version::synergia_version_year);
    version.append(".");
    version.append(std::to_string(synergia_version::synergia_version_month));
    version.append(".");
    version.append(std::to_string(synergia_version::synergia_version_day));
    version.append("-");
    version.append(synergia_version::synergia_git_hash);

    file->setSoftware("synergia3", version);

#else
    file.emplace(get_filename(), Hdf5_file::Flag::truncate, comm);
#endif
}

Diagnostics_io::Diagnostics_io()
    : file()
    , single_file(true)
    , iteration_count(0)
    , filename()
    , filename_base()
    , filename_suffix()
{}

Diagnostics_io::~Diagnostics_io()
{
    if (file.has_value()) { close_file(); }
}

io_device&
Diagnostics_io::get_io_device()
{
#ifdef SYNERGIA_HAVE_OPENPMD
#else
    open_file();
#endif
    return file.value();
}

std::string
Diagnostics_io::get_filename()
{
    std::stringstream sstream;

    sstream << filename_base;

    if (!single_file) {
        sstream << "_";
        sstream << std::setw(5);
        sstream << std::setfill('0');
        sstream << iteration_count;
    }

    sstream << filename_suffix;

    return sstream.str();
}

void
Diagnostics_io::open_file()
{
    if (!file.has_value()) {
#ifdef SYNERGIA_HAVE_OPENPMD
        file.emplace(filename, openPMD::Access::CREATE, MPI_Comm(*comm));
#else
        file.emplace(get_filename(), Hdf5_file::Flag::truncate, comm);
#endif
    }
    return;
}

void
Diagnostics_io::flush_file()
{
    if (file.has_value() && (iteration_count % flush_period == 0))
        file.value().flush();
}

void
Diagnostics_io::close_file()
{
    if (file.has_value()) {
#ifdef SYNERGIA_HAVE_OPENPMD
#else
        file.value().close();
#endif
        file.reset();
    } else {
        throw std::runtime_error(
            "Diagnostics_io: Close called with no active file!");
    }
}

void
Diagnostics_io::finish_write()
{
    if (single_file)
        flush_file();
    else
        close_file();

    increment_count();
}

int
Diagnostics_io::get_root_rank()
{

#ifdef SYNERGIA_HAVE_OPENPMD
    return -1;
#else
    if (file.has_value()) {
        return file.value().get_root_rank();
    } else {
        throw std::runtime_error("Diagnostics_io: Attempting to get root rank "
                                 "when no HDF5 file is present");
    }
#endif
}
