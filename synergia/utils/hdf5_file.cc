#include <thread>

#include "hdf5_file.h"
#include "hdf5_misc.h"

Hdf5_file::Hdf5_file(std::string const& file_name, Flag flag, Commxx const& c)
    : comm(std::make_shared<Commxx>(c))
    , file_name(file_name)
    , h5file()
    , root_rank(c.size() - 1)
    , is_open(false)
    , current_flag(flag)
#ifdef USE_PARALLEL_HDF5
    , has_file(true)
#else
    , has_file(c.rank() == root_rank)
#endif
{
    // turn off the automatic error printing
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    // open file
    open(flag);
}

Hdf5_file::Hdf5_file(std::string const& file_name,
                     Flag flag,
                     std::shared_ptr<Commxx> const& c)
    : comm(c)
    , file_name(file_name)
    , h5file()
    , root_rank(c->size() - 1)
    , is_open(false)
    , current_flag(flag)
#ifdef USE_PARALLEL_HDF5
    , has_file(true)
#else
    , has_file(c->rank() == root_rank)
#endif
{
    // turn off the automatic error printing
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    // open file
    open(flag);
}

void
Hdf5_file::open(Flag flag)
{
    // open file on this rank?
    if (!has_file) return;

    // already opened
    if (is_open) close();

    int attempts = 0;
    bool fail = true;

    while ((attempts < 5) && fail) {
        try {
            Hdf5_handler plist_id = H5Pcreate(H5P_FILE_ACCESS);

#ifdef USE_PARALLEL_HDF5
            MPI_Info info = MPI_INFO_NULL;
            H5Pset_fapl_mpio(plist_id, *comm, info);
#endif

            if (flag == Hdf5_file::Flag::truncate) {
                // create
                h5file = H5Fcreate(
                    file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
                current_flag = Hdf5_file::Flag::read_write;
            } else {
                // open
                h5file = H5Fopen(
                    file_name.c_str(), flag_to_h5_flags(flag), plist_id);
            }

            fail = false;
        }
        catch (Hdf5_exception& e) {
            ++attempts;
            fail = true;

            std::cout << e.what() << "\n";
            std::cout << "caught hdf5 open file error, attempts number="
                      << attempts << " on rank=" << Commxx::world_rank()
                      << std::endl;

            std::this_thread::sleep_for(std::chrono::seconds(3));
        }
    }

    is_open = true;
}

void
Hdf5_file::close()
{
    if (is_open) {
        // delete all seq writers before closing the file,
        // otherwise reopening the file would cause error
        seq_writers.clear();
        h5file.close();
        is_open = false;
    }
}

void
Hdf5_file::flush() const
{
    if (is_open) H5Fflush(h5file.hid, H5F_SCOPE_GLOBAL);
}
