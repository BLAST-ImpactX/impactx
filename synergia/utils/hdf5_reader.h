#ifndef HDF5_READER_H
#define HDF5_READER_H

#include "synergia/utils/commxx.h"
#include "synergia/utils/hdf5_misc.h"

class Hdf5_reader {
  public:
    Hdf5_reader() = delete;

    // if T is a scalar type, raed in a single value
    // if T is kokkos::view, read in the array data
    //
    // notes:
    // 1. the array container must be in c ordering
    //
    // 2. the read method does a dimension check. e.g.,
    //    a 3d dataset can be read back only via a
    //    3d container (kokkos::view) as the template
    //    parameter
    //
    // 3. the returned container will have a shape
    //    suitable for holding the entire dataset
    //
    // 4. no collective read allowed for now. the same
    //    data object containing the entire array will
    //    be returned on every rank
    //
    template <class T>
    static T
    read(Hdf5_handler const& file,
         std::string const& name,
         Commxx const& comm,
         int root_rank)
    {
        T retval;
        auto di = syn::extract_data_info(retval);

        // collect dim0
        auto all_dim0 = syn::collect_dims(di.dims, false, comm, root_rank);

        read_verify_and_set_data_dims(
            file, name, all_dim0, di, false, comm, root_rank);

        syn::resize_data_obj(retval, di);

        read_impl(file, name, all_dim0, di, false, comm, root_rank);

        return retval;
    }

    template <class T>
    static T
    read(Hdf5_handler const& file,
         std::string const& name,
         size_t len,
         Commxx const& comm,
         int root_rank)
    {
        T retval;
        auto di = syn::extract_data_info(retval);

        // collective read from an array into a
        // scalar doesnt make much sense
        if (di.dims.size() == 0)
            throw std::runtime_error(
                "hdf5_reader::read() read from an array into a scalar");

        // set the first dim, rests are still 0s
        di.dims[0] = len;

        // collect
        auto all_dim0 = syn::collect_dims(di.dims, true, comm, root_rank);

        read_verify_and_set_data_dims(
            file, name, all_dim0, di, true, comm, root_rank);

        syn::resize_data_obj(retval, di);

        read_impl(file, name, all_dim0, di, true, comm, root_rank);

        return retval;
    }

    // read 1d array
    template <class T>
    static void
    read(Hdf5_handler const& file,
         std::string const& name,
         T* data,
         Commxx const& comm,
         int root_rank)
    {
        throw std::runtime_error("hdf5_reader::read_3() not implemented");
    }

    // collective read a 1d array
    template <class T>
    static void
    read(Hdf5_handler const& file,
         std::string const& name,
         T* data,
         size_t len,
         Commxx const& comm,
         int root_rank)
    {
        syn::data_info_t di{
            data, {len}, hdf5_atomic_data_type<T>(), sizeof(T), len};

        auto all_dim0 = syn::collect_dims(di.dims, true, comm, root_rank);

        read_verify_and_set_data_dims(
            file, name, all_dim0, di, true, comm, root_rank);

        read_impl(file, name, all_dim0, di, true, comm, root_rank);
    }

    static std::vector<hsize_t>
    get_dims(Hdf5_handler const& file,
             std::string const& name,
             Commxx const& comm,
             int root_rank)
    {
        int mpi_size = comm.size();
        int mpi_rank = comm.rank();

        if (!has_dataset(file, name, comm, root_rank))
            throw std::runtime_error("get_dims() dataset does not exist: " +
                                     name);

        int rank = 0;
        std::vector<hsize_t> dims = {};

#ifndef USE_PARALLEL_HDF5
        // get and verify done only on the root_rank
        if (mpi_rank == root_rank)
#endif
        {
            if (!file.valid()) throw std::runtime_error("invalid file handler");

            Hdf5_handler dset = H5Dopen(file, name.c_str(), H5P_DEFAULT);
            Hdf5_handler filespace = H5Dget_space(dset);

            // get and check the data rank
            rank = H5Sget_simple_extent_ndims(filespace);
            if (rank < 0) throw Hdf5_exception("error at getting rank");

            // get, check, and broadcast the dims for non-scalar
            if (rank) {
                dims.resize(rank);
                herr_t res =
                    H5Sget_simple_extent_dims(filespace, dims.data(), NULL);
                if (res < 0) throw Hdf5_exception("error at getting dims");
            }
        }

#ifndef USE_PARALLEL_HDF5
        // boradcast data rank
        MPI_Bcast(&rank, 1, MPI_INT, root_rank, comm);

        // broadcast the dims to all ranks
        if (rank) {
            dims.resize(rank);
            MPI_Bcast(dims.data(), dims.size(), MPI_UINT64_T, root_rank, comm);
        }
#endif

        return dims;
    }

    static void
    read_verify_and_set_data_dims(Hdf5_handler const& file,
                                  std::string const& name,
                                  std::vector<hsize_t>& all_dims_0,
                                  syn::data_info_t& di,
                                  bool collective,
                                  Commxx const& comm,
                                  int root_rank)
    {
        int mpi_size = comm.size();
        int mpi_rank = comm.rank();

        if (!has_dataset(file, name, comm, root_rank))
            throw std::runtime_error("hdf5_read() dataset does not exist: " +
                                     name);

        // data rank
        auto data_rank = di.dims.size();

        // offsets for each rank (offsets of dim0 in the combined array)
        std::vector<hsize_t> offsets(mpi_size, 0);
        for (int r = 0; r < mpi_size - 1; ++r)
            offsets[r + 1] = offsets[r] + all_dims_0[r];

        // dim0 of the combined array
        hsize_t dim0 = offsets[mpi_size - 1] + all_dims_0[mpi_size - 1];

        // get and verify done only on the root_rank
        int error = 0;

#ifndef USE_PARALLEL_HDF5
        if (mpi_rank == root_rank)
#endif
        {
            if (!file.valid()) {
                error = 1;
                MPI_Bcast(&error, 1, MPI_INT, root_rank, comm);
                throw std::runtime_error("invalid file handler");
            }

            Hdf5_handler dset = H5Dopen(file, name.c_str(), H5P_DEFAULT);
            Hdf5_handler filespace = H5Dget_space(dset);

            // get and check the datatype
            Hdf5_handler type = H5Dget_type(dset);
            if (H5Tequal(type, di.atomic_type) <= 0) {
                error = 1;
                MPI_Bcast(&error, 1, MPI_INT, root_rank, comm);
                throw std::runtime_error("wrong atomic data type");
            }

            // get and check the data rank
            int rank = H5Sget_simple_extent_ndims(filespace);
            if (rank < 0) {
                error = 1;
                MPI_Bcast(&error, 1, MPI_INT, root_rank, comm);
                throw Hdf5_exception("error at getting rank");
            }

            if (rank != data_rank) {
                error = 1;
                MPI_Bcast(&error, 1, MPI_INT, root_rank, comm);
                throw std::runtime_error("Hdf5_reader: wrong rank");
            }

            // get, check, and broadcast the dims for non-scalar
            if (rank) {
                std::vector<hsize_t> dims(rank);
                herr_t res =
                    H5Sget_simple_extent_dims(filespace, dims.data(), NULL);
                if (res < 0) {
                    error = 1;
                    MPI_Bcast(&error, 1, MPI_INT, root_rank, comm);
                    throw Hdf5_exception("error at getting dims");
                }

                // only the collective read of non-scalar needs the dim check
                if (collective && (dims[0] != dim0)) {
                    error = 1;
                    MPI_Bcast(&error, 1, MPI_INT, root_rank, comm);
                    throw std::runtime_error(
                        "Hdf5_reader: combined dim[0] is inconsistent with "
                        "the dataset from file");
                }

                di.dims = dims;
            }
        }

#ifndef USE_PARALLEL_HDF5
        MPI_Bcast(&error, 1, MPI_INT, root_rank, comm);
        if (error) throw std::runtime_error("err at root_rank");

        // broadcast the dims to all ranks
        if (data_rank)
            MPI_Bcast(di.dims.data(), data_rank, MPI_UINT64_T, root_rank, comm);
#endif

        // in collective vector read, all_dim0 is init with everyone's read
        // size. after read_verify, di.dims is populated with the full array
        // shape. so the final set is needed to get a proper array shape for
        // each rank
        //
        // collective scalar read is prohibited
        //
        // in non-collective vector read, all_dim0 is inited with all 0, and
        // di.dims are set to the dims of the full array. In this case, set the
        // root rank of all_dims0 to the full dim0, and leave all other ranks
        // to 0. Later in read_impl it doesnt need to read the file again in
        // order to figure out the full array size (so when the size is 0, the
        // read should be skipped)
        //
        // in non-collective scalar read, all_dim0 already init with 1 in the
        // root rank, and 0 in others. so no need to do any more setups.
        //
        if (collective)
            di.dims[0] = all_dims_0[mpi_rank];
        else if (data_rank)
            all_dims_0[root_rank] = di.dims[0];
    }

    static void
    read_impl(Hdf5_handler const& file,
              std::string const& name,
              std::vector<hsize_t> const& all_dims_0,
              syn::data_info_t const& di,
              bool collective,
              Commxx const& comm,
              int root_rank)
    {
        int mpi_size = comm.size();
        int mpi_rank = comm.rank();

        // data rank
        auto data_rank = di.dims.size();

        // offsets for each rank (offsets of dim0 in the combined array)
        std::vector<hsize_t> offsets(mpi_size, 0);
        for (int r = 0; r < mpi_size - 1; ++r)
            offsets[r + 1] = offsets[r] + all_dims_0[r];

        // total dim0
        hsize_t dim0 = offsets[mpi_size - 1] + all_dims_0[mpi_size - 1];

        // do not read if the file data is empty
        if (data_rank && dim0 == 0) return;

#ifdef USE_PARALLEL_HDF5

        if (!file.valid()) throw std::runtime_error("invalid file handler");

        // open dataset
        Hdf5_handler dset = H5Dopen(file, name.c_str(), H5P_DEFAULT);
        Hdf5_handler fspace = H5Dget_space(dset);
        Hdf5_handler mspace = H5Screate_simple(data_rank, di.dims.data(), NULL);

        if (data_rank) {
            // select the slab
            auto offset = std::vector<hsize_t>(data_rank, 0);
            offset[0] = offsets[mpi_rank];

            auto res = H5Sselect_hyperslab(fspace,
                                           H5S_SELECT_SET,
                                           offset.data(),
                                           NULL,
                                           di.dims.data(),
                                           NULL);

            if (res < 0) throw Hdf5_exception("error at select hyperslab");
        }

        // collective write
        Hdf5_handler plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        // read
        auto res = H5Dread(
            dset, di.atomic_type, mspace, fspace, plist_id, (void*)di.ptr);
        if (res < 0) throw Hdf5_exception("error read");
            // if (res < 0) H5Eprint(H5E_DEFAULT, stdout);

#else

        if (collective) {
            // collective read
            // read a portion and point-to-point send to each rank
            if (mpi_rank == root_rank) {
                if (!file.valid())
                    throw std::runtime_error("invalid file handler");

                // open dataset
                Hdf5_handler dset = H5Dopen(file, name.c_str(), H5P_DEFAULT);
                Hdf5_handler fspace = H5Dget_space(dset);

                // local dims
                auto dimsm = di.dims;

                for (int r = 0; r < mpi_size; ++r) {
                    // set the read count
                    if (dimsm.size()) dimsm[0] = all_dims_0[r];

                    // memspace
                    Hdf5_handler mspace =
                        H5Screate_simple(data_rank, dimsm.data(), NULL);

                    if (data_rank) {
                        // select the slab
                        auto offset = std::vector<hsize_t>(data_rank, 0);
                        offset[0] = offsets[r];

                        auto res = H5Sselect_hyperslab(fspace,
                                                       H5S_SELECT_SET,
                                                       offset.data(),
                                                       NULL,
                                                       dimsm.data(),
                                                       NULL);

                        if (res < 0)
                            throw Hdf5_exception("error at select hyperslab");
                    }

                    if (r == mpi_rank) {
                        auto res = H5Dread(dset,
                                           di.atomic_type,
                                           mspace,
                                           fspace,
                                           H5P_DEFAULT,
                                           (void*)di.ptr);
                        if (res < 0) throw Hdf5_exception("error read");
                    } else {
                        size_t sz = di.atomic_data_size;
                        for (auto d : dimsm)
                            sz *= d;

                        if (sz) {
                            std::vector<uint8_t> buf(sz);

                            // read
                            auto res = H5Dread(dset,
                                               di.atomic_type,
                                               mspace,
                                               fspace,
                                               H5P_DEFAULT,
                                               (void*)buf.data());

                            if (res < 0) throw Hdf5_exception("error read");

                            // send to rank r
                            MPI_Send((void*)buf.data(),
                                     buf.size(),
                                     MPI_BYTE,
                                     r,
                                     0,
                                     comm);
                        }
                    }
                }
            } else {
                // local dims(counts)
                auto dimsm = di.dims;
                if (dimsm.size()) dimsm[0] = all_dims_0[mpi_rank];

                size_t sz = di.atomic_data_size;
                for (auto d : dimsm)
                    sz *= d;

                if (sz) {
                    MPI_Status status;
                    MPI_Recv((void*)di.ptr,
                             sz,
                             MPI_BYTE,
                             root_rank,
                             0,
                             comm,
                             &status);
                }
            }
        } else {
            // non-collective read
            // read entire array, and bcast to all ranks
            if (mpi_rank == root_rank) {
                if (!file.valid())
                    throw std::runtime_error("invalid file handler");

                // open dataset
                Hdf5_handler dset = H5Dopen(file, name.c_str(), H5P_DEFAULT);
                Hdf5_handler fspace = H5Dget_space(dset);

                // memspace
                Hdf5_handler mspace =
                    H5Screate_simple(data_rank, di.dims.data(), NULL);

                // read
                auto res = H5Dread(dset,
                                   di.atomic_type,
                                   mspace,
                                   fspace,
                                   H5P_DEFAULT,
                                   (void*)di.ptr);
                if (res < 0) throw Hdf5_exception("error read");
            }

            // buf size
            size_t sz = di.atomic_data_size;
            for (auto d : di.dims)
                sz *= d;

            // brodcast
            if (sz) MPI_Bcast((void*)di.ptr, sz, MPI_BYTE, root_rank, comm);
        }
#endif
    }

    static void
    bcast_vec_str_send(std::vector<std::string> const& vs,
                       Commxx const& comm,
                       int root)
    {
        size_t total = 0;
        std::vector<size_t> lens;

        for (auto const& n : vs) {
            total += n.size();
            lens.push_back(n.size());
        }

        lens.push_back(total);
        size_t l = lens.size();

        size_t off = 0;
        std::vector<char> buf(total);

        for (auto const& n : vs) {
            strcpy(&buf[off], n.c_str());
            off += n.size();
        }

        MPI_Bcast(&l, 1, MPI_UINT64_T, root, comm);
        MPI_Bcast(&lens[0], lens.size(), MPI_UINT64_T, root, comm);
        if (lens.back()) MPI_Bcast(&buf[0], buf.size(), MPI_BYTE, root, comm);
    }

    static void
    bcast_vec_str_recv(std::vector<std::string>& vs,
                       Commxx const& comm,
                       int root)
    {
        size_t l = 0;
        MPI_Bcast(&l, 1, MPI_UINT64_T, root, comm);

        std::vector<size_t> lens(l);
        MPI_Bcast(&lens[0], lens.size(), MPI_UINT64_T, root, comm);

        if (lens.back()) {
            std::vector<char> buf(lens.back());
            MPI_Bcast(&buf[0], buf.size(), MPI_BYTE, root, comm);

            size_t off = 0;
            for (int i = 0; i < lens.size() - 1; ++i) {
                vs.push_back(std::string(&buf[off], lens[i]));
                off += lens[i];
            }
        }
    }

    static std::vector<std::string>
    get_dataset_names(Hdf5_handler const& file,
                      Commxx const& comm,
                      int root_rank)
    {

        std::vector<std::string> names;

        // declaration of cb is reused later, so define it here.
        auto cb =
#if H5_VERSION_GE(1, 12, 0)
            [](hid_t oid, const char* name, const H5O_info2_t* info, void* op) {
                if (info->type == H5O_TYPE_DATASET)
                    ((std::vector<std::string>*)op)->push_back(name);
                return 0;
            };
#else
            [](hid_t oid, const char* name, const H5O_info_t* info, void* op) {
                if (info->type == H5O_TYPE_DATASET)
                    ((std::vector<std::string>*)op)->push_back(name);
                return 0;
            };

#endif // H5_VERSION_GE(1,12,0)

#ifdef USE_PARALLEL_HDF5

#if H5_VERSION_GE(1, 12, 0)
        auto res = H5Ovisit3(file,
                             H5_INDEX_NAME,
                             H5_ITER_NATIVE,
                             cb,
                             (void*)&names,
                             H5O_INFO_BASIC);

#else
        auto res =
            H5Ovisit(file, H5_INDEX_NAME, H5_ITER_NATIVE, cb, (void*)&names);

#endif // H5_VERSION_GE(1,12,0)

#else

        if (comm.rank() == root_rank) {

#if H5_VERSION_GE(1, 12, 0)
            auto res = H5Ovisit3(file,
                                 H5_INDEX_NAME,
                                 H5_ITER_NATIVE,
                                 cb,
                                 (void*)&names,
                                 H5O_INFO_BASIC);

#else

            auto res = H5Ovisit(
                file, H5_INDEX_NAME, H5_ITER_NATIVE, cb, (void*)&names);
#endif // H5_VERSION_GE(1,12,0)

            bcast_vec_str_send(names, comm, root_rank);
        } else {
            bcast_vec_str_recv(names, comm, root_rank);
        }

#endif

        return names;
    }

    static bool
    has_dataset(Hdf5_handler const& file,
                std::string const& name,
                Commxx const& comm,
                int root_rank)
    {
        int has = 0;

#ifndef USE_PARALLEL_HDF5
        if (comm.rank() == root_rank)
#endif
        {
            // if the name exists
            if (H5Oexists_by_name(file, name.c_str(), H5P_DEFAULT) > 0) {
                // and it is a dataset

#if H5_VERSION_GE(1, 12, 0)
                H5O_info2_t info;
                auto res = H5Oget_info_by_name3(
                    file, name.c_str(), &info, H5O_INFO_BASIC, H5P_DEFAULT);
#else
                H5O_info_t info;
                auto res =
                    H5Oget_info_by_name(file, name.c_str(), &info, H5P_DEFAULT);

#endif // H5_VERSION_GE(1,12,0)
                if (res < 0)
                    throw std::runtime_error("error when getting obj info");

                has = (info.type == H5O_TYPE_DATASET);
            }
        }

#ifndef USE_PARALLEL_HDF5
        MPI_Bcast(&has, 1, MPI_INT, root_rank, comm);
#endif

        return has;
    }
};

#endif
