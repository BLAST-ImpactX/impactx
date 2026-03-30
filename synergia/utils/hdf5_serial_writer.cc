#include "hdf5_serial_writer.h"

void
Hdf5_serial_writer::do_setup(std::vector<hsize_t> const& data_dims)
{
    std::vector<hsize_t> chunk_dims(data_rank + 1);

    dims.resize(data_rank + 1);
    max_dims.resize(data_rank + 1);
    size.resize(data_rank + 1);
    offset.resize(data_rank + 1);

    for (int i = 0; i < data_rank; ++i) {
        dims[i] = data_dims.at(i);
        max_dims[i] = data_dims.at(i);
        size[i] = data_dims.at(i);
        chunk_dims[i] = data_dims.at(i);
        offset[i] = 0;
    }

    // set the last dim (the extensible dim)
    dims[data_rank] = 1;
    max_dims[data_rank] = H5S_UNLIMITED;

    const size_t good_chunk_size = 8192; // pulled out of air
    chunk_dims[data_rank] =
        (data_size < good_chunk_size) ? good_chunk_size / data_size : 1;

    try {
        // try to open the dataset
        dataset = H5Dopen(file_ptr, name.c_str(), H5P_DEFAULT);
        Hdf5_handler dataspace = H5Dget_space(dataset);

        // check the data_rank
        if (H5Sget_simple_extent_ndims(dataspace) != data_rank + 1) {
            throw std::runtime_error(
                "Hdf5_serial_writer::resumed data has wrong rank");
        }

        std::vector<hsize_t> fdims(data_rank + 1);
        herr_t res = H5Sget_simple_extent_dims(dataspace, &fdims[0], NULL);
        if (res < 0)
            throw Hdf5_exception("Error getting dims of the dataspace");

        // check the dims of each rank
        for (int i = 0; i < data_rank; ++i) {
            if (dims[i] != fdims[i])
                throw std::runtime_error("inconsistent data dimensions");
        }

        size[data_rank] = fdims[data_rank];
        offset[data_rank] = fdims[data_rank];
    }
    catch (Hdf5_exception& e) {
        // dataset not exist, create a new one
        size[data_rank] = 0;
        offset[data_rank] = 0;

        if (data_size == 0) {
            throw std::runtime_error(
                "Hdf5_serial_writer: zero data size encountered");
        }

        Hdf5_handler cparms = H5Pcreate(H5P_DATASET_CREATE);
        herr_t res = H5Pset_chunk(cparms, data_rank + 1, &chunk_dims[0]);
        if (res < 0) throw Hdf5_exception();

        Hdf5_handler dataspace =
            H5Screate_simple(data_rank + 1, &dims[0], &max_dims[0]);
        dataset = H5Dcreate(file_ptr,
                            name.c_str(),
                            atomic_type,
                            dataspace,
                            H5P_DEFAULT,
                            cparms,
                            H5P_DEFAULT);
    }

    have_setup = true;
}

void
Hdf5_serial_writer::do_append(void* ptr)
{
    // create dataspace for current data block (it looks like max_dims can be
    // null)
    Hdf5_handler dataspace =
        H5Screate_simple(data_rank + 1, &dims[0], &max_dims[0]);

    // increment and extend the dataset to the new size (last_dim+1)
    ++size[data_rank];
    herr_t res = H5Dextend(dataset, &size[0]);
    if (res < 0) throw Hdf5_exception();

    // select the slab to write
    Hdf5_handler filespace = H5Dget_space(dataset);
    res = H5Sselect_hyperslab(
        filespace, H5S_SELECT_SET, &offset[0], NULL, &dims[0], NULL);
    if (res < 0) throw Hdf5_exception();

    // write
    res =
        H5Dwrite(dataset, atomic_type, dataspace, filespace, H5P_DEFAULT, ptr);
    if (res < 0) throw Hdf5_exception();

    // increment the offset
    ++offset[data_rank];
}
