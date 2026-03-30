#ifndef HDF5_SERIAL_WRITER_H_
#define HDF5_SERIAL_WRITER_H_

#include <string>
#include <vector>

#include <hdf5.h>

#include "synergia/utils/cereal.h"
#include "synergia/utils/hdf5_misc.h"
#include "synergia/utils/kokkos_views.h"

class Hdf5_serial_writer {
  private:
    int data_rank;
    std::string name;
    std::vector<hsize_t> dims, max_dims, size, offset;

    hid_t file_ptr;
    Hdf5_handler dataset;
    Hdf5_handler atomic_type;

    bool have_setup;
    size_t data_size;

  private:
    void do_setup(std::vector<hsize_t> const& data_dims);
    void do_append(void* ptr);

    // default for karray types
    template <typename T>
    std::vector<hsize_t>
    pre_setup(T const& data)
    {
        data_rank = T::Rank;
        atomic_type = hdf5_atomic_data_type<double>();
        data_size = sizeof(double) * data.size();

        std::vector<hsize_t> data_dims(data_rank);
        for (int i = 0; i < data_rank; ++i)
            data_dims[i] = data.extent(i);
        return data_dims;
    }

    // default for karray types
    template <typename T>
    void*
    get_data_ptr(T const& data)
    {
        return (void*)data.data();
    }

  public:
    Hdf5_serial_writer(hid_t file_ptr, std::string const& name)
        : data_rank(0)
        , name(name)
        , dims()
        , max_dims()
        , size()
        , offset()
        , file_ptr(file_ptr)
        , dataset()
        , atomic_type()
        , have_setup(false)
        , data_size(0)
    {}

    template <typename T>
    void
    append(T const& data)
    {
        if (!have_setup) do_setup(pre_setup(data));
        do_append(get_data_ptr(data));
    }
};

template <>
inline std::vector<hsize_t>
Hdf5_serial_writer::pre_setup<int>(int const& data)
{
    data_rank = 0;
    atomic_type = hdf5_atomic_data_type<int>();
    data_size = sizeof(int);

    // dummy variable -- length really should
    // be 0, but that would not compile
    return std::vector<hsize_t>(1);
}

template <>
inline std::vector<hsize_t>
Hdf5_serial_writer::pre_setup<double>(double const& data)
{
    data_rank = 0;
    atomic_type = hdf5_atomic_data_type<double>();
    data_size = sizeof(double);

    // dummy variable -- length really should
    // be 0, but that would not compile
    return std::vector<hsize_t>(1);
}

template <>
inline void*
Hdf5_serial_writer::get_data_ptr<int>(int const& data)
{
    return (void*)&data;
}

template <>
inline void*
Hdf5_serial_writer::get_data_ptr<double>(double const& data)
{
    return (void*)&data;
}

#endif /* HDF5_SERIAL_WRITER_H_ */
