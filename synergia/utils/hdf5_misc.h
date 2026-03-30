#ifndef HDF5_MISC_H_
#define HDF5_MISC_H_

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <Kokkos_Core.hpp>
#include <hdf5.h>

namespace storage_order {
    constexpr int c = 0;
    constexpr int fortran = 1;

    constexpr int right = c;
    constexpr int left = fortran;

    constexpr int row = c;
    constexpr int col = fortran;

    constexpr int hdf5_default = c;
}

// hdf5 exceptions
struct Hdf5_exception : public std::exception {
    Hdf5_exception(std::string const& msg = "")
        : hdf5_msg(), user_msg(msg), what_msg()
    {
        std::stringstream buf;

        // H5Eprint(H5E_DEFAULT, stderr);
        H5Ewalk(H5E_DEFAULT, H5E_WALK_UPWARD, &Hdf5_exception::err_walk, &buf);

        hdf5_msg = buf.str();

        what_msg = std::string("\n") + "===================================\n" +
                   "USER MESSAGE:\n  " + user_msg + "\n\n" + "HDF5 MESSAGE:\n" +
                   hdf5_msg + "===================================\n";
    }

    ~Hdf5_exception() throw() {}

    const char*
    what() const noexcept override
    {
        return what_msg.c_str();
    }

  private:
    static herr_t
    err_walk(unsigned int n, H5E_error_t const* err_desc, void* client_data)
    {
        std::stringstream* ss = (std::stringstream*)client_data;

        /* Check arguments */
        if (!client_data) return 0;

        /* Get descriptions for the major and minor error numbers */
        const char* maj_str = H5Eget_major(err_desc->maj_num);
        const char* min_str = H5Eget_minor(err_desc->min_num);

        /* Print error message */
        (*ss) << "  #" << std::setfill('0') << std::setw(4) << n << ": "
              << err_desc->file_name << " line " << err_desc->line << " in "
              << err_desc->func_name << "(): " << err_desc->desc << "\n"
              << "    major(" << std::setfill('0') << std::setw(3)
              << (int)err_desc->maj_num << "): " << maj_str << "\n"
              << "    minor(" << std::setfill('0') << std::setw(3)
              << (int)err_desc->min_num << "): " << min_str << "\n";

        return 0;
    }

    std::string hdf5_msg;
    std::string user_msg;
    std::string what_msg;
};

// handles the closing of resources in the RAII way
struct Hdf5_handler {
    Hdf5_handler(hid_t handler = 0) : hid(handler)
    {
        if (hid < 0) throw Hdf5_exception("Bad HDF5 Handler");
    }

    Hdf5_handler(Hdf5_handler&& o) noexcept : hid(o.hid) { o.hid = 0; }

    Hdf5_handler&
    operator=(hid_t handler)
    {
        if (handler < 0) throw Hdf5_exception("Bad HDF5 Handler");
        hid = handler;
        return *this;
    }

    // no copy and copy assignment
    Hdf5_handler(Hdf5_handler const&) = delete;
    Hdf5_handler& operator=(Hdf5_handler const&) = delete;

    ~Hdf5_handler() { close(); }

    void
    close()
    {
        if (hid == 0) return;

        switch (H5Iget_type(hid)) {
            case H5I_FILE:
                H5Fclose(hid);
                break;
            case H5I_GROUP:
                H5Gclose(hid);
                break;
            case H5I_DATATYPE:
                H5Tclose(hid);
                break;
            case H5I_DATASPACE:
                H5Sclose(hid);
                break;
            case H5I_DATASET:
                H5Dclose(hid);
                break;
            case H5I_ATTR:
                H5Aclose(hid);
                break;
            case H5I_GENPROP_LST:
                H5Pclose(hid);
                break;
            default:
                throw std::runtime_error(
                    "close of unhandled Hdf5 handler type");
        }

        hid = 0;
    }

    operator hid_t() const { return hid; }

    bool
    valid() const
    {
        return hid != 0;
    }

    hid_t hid;
};

// The generic (T) version of h5_atomic_data_type is undefined.
// Only versions with specializations will compile.
template <typename T>
inline hid_t
hdf5_atomic_data_type()
{
    return hdf5_atomic_data_type<typename T::value_type>();
}
//{ return hdf5_traits<T>::value_type; }

template <>
inline hid_t
hdf5_atomic_data_type<uint8_t>()
{
    return H5Tcopy(H5T_NATIVE_UCHAR);
}

template <>
inline hid_t
hdf5_atomic_data_type<int>()
{
    return H5Tcopy(H5T_NATIVE_INT);
}

template <>
inline hid_t
hdf5_atomic_data_type<double>()
{
    return H5Tcopy(H5T_NATIVE_DOUBLE);
}

class Commxx;

namespace syn {
    struct data_info_t {
        void const* ptr;
        std::vector<hsize_t> dims;
        Hdf5_handler atomic_type;
        size_t atomic_data_size;
        size_t size;
    };

    template <class T>
    data_info_t
    extract_data_info(T const& t)
    {
        if constexpr (std::is_arithmetic_v<T>) {
            return {&t, {}, hdf5_atomic_data_type<T>(), sizeof(T), 1};
        } else if constexpr (Kokkos::is_view<T>::value) {
            size_t size = 1;
            std::vector<hsize_t> dims(T::rank);

            for (int i = 0; i < T::rank; ++i) {
                dims[i] = t.extent(i);
                size *= t.extent(i);
            }

            return {t.data(),
                    dims,
                    hdf5_atomic_data_type<typename T::value_type>(),
                    sizeof(typename T::value_type),
                    size};
        };
    }

    template <class T>
    void
    resize_data_obj(T& t, data_info_t& di)
    {
        if constexpr (std::is_arithmetic_v<T>) {
            if (di.dims.size())
                throw std::runtime_error(
                    "resize_data_obj: non-zero rank resizing scalar");
        } else if constexpr (Kokkos::is_view<T>::value) {
            if (di.dims.size() != T::rank)
                throw std::runtime_error(
                    "resize_data_obj: inconsistent data rank");

            switch (di.dims.size()) {
                case 1:
                    Kokkos::resize(t, di.dims[0]);
                    break;
                case 2:
                    Kokkos::resize(t, di.dims[0], di.dims[1]);
                    break;
                case 3:
                    Kokkos::resize(t, di.dims[0], di.dims[1], di.dims[2]);
                    break;
                default:
                    throw std::runtime_error("unsupported dims");
            }

            di.ptr = t.data();
        };
    }

    // collects dims of the data object thats about to be written/read
    // from all mpi ranks, does the sanity check, and produces an array
    // of all_dim0[mpi_rank] containing the portion of data responsible
    // for the particular mpi_rank.
    //
    // examples (root_rank = 0):
    //
    // no | r0        | r1        | col   | all_dim_0 | notes
    //----------------------------------------------------------------------
    // 1  | {}        | {}        | true  | error     | scalars in col write
    //    |           |           |       |           | must be promoted
    //----+-----------+-----------+-------+-----------+----------------------
    // 2  | {}        | {}        | false | [1, 0]    | write from r0 only
    //----+-----------+-----------+-------+-----------+----------------------
    // 3  | {3}       | {4}       | true  | [3, 4]    |
    //----+-----------+-----------+-------+-----------+----------------------
    // 4  | {3}       | {0}       | true  | [3, 0]    | it is allowed
    //----+-----------+-----------+-------+-----------+----------------------
    // 5  | {0}       | {0}       | true  | [0, 0]    | also allowed
    //----+-----------+-----------+-------+-----------+----------------------
    // 6  | {3}       | {4}       | false | error     | dims must be the same
    //----+-----------+-----------+-------+-----------+----------------------
    // 7  | {4}       | {4}       | false | [4, 0]    |
    //----+-----------+-----------+-------+-----------+----------------------
    // 8  | {0}       | {0}       | false | [0, 0]    |
    //----+-----------+-----------+-------+-----------+----------------------
    // 9  | {3, 3}    | {4, 4}    | true  | error     | only 1st dim can be diff
    //----+-----------+-----------+-------+-----------+----------------------
    // 10 | {3, 5}    | {4, 5}    | true  | [3, 4]    |
    //----+-----------+-----------+-------+-----------+----------------------
    // 11 | {3, 5}    | {3, 5}    | true  | [3, 3]    |
    //----+-----------+-----------+-------+-----------+----------------------
    // 12 | {3, 3}    | {4, 3}    | false | error     | dims must be the same
    //----+-----------+-----------+-------+-----------+----------------------
    // 13 | {3, 4}    | {3, 4}    | false | [3, 0]    |
    //----+-----------+-----------+-------+-----------+----------------------
    // 14 | {3, 4, 5} | {6, 4, 4} | true  | error     |
    //----+-----------+-----------+-------+-----------+----------------------
    // 15 | {3, 4, 5} | {6, 4, 5} | true  | [3, 6]    |
    //----+-----------+-----------+-------+-----------+----------------------
    // 16 | {3, 4, 5} | {6, 4, 5} | false | error     |
    //----+-----------+-----------+-------+-----------+----------------------
    // 17 | {3, 4, 5} | {3, 4, 5} | false | [3, 0]    |
    //----+-----------+-----------+-------+-----------+----------------------
    // 18 | {3}       | {3, 4}    | true  | error     | data rank must be same
    //----+-----------+-----------+-------+-----------+----------------------
    // 19 | {3}       | {3, 4}    | false | error     | data rank must be same
    //----+-----------+-----------+-------+-----------+----------------------
    //
    std::vector<hsize_t> collect_dims(std::vector<hsize_t> const& dims,
                                      bool collective,
                                      Commxx const& comm,
                                      int root_rank);
}

#endif /* HDF5_MISC_H_ */
