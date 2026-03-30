#ifndef HDF5_FILE_H_
#define HDF5_FILE_H_

#include <memory>
#include <string>

#include "synergia/utils/hdf5_misc.h"
#include "synergia/utils/hdf5_reader.h"
#include "synergia/utils/hdf5_seq_writer.h"
#include "synergia/utils/hdf5_writer.h"

#include "synergia/utils/cereal.h"
#include "synergia/utils/commxx.h"

class Hdf5_file {
  public:
    enum class Flag { truncate, read_write, read_only };
    // enum Atomic_type { double_type, int_type };

    // set the type based on cmake config
#ifdef USE_PARALLEL_HDF5
    using use_parallel = std::true_type;
#else
    using use_parallel = std::false_type;
#endif

  private:
    std::shared_ptr<Commxx> comm;
    std::string file_name;
    Hdf5_handler h5file;
    int root_rank;
    bool is_open;
    Flag current_flag;
    bool has_file;

    std::map<std::string, Hdf5_seq_writer> seq_writers;

    static unsigned int
    flag_to_h5_flags(Flag flag)
    {
        if (flag == Hdf5_file::Flag::truncate) return H5F_ACC_TRUNC;
        if (flag == Hdf5_file::Flag::read_write) return H5F_ACC_RDWR;
        if (flag == Hdf5_file::Flag::read_only) return H5F_ACC_RDONLY;

        return 0;
    }

  public:
    // since C++17 the shared_from_this() called from an unmanaged
    // shared_ptr would throw the bad_weak_ptr exception, instead
    // of having an undefined behavior. So if we move to C++17 we
    // would be able to reduce the constructor to a single one by
    // testing if we are able to call the shared_from_this(),
    //    ...
    //    try { comm = c.shared_from_this();}
    //    catch(...) { comm = std::make_shared<const Commxx>(c); }
    //    ...
    //
    // for now two separate constructors are both provided
    //
    Hdf5_file(std::string const& file_name,
              Flag flag,
              std::shared_ptr<Commxx> const& comm);

    Hdf5_file(std::string const& file_name,
              Flag flag,
              Commxx const& comm = Commxx());

    ~Hdf5_file() { close(); }

    // this is a resource handler, so disable the copy and assignment
    Hdf5_file(Hdf5_file const&) = delete;
    Hdf5_file& operator=(Hdf5_file const&) = delete;

    // move is still allowed
    Hdf5_file(Hdf5_file&&) noexcept = default;

    void open(Flag flag);
    void close();
    void flush() const;

    std::string const&
    get_filename() const
    {
        return file_name;
    }

#if 0
    std::vector<std::string> get_member_names();
    Atomic_type get_atomic_type(std::string const& name);
    std::vector<int > get_dims(std::string const& name);
#endif

    hid_t
    get_h5file()
    {
        return h5file.hid;
    }

    int
    get_root_rank() const
    {
        return root_rank;
    }

    // gather on the first dimension. all other dimensions must be of the same
    // extents calling from 4 ranks:
    //   write_collective("ds", pz) -> "ds" : [pz, pz, ...]
    //   write_collective("part", part[0:1][0:6]) -> "part" : part[0:3][0:6]
    template <typename T>
    void
    write_collective(std::string const& name, T const& data) const
    {
        write(name, data, true);
    }

    // no gather, only the root rank will execute the write
    // calling from 4 ranks:
    //   write_single("ds", pz) -> "ds" : pz
    //   write_single("part", part[0:1][0:6]) -> "part" : part[0:1][0:6]
    template <typename T>
    void
    write_single(std::string const& name, T const& data) const
    {
        write(name, data, false);
    }

    template <typename T>
    void
    write(std::string const& name, T const& data, bool collective = false) const
    {
        Hdf5_writer::write(h5file, name, data, collective, *comm, root_rank);
    }

    template <typename T>
    void
    write(std::string const& name,
          T const* data,
          size_t len,
          bool collective = false) const
    {
        Hdf5_writer::write(
            h5file, name, data, len, collective, *comm, root_rank);
    }

    // same as write_single(), except this will do append instead of overwrite
    template <typename T>
    void
    append_single(std::string const& name, T const& data)
    {
        append(name, data, false);
    }

    template <typename T>
    void
    append_collective(std::string const& name, T const& data)
    {
        append(name, data, true);
    }

    // same as write(), except this will do append instead of overwrite
    template <typename T>
    void
    append(std::string const& name, T const& data, bool collective = false)
    {
        auto w = seq_writers.find(name);
        if (w == seq_writers.end()) {
            w = seq_writers
                    .emplace(name,
                             Hdf5_seq_writer(h5file, name, *comm, root_rank))
                    .first;
        }

        w->second.append(data, collective);
    }

    // read the dataset to all ranks
    template <typename T>
    T
    read(std::string const& name) const
    {
        return Hdf5_reader::read<T>(h5file, name, *comm, root_rank);
    }

    // does a collective read on the dataset
    // rank r will get a container that has the len extent in the first dim
    template <typename T>
    T
    read(std::string const& name, size_t len) const
    {
        return Hdf5_reader::read<T>(h5file, name, len, *comm, root_rank);
    }

#if 0
    // a single read from a 1d dataset into the memory
    // memory must be larger or the same as the dataset size
    template<typename T>
    void read(std::string const& name, T * const data) const
    { Hdf5_reader::read<T>(h5file, name, data, *comm, root_rank); }
#endif

    // does a collective read of the named dataset into the raw array
    template <typename T>
    void
    read(std::string const& name, T* const data, size_t len) const
    {
        Hdf5_reader::read<T>(h5file, name, data, len, *comm, root_rank);
    }

    std::vector<hsize_t>
    get_dims(std::string const& name) const
    {
        return Hdf5_reader::get_dims(h5file, name, *comm, root_rank);
    }

    std::vector<std::string>
    get_dataset_names() const
    {
        return Hdf5_reader::get_dataset_names(h5file, *comm, root_rank);
    }

    bool
    has_dataset(std::string const& name) const
    {
        return Hdf5_reader::has_dataset(h5file, name, *comm, root_rank);
    }

  private:
    friend class cereal::access;

    Hdf5_file()
        : comm()
        , file_name()
        , h5file()
        , root_rank(0)
        , is_open(false)
        , current_flag(Hdf5_file::Flag::read_only)
        , has_file(false)
    {}

    template <class Archive>
    void
    save(Archive& ar) const
    {
        ar(CEREAL_NVP(comm));
        ar(CEREAL_NVP(file_name));
        ar(CEREAL_NVP(root_rank));
        ar(CEREAL_NVP(is_open));
        ar(CEREAL_NVP(current_flag));
        ar(CEREAL_NVP(has_file));

        if (is_open) {
            flush();
            // TODO: copy_to_serialization_directory(file_name);
        }
    }

    template <class Archive>
    void
    load(Archive& ar)
    {
        ar(CEREAL_NVP(comm));
        ar(CEREAL_NVP(file_name));
        ar(CEREAL_NVP(root_rank));
        ar(CEREAL_NVP(is_open));
        ar(CEREAL_NVP(current_flag));
        ar(CEREAL_NVP(has_file));

        if (is_open) {
            // TODO: copy_from_serialization_directory(file_name)
            is_open = false;
            open(current_flag);
        }
    }
};

#endif /* HDF5_FILE_H_ */
