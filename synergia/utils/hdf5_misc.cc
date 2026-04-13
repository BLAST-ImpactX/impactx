#include <cassert>

#include "synergia/utils/commxx.h"
#include "synergia/utils/hdf5_misc.h"

namespace {
    bool
    same_nth_dim(std::vector<hsize_t> const& dims,
                 size_t dim,
                 size_t ndim,
                 size_t mpi_size)
    {
        assert(mpi_size > 0);
        assert(dims.size() == ndim * mpi_size);

        auto val = dims[dim];

        for (int r = 0; r < mpi_size; ++r)
            if (dims[r * ndim + dim] != val) return false;

        return true;
    }

    bool
    same_data_ranks(std::vector<size_t> const& ranks)
    {
        auto it = std::adjacent_find(
            ranks.begin(), ranks.end(), std::not_equal_to<>());
        return it == ranks.end();
    }
}

std::vector<hsize_t>
syn::collect_dims(std::vector<hsize_t> const& dims,
                  bool collective,
                  Commxx const& comm,
                  int root_rank)
{
    const int mpi_size = comm.size();
    const int mpi_rank = comm.rank();

    // parameter check
    if (collective && !dims.size()) {
        throw std::runtime_error(
            "Hdf5_writer: collective read/write on a scalar, "
            "should have been promoted to 1d array");
    }

    // create a local copy of the dims, always do a promotion
    // on scalar so we don't have to branch the code
    const auto local_dims = dims.size() ? dims : std::vector<hsize_t>{1};
    const auto data_rank = local_dims.size();

    // collect all data_ranks
    std::vector<size_t> all_ranks(mpi_size);
    MPI_Allgather(
        &data_rank, 1, MPI_UINT64_T, all_ranks.data(), 1, MPI_UINT64_T, comm);

    // data rank check -- all data_ranks should be the same
    if (!same_data_ranks(all_ranks)) {
        throw std::runtime_error("Hdf5_writer: inconsistent data ranks");
    }

    std::vector<hsize_t> all_dims(mpi_size * data_rank);
    std::vector<hsize_t> all_dim0(mpi_size, 0);

    // dimension check on non-scalars
    MPI_Allgather(local_dims.data(),
                  data_rank,
                  MPI_UINT64_T,
                  all_dims.data(),
                  data_rank,
                  MPI_UINT64_T,
                  comm);

    // for collective write: higher (>0) order dimensions must be the
    //   same across all ranks
    // for single write: all dimensions must be the same
    size_t d0 = collective ? 1 : 0;
    for (int d = d0; d < data_rank; ++d) {
        if (!same_nth_dim(all_dims, d, data_rank, mpi_size))
            throw std::runtime_error(
                "Hdf5_writer: inconsistent data dimensions");
    }

    if (collective) {
        // everyone has its own share
        for (int r = 0; r < mpi_size; ++r)
            all_dim0[r] = all_dims[r * data_rank];
    } else {
        // only the root rank gets a non-zero share
        all_dim0[root_rank] = local_dims[0];
    }

    return all_dim0;
}
