#include "synergia/simulation/checkpoint.h"

#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/propagator.h"

namespace syn {
    void checkpoint_save_as_json(std::string const& prop_str,
                                 std::string const& sims_str,
                                 std::vector<int> const& displs,
                                 std::vector<int> const& lens);

    std::pair<std::string, std::string> checkpoint_load_json(
        std::vector<char> const& buf,
        int rank);
}

void
syn::checkpoint_save(Propagator const& prop, Bunch_simulator const& sim)
{
    std::string prop_str = prop.dump();
    std::string sim_str = sim.dump();

    // collect sim_str to the root rank
    auto const& comm = sim.get_comm();

    const int root = 0;
    const int mpi_size = comm.size();
    const int mpi_rank = comm.rank();

    int len = sim_str.size();
    std::vector<int> lens(mpi_size, 0);

    // gather string size
    MPI_Gather(&len, 1, MPI_INT, lens.data(), 1, MPI_INT, root, comm);

    // accumulate sizes
    std::vector<int> displs(mpi_size, 0);
    for (int i = 1; i < mpi_size; ++i)
        displs[i] = displs[i - 1] + lens[i - 1];

    // recv string buffer
    int total_len = displs[mpi_size - 1] + lens[mpi_size - 1];
    std::string sims_str(total_len, ' ');

    // gather strings
    MPI_Gatherv(sim_str.data(),
                sim_str.size(),
                MPI_CHAR,
                (void*)sims_str.data(),
                lens.data(),
                displs.data(),
                MPI_CHAR,
                root,
                comm);

    // extract each string and parse into a JSON object
    if (mpi_rank == root)
        checkpoint_save_as_json(prop_str, sims_str, displs, lens);
}

std::pair<Propagator, Bunch_simulator>
syn::checkpoint_load()
{
    const int root = 0;
    const int mpi_size = Commxx::world_size();
    const int mpi_rank = Commxx::world_rank();

    size_t len;
    std::vector<char> buf;

    if (mpi_rank == root) {
        // read states from file
        std::ifstream file("cp_state.json");
        if (!file.good())
            throw std::runtime_error("Error at openning checkpointing file");

        file.seekg(0, std::ios::end);
        len = file.tellg();

        // broadcast the buffer length
        MPI_Bcast(&len, 1, MPI_UINT64_T, root, Commxx::World);

        // read the buffer
        buf.resize(len);
        file.seekg(0);
        file.read(&buf[0], len);

        // broadcast the buffer
        MPI_Bcast(&buf[0], len, MPI_BYTE, root, Commxx::World);
    } else {
        // receive the buffer length
        MPI_Bcast(&len, 1, MPI_UINT64_T, root, Commxx::World);

        // prepare buffer
        buf.resize(len);

        // receive the buffer
        MPI_Bcast(&buf[0], len, MPI_BYTE, root, Commxx::World);
    }

    // parse the json object
    auto cp = syn::checkpoint_load_json(buf, mpi_rank);

    // recreate the objects
    return std::make_pair(Propagator::load_from_string(cp.first),
                          Bunch_simulator::load_from_string(cp.second));
}

void
syn::resume()
{
    Logger simlog(0, LoggerV::INFO_STEP);

    auto cp = syn::checkpoint_load();

    Propagator& propagator = cp.first;
    Bunch_simulator& sim = cp.second;

    propagator.propagate(sim, simlog);
}
