#include <numeric>
#include <random>
#include <sstream>

#include "synergia/bunch/populate_global.h"
#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/independent_operation.h"
#include "synergia/simulation/operator.h"

namespace impl {

    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_int_distribution<> dis(0, 15);
    static std::uniform_int_distribution<> dis2(8, 11);

    std::string
    generate_uuid_v4()
    {
        std::stringstream ss;
        int i;
        ss << std::hex;
        for (i = 0; i < 8; i++) {
            ss << dis(gen);
        }
        ss << "-";
        for (i = 0; i < 4; i++) {
            ss << dis(gen);
        }
        ss << "-4";
        for (i = 0; i < 3; i++) {
            ss << dis(gen);
        }
        ss << "-";
        ss << dis2(gen);
        for (i = 0; i < 3; i++) {
            ss << dis(gen);
        }
        ss << "-";
        for (i = 0; i < 12; i++) {
            ss << dis(gen);
        };
        return ss.str();
    }

    void
    divide_bunches(int size,
                   size_t num_bunches_pri,
                   size_t num_bunches_sec,
                   std::vector<int>& p_ranks,
                   std::vector<int>& s_ranks)
    {
        assert(size > 0);

        const size_t num_bunches = num_bunches_pri + num_bunches_sec;
        size_t st_start_rank = 0;

        if (num_bunches == 0) {
            // no bunch at all
            p_ranks.resize(0);
            s_ranks.resize(0);
            st_start_rank = 0;
        } else if (size == 1) {
            // all bunches on a single rank
            p_ranks.resize(num_bunches_pri ? 1 : 0);
            s_ranks.resize(num_bunches_sec ? 1 : 0);
            st_start_rank = 0;
        } else if (size < num_bunches) {
            // multiple ranks, each rank stores 1 train at max
            if (num_bunches % size != 0) {
                throw std::runtime_error(
                    "Bunch_simulator::create_bunch_train_simulator() "
                    "the number of bunches must be divisible by the number of "
                    "ranks");
            }

            int const bunch_per_rank = num_bunches / size;

            if (num_bunches_pri % bunch_per_rank != 0 ||
                num_bunches_sec % bunch_per_rank != 0) {
                throw std::runtime_error(
                    "Bunch_simulator::create_bunch_train_simulator() "
                    "the number of bunches in primary or secondary train must "
                    "be "
                    "divisible by the number of bunches per rank");
            }

            p_ranks.resize(num_bunches_pri / bunch_per_rank);
            s_ranks.resize(num_bunches_sec / bunch_per_rank);
            st_start_rank = p_ranks.size();
        } else {
            // now size >= num_bunches, one bunch (or a fraciton of a bunch) per
            // rank
            if (size % num_bunches != 0) {
                throw std::runtime_error(
                    "Bunch_simulator::create_bunch_train_simulator() "
                    "the number of ranks must be divisible by the number of "
                    "bunches");
            }

            int const rank_per_bunch = size / num_bunches;

            p_ranks.resize(rank_per_bunch * num_bunches_pri);
            s_ranks.resize(rank_per_bunch * num_bunches_sec);
            st_start_rank = p_ranks.size();
        }

        std::iota(p_ranks.begin(), p_ranks.end(), 0);
        std::iota(s_ranks.begin(), s_ranks.end(), st_start_rank);

        return;
    }

}

Bunch_simulator
Bunch_simulator::create_empty_bunch_simulator()
{
    return construct(Reference_particle(),
                     Reference_particle(),
                     0,   // num_part
                     0,   // num_spectator
                     1.0, // num_real_particle
                     0,   // primary bunch
                     0,   // secondary bunch
                     1.0,
                     1.0, // spacing
                     Commxx());
}

Bunch_simulator
Bunch_simulator::create_single_bunch_simulator(Reference_particle const& ref,
                                               size_t num_particles,
                                               double num_real_particles,
                                               Commxx const& comm,
                                               size_t num_spectators)
{
    return construct(ref,
                     ref,
                     num_particles,
                     num_spectators,
                     num_real_particles,
                     1,
                     0, // num_bunches
                     1.0,
                     1.0, // spacing
                     comm);
}

Bunch_simulator
Bunch_simulator::create_bunch_train_simulator(Reference_particle const& ref,
                                              size_t num_particles,
                                              double num_real_particles,
                                              size_t num_bunches,
                                              double spacing,
                                              Commxx const& comm,
                                              size_t num_spectators)
{
    return construct(ref,
                     ref,
                     num_particles,
                     num_spectators,
                     num_real_particles,
                     num_bunches,
                     0, // num_bunches
                     spacing,
                     1.0, // spacing
                     comm);
}

Bunch_simulator
Bunch_simulator::create_two_trains_simulator(Reference_particle const& ref_pri,
                                             Reference_particle const& ref_sec,
                                             size_t num_particles,
                                             double num_real_particles,
                                             size_t num_bunches_pri,
                                             size_t num_bunches_sec,
                                             double spacing_pri,
                                             double spacing_sec,
                                             Commxx const& comm,
                                             size_t num_spectators)
{
    return construct(ref_pri,
                     ref_sec,
                     num_particles,
                     num_spectators,
                     num_real_particles,
                     num_bunches_pri,
                     num_bunches_sec,
                     spacing_pri,
                     spacing_sec,
                     comm);
}

Bunch_simulator
Bunch_simulator::construct(Reference_particle const& ref_pri,
                           Reference_particle const& ref_sec,
                           size_t num_part,
                           size_t num_spec,
                           double num_real_part,
                           size_t num_bunches_pri,
                           size_t num_bunches_sec,
                           double spacing_pri,
                           double spacing_sec,
                           Commxx const& comm)
{
    auto comm_ptr = std::make_shared<Commxx>(comm);

    std::vector<int> p_ranks;
    std::vector<int> s_ranks;

    impl::divide_bunches(
        comm.size(), num_bunches_pri, num_bunches_sec, p_ranks, s_ranks);

    auto comm_pri = comm_ptr->group(p_ranks);
    auto comm_sec = comm_ptr->group(s_ranks);

    return Bunch_simulator(Bunch_train(ref_pri,
                                       num_bunches_pri,
                                       num_part,
                                       num_real_part,
                                       spacing_pri,
                                       comm_pri,
                                       num_spec,
                                       0),
                           Bunch_train(ref_sec,
                                       num_bunches_sec,
                                       num_part,
                                       num_real_part,
                                       spacing_sec,
                                       comm_sec,
                                       num_spec,
                                       1),
                           comm_ptr);

    // A more general approach to the above logic could be like:
    //
    // bunch_ranks[ti][bi][ri]:
    //
    //   stores the MPI rank in the simulator communicator of given
    //   bunch index and rank inex
    //
    //   ti - index of the train (0 or 1)
    //   bi - index of the bunch (0 to num_bunches of the train)
    //   ri - index from 0 to number of ranks of this bunch
    //
    // bunch_idx_map[ti][bi]
    //
    //   stores the array index in the local bunch array of the
    //   given train/bunch. -1 if the bunch is not present in the
    //   current rank
    //
    //   ti: index of the train, 0 or 1
    //   bi: index of the bunch, 0 to num_bunches in the train
    //
    // E.g.:
    //
    //   1. num_bunches_pri = 1, num_bunches_sec = 1
    //      total 2 bunches, 8 ranks (n_ranks > n_bunches)
    //
    //      bunch_ranks = {
    //          { {0, 1, 2, 3} },
    //          { {4, 5, 6, 7} }
    //      }
    //
    //      First bunch spans across 4 MPI processors of rank (0, 1,
    //      2, and 3), in the context of the bunch simulator
    //      communicator. The second bunch on rank (4, 5, 6, and 7).
    //
    //      bunch_idx_map = { {0}, {-1} } on rank 0
    //      bunch_idx_map = { {0}, {-1} } on rank 1
    //      bunch_idx_map = { {0}, {-1} } on rank 2
    //      bunch_idx_map = { {0}, {-1} } on rank 3
    //      bunch_idx_map = { {-1}, {0} } on rank 4
    //      bunch_idx_map = { {-1}, {0} } on rank 5
    //      bunch_idx_map = { {-1}, {0} } on rank 6
    //      bunch_idx_map = { {-1}, {0} } on rank 7
    //
    //
    //   2. num_bunches_pri = 2, num_bunches_sec = 2
    //      4 bunches, 4 ranks (n_ranks = n_bunches)
    //
    //      bunch_ranks = {
    //          { {0}, {1} },
    //          { {2}, {3} }
    //      }
    //
    //      Each of the bunch takes over a single rank from 0 to 3.
    //
    //      bunch_idx_map = { { 0, -1}, {-1, -1} } on rank 0
    //      bunch_idx_map = { {-1,  0}, {-1, -1} } on rank 1
    //      bunch_idx_map = { {-1, -1}, { 0, -1} } on rank 2
    //      bunch_idx_map = { {-1, -1}, {-1,  0} } on rank 3
    //
    //
    //   3. num_bunches_pri = 2, num_bunches_sec = 2
    //      4 bunches, 2 ranks (n_ranks < n_bunches)
    //
    //      bunch_ranks = {
    //          { {0}, {0} },
    //          { {1], {1} }
    //      }
    //
    //      The first two bunches reside on rank 0, and the next two
    //      on rank 1.
    //
    //      bunch_idx_map = { { 0,  1}, {-1, -1} } on rank 0
    //      bunch_idx_map = { {-1, -1}, { 0,  1} } on rank 1
    //
    //   4. num_bunches_pri = 4, num_bunches_sec = 2
    //      6 bunches, 3 ranks (n_ranks < n_bunches)
    //
    //      bunch_ranks = {
    //          { {0}, {0}, {1}, {1} },
    //          { {2}, {2} }
    //      }
    //
    //      bunch_idx_map = { { 0,  1, -1, -1}, {-1, -1} } on rank 0
    //      bunch_idx_map = { {-1, -1,  0,  1}, {-1, -1} } on rank 1
    //      bunch_idx_map = { {-1, -1, -1, -1}, { 0,  1} } on rank 2
    //
}

Bunch_simulator::Bunch_simulator(Bunch_train&& pt,
                                 Bunch_train&& st,
                                 std::shared_ptr<Commxx> const& comm)
    : uuid(impl::generate_uuid_v4())
    , comm(std::move(comm))
    , trains{std::move(pt), std::move(st)}
    , diags_step_period()
    , diags_turn_listed()
    , diags_element()
    , prop_actions()
    , prop_actions_step_end()
    , prop_actions_turn_end()
{}

int
Bunch_simulator::get_bunch_array_idx(int train, int bunch) const
{
    if (train > 1) return -1;
    if (bunch >= trains[train].get_num_bunches()) return -1;

    return trains[train].get_array_idx_of_bunch(bunch);
}

bool
Bunch_simulator::has_local_bunch(size_t train, size_t bunch) const
{
    return get_bunch_array_idx(train, bunch) != -1;
}

Bunch&
Bunch_simulator::get_bunch(size_t train, size_t bunch)
{
    auto idx = get_bunch_array_idx(train, bunch);
    if (idx == -1) throw std::runtime_error("bunch not avaialble on the rank");
    return trains[train][idx];
}

Bunch const&
Bunch_simulator::get_bunch(size_t train, size_t bunch) const
{
    auto idx = get_bunch_array_idx(train, bunch);
    if (idx == -1) throw std::runtime_error("bunch not avaialble on the rank");
    return trains[train][idx];
}

std::vector<int>
Bunch_simulator::get_bunch_ranks(size_t train, size_t bunch) const
{
    int rank_per_bunch =
        std::ceil(1.0 * comm->size() /
                  (trains[0].get_num_bunches() + trains[1].get_num_bunches()));

    std::vector<int> ranks(rank_per_bunch);
    std::iota(ranks.begin(), ranks.end(), bunch * rank_per_bunch);

    return ranks;
}

void
Bunch_simulator::diag_action_step_and_turn(int turn_num, int step_num)
{
    for (auto const& dt : diags_step_period) {
        if (dt.trigger(turn_num, step_num)) {
            trains[dt.train][dt.bunch].diag_update_and_write(dt.diag_id);
        }
    }

    for (auto const& dt : diags_turn_listed) {
        if (dt.trigger(turn_num, step_num)) {
            trains[dt.train][dt.bunch].diag_update_and_write(dt.diag_id);
        }
    }
}

void
Bunch_simulator::diag_action_element(Lattice_element const& element)
{
    for (auto const& dt : diags_element) {
        if (dt.trigger(current_turn(), element)) {
            trains[dt.train][dt.bunch].diag_update_and_write(dt.diag_id);
        }
    }
}

void
Bunch_simulator::diag_action_operator(Operator const& opr)
{}

void
Bunch_simulator::diag_action_operation(Independent_operation const& opn)
{}

void
Bunch_simulator::reg_prop_action_step_end(action_step_t fun)
{
    prop_actions_step_end.push_back(fun);
}

void
Bunch_simulator::reg_prop_action_turn_end(action_turn_t fun)
{
    prop_actions_turn_end.push_back(fun);
}

void
Bunch_simulator::reg_prop_action_step_end(action_data_step_t fun, void* data)
{
    using namespace std::placeholders;
    auto fun2 = std::bind(fun, _1, _2, _3, _4, data);
    prop_actions_step_end.push_back(fun2);
}

void
Bunch_simulator::reg_prop_action_turn_end(action_data_turn_t fun, void* data)
{
    using namespace std::placeholders;
    auto fun2 = std::bind(fun, _1, _2, _3, data);
    prop_actions_turn_end.push_back(fun2);
}

void
Bunch_simulator::prop_action_first(Lattice& lattice)
{
    if (prop_actions) prop_actions->first(*this, lattice);
}

void
Bunch_simulator::prop_action_step_end(Lattice& lattice, int turn, int step)
{
    if (prop_actions) prop_actions->step_end(*this, lattice, turn, step);

    for (auto const& action : prop_actions_step_end)
        action(*this, lattice, turn, step);
}

void
Bunch_simulator::prop_action_turn_end(Lattice& lattice, int turn)
{
    if (prop_actions) prop_actions->turn_end(*this, lattice, turn);

    for (auto const& action : prop_actions_turn_end)
        action(*this, lattice, turn);
}

void
Bunch_simulator::set_lattice_reference_particle(Reference_particle const& ref)
{
    for (auto& train : get_trains()) {
        for (auto& bunch : train.get_bunches()) {
            bunch.set_design_reference_particle(ref);
        }
    }
}

namespace {
    std::vector<Bunch const*>
    get_bunch_ptrs(std::array<Bunch_train, 2> const& trains)
    {
        std::vector<Bunch const*> bunches;

        for (auto const& t : trains)
            for (auto const& b : t.get_bunches())
                bunches.push_back(&b);

        return bunches;
    }

    std::vector<Bunch*>
    get_bunch_ptrs(std::array<Bunch_train, 2>& trains)
    {
        std::vector<Bunch*> bunches;

        for (auto& t : trains)
            for (auto& b : t.get_bunches())
                bunches.push_back(&b);

        return bunches;
    }
}

void
Bunch_simulator::save_checkpoint_particles(std::string const& fname)
{
    Hdf5_file file(fname, Hdf5_file::Flag::truncate, *comm);
    auto bunches = get_bunch_ptrs(trains);

    for (int i = 0; i < bunches.size(); ++i)
        bunches[i]->save_checkpoint_particles(file, i);
}

void
Bunch_simulator::load_checkpoint_particles(std::string const& fname)
{
    Hdf5_file file(fname, Hdf5_file::Flag::read_only, *comm);
    auto bunches = get_bunch_ptrs(trains);

    for (int i = 0; i < bunches.size(); ++i)
        bunches[i]->load_checkpoint_particles(file, i);
}

void
Bunch_simulator::populate_6d(uint64_t seed,
                             const_karray1d means,
                             const_karray2d_row covariances)
{
    karray1d limits("limits", 6);
    for (int i = 0; i < 6; ++i)
        limits[i] = 0.0;

    populate_6d_truncated(seed, means, covariances, limits);
}

void
Bunch_simulator::populate_6d_truncated(uint64_t seed,
                                       const_karray1d means,
                                       const_karray2d_row covariances,
                                       const_karray1d limits)
{
    int train_idx = 0;

    for (auto& train : get_trains()) {
        for (auto& bunch : train.get_bunches()) {
            // assign unique particle ids across
            // the bunch simulator,
            bunch.assign_particle_ids(train_idx);

            // coordinated populate depends on the
            // global particle id
            populate_global_6d_truncated(
                seed, bunch, means, covariances, limits);
        }

        ++train_idx;
    }
}
