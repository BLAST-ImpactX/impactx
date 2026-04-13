#ifndef SIMULATION_BUNCH_SIMULATOR_H
#define SIMULATION_BUNCH_SIMULATOR_H

#include "synergia/bunch/bunch_train.h"
#include "synergia/lattice/lattice.h"
#include "synergia/simulation/propagate_actions.h"

#include <cereal/types/utility.hpp> // std::pair
#include <cereal/types/vector.hpp>

class Operator;
class Independent_operation;

class Diagnostics_loss;
class Bunch_simulator;

class Bunch_simulator {

  public:
    static const int PRE_TURN = -1;
    static const int FINAL_STEP = -1;

  private:
#if 0
    // bool trigger(int turn, int step)
    using trigger_step_t = std::function<bool(int, int)>;

    // bool trigger(Lattice_element const& ele)
    using trigger_ele_t  = std::function<bool(Lattice_element const&)>;

    // bool trigger(Operator const& opr)
    using trigger_opr_t  = std::function<bool(Operator const&)>;

    // bool trigger(Operation const& opn)
    using trigger_opn_t  = std::function<bool(Independent_operation const&)>;
#endif

    struct trigger_step_period {
        int turn_period, step_period;

        bool
        operator()(int turn, int step) const
        {
            // always trigger at before_start
            if (turn == PRE_TURN && step == FINAL_STEP) return true;

            // only the turn number matters
            if (step_period < 0)
                return (turn % turn_period == 0) && (step == FINAL_STEP);

            // at both the right turn and step
            return (turn % turn_period == 0) && (step % step_period == 0);
        }

        template <class AR>
        void
        serialize(AR& ar)
        {
            ar(turn_period, step_period);
        }
    };

    struct trigger_turn_listed {
        std::vector<int> numbers;

        bool
        operator()(int turn, int step) const
        {
            return (step == FINAL_STEP) &&
                   std::find(numbers.begin(), numbers.end(), turn) !=
                       numbers.end();
        }

        template <class AR>
        void
        serialize(AR& ar)
        {
            ar(numbers);
        }
    };

    struct trigger_element {
        std::string name;
        int turn_period;

        bool
        operator()(int turn, Lattice_element const& ele) const
        {
            return (turn % turn_period == 0) && (ele.get_name() == name);
        }

        template <class AR>
        void
        serialize(AR& ar)
        {
            ar(name, turn_period);
        }
    };

    template <typename TriggerT>
    struct diag_tuple_t {
        int train;
        int bunch;
        int diag_id;
        TriggerT trigger;

        template <class AR>
        void
        serialize(AR& ar)
        {
            ar(train, bunch, diag_id, trigger);
        }
    };

    using dt_step_period = diag_tuple_t<trigger_step_period>;
    using dt_turn_listed = diag_tuple_t<trigger_turn_listed>;
    using dt_element = diag_tuple_t<trigger_element>;

    // void action(Bunch_simulator&, Lattice&, int turn, int step, void* data)
    using action_step_t =
        std::function<void(Bunch_simulator&, Lattice&, int, int)>;
    using action_data_step_t =
        std::function<void(Bunch_simulator&, Lattice&, int, int, void*)>;

    // void action(Bunch_simulator&, Lattice&, int turn, void* data)
    using action_turn_t = std::function<void(Bunch_simulator&, Lattice&, int)>;
    using action_data_turn_t =
        std::function<void(Bunch_simulator&, Lattice&, int, void*)>;

  private:
    // constructor
    Bunch_simulator(Bunch_train&& pt,
                    Bunch_train&& st,
                    std::shared_ptr<Commxx> const& comm);

    static Bunch_simulator construct(Reference_particle const& ref_pri,
                                     Reference_particle const& ref_sec,
                                     size_t num_part,
                                     size_t num_spec,
                                     double num_real_part,
                                     size_t num_bunches_pri,
                                     size_t num_bunches_sec,
                                     double spacing_pri,
                                     double spacing_sec,
                                     Commxx const& comm);

    int get_bunch_array_idx(int train, int bunch) const;

  public:
    Bunch_simulator(Bunch_simulator const&) = delete;
    Bunch_simulator(Bunch_simulator&&) = default;

    // factory methods
    static Bunch_simulator create_empty_bunch_simulator();

    static Bunch_simulator create_single_bunch_simulator(
        Reference_particle const& ref,
        size_t num_particles,
        double num_real_particles,
        Commxx const& comm = Commxx(),
        size_t num_spectators = 0);

    static Bunch_simulator create_bunch_train_simulator(
        Reference_particle const& ref,
        size_t num_particles,
        double num_real_particles,
        size_t num_bunches,
        double spacing,
        Commxx const& comm = Commxx(),
        size_t num_spectators = 0);

    static Bunch_simulator create_two_trains_simulator(
        Reference_particle const& ref_pri,
        Reference_particle const& ref_sec,
        size_t num_particles,
        double num_real_particles,
        size_t num_bunches_pri,
        size_t num_bunches_sec,
        double spacing_pri,
        double spacing_sec,
        Commxx const& comm = Commxx(),
        size_t num_spectators = 0);

    // bunch simulator id
    std::string const&
    id() const
    {
        return uuid;
    }

    // register a per-turn diagnostics (with an optional turn period)
    // the train and bunch are indexed based on actual number of
    // bunches per train
    template <class Diag>
    Diagnostics_handler
    reg_diag_period(Diag const& diag,
                    int turn_period,
                    int step_period,
                    int bunch,
                    int train)
    {
        if (turn_period <= 0 || step_period == 0)
            throw std::runtime_error("reg_diag_per_turn/step(): "
                                     "the turn or step period cannot be <=0");

        int bunch_idx = get_bunch_array_idx(train, bunch);
        if (bunch_idx == -1) return Diagnostics_handler();

        auto handler = trains[train][bunch_idx].add_diagnostics(diag);

        dt_step_period dt{train,
                          bunch_idx,
                          handler.second,
                          trigger_step_period{turn_period, step_period}};

        diags_step_period.push_back(dt);
        return handler.first;
    }

    template <class Diag>
    Diagnostics_handler
    reg_diag_per_turn(Diag const& diag,
                      int period = 1,
                      int bunch = 0,
                      int train = 0)
    {
        return reg_diag_period(diag, period, -1, bunch, train);
    }

    template <class Diag>
    Diagnostics_handler
    reg_diag_per_step(Diag const& diag,
                      int period = 1,
                      int bunch = 0,
                      int train = 0)
    {
        return reg_diag_period(diag, 1, period, bunch, train);
    }

    // register a diagnsotics at specified number of turns
    template <class Diag>
    Diagnostics_handler
    reg_diag_turn_listed(Diag const& diag,
                         std::vector<int> const& numbers = {},
                         int bunch = 0,
                         int train = 0)
    {
        int bunch_idx = get_bunch_array_idx(train, bunch);
        if (bunch_idx == -1) return Diagnostics_handler();

        auto handler = trains[train][bunch_idx].add_diagnostics(diag);

        dt_turn_listed dt{
            train, bunch_idx, handler.second, trigger_turn_listed{numbers}};

        diags_turn_listed.push_back(dt);
        return handler.first;
    }

    // register a diagnostics at specific element and
    // the period of turns
    // * this is the former forced diagnostics
    template <class Diag>
    Diagnostics_handler
    reg_diag_at_element(Diag const& diag,
                        Lattice_element const& ele,
                        int turn_period = 1,
                        int bunch = 0,
                        int train = 0)
    {
        if (turn_period <= 0)
            throw std::runtime_error("reg_diag_at_element(): "
                                     "the turn period cannot be <=0");

        int bunch_idx = get_bunch_array_idx(train, bunch);
        if (bunch_idx == -1) return Diagnostics_handler();

        auto handler = trains[train][bunch_idx].add_diagnostics(diag);

        dt_element dt{train,
                      bunch_idx,
                      handler.second,
                      trigger_element{ele.get_name(), turn_period}};

        diags_element.push_back(dt);
        return handler.first;
    }

    // diag loss
    void
    reg_diag_loss_aperture(std::string const& filename,
                           int bunch = 0,
                           int train = 0)
    {
        get_bunch(train, bunch).set_diag_loss_aperture(filename);
    }

    void
    reg_diag_loss_zcut(std::string const& filename,
                       int bunch = 0,
                       int train = 0)
    {
        get_bunch(train, bunch).set_diag_loss_zcut(filename);
    }

#if 0
    // diag per operator
    void reg_diag_at_operator(
            std::string const& name, Diagnostics & diag,
            trigger_opr_t, int bunch = 0, int train = 0 )
    { }

    // diag per operation
    void reg_diag_at_operation(
            std::string const& name, Diagnostics & diag,
            trigger_opn_t, int bunch = 0, int train = 0 )
    { }
#endif

    // register propagation actions
    template <class PA>
    void
    reg_prop_actions(PA const& pa)
    {
        prop_actions = std::make_shared<PA>(pa);
    }

    // register individual propagation actions
    void reg_prop_action_step_end(action_step_t fun);
    void reg_prop_action_turn_end(action_turn_t fun);

    void reg_prop_action_step_end(action_data_step_t fun, void* data = nullptr);
    void reg_prop_action_turn_end(action_data_turn_t fun, void* data = nullptr);

    // diag actions
    void diag_action_step_and_turn(int turn_num, int step_num);
    void diag_action_element(Lattice_element const& element);
    void diag_action_operator(Operator const& opr);
    void diag_action_operation(Independent_operation const& opn);

    // propagate actions
    void prop_action_first(Lattice& lattice);
    void prop_action_step_end(Lattice& lattice, int turn, int step);
    void prop_action_turn_end(Lattice& lattice, int turn);

    // set reference particle
    void set_lattice_reference_particle(Reference_particle const& ref);

    // accessors
    Commxx const&
    get_comm() const
    {
        return *comm;
    }

    std::array<Bunch_train, 2>&
    get_trains()
    {
        return trains;
    }
    std::array<Bunch_train, 2> const&
    get_trains() const
    {
        return trains;
    }

    Bunch_train&
    operator[](size_t idx)
    {
        return trains[idx];
    }
    Bunch_train const&
    operator[](size_t idx) const
    {
        return trains[idx];
    }

    // retrive the bunch according to its index in the train,
    // not the index in the array
    Bunch& get_bunch(size_t train = 0, size_t bunch = 0);
    Bunch const& get_bunch(size_t train = 0, size_t bunch = 0) const;

    // test. returns false when train and bunch out of bounds
    bool has_local_bunch(size_t train, size_t bunch) const;

    // which ranks are a given bunch on
    std::vector<int> get_bunch_ranks(size_t train, size_t bunch) const;

    // turns
    void
    inc_turn()
    {
        ++curr_turn;
    }
    void
    set_max_turns(int turns)
    {
        num_turns = turns;
    }

    int
    current_turn() const
    {
        return curr_turn;
    }
    int
    max_turns() const
    {
        return num_turns;
    }

    // longitudinal boundary
    void
    set_longitudinal_boundary(LongitudinalBoundary lb, double param = 0.0)
    {
        trains[0].set_longitudinal_boundary(lb, param);
        trains[1].set_longitudinal_boundary(lb, param);
    }

    // coordinated populate
    // the populated bunch simulator will be deterministic
    // regardless of the number of ranks in the simulation
    void populate_6d(uint64_t seed,
                     const_karray1d means,
                     const_karray2d_row covariances);

    void populate_6d_truncated(uint64_t seed,
                               const_karray1d means,
                               const_karray2d_row covariances,
                               const_karray1d limits);

    // serialization helper
    void save_checkpoint_particles(std::string const& fname);
    void load_checkpoint_particles(std::string const& fname);

    std::string
    dump() const
    {
        std::stringstream ss;
        {
            cereal::JSONOutputArchive ar(ss);
            ar(*this);
        }
        return ss.str();
    }

    static Bunch_simulator
    load_from_string(std::string const& str)
    {
        std::stringstream ss(str);
        cereal::JSONInputArchive ar(ss);

        auto bs = create_empty_bunch_simulator();
        ar(bs);

        return bs;
    }

  private:
    std::string uuid;

    int curr_turn = 0;  // current progress in turns
    int num_turns = -1; // total number of turns (-1 no limit)

    std::shared_ptr<Commxx> comm;
    std::array<Bunch_train, 2> trains;

    // diagnostics action trigger conditions
    std::vector<dt_step_period> diags_step_period;
    std::vector<dt_turn_listed> diags_turn_listed;
    std::vector<dt_element> diags_element;

    // it would survive the checkpoint load
    std::shared_ptr<Propagate_actions> prop_actions;

    // non-persistent propagate actions --
    // reg again after checkpoint load
    std::vector<action_step_t> prop_actions_step_end;
    std::vector<action_turn_t> prop_actions_turn_end;

  private:
    friend class cereal::access;

    template <class AR>
    void
    serialize(AR& ar)
    {
        ar(CEREAL_NVP(uuid));

        ar(CEREAL_NVP(curr_turn));
        ar(CEREAL_NVP(num_turns));

        ar(CEREAL_NVP(comm));
        ar(CEREAL_NVP(trains));

        ar(CEREAL_NVP(diags_step_period));
        ar(CEREAL_NVP(diags_turn_listed));
        ar(CEREAL_NVP(diags_element));

        ar(CEREAL_NVP(prop_actions));

        // save/load particles with parallel hdf5
        if (AR::is_saving::value) {
            // save particles
            std::stringstream ss;
            ss << "bunch_simulator.h5";

            std::string particle_fname = ss.str();
            ar(CEREAL_NVP(particle_fname));

            save_checkpoint_particles(particle_fname);
        } else {
            // load particle
            std::string particle_fname;
            ar(CEREAL_NVP(particle_fname));

            load_checkpoint_particles(particle_fname);
        }
    }
};

template <>
inline void
Bunch_simulator::reg_prop_actions(std::shared_ptr<Propagate_actions> const& pa)
{
    prop_actions = pa;
}

#endif
