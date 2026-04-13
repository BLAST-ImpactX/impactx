#ifndef IMPEDANCE_H_
#define IMPEDANCE_H_

#include "wake_field.h"

#include "synergia/simulation/collective_operator.h"
#include "synergia/simulation/implemented_collective_options.h"

class Impedance;

struct Bunch_props {
  Bunch_props(int num_bunches, int max_turns)
    : num_bunches(num_bunches)
    , max_turns(max_turns)
    , registered_turns(0)
    , write_pos(0)
    , xmean("xmean", num_bunches * max_turns)
    , ymean("ymean", num_bunches * max_turns)
    , zmean("zmean", num_bunches * max_turns)
    , realnum("realnum", num_bunches * max_turns)
    , bucket_index("bucket", num_bunches * max_turns)
  {}

  // done writing a turn. increment the pointers
  void
  increment_registered_turns()
  {
    if (++registered_turns >= max_turns) registered_turns = max_turns;
    if (++write_pos >= max_turns) write_pos = 0;
  }

  // write position for current turn and the specified bunch
  KOKKOS_INLINE_FUNCTION
  int
  get_write_index(int bunch) const
  {
    return write_pos * num_bunches + bunch;
  }

  // turn: 0 is current, -1 is prev, -2 ...
  KOKKOS_INLINE_FUNCTION
  int
  get_read_index(int turn, int bunch) const
  {
    // read_pos is one behind the write_pos
    int pos = write_pos + turn - 1;
    if (pos < 0) pos += max_turns;
    return pos * num_bunches + bunch;
  }

  int num_bunches;
  int max_turns;

  int registered_turns;
  int write_pos;

  karray1d_dev xmean;
  karray1d_dev ymean;
  karray1d_dev zmean;
  karray1d_dev realnum;
  Kokkos::View<int*> bucket_index;
};

struct Bunch_params {
  double z_mean;
  double z_left;
  double cell_size_z;
  double N_factor;
  int bucket;
};

class Impedance : public Collective_operator {
private:
  const Impedance_options opts;
  std::string bunch_sim_id;

  // bunch properties over turns
  Bunch_props bps;

  // z_grid*3, in the fortran order for
  // zdensity, xmom, ymom
  karray1d_dev zbinning;
  karray1d_hst h_zbinning;

  // buffer for wake fields
  // z_grid*5, in the fortran order for
  // xwake_leading, xwake_trailing,
  // ywake_leading, ywake_trailing,
  // zwake0
  karray1d_dev wakes;
  karray1d_hst h_wakes;

  Wake_field wake_field;

private:
  void apply_impl(Bunch_simulator& simulator,
                  double time_step,
                  Logger& logger) override;

  void apply_bunch(Bunch& bunch, double time_step, Logger& logger);

  void construct_workspaces(Bunch_simulator const& sim);

  void store_bunches_data(Bunch_simulator const& sim);

  Bunch_params calculate_moments_and_partitions(Bunch const& bunch);

  void calculate_kicks(Bunch const& bunch, Bunch_params const& bp);

  void apply_impedance_kick(Bunch& bunch,
                            Bunch_params const& bp,
                            double wake_factor);

public:
  Impedance(Impedance_options const& ops);
};

#endif /* IMPEDANCE_H_ */
