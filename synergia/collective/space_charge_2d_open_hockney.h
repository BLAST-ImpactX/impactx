#ifndef SPACE_CHARGE_2D_OPEN_HOCKNEY_H_
#define SPACE_CHARGE_2D_OPEN_HOCKNEY_H_

#include "synergia/simulation/collective_operator.h"
#include "synergia/simulation/implemented_collective_options.h"

#include "synergia/utils/distributed_fft2d.h"

#include "synergia/collective/rectangular_grid.h"
#include "synergia/collective/rectangular_grid_domain.h"

class Space_charge_2d_open_hockney : public Collective_operator {

private:
  const Space_charge_2d_open_hockney_options options;

  // cached bunch simulator id
  // if the id is changed, the workspace needs to be reconstructed
  std::string bunch_sim_id;

  Rectangular_grid_domain domain;
  Rectangular_grid_domain doubled_domain;

  karray2d_dev particle_bin;

  std::array<std::vector<Distributed_fft2d>, 2> ffts;

  karray1d_dev rho2;
  karray1d_dev phi2;
  karray1d_dev g2;

  karray1d_hst h_rho2;
  karray1d_hst h_phi2;

private:
  void apply_impl(Bunch_simulator& simulator,
                  double time_step,
                  Logger& logger) override;

  void apply_bunch(Bunch& bunch,
                   Distributed_fft2d& fft,
                   double time_step,
                   Logger& logger);

  void construct_workspaces(Bunch_simulator const& sim);

  void update_domain(Bunch const& bunch);

  void get_local_charge_density(Bunch const& bunch);

  void get_global_charge_density(Bunch const& bunch);

  void get_green_fn2_pointlike();

  void get_local_force2(Distributed_fft2d& fft);

  void get_global_force2(Commxx const& comm);

  void apply_kick(Bunch& bunch, double fn_norm, double time_step);

  double get_normalization_force(Bunch const& bunch,
                                 Distributed_fft2d const& fft);

public:
  Space_charge_2d_open_hockney(Space_charge_2d_open_hockney_options const& ops);
};

#endif /* SPACE_CHARGE_2D_OPEN_HOCKNEY_H_ */
