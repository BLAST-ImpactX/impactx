#ifndef SPACE_CHARGE_RECTANGULAR_H_
#define SPACE_CHARGE_RECTANGULAR_H_

#include "synergia/simulation/collective_operator.h"
#include "synergia/simulation/implemented_collective_options.h"

#include "synergia/collective/rectangular_grid_domain.h"
#include "synergia/utils/distributed_fft3d_rect.h"

class Space_charge_rectangular : public Collective_operator {

private:
  const Space_charge_rectangular_options options;

  std::string bunch_sim_id;

  Rectangular_grid_domain domain;

  std::array<std::vector<Distributed_fft3d_rect>, 2> ffts;

  karray1d_dev rho;
  karray1d_dev phi;
  karray1d_dev phihat;

  karray1d_hst h_rho;
  karray1d_hst h_phi;

  karray1d_dev enx;
  karray1d_dev eny;
  karray1d_dev enz;

private:
  void apply_impl(Bunch_simulator& simulator, double time_step, Logger& logger);

  void apply_bunch(Bunch& bunch,
                   Distributed_fft3d_rect& fft,
                   double time_step,
                   Logger& logger);

  void construct_workspaces(Bunch_simulator const& sim);

  void update_domain(Bunch const& bunch);

  void get_local_charge_density(Bunch const& bunch);

  void get_global_charge_density(Bunch const& bunch);

  void get_local_phi(Distributed_fft3d_rect& fft, double gamma);

  void get_global_phi(Distributed_fft3d_rect const& fft);

  void extract_force();
  double get_normalization_force();

  void apply_kick(Bunch& bunch, double fn_norm, double time_step);

public:
  Space_charge_rectangular(Space_charge_rectangular_options const& ops);
};

#endif /* SPACE_CHARGE_RECTANGULAR_H_ */
