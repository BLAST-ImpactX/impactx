#ifndef SPACE_CHARGE_3D_OPEN_HOCKNEY_H_
#define SPACE_CHARGE_3D_OPEN_HOCKNEY_H_

#include "synergia/simulation/collective_operator.h"
#include "synergia/simulation/implemented_collective_options.h"

#include "synergia/collective/rectangular_grid.h"
#include "synergia/collective/rectangular_grid_domain.h"

#include "synergia/utils/distributed_fft3d.h"

class Space_charge_3d_open_hockney;

/// Note: internal grid is stored in [z][y][x] order, but
/// grid shape expects [x][y][z] order.
class Space_charge_3d_open_hockney : public Collective_operator {
  private:
    const Space_charge_3d_open_hockney_options options;

    std::string bunch_sim_id;

    Rectangular_grid_domain domain;
    Rectangular_grid_domain doubled_domain;
    bool use_fixed_domain;

    std::array<std::vector<Distributed_fft3d>, 2> ffts;

    karray1d_dev rho2;
    karray1d_dev phi2;
    karray1d_dev g2;

    karray1d_hst h_rho2;
    karray1d_hst h_phi2;

    karray1d_dev enx;
    karray1d_dev eny;
    karray1d_dev enz;

  private:
    void apply_impl(Bunch_simulator& simulator,
                    double time_step,
                    Logger& logger);

    void apply_bunch(Bunch& bunch,
                     Distributed_fft3d& fft,
                     double time_step,
                     Logger& logger);

    void construct_workspaces(Bunch_simulator const& sim);

    void update_domain(Bunch const& bunch);

    void get_local_charge_density(Bunch const& bunch);

    void get_global_charge_density(Bunch const& bunch);

    void apply_kick(Bunch& bunch, double fn_norm, double time_step);

    void get_green_fn2_pointlike();
    void get_green_fn2_linear();

    void get_local_phi2(Distributed_fft3d& fft);

    void get_global_phi2(Distributed_fft3d const& fft);

    void get_force();

    double get_normalization_force(Distributed_fft3d const& fft);

  public:
    Space_charge_3d_open_hockney(
        Space_charge_3d_open_hockney_options const& ops);
};

#endif /* SPACE_CHARGE_3D_OPEN_HOCKNEY_H_ */
