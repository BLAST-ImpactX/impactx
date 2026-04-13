#ifndef SPACE_CHARGE_3D_FD_H_
#define SPACE_CHARGE_3D_FD_H_

#include "synergia/simulation/collective_operator.h"
#include "synergia/simulation/implemented_collective_options.h"

#include "rectangular_grid_domain.h"

#include "space_charge_3d_fd_impl.h"

class Space_charge_3d_fd;

/// Note: internal grid is stored in [z][y][x] order, but
/// grid shape expects [x][y][z] order.
class Space_charge_3d_fd : public Collective_operator {
  private:
    const Space_charge_3d_fd_options
        options; /* Options to initialize sc-3d-fd */
    std::string bunch_sim_id;
    Rectangular_grid_domain domain;
    bool use_fixed_domain;
    bool allocated = false;
    bool initialized = false;
    double scale_x_threshold =
        15; /*! rebuild preconditioner if the size
              of the domain along x changes by this factor */
    double scale_y_threshold =
        15; /*! rebuild preconditioner if the size
              of the domain along y changes by this factor */
    double scale_z_threshold =
        2.5; /*! rebuild preconditioner if the size
              of the domain along z changes by this factor */

    GlobalCtx gctx;
    SubcommCtx sctx;
    LocalCtx lctx;

  private:
    void set_fixed_domain(std::array<double, 3> offset,
                          std::array<double, 3> size);

    void apply_impl(Bunch_simulator& simulator,
                    double time_step,
                    Logger& logger);

    void get_local_charge_density(const Bunch& bunch);

    PetscErrorCode apply_bunch(Bunch& bunch, double time_step, Logger& logger);

    void get_force();

    void apply_kick(Bunch& bunch, double time_step);

    PetscErrorCode allocate_sc3d_fd(const Bunch& bunch);

    PetscErrorCode init_solver_sc3d_fd(Logger& logger);

    PetscErrorCode destroy_sc3d_fd();

    PetscErrorCode update_domain(Bunch const& bunch);

  public:
    Space_charge_3d_fd(Space_charge_3d_fd_options const& ops);

    ~Space_charge_3d_fd();
};

#endif /* SPACE_CHARGE_3D_FD_H_ */
