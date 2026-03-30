#ifndef SPACE_CHARGE_2D_KV_H_
#define SPACE_CHARGE_2D_KV_H_

#include "synergia/simulation/collective_operator.h"
#include "synergia/simulation/implemented_collective_options.h"

class Space_charge_2d_kv : public Collective_operator {

private:
  const Space_charge_2d_kv_options opts;

  void apply_impl(Bunch_simulator& simulator,
                  double time_step,
                  Logger& logger) override;

  void apply_bunch(Bunch& bunch, double time_step, Logger& logger);

public:
  Space_charge_2d_kv(Space_charge_2d_kv_options const& opts);
};

#endif /* SPACE_CHARGE_2D_KV_H_ */
