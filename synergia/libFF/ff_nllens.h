#ifndef FF_NONLINEARLENS_H
#define FF_NONLINEARLENS_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/libFF/ff_patterned_propagator.h"
#include "synergia/utils/simple_timer.h"

namespace FF_nllens {
  template <class BunchT>
  void
  apply(Lattice_element_slice const& slice, BunchT& bunch)
  {
    auto const& elem = slice.get_lattice_element();

    const double l = elem.get_double_attribute("l", 0.0);
    const double knll = elem.get_double_attribute("knll");
    const double cnll = elem.get_double_attribute("cnll");

    if (!close_to_zero(l))
      throw std::runtime_error("nonlinear lens has non-zero length");

    const double icnll = 1.0 / cnll;
    const double kick = -knll * icnll;
    const double k[2] = {icnll, kick};

    using gsv_t = typename std::conditional<
      std::is_floating_point<typename BunchT::part_t>::value,
      GSVec<double>,
      typename BunchT::gsv_t>::type;

    using pp = FF_patterned_propagator<BunchT,
                                       gsv_t,
                                       FF_algorithm::nllens_unit<gsv_t>,
                                       FF_algorithm::nllens_unit<double>>;

    Reference_particle& ref = bunch.get_design_reference_particle();
    pp::get_reference_cdt_zero(ref, k);

    double xp = ref.get_state()[Bunch::xp];
    double yp = ref.get_state()[Bunch::yp];

    if (std::isnan(xp) || std::isnan(yp)) {
      throw std::runtime_error(
        "the bunch reference particle hits the sigularity region "
        "when propagating through a nonlinearlens");
    }

    pp::apply_thin_kick(bunch, ParticleGroup::regular, k);
    pp::apply_thin_kick(bunch, ParticleGroup::spectator, k);
  }
}

#endif // FF_NONLINEARLENS_H
