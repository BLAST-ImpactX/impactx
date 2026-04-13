#ifndef FF_ELENS_H
#define FF_ELENS_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/libFF/ff_patterned_propagator.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/foundation/physical_constants.h"

namespace FF_elens {
  template <class BunchT>
  void
  apply(Lattice_element_slice const& slice, BunchT& bunch)
  {
    scoped_simple_timer timer("libFF_elens");

    double length = slice.get_right() - slice.get_left();

    auto const& elem = slice.get_lattice_element();

    // hk and vk are the hk/vk under lattice reference momentum
    double lengh = elem.get_double_attribute("l");
    double current = elem.get_double_attribute("current");
    double eenergy = elem.get_double_attribute("eenergy") * 0.001;
    double radius = elem.get_double_attribute("radius");

    bool gaussian = elem.get_double_attribute("gaussian", 0.0) != 0.0;
    bool uniform = elem.get_double_attribute("uniform", 0.0) != 0.0;

    if (!(uniform || gaussian)) {
      throw std::runtime_error(
        "elens must set either gaussian or uniform attribute");
    }

    if (gaussian && uniform) {
      throw std::runtime_error(
        "elens must not set both gaussian and uniform attributes");
    }

    const double current_over_e = current / pconstants::e;

    // references
    Reference_particle& ref_l = bunch.get_design_reference_particle();
    Reference_particle const& ref_b = bunch.get_reference_particle();

    const double gamma_e = (eenergy + pconstants::me) / pconstants::me;
    const double beta_e = std::sqrt(1.0 - 1.0 / (gamma_e * gamma_e));
    const double gamma_b = ref_b.get_gamma();
    const double beta_b = ref_b.get_beta();

    const double pref_l = ref_l.get_momentum();
    const double pref_b = ref_b.get_momentum();

    const double m_l = ref_l.get_mass();
    const double m_b = ref_b.get_mass();

    const double k[6] = {
      beta_b, gamma_b, beta_e, current_over_e, length, radius};

    using gsv_t = typename std::conditional<
      std::is_floating_point<typename BunchT::part_t>::value,
      GSVec<double>,
      typename BunchT::gsv_t>::type;

    if (gaussian) {
      using pp =
        FF_patterned_propagator<BunchT,
                                gsv_t,
                                FF_algorithm::elens_kick_gaussian<gsv_t>,
                                FF_algorithm::elens_kick_gaussian<double>>;

      if (close_to_zero(length)) {
        pp::get_reference_cdt_zero(ref_l, k);

        pp::apply_thin_kick(bunch, ParticleGroup::regular, k);
        pp::apply_thin_kick(bunch, ParticleGroup::spectator, k);
      } else {
        double ref_cdt = pp::get_reference_cdt_simple(ref_l, length, k);

        pp::apply_simple_kick(
          bunch, ParticleGroup::regular, pref_b, m_b, ref_cdt, length, k);

        pp::apply_simple_kick(
          bunch, ParticleGroup::spectator, pref_b, m_b, ref_cdt, length, k);
      }
    } else {
      using pp =
        FF_patterned_propagator<BunchT,
                                gsv_t,
                                FF_algorithm::elens_kick_uniform<gsv_t>,
                                FF_algorithm::elens_kick_uniform<double>>;

      if (close_to_zero(length)) {
        pp::get_reference_cdt_zero(ref_l, k);

        pp::apply_thin_kick(bunch, ParticleGroup::regular, k);
        pp::apply_thin_kick(bunch, ParticleGroup::spectator, k);
      } else {
        double ref_cdt = pp::get_reference_cdt_simple(ref_l, length, k);

        pp::apply_simple_kick(
          bunch, ParticleGroup::regular, pref_b, m_b, ref_cdt, length, k);

        pp::apply_simple_kick(
          bunch, ParticleGroup::spectator, pref_b, m_b, ref_cdt, length, k);
      }
    }

    // increment the trajectory
    bunch.get_reference_particle().increment_trajectory(length);
    // absolute time
    double const velocity = bunch.get_reference_particle().get_beta()*pconstants::c;
    bunch.get_reference_particle().increment_bunch_abs_time(length/velocity);

  }
}

#endif // FF_HKICKER_H
