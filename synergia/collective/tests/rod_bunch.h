#include "synergia/collective/space_charge_3d_open_hockney.h"
#include "synergia/foundation/physical_constants.h"

const int charge = pconstants::proton_charge;
const double mass = pconstants::mp;
const double real_num = 1.7e11;
const int total_num = 10000;
const double total_energy = 125.0;

// 1 + (odd number) * 8 + 2 to expand the domein on both sides
const int rod_num_particles = 1 + 250001 * 8;

const double rod_lowgamma = 61.0 / 60.0;
const double rod_highgamma = 61.0 / 11.0;
const double rod_real_num = 5.0e9;
const double rod_length = 0.1;

// radius of rod particles
const double rod_radius = 1.0e-6;

// radius of the probe particle
const double rod_probe = 1.0e-3;

struct Rod_bunch_fixture_lowgamma {
  Rod_bunch_fixture_lowgamma()
    : bsim(Bunch_simulator::create_single_bunch_simulator(
        Reference_particle(charge, Four_momentum(mass, mass * rod_lowgamma)),
        rod_num_particles,
        rod_real_num))
  {
    auto& bunch = bsim.get_bunch();

    // bunch.set_longitudinal_boundary(LongitudinalBoundary::periodic,
    // rod_length);
    auto local_particles = bunch.get_host_particles();

    // a ring of 8 particles around each longitudinal location
    int num_longitudinal = (rod_num_particles - 1) / 8;
    double dz = rod_length / (num_longitudinal - 1);
    double r2o2 = std::sqrt(2.0) / 2.0;
    double z = -rod_length / 2.0;

    auto const& ref = bunch.get_reference_particle();
    double rod_beta = ref.get_beta();

    for (int i = 1; i < rod_num_particles; i += 8, z += dz) {
      local_particles(i, Bunch::x) = rod_radius;
      local_particles(i, Bunch::y) = 0.0;

      local_particles(i + 1, Bunch::x) = rod_radius * r2o2;
      local_particles(i + 1, Bunch::y) = rod_radius * r2o2;

      local_particles(i + 2, Bunch::x) = 0.0;
      local_particles(i + 2, Bunch::y) = rod_radius;

      local_particles(i + 3, Bunch::x) = -rod_radius * r2o2;
      local_particles(i + 3, Bunch::y) = rod_radius * r2o2;

      local_particles(i + 4, Bunch::x) = -rod_radius;
      local_particles(i + 4, Bunch::y) = 0.0;

      local_particles(i + 5, Bunch::x) = -rod_radius * r2o2;
      local_particles(i + 5, Bunch::y) = -rod_radius * r2o2;

      local_particles(i + 6, Bunch::x) = 0.0;
      local_particles(i + 6, Bunch::y) = -rod_radius;

      local_particles(i + 7, Bunch::x) = rod_radius * r2o2;
      local_particles(i + 7, Bunch::y) = -rod_radius * r2o2;

      for (int j = i; j < i + 8; ++j) {
        local_particles(j, Bunch::cdt) = z / rod_beta;
        local_particles(j, Bunch::xp) = 0.0;
        local_particles(j, Bunch::yp) = 0.0;
        local_particles(j, Bunch::dpop) = 0.0;
        local_particles(j, Bunch::id) = j;
      }
    }

    // when the probe is too far, it falls outside the domain and does
    // not get any sc kicks.
    local_particles(0, Bunch::x) = 80 * rod_radius;
    local_particles(0, Bunch::y) = 0.0;
    local_particles(0, Bunch::xp) = 0.0;
    local_particles(0, Bunch::yp) = 0.0;
    local_particles(0, Bunch::cdt) = 0.0;
    local_particles(0, Bunch::dpop) = 0.0;
    local_particles(0, Bunch::id) = 0.0;

    // check in particles
    bunch.checkin_particles();
  }

  Bunch_simulator bsim;
};

struct Rod_bunch_fixture_lowgamma_offaxis: public Rod_bunch_fixture_lowgamma {
 Rod_bunch_fixture_lowgamma_offaxis() : Rod_bunch_fixture_lowgamma()
  {
    auto& bunch = bsim.get_bunch();
    bunch.checkout_particles();
    auto local_particles = bunch.get_host_particles();
    // Rotate particle 0 30 degrees
    double s30 = 0.5; // sin 30degrees
    double c30 = std::sqrt(3.0)/2.0; // cos 30 degrees
    double lpx = local_particles(0, 0);
    double lpy = local_particles(0, 2);
    local_particles(0, 0) = c30*lpx - s30*lpy;
    local_particles(0, 2) = s30*lpx + c30*lpy;
    bunch.checkout_particles();
  }
};

struct Rod_bunch_fixture_highgamma {
  Rod_bunch_fixture_highgamma()
    : bsim(Bunch_simulator::create_single_bunch_simulator(
        Reference_particle(charge, Four_momentum(mass, mass * rod_highgamma)),
        rod_num_particles,
        rod_real_num))
  {
    auto& bunch = bsim.get_bunch();

    // bunch.set_longitudinal_boundary(LongitudinalBoundary::periodic,
    // rod_length);
    auto local_particles = bunch.get_host_particles();

    // a ring of 8 particles around each longitudinal location
    int num_longitudinal = (rod_num_particles - 1) / 8;
    double dz = rod_length / (num_longitudinal - 1);
    double r2o2 = std::sqrt(2.0) / 2.0;
    double z = -rod_length / 2.0;

    auto const& ref = bunch.get_reference_particle();
    double rod_beta = ref.get_beta();

    for (int i = 1; i < rod_num_particles; i += 8, z += dz) {
      local_particles(i, Bunch::x) = rod_radius;
      local_particles(i, Bunch::y) = 0.0;

      local_particles(i + 1, Bunch::x) = rod_radius * r2o2;
      local_particles(i + 1, Bunch::y) = rod_radius * r2o2;

      local_particles(i + 2, Bunch::x) = 0.0;
      local_particles(i + 2, Bunch::y) = rod_radius;

      local_particles(i + 3, Bunch::x) = -rod_radius * r2o2;
      local_particles(i + 3, Bunch::y) = rod_radius * r2o2;

      local_particles(i + 4, Bunch::x) = -rod_radius;
      local_particles(i + 4, Bunch::y) = 0.0;

      local_particles(i + 5, Bunch::x) = -rod_radius * r2o2;
      local_particles(i + 5, Bunch::y) = -rod_radius * r2o2;

      local_particles(i + 6, Bunch::x) = 0.0;
      local_particles(i + 6, Bunch::y) = -rod_radius;

      local_particles(i + 7, Bunch::x) = rod_radius * r2o2;
      local_particles(i + 7, Bunch::y) = -rod_radius * r2o2;

      for (int j = i; j < i + 8; ++j) {
        local_particles(j, Bunch::cdt) = z / rod_beta;
        local_particles(j, Bunch::xp) = 0.0;
        local_particles(j, Bunch::yp) = 0.0;
        local_particles(j, Bunch::dpop) = 0.0;
        local_particles(j, Bunch::id) = j;
      }
    }

    // when the probe is too far, it falls outside the domain and does
    // not get any sc kicks.
    local_particles(0, Bunch::x) = 80 * rod_radius;
    local_particles(0, Bunch::y) = 0.0;
    local_particles(0, Bunch::xp) = 0.0;
    local_particles(0, Bunch::yp) = 0.0;
    local_particles(0, Bunch::cdt) = 0.0;
    local_particles(0, Bunch::dpop) = 0.0;
    local_particles(0, Bunch::id) = 0.0;

    // check in particles
    bunch.checkin_particles();
  }

  Bunch_simulator bsim;
};

struct Rod_bunch_fixture_highgamma_offaxis: public Rod_bunch_fixture_highgamma {
 Rod_bunch_fixture_highgamma_offaxis() : Rod_bunch_fixture_highgamma()
  {
    auto& bunch = bsim.get_bunch();
    bunch.checkout_particles();
    auto local_particles = bunch.get_host_particles();
    // Rotate particle 0 30 degrees
    double s30 = 0.5; // sin 30degrees
    double c30 = std::sqrt(3.0)/2.0; // cos 30 degrees
    double lpx = local_particles(0, 0);
    double lpy = local_particles(0, 2);
    local_particles(0, 0) = c30*lpx - s30*lpy;
    local_particles(0, 2) = s30*lpx + c30*lpy;
    bunch.checkout_particles();
  }
};
