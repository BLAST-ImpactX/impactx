#ifndef REFERENCE_PARTICLE_H_
#define REFERENCE_PARTICLE_H_

//#include "synergia/utils/kokkos_views.h"
//#include "synergia/utils/multi_array_serialization.h"
#include "synergia/foundation/four_momentum.h"
#include <cereal/types/array.hpp>

/// Reference_particle stores the four momentum of the reference frame
/// with respect to  the lab frame (defined to be along the axis of the
/// accelerator) as well as the six-dimensional state vector of the the
/// reference particle in the reference frame. Reference particle
/// also keeps track of the total path length of the reference particle
/// trajectory.
class Reference_particle {
private:
  int charge;
  Four_momentum four_momentum;
  std::array<double, 6> state;
  int repetition;
  double s;
  double s_n;
  double abs_time;
  double abs_offset;

public:
  /// Default constructor for internal use only
  Reference_particle();

  /// Construct a Reference_particle with a given mass and total energy.
  /// @param mass in GeV/c^2
  /// @param charge in units of e
  /// @param total_energy in GeV in the lab frame
  Reference_particle(int charge, double mass, double total_energy);

  /// Construct a Reference_particle with a given four momentum.
  /// @param charge in units of e
  /// @param four_momentum in the lab frame
  Reference_particle(int charge, Four_momentum const& four_momentum);

  /// Construct a Reference_particle with a given four momentum and state
  /// in the reference frame.
  /// @param charge in units of e
  /// @param four_momentum in the lab frame
  /// @param state is a six-dimensional state vector
  Reference_particle(int charge,
                     Four_momentum const& four_momentum,
                     std::array<double, 6> const& state);

  /// Construct a Reference_particle from the Lsexpr representation
  /// @param lsexpr representation
  Reference_particle(Lsexpr const& lsexpr);

  /// Extract an Lsexpr representation of the Reference_particle
  Lsexpr as_lsexpr() const;

  /// MadX beam command
  std::string as_madx() const;

  /// Set the four momentum.
  /// @param four_momentum in the lab frame
  void set_four_momentum(Four_momentum const& four_momentum);

  /// Set the state vector in the reference frame.
  /// @param state is a six-dimensional state vector
  void set_state(std::array<double, 6> const& state);

  /// Set the state vector in the reference frame.
  /// @param x
  /// @param xp
  /// @param y
  /// @param yp
  /// @param cdt
  /// @param dpop
  void set_state(double x,
                 double xp,
                 double y,
                 double yp,
                 double cdt,
                 double dpop);

  void set_state_x(double x);

  void set_state_xp(double xp);

  void set_state_y(double y);

  void set_state_yp(double yp);

  void set_state_cdt(double cdt);

  void set_state_dpop(double dpop);

  /// Set the total energy.
  /// @param total_energy in GeV in the lab frame
  void set_total_energy(double total_energy);

  /// Increment the trajectory length.
  /// @param length in m
  void increment_trajectory(double length);

  /// Start a new repetition
  void start_repetition();

  /// Manually set trajectory parameters
  /// @param repetition starting at 0
  /// @param repetition_length in m
  /// @param s in m
  void set_trajectory(int repetition, double repetition_length, double s);

  /// Return the Reference_particle charge in units of e
  int get_charge() const;

  /// Return the Reference_particle mass in units of GeV/c
  double get_mass() const;

  /// Get the four momentum in the lab frame.
  Four_momentum const& get_four_momentum() const;

  /// Get the six-dimensional state vector in the reference frame.
  std::array<double, 6> const& get_state() const;

  /// Get the relativistic beta in the lab frame.
  double get_beta() const;

  /// Get the relativistic gamma in the lab frame.
  double get_gamma() const;

  /// Get the momentum in GeV/c in the lab frame.
  double get_momentum() const;

  /// Get the total energy in GeV in the lab frame.
  double get_total_energy() const;

  /// Get the total path length in m of the reference
  /// particle trajectory
  double get_s() const;

  /// Get the distance traveled in m since the beginning
  /// of the current repetition.
  double get_s_n() const;

  /// Get the number of repetition.
  int get_repetition() const;

  /// Get the repetition length in m.
  double get_repetition_length() const;

  // get the absolute time
  double get_bunch_abs_time() const;

  // get the time offset
  double get_bunch_abs_offset() const;

  void
  set_bunch_abs_time(double const& t);

  void
  increment_bunch_abs_time(double const& dt);

  void
  set_bunch_abs_offset(double const& toffset);

  /// Check equality to the given tolerance
  /// @param reference_particle another Reference_particle
  /// @param tolerance fractional accuracy
  bool equal(Reference_particle const& reference_particle,
             double tolerance) const;

  /// Serialization support
  template <class Archive>
  void
  serialize(Archive& ar)
  {
    ar(CEREAL_NVP(charge));
    ar(CEREAL_NVP(four_momentum));
    ar(CEREAL_NVP(state));
    ar(CEREAL_NVP(repetition));
    ar(CEREAL_NVP(s));
    ar(CEREAL_NVP(s_n));
    ar(CEREAL_NVP(abs_time));
    ar(CEREAL_NVP(abs_offset));
  }
};

#endif /* REFERENCE_PARTICLE_H_ */
