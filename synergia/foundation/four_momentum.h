#ifndef FOUR_MOMENTUM_H_
#define FOUR_MOMENTUM_H_

#include "synergia/utils/cereal.h"
#include "synergia/utils/lsexpr.h"

/// Four_momentum provides conversion between various relativistic kinematic
/// parameters.
class Four_momentum {
  private:
    double mass, energy, momentum, gamma, beta;
    void update_from_gamma();

  public:
    /// Default constructor for internal use only
    Four_momentum();

    /// Construct a Four_momentum in the rest frame
    /// @param mass in GeV/c^2
    Four_momentum(double mass);

    /// Construct a Four_momentum with the given total energy
    /// @param mass in GeV/c^2
    /// @param total_energy in GeV
    Four_momentum(double mass, double total_energy);

    /// Construct a Four_momentum from the Lsexpr representation
    /// @param lsexpr representation
    Four_momentum(Lsexpr const& lsexpr);

    /// Extract an Lsexpr representation of the Four_momentum
    Lsexpr as_lsexpr() const;

    /// Set the total energy
    /// @param total_energy in GeV
    void set_total_energy(double total_energy);

    /// Set the kinetic energy
    /// @param kinetic_energy in GeV
    void set_kinetic_energy(double kinetic_energy);

    /// Set the momentum
    /// @param momentum in GeV/c
    void set_momentum(double momentum);

    /// Set the relativistic gamma factor
    /// @param gamma unitless
    void set_gamma(double gamma);

    /// Set the relativistic beta factor
    /// @param beta unitless
    void set_beta(double beta);

    /// Get the mass in GeV/c^2
    double get_mass() const;

    /// Get the total energy in GeV
    double get_total_energy() const;

    /// Get the kinetic energy in GeV
    double get_kinetic_energy() const;

    /// Get momentum in GeV/c
    double get_momentum() const;

    /// Get the relativistic gamma factor
    double get_gamma() const;

    /// Get the relativistic beta factor
    double get_beta() const;

    /// Check equality to the given tolerance
    /// @param four_momentum another Four_momentum
    /// @param tolerance fractional accuracy for beta and gamma
    bool equal(Four_momentum const& four_momentum, double tolerance) const;

    /// Serialization support
    template <class Archive>
    void
    serialize(Archive& ar)
    {
        ar(CEREAL_NVP(mass));
        ar(CEREAL_NVP(energy));
        ar(CEREAL_NVP(momentum));
        ar(CEREAL_NVP(gamma));
        ar(CEREAL_NVP(beta));
    }
};

#endif /* FOUR_MOMENTUM_H_ */
