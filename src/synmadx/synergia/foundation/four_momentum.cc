#include <cmath>
#include <stdexcept>

#include "synergia/utils/floating_point.h"

#include "four_momentum.h"

void
Four_momentum::update_from_gamma()
{
    if (gamma < 1.0) {
        throw std::range_error("Four_momentum: gamma not >= 1.0");
    }
    energy = gamma * mass;
    beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    momentum = gamma * beta * mass;
}

Four_momentum::Four_momentum()
{
    mass = 0.0;
    gamma = 1.0;
}

Four_momentum::Four_momentum(double mass)
{
    this->mass = mass;
    gamma = 1.0;
    update_from_gamma();
}

Four_momentum::Four_momentum(double mass, double total_energy)
{
    this->mass = mass;
    set_total_energy(total_energy);
}

Four_momentum::Four_momentum(Lsexpr const& lsexpr)
{
    bool found_mass(false), found_gamma(false);
    for (Lsexpr::const_iterator_t it = lsexpr.begin(); it != lsexpr.end();
         ++it) {
        if (it->is_labeled()) {
            if (it->get_label() == "mass") {
                mass = it->get_double();
                found_mass = true;
            } else if (it->get_label() == "gamma") {
                gamma = it->get_double();
                found_gamma = true;
            } else {
                throw std::runtime_error(
                    "Four_momentum: Lsexpr unknown label '" + it->get_label() +
                    "'");
            }
        }
    }
    if (!found_mass) {
        throw std::runtime_error("Four_momentum: Lsexpr missing mass");
    }
    if (!found_gamma) {
        throw std::runtime_error("Four_momentum: Lsexpr missing gamma");
    }
    update_from_gamma();
}

Lsexpr
Four_momentum::as_lsexpr() const
{
    Lsexpr retval;
    retval.set_label("Four_momentum");
    retval.push_back(Lsexpr(mass, "mass"));
    retval.push_back(Lsexpr(gamma, "gamma"));
    return retval;
}
void
Four_momentum::set_total_energy(double total_energy)
{
    gamma = total_energy / mass;
    update_from_gamma();
}

void
Four_momentum::set_kinetic_energy(double kinetic_energy)
{
    gamma = (mass + kinetic_energy) / mass;
    update_from_gamma();
}

void
Four_momentum::set_momentum(double momentum)
{
    double r2 = momentum * momentum / (mass * mass);
    beta = sqrt(r2 / (1 + r2));
    set_beta(beta);
}

void
Four_momentum::set_gamma(double gamma)
{
    this->gamma = gamma;
    update_from_gamma();
}

void
Four_momentum::set_beta(double beta)
{
    if ((beta < 0.0) || (beta >= 1.0)) {
        throw std::range_error(
            "Four_momentum: beta not in range 0.0 <= beta < 1.0");
    }
    gamma = 1.0 / sqrt(1.0 - beta * beta);
    update_from_gamma();
}

double
Four_momentum::get_mass() const
{
    return mass;
}

double
Four_momentum::get_total_energy() const
{
    return energy;
}

double
Four_momentum::get_kinetic_energy() const
{
    return energy - mass;
}

double
Four_momentum::get_momentum() const
{
    return momentum;
}

double
Four_momentum::get_gamma() const
{
    return gamma;
}

double
Four_momentum::get_beta() const
{
    return beta;
}

bool
Four_momentum::equal(Four_momentum const& four_momentum, double tolerance) const
{
    if (!floating_point_equal(mass, four_momentum.get_mass(), tolerance)) {
        return false;
    }
    if (!floating_point_equal(gamma, four_momentum.get_gamma(), tolerance)) {
        return false;
    }
    if (!floating_point_equal(beta, four_momentum.get_beta(), tolerance)) {
        return false;
    }
    return true;
}
