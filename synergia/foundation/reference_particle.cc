#include "synergia/utils/floating_point.h"

#include "physical_constants.h"
#include "reference_particle.h"

Reference_particle::Reference_particle()
    : charge(0)
    , four_momentum(1.0, 1.0)
    , state{0, 0, 0, 0, 0, 0}
    , repetition(0)
    , s(0)
    , s_n(0)
    , abs_time(0)
    , abs_offset(0)
{}

Reference_particle::Reference_particle(int charge,
                                       double mass,
                                       double total_energy)
    : charge(charge)
    , four_momentum(mass, total_energy)
    , state{0, 0, 0, 0, 0, 0}
    , repetition(0)
    , s(0)
    , s_n(0)
    , abs_time(0)
    , abs_offset(0)
{}

Reference_particle::Reference_particle(int charge,
                                       Four_momentum const& four_momentum_in)
    : charge(charge)
    , four_momentum(four_momentum_in)
    , state{0, 0, 0, 0, 0, 0}
    , repetition(0)
    , s(0)
    , s_n(0)
    , abs_time(0)
    , abs_offset(0)
{}

Reference_particle::Reference_particle(int charge,
                                       Four_momentum const& four_momentum_in,
                                       std::array<double, 6> const& state)
    : charge(charge)
    , four_momentum(four_momentum_in)
    , state(state)
    , repetition(0)
    , s(0)
    , s_n(0)
    , abs_time(0)
    , abs_offset(0)
{}

Reference_particle::Reference_particle(Lsexpr const& lsexpr)
    : state{0, 0, 0, 0, 0, 0}, repetition(0), s(0), s_n(0)
{
    bool found_charge(false), found_four_momentum(false);
    for (Lsexpr::const_iterator_t it = lsexpr.begin(); it != lsexpr.end();
         ++it) {
        if (it->is_labeled()) {
            if (it->get_label() == "charge") {
                charge = it->get_int();
                found_charge = true;
            } else if (it->get_label() == "four_momentum") {
                four_momentum = Four_momentum(*it);
                found_four_momentum = true;
            } else if (it->get_label() == "repetition") {
                repetition = it->get_int();
            } else if (it->get_label() == "s") {
                s = it->get_double();
            } else if (it->get_label() == "s_n") {
                s_n = it->get_double();
            } else if (it->get_label() == "state") {
                std::vector<double> state_vector(it->get_double_vector());
                if (state_vector.size() != 6) {
                    throw std::runtime_error("Reference_particle from lsexpr: "
                                             "state vector must have length 6");
                }
                std::copy(
                    state_vector.begin(), state_vector.end(), state.begin());
            } else {
                throw std::runtime_error(
                    "Reference_particle: Lsexpr unknown label \"" +
                    it->get_label() + "\"");
            }
        }
    }
    if (!found_charge) {
        throw std::runtime_error("Reference_particle: Lsexpr missing charge");
    }
    if (!found_four_momentum) {
        throw std::runtime_error(
            "Reference_particle: Lsexpr missing four_momentum");
    }
}

Lsexpr
Reference_particle::as_lsexpr() const
{
    Lsexpr retval;

    Lsexpr four_momentum_lsexpr(four_momentum.as_lsexpr());
    four_momentum_lsexpr.set_label("four_momentum");
    retval.push_back(four_momentum_lsexpr);

    retval.push_back(Lsexpr(charge, "charge"));

    if (repetition != 0) { retval.push_back(Lsexpr(repetition, "repetition")); }

    if (s != 0) { retval.push_back(Lsexpr(s, "s")); }

    if (s_n != 0) { retval.push_back(Lsexpr(s_n, "s_n")); }

    if ((state[0] != 0) || (state[1] != 0) || (state[2] != 0) ||
        (state[3] != 0) || (state[4] != 0) || (state[5] != 0)) {
        std::vector<double> state_vector(6);
        std::copy(state.begin(), state.end(), state_vector.begin());
        retval.push_back(Lsexpr(state_vector, "state"));
    }

    return retval;
}

std::string
Reference_particle::as_madx() const
{
    std::stringstream ss;
    ss << "beam, pc=" << get_momentum();

    // particle type
    double mass = get_mass();
    int charge = get_charge();

    if (mass == pconstants::mp && charge == pconstants::proton_charge)
        ss << ", particle=proton";
    else if (mass == pconstants::mp && charge == pconstants::antiproton_charge)
        ss << ", particle=antiproton";
    else if (mass == pconstants::me && charge == pconstants::electron_charge)
        ss << ", particle=electron";
    else if (mass == pconstants::me && charge == pconstants::positron_charge)
        ss << ", particle=positron";
    else if (mass == pconstants::mmu && charge == pconstants::muon_charge)
        ss << ", particle=negmuon";
    else if (mass == pconstants::mmu && charge == pconstants::antimuon_charge)
        ss << ", particle=posmuon";
    else
        ss << ", mass=" << mass << ", charge=" << charge;

    ss << ";";
    return ss.str();
}

void
Reference_particle::set_four_momentum(Four_momentum const& four_momentum)
{
    this->four_momentum = four_momentum;
}

void
Reference_particle::set_state(std::array<double, 6> const& state)
{
    this->state = state;
}

void
Reference_particle::set_state(double x,
                              double xp,
                              double y,
                              double yp,
                              double cdt,
                              double dpop)
{
    this->state[0] = x;
    this->state[1] = xp;
    this->state[2] = y;
    this->state[3] = yp;
    this->state[4] = cdt;
    this->state[5] = dpop;
}

void
Reference_particle::set_state_x(double x)
{
    this->state[0] = x;
}

void
Reference_particle::set_state_xp(double xp)
{
    this->state[1] = xp;
}

void
Reference_particle::set_state_y(double y)
{
    this->state[2] = y;
}

void
Reference_particle::set_state_yp(double yp)
{
    this->state[3] = yp;
}

void
Reference_particle::set_state_cdt(double cdt)
{
    this->state[4] = cdt;
}

void
Reference_particle::set_state_dpop(double dpop)
{
    this->state[5] = dpop;
}

void
Reference_particle::set_total_energy(double total_energy)
{
    four_momentum.set_total_energy(total_energy);
}

void
Reference_particle::increment_trajectory(double length)
{
    s_n += length;
}

void
Reference_particle::start_repetition()
{
    if (s == 0.0) { s = s_n; }
    if (s_n > 0.0) { repetition += 1; }
    s_n = 0.0;
}

void
Reference_particle::set_trajectory(int repetition,
                                   double repetition_length,
                                   double s)
{
    this->repetition = repetition;
    this->s = repetition_length;
    this->s_n = s;
}

int
Reference_particle::get_charge() const
{
    return charge;
}

double
Reference_particle::get_mass() const
{
    return four_momentum.get_mass();
}

Four_momentum const&
Reference_particle::get_four_momentum() const
{
    return four_momentum;
}

std::array<double, 6> const&
Reference_particle::get_state() const
{
    return state;
}

double
Reference_particle::get_gamma() const
{
    return four_momentum.get_gamma();
}

double
Reference_particle::get_beta() const
{
    return four_momentum.get_beta();
}

double
Reference_particle::get_momentum() const
{
    return four_momentum.get_momentum();
}

double
Reference_particle::get_total_energy() const
{
    return four_momentum.get_total_energy();
}

double
Reference_particle::get_s() const
{
    return repetition * s + s_n;
}

double
Reference_particle::get_s_n() const
{
    return s_n;
}

int
Reference_particle::get_repetition() const
{
    return repetition;
}

double
Reference_particle::get_repetition_length() const
{
    return s;
}

double
Reference_particle::get_bunch_abs_time() const
{
    return abs_time;
}

double
Reference_particle::get_bunch_abs_offset() const
{
    return abs_offset;
}

void
Reference_particle::set_bunch_abs_time(double const& t)
{
    abs_time = t;
}

void
Reference_particle::increment_bunch_abs_time(double const& incr)
{
    abs_time += incr;
}

void
Reference_particle::set_bunch_abs_offset(double const& offset)
{
    abs_offset = offset;
}

bool
Reference_particle::equal(Reference_particle const& reference_particle,
                          double tolerance) const
{
    if (charge != reference_particle.get_charge()) { return false; }
    if (!four_momentum.equal(reference_particle.get_four_momentum(),
                             tolerance)) {
        return false;
    }
    for (int i = 0; i < 6; ++i) {
        if (!floating_point_equal(
                state[i], reference_particle.get_state()[i], tolerance)) {
            return false;
        }
    }
    if (!floating_point_equal(
            abs_time, reference_particle.get_bunch_abs_time(), tolerance) ||
        !floating_point_equal(
            abs_offset, reference_particle.get_bunch_abs_offset(), tolerance)) {
        return false;
    }

    return true;
}
