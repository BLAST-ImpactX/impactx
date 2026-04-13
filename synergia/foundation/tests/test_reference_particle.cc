#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "catch2/matchers/catch_matchers.hpp"
#include "synergia/foundation/reference_particle.h"
#include "synergia/utils/cereal_files.h"

const double tolerance = 1.0e-13;
const double mass = 100.0;
const double total_energy = 125.0;
const int charge = 1;
const double step_length = 1.234;
const int steps = 17;
const int turns = 7;

TEST_CASE("construct")
{
    REQUIRE_NOTHROW(Reference_particle(charge, mass, total_energy));
}

TEST_CASE("construct2")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    REQUIRE_NOTHROW(Reference_particle(charge, four_momentum));
}

TEST_CASE("construct2_lsexpr")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle original(charge, four_momentum);

    Lsexpr original_as_lsexpr(original.as_lsexpr());
    Reference_particle from_lsexpr(original_as_lsexpr);

    CHECK(from_lsexpr.get_charge() == charge);
    REQUIRE_THAT(from_lsexpr.get_four_momentum().get_total_energy(),
                 Catch::Matchers::WithinAbs(total_energy, tolerance));

    for (int i = 0; i < 6; ++i)
        REQUIRE_THAT(from_lsexpr.get_state()[i],
                     Catch::Matchers::WithinAbs(0.0, tolerance));

    CHECK(from_lsexpr.get_repetition() == 0);

    REQUIRE_THAT(from_lsexpr.get_repetition_length(),
                 Catch::Matchers::WithinAbs(0.0, tolerance));
    REQUIRE_THAT(from_lsexpr.get_s_n(),
                 Catch::Matchers::WithinAbs(0.0, tolerance));
}

TEST_CASE("construct3")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);

    std::array<double, 6> state;
    for (int i = 0; i < 6; ++i)
        state[i] = 1.1 * i;

    REQUIRE_NOTHROW(Reference_particle(charge, four_momentum, state));
}

TEST_CASE("construct3_lsexpr")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);

    std::array<double, 6> state;
    for (int i = 0; i < 6; ++i)
        state[i] = 1.1 * i;

    Reference_particle original(charge, four_momentum, state);

    const double partial_s = 3 * step_length;
    original.set_trajectory(turns, steps * step_length, partial_s);

    Lsexpr original_as_lsexpr(original.as_lsexpr());
    Reference_particle from_lsexpr(original_as_lsexpr);

    CHECK(from_lsexpr.get_charge() == charge);

    REQUIRE_THAT(from_lsexpr.get_four_momentum().get_total_energy(),
                 Catch::Matchers::WithinAbs(total_energy, tolerance));

    for (int i = 0; i < 6; ++i)
        REQUIRE_THAT(
            from_lsexpr.get_state()[i],
            Catch::Matchers::WithinAbs(original.get_state()[i], tolerance));

    CHECK(from_lsexpr.get_repetition() == turns);
    REQUIRE_THAT(from_lsexpr.get_repetition_length(),
                 Catch::Matchers::WithinAbs(steps * step_length, tolerance));

    REQUIRE_THAT(from_lsexpr.get_s_n(),
                 Catch::Matchers::WithinAbs(partial_s, tolerance));
}

TEST_CASE("get_charge")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    CHECK(reference_particle.get_charge() == charge);
}

TEST_CASE("get_four_momentum")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    REQUIRE_THAT(reference_particle.get_four_momentum().get_total_energy(),
                 Catch::Matchers::WithinAbs(four_momentum.get_total_energy(),
                                            tolerance));
}

TEST_CASE("get_state")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, total_energy);

    for (int i = 0; i < 6; ++i)
        REQUIRE_THAT(reference_particle.get_state()[i],
                     Catch::Matchers::WithinAbs(0.0, tolerance));
}

TEST_CASE("get_beta")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    REQUIRE_THAT(
        reference_particle.get_beta(),
        Catch::Matchers::WithinAbs(four_momentum.get_beta(), tolerance));
}

TEST_CASE("get_gamma")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    REQUIRE_THAT(
        reference_particle.get_gamma(),
        Catch::Matchers::WithinAbs(four_momentum.get_gamma(), tolerance));
}

TEST_CASE("get_momentum")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    REQUIRE_THAT(
        reference_particle.get_momentum(),
        Catch::Matchers::WithinAbs(four_momentum.get_momentum(), tolerance));
}

TEST_CASE("get_total_energy")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    REQUIRE_THAT(reference_particle.get_total_energy(),
                 Catch::Matchers::WithinAbs(four_momentum.get_total_energy(),
                                            tolerance));
}

TEST_CASE("increment_trajectory")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);

    for (int step = 0; step < steps; ++step)
        reference_particle.increment_trajectory(step_length);

    REQUIRE_THAT(reference_particle.get_s_n(),
                 Catch::Matchers::WithinAbs(steps * step_length, tolerance));

    REQUIRE_THAT(reference_particle.get_s(),
                 Catch::Matchers::WithinAbs(steps * step_length, tolerance));
}

TEST_CASE("start_repetition")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);

    Reference_particle reference_particle(charge, four_momentum);
    reference_particle.start_repetition();

    CHECK(reference_particle.get_repetition() == 0);
    CHECK(reference_particle.get_s_n() == 0);

    for (int step = 0; step < steps; ++step) {
        reference_particle.increment_trajectory(step_length);
    }

    reference_particle.start_repetition();

    CHECK(reference_particle.get_repetition() == 1);
    CHECK(reference_particle.get_s_n() == 0.0);
    REQUIRE_THAT(
        reference_particle.get_repetition_length(),
        Catch::Matchers::WithinAbsMatcher(steps * step_length, tolerance));
}

TEST_CASE("set_trajectory")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);

    const double partial_s = 3 * step_length;
    reference_particle.set_trajectory(turns, steps * step_length, partial_s);

    CHECK(reference_particle.get_repetition() == turns);
    REQUIRE_THAT(
        reference_particle.get_repetition_length(),
        Catch::Matchers::WithinAbsMatcher(steps * step_length, tolerance));
    REQUIRE_THAT(reference_particle.get_s_n(),
                 Catch::Matchers::WithinAbsMatcher(partial_s, tolerance));
}

TEST_CASE("set_four_momentum")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);

    double new_total_energy = total_energy * 1.1;
    Four_momentum new_four_momentum(mass);
    new_four_momentum.set_total_energy(new_total_energy);
    reference_particle.set_four_momentum(new_four_momentum);

    REQUIRE_THAT(
        reference_particle.get_four_momentum().get_total_energy(),
        Catch::Matchers::WithinAbsMatcher(new_total_energy, tolerance));
}

TEST_CASE("set_state")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);

    std::array<double, 6> new_state;
    for (int i = 0; i < 6; ++i)
        new_state[i] = 1.1 * i;

    reference_particle.set_state(new_state);

    for (int i = 0; i < 6; ++i) {
        REQUIRE_THAT(
            reference_particle.get_state()[i],
            Catch::Matchers::WithinAbsMatcher(new_state[i], tolerance));
    }
}

void
kick_reference_particle_state(Reference_particle& rp)
{
    auto new_state = rp.get_state();
    new_state[1] += 1.0e-2;
    new_state[3] += 1.0e-2;

    rp.set_state(new_state);
}

TEST_CASE("set_state_in_function")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);

    std::array<double, 6> start_state;
    for (int i = 0; i < 6; ++i)
        start_state[i] = 1.1 * i;

    reference_particle.set_state(start_state);
    kick_reference_particle_state(reference_particle);

    for (int i = 0; i < 6; ++i) {
        if ((i != 1) && (i != 3)) {
            REQUIRE_THAT(
                reference_particle.get_state()[i],
                Catch::Matchers::WithinAbsMatcher(start_state[i], tolerance));
        } else {
            REQUIRE_THAT(reference_particle.get_state()[i],
                         Catch::Matchers::WithinAbsMatcher(
                             start_state[i] + 1.0e-2, tolerance));
        }
    }
}

TEST_CASE("set_total_energy")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);

    double new_total_energy = total_energy * 1.1;
    reference_particle.set_total_energy(new_total_energy);

    REQUIRE_THAT(
        reference_particle.get_four_momentum().get_total_energy(),
        Catch::Matchers::WithinAbsMatcher(new_total_energy, tolerance));
}

TEST_CASE("set_abs_time")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    reference_particle.set_bunch_abs_time(1 / 4096.0); // Exactly representable
    CHECK(reference_particle.get_bunch_abs_time() == (1 / 4096.0));
}

TEST_CASE("incr_abs_time")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    reference_particle.set_bunch_abs_time(1 / 4096.0); // Exactly representable
    reference_particle.increment_bunch_abs_time(
        1 / 32768.0); // exactly representable
    CHECK(reference_particle.get_bunch_abs_time() ==
          (9.0 / 32768.0)); // 1/4096 + 1/32768 = 9/32768
}

TEST_CASE("set_abs_offset")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    reference_particle.set_bunch_abs_offset(1 /
                                            32768.0); // Exactly representable
    CHECK(reference_particle.get_bunch_abs_offset() == (1 / 32768.0));
}

TEST_CASE("copy")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);

    std::array<double, 6> state;
    for (int i = 0; i < 6; ++i)
        state[i] = 1.1 * i;

    Reference_particle original_reference_particle(
        charge, four_momentum, state);
    original_reference_particle.set_bunch_abs_time(
        1 / 2048.0); // exactly representable
    original_reference_particle.set_bunch_abs_offset(
        1 / 32768.0); // exactly representable
    std::cout << "original reference particle abs time: "
              << original_reference_particle.get_bunch_abs_time() << std::endl;
    std::cout << "original reference particle abs offset: "
              << original_reference_particle.get_bunch_abs_offset()
              << std::endl;
    Reference_particle reference_particle(original_reference_particle);
    std::cout << "copied reference particle abs time: "
              << reference_particle.get_bunch_abs_time() << std::endl;
    std::cout << "copied reference particle abs offset: "
              << reference_particle.get_bunch_abs_offset() << std::endl;

    REQUIRE_THAT(reference_particle.get_four_momentum().get_total_energy(),
                 Catch::Matchers::WithinAbsMatcher(total_energy, tolerance));

    for (int i = 0; i < 6; ++i) {
        REQUIRE_THAT(reference_particle.get_state()[i],
                     Catch::Matchers::WithinAbsMatcher(state[i], tolerance));
    }
    CHECK(reference_particle.get_bunch_abs_time() == (1 / 2048.0));
    CHECK(reference_particle.get_bunch_abs_offset() == (1 / 32768.0));
}

TEST_CASE("copy2")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle original_reference_particle(charge, four_momentum);
    original_reference_particle.set_bunch_abs_time(
        1.0 / 1024.0); // exactly representable
    Reference_particle reference_particle(original_reference_particle);

    double new_total_energy = total_energy * 1.1;
    Four_momentum new_four_momentum(reference_particle.get_four_momentum());
    new_four_momentum.set_total_energy(new_total_energy);
    original_reference_particle.set_four_momentum(new_four_momentum);

    std::array<double, 6> new_state;
    for (int i = 0; i < 6; ++i)
        new_state[i] = 1.1 * i;

    original_reference_particle.set_state(new_state);

    REQUIRE_THAT(
        original_reference_particle.get_four_momentum().get_total_energy(),
        Catch::Matchers::WithinAbsMatcher(new_total_energy, tolerance));

    REQUIRE_THAT(reference_particle.get_four_momentum().get_total_energy(),
                 Catch::Matchers::WithinAbsMatcher(total_energy, tolerance));

    for (int i = 0; i < 6; ++i) {
        REQUIRE_THAT(
            original_reference_particle.get_state()[i],
            Catch::Matchers::WithinAbsMatcher(new_state[i], tolerance));
        REQUIRE_THAT(reference_particle.get_state()[i],
                     Catch::Matchers::WithinAbsMatcher(0.0, tolerance));
    }

    std::cout << "reference_particle abs_time before: "
              << reference_particle.get_bunch_abs_time() << std::endl;
    reference_particle.increment_bunch_abs_time(1 / 32768.0);
    std::cout << "reference_particle abs_time after: "
              << reference_particle.get_bunch_abs_time() << std::endl;
    CHECK(reference_particle.get_bunch_abs_time() ==
          (33.0 / 32768.0)); // exactly representable
}

TEST_CASE("get_s")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);

    for (int turn = 0; turn < turns; ++turn) {
        reference_particle.start_repetition();

        for (int step = 0; step < steps; ++step) {
            reference_particle.increment_trajectory(step_length);
        }
    }

    REQUIRE_THAT(reference_particle.get_s(),
                 Catch::Matchers::WithinAbsMatcher(turns * steps * step_length,
                                                   tolerance));
}

const double equal_tolerance = 1.0e-12;

TEST_CASE("equal")
{
    Four_momentum four_momentum(mass);
    Reference_particle reference_particle1(charge, four_momentum);

    std::array<double, 6> new_state;
    for (int i = 0; i < 6; ++i)
        new_state[i] = 1.1 * i;

    reference_particle1.set_state(new_state);
    Reference_particle reference_particle2(charge, four_momentum);
    reference_particle2.set_state(new_state);

    CHECK(reference_particle1.equal(reference_particle2, equal_tolerance));
}

TEST_CASE("equal_different_four_momentum")
{
    Four_momentum four_momentum1(mass, total_energy);
    Reference_particle reference_particle1(charge, four_momentum1);

    Four_momentum four_momentum2(mass, total_energy * 1.1);
    Reference_particle reference_particle2(charge, four_momentum2);

    CHECK(!reference_particle1.equal(reference_particle2, equal_tolerance));
}

TEST_CASE("equal_different_state")
{
    Four_momentum four_momentum(mass);
    Reference_particle reference_particle1(charge, four_momentum);

    std::array<double, 6> state1;
    for (int i = 0; i < 6; ++i)
        state1[i] = 1.1 * i;

    reference_particle1.set_state(state1);

    Reference_particle reference_particle2(charge, four_momentum);

    std::array<double, 6> state2;
    for (int i = 0; i < 6; ++i)
        state2[i] = state1[i];
    state2[1] *= 1.1;

    reference_particle2.set_state(state2);

    CHECK(!reference_particle1.equal(reference_particle2, equal_tolerance));
}

TEST_CASE("equal_different_charge")
{
    Four_momentum four_momentum(mass);
    Reference_particle reference_particle1(charge, four_momentum);

    std::array<double, 6> new_state;
    for (int i = 0; i < 6; ++i)
        new_state[i] = 1.1 * i;

    int charge2 = -1;
    reference_particle1.set_state(new_state);

    Reference_particle reference_particle2(charge2, four_momentum);
    reference_particle2.set_state(new_state);

    CHECK(!reference_particle1.equal(reference_particle2, equal_tolerance));
}

TEST_CASE("equal_same_abs_time")
{
    Four_momentum four_momentum(mass);
    Reference_particle reference_particle1(charge, four_momentum);

    std::array<double, 6> new_state;
    for (int i = 0; i < 6; ++i)
        new_state[i] = 1.1 * i;
    reference_particle1.set_bunch_abs_time(1.0 / 1024.0);

    Reference_particle reference_particle2(reference_particle1);

    CHECK(reference_particle1.equal(reference_particle2, equal_tolerance));
}

TEST_CASE("equal_different_abs_time")
{
    Four_momentum four_momentum(mass);
    Reference_particle reference_particle1(charge, four_momentum);

    std::array<double, 6> new_state;
    for (int i = 0; i < 6; ++i)
        new_state[i] = 1.1 * i;
    reference_particle1.set_bunch_abs_time(1.0 / 1024.0);

    Reference_particle reference_particle2(reference_particle1);
    reference_particle2.set_bunch_abs_time(3 / 2048.0);

    CHECK(!reference_particle1.equal(reference_particle2, equal_tolerance));
}

TEST_CASE("test_serialize_xml")
{
    Reference_particle reference_particle(charge, mass, total_energy);
    xml_save<Reference_particle>(reference_particle, "reference_particle.xml");

    Reference_particle loaded;
    xml_load<Reference_particle>(loaded, "reference_particle.xml");
    CHECK(reference_particle.equal(loaded, tolerance));
}

TEST_CASE("test_serialize_json")
{
    Reference_particle reference_particle(charge, mass, total_energy);
    json_save<Reference_particle>(reference_particle,
                                  "reference_particle.json");

    Reference_particle loaded;
    json_load<Reference_particle>(loaded, "reference_particle.json");
    CHECK(reference_particle.equal(loaded, tolerance));
}
