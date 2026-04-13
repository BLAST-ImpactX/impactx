
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "synergia/foundation/four_momentum.h"

#include "synergia/utils/cereal_files.h"

const double tolerance = 1.0e-15;

// some kinematics with simple decimal representations
const double mass = 100.0;
const double my_gamma = 1.25; // 5/4
const double beta = 0.6;      // = sqrt(1-1/my_gamma^2) = 3/5
const double total_energy = 125.0;

TEST_CASE("construct")
{
    REQUIRE_NOTHROW(Four_momentum(mass));
}

TEST_CASE("construct2")
{
    REQUIRE_NOTHROW(Four_momentum(mass, total_energy));
}

TEST_CASE("lsexpr")
{
    Four_momentum original(mass, total_energy);
    Lsexpr original_as_lsexpr(original.as_lsexpr());
    Four_momentum from_lsexpr(original_as_lsexpr);

    REQUIRE_THAT(from_lsexpr.get_mass(),
                 Catch::Matchers::WithinAbs(mass, tolerance));
    REQUIRE_THAT(from_lsexpr.get_total_energy(),
                 Catch::Matchers::WithinAbs(total_energy, tolerance));
}

TEST_CASE("get_mass")
{
    Four_momentum four_momentum(mass);
    REQUIRE_THAT(four_momentum.get_mass(),
                 Catch::Matchers::WithinAbs(mass, tolerance));
}

TEST_CASE("get_total_energy")
{
    Four_momentum four_momentum(mass);
    REQUIRE_THAT(four_momentum.get_total_energy(),
                 Catch::Matchers::WithinAbs(mass, tolerance));
}

TEST_CASE("set_and_get_total_energy")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    REQUIRE_THAT(four_momentum.get_total_energy(),
                 Catch::Matchers::WithinAbs(total_energy, tolerance));
}

TEST_CASE("get_kinetic_energy")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);

    REQUIRE_THAT(four_momentum.get_kinetic_energy(),
                 Catch::Matchers::WithinAbs(total_energy - mass, tolerance));
}

TEST_CASE("get_momentum")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);

    REQUIRE_THAT(four_momentum.get_momentum(),
                 Catch::Matchers::WithinAbs(my_gamma * beta * mass, tolerance));
}

TEST_CASE("get_gamma")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    REQUIRE_THAT(four_momentum.get_gamma(),
                 Catch::Matchers::WithinAbs(my_gamma, tolerance));
}

TEST_CASE("get_beta")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);

    REQUIRE_THAT(four_momentum.get_beta(),
                 Catch::Matchers::WithinAbs(beta, tolerance));
}

TEST_CASE("set_kinetic_energy")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_kinetic_energy(total_energy - mass);
    REQUIRE_THAT(four_momentum.get_total_energy(),
                 Catch::Matchers::WithinAbs(total_energy, tolerance));
}

TEST_CASE("set_momentum")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_momentum(my_gamma * beta * mass);
    REQUIRE_THAT(four_momentum.get_total_energy(),
                 Catch::Matchers::WithinAbs(total_energy, tolerance));
}

TEST_CASE("set_gamma")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_gamma(my_gamma);
    REQUIRE_THAT(four_momentum.get_total_energy(),
                 Catch::Matchers::WithinAbs(total_energy, tolerance));
}

TEST_CASE("set_gamma_invalid")
{
    Four_momentum four_momentum(mass);
    const double too_small = 0.8;
    CHECK_THROWS(four_momentum.set_gamma(too_small));
}

TEST_CASE("set_beta")
{
    Four_momentum four_momentum(mass);
    four_momentum.set_beta(beta);
    REQUIRE_THAT(four_momentum.get_total_energy(),
                 Catch::Matchers::WithinAbs(total_energy, tolerance));
}

TEST_CASE("set_beta_invalid")
{
    Four_momentum four_momentum(mass);
    const double too_large = 4.0;
    CHECK_THROWS(four_momentum.set_beta(too_large));
}

TEST_CASE("set_beta_invalid2")
{
    Four_momentum four_momentum(mass);
    const double too_small = -0.5;
    CHECK_THROWS(four_momentum.set_beta(too_small));
}

// tolerance for the equal member function
const double equal_tolerance = 1.0e-12;

TEST_CASE("equal_equal")
{
    Four_momentum four_momentum1(mass, total_energy);
    Four_momentum four_momentum2(mass, total_energy);
    CHECK(four_momentum1.equal(four_momentum2, equal_tolerance));
}

TEST_CASE("equal_different_masses")
{
    Four_momentum four_momentum1(mass, total_energy);
    Four_momentum four_momentum2(mass * 1.01, total_energy);
    CHECK(not four_momentum1.equal(four_momentum2, equal_tolerance));
}

TEST_CASE("equal_different_energies")
{
    Four_momentum four_momentum1(mass, total_energy);
    Four_momentum four_momentum2(mass, total_energy * 1.01);
    CHECK(not four_momentum1.equal(four_momentum2, equal_tolerance));
}

TEST_CASE("equal_different_small_beta")
{
    Four_momentum four_momentum1(mass);
    const double beta = 1.0e-6;
    four_momentum1.set_beta(beta);
    Four_momentum four_momentum2(mass);
    four_momentum2.set_beta(beta * 1.01);
    CHECK(not four_momentum1.equal(four_momentum2, equal_tolerance));
}

TEST_CASE("equal_different_large_gamma")
{
    Four_momentum four_momentum1(mass);
    const double gamma = 1.0e10;
    four_momentum1.set_gamma(gamma);
    Four_momentum four_momentum2(mass);
    four_momentum2.set_gamma(gamma * 1.01);
    CHECK(not four_momentum1.equal(four_momentum2, equal_tolerance));
}

TEST_CASE("test_serialize_xml")
{
    Four_momentum four_momentum(mass, total_energy);
    xml_save<Four_momentum>(four_momentum, "four_momentum.xml");

    Four_momentum loaded;
    xml_load<Four_momentum>(loaded, "four_momentum.xml");
    CHECK(four_momentum.equal(loaded, tolerance));
}

TEST_CASE("test_serialize_json")
{
    Four_momentum four_momentum(mass, total_energy);
    json_save<Four_momentum>(four_momentum, "four_momentum.json");

    Four_momentum loaded;
    json_load<Four_momentum>(loaded, "four_momentum.json");
    CHECK(four_momentum.equal(loaded, tolerance));
}
