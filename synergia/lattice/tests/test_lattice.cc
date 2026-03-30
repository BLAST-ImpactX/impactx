#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "catch2/matchers/catch_matchers.hpp"
#include "synergia/lattice/lattice.h"
#include "synergia/lattice/madx_reader.h"

#include "synergia/utils/cereal_files.h"

const std::string name("foo");
const int charge = -1;
const double mass = 100.0;
const double total_energy = 125.0;
const double tolerance = 1.0e-12;

TEST_CASE("construct_lattice")
{
    Lattice lattice(name);
}

TEST_CASE("get_name_lattice")
{
    Lattice lattice(name);
    CHECK(lattice.get_name() == name);
}

TEST_CASE("set_reference_particle")
{
    Lattice lattice(name);
    Reference_particle reference_particle(charge, mass, total_energy);
    lattice.set_reference_particle(reference_particle);
}

TEST_CASE("get_reference_particle")
{
    Lattice lattice(name);
    Reference_particle reference_particle(charge, mass, total_energy);
    lattice.set_reference_particle(reference_particle);
    CHECK(
        lattice.get_reference_particle().equal(reference_particle, tolerance));
}

const double quad_length = 0.2;
const double drift_length = 3.0;
const double bend_length = 4.0;

TEST_CASE("append_fodo")
{
    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element d("quadrupole", "d");
    d.set_double_attribute("l", quad_length);

    Lattice lattice(name);
    lattice.append(f);
    lattice.append(o);
    lattice.append(d);
    lattice.append(o);

    auto it = (lattice.get_elements()).begin();
    CHECK(it != lattice.get_elements().end());
    CHECK(it->get_name() == "f");
    CHECK(it->get_type_name() == "quadrupole");
    REQUIRE_THAT(it->get_double_attribute("l"),
                 Catch::Matchers::WithinAbs(quad_length, tolerance));
    ++it;

    CHECK(it->get_name() == "o");
    CHECK(it->get_type_name() == "drift");
    REQUIRE_THAT(it->get_double_attribute("l"),
                 Catch::Matchers::WithinAbs(drift_length, tolerance));
    ++it;

    CHECK(it->get_name() == "d");
    CHECK(it->get_type_name() == "quadrupole");
    REQUIRE_THAT(it->get_double_attribute("l"),
                 Catch::Matchers::WithinAbs(quad_length, tolerance));
    ++it;

    CHECK(it->get_name() == "o");
    CHECK(it->get_type_name() == "drift");
    REQUIRE_THAT(it->get_double_attribute("l"),
                 Catch::Matchers::WithinAbs(drift_length, tolerance));
}

#if 0
TEST_CASE(derive_internal_attributes)
{
    Lattice_element b("rbend", "b");
    b.set_double_attribute("l", bend_length);
    double bend_angle = 2 * mconstants::pi / 24;
    b.set_double_attribute("angle", bend_angle);

    Lattice lattice(name);
    lattice.append(b);
    lattice.set_defaults();
    lattice.derive_internal_attributes();

    double arc_length = bend_angle * bend_length / (2
            * std::sin(bend_angle / 2));
    CHECK(lattice.get_length() == Approx(arc_length).margin(tolerance));
}
#endif

TEST_CASE("get_length")
{
    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element d("quadrupole", "d");
    d.set_double_attribute("l", quad_length);

    Lattice lattice(name);
    lattice.append(f);
    lattice.append(o);
    lattice.append(d);
    lattice.append(o);

    REQUIRE_THAT(lattice.get_length(),
                 Catch::Matchers::WithinAbs(2 * quad_length + 2 * drift_length,
                                            tolerance));
}

TEST_CASE("get_total_angle1")
{
    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element d("quadrupole", "d");
    d.set_double_attribute("l", quad_length);

    Lattice lattice(name);
    lattice.append(f);
    lattice.append(o);
    lattice.append(d);
    lattice.append(o);

    REQUIRE_THAT(lattice.get_total_angle(),
                 Catch::Matchers::WithinAbs(0.0, tolerance));
}

TEST_CASE("get_total_angle2")
{
    const double pi = 3.141592653589793238462643;
    const int n_cells = 8;
    double bend_angle = 2 * pi / (2 * n_cells);

    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element b("sbend", "b");
    b.set_double_attribute("l", bend_length);
    b.set_double_attribute("angle", bend_angle);
    Lattice_element d("quadrupole", "d");
    d.set_double_attribute("l", quad_length);

    Lattice lattice(name);
    for (int i = 0; i < n_cells; ++i) {
        lattice.append(f);
        lattice.append(o);
        lattice.append(b);
        lattice.append(o);
        lattice.append(d);
        lattice.append(o);
        lattice.append(b);
        lattice.append(o);
    }

    REQUIRE_THAT(lattice.get_total_angle(),
                 Catch::Matchers::WithinAbs(2 * pi, tolerance));
}

TEST_CASE("lattice_from_lattice_element")
{
    Lattice_element f("drift", "foo");
    Lattice lattice(name);
    lattice.append(f);
    for (auto it = lattice.get_elements().begin();
         it != lattice.get_elements().end();
         ++it) {
        CHECK(it->has_lattice());
        CHECK(&(it->get_lattice()) == &lattice);
    }
}

TEST_CASE("copy_lattice")
{
    Lattice_element f("drift", "foo");
    const double foo_length = 1.0;
    f.set_double_attribute("l", foo_length);
    Lattice lattice(name);
    lattice.append(f);
    Reference_particle reference_particle(charge, mass, total_energy);
    lattice.set_reference_particle(reference_particle);

    Lattice copied_lattice(lattice);

    CHECK(lattice.get_elements().size() ==
          copied_lattice.get_elements().size());
    CHECK(lattice.get_elements().begin()->get_name() ==
          copied_lattice.get_elements().begin()->get_name());

    REQUIRE_THAT(lattice.get_elements().begin()->get_length(),
                 Catch::Matchers::WithinAbs(foo_length, tolerance));
    REQUIRE_THAT(copied_lattice.get_elements().begin()->get_length(),
                 Catch::Matchers::WithinAbs(foo_length, tolerance));

    for (auto it = copied_lattice.get_elements().begin();
         it != copied_lattice.get_elements().end();
         ++it) {
        CHECK(it->has_lattice());
        CHECK(&(it->get_lattice()) == &copied_lattice);
    }

    const double new_length = 2.0;
    copied_lattice.get_elements().begin()->set_double_attribute("l",
                                                                new_length);

    REQUIRE_THAT(lattice.get_elements().begin()->get_length(),
                 Catch::Matchers::WithinAbs(foo_length, tolerance));
    REQUIRE_THAT(copied_lattice.get_elements().begin()->get_length(),
                 Catch::Matchers::WithinAbs(new_length, tolerance));

    const double new_energy = 2 * total_energy;
    Reference_particle new_reference_particle(charge, mass, new_energy);
    copied_lattice.set_reference_particle(new_reference_particle);
    REQUIRE_THAT(lattice.get_reference_particle().get_total_energy(),
                 Catch::Matchers::WithinAbs(total_energy, tolerance));
    REQUIRE_THAT(copied_lattice.get_reference_particle().get_total_energy(),
                 Catch::Matchers::WithinAbs(new_energy, tolerance));
}

TEST_CASE("copy_lattice2")
{
    Lattice lattice(name);
    Lattice copied_lattice(lattice);

    Reference_particle reference_particle(charge, mass, total_energy);
    copied_lattice.set_reference_particle(reference_particle);
}

TEST_CASE("copy_lattice_from_lattice_sptr")
{
    Lattice_element f("drift", "foo");
    const double foo_length = 1.0;
    f.set_double_attribute("l", foo_length);

    Lattice lattice(name);
    lattice.append(f);

    Reference_particle reference_particle(charge, mass, total_energy);
    lattice.set_reference_particle(reference_particle);

    Lattice copied_lattice = lattice;

    CHECK(lattice.get_elements().size() ==
          copied_lattice.get_elements().size());
    CHECK(lattice.get_elements().begin()->get_name() ==
          copied_lattice.get_elements().begin()->get_name());

    REQUIRE_THAT(lattice.get_elements().begin()->get_length(),
                 Catch::Matchers::WithinAbs(foo_length, tolerance));
    REQUIRE_THAT(copied_lattice.get_elements().begin()->get_length(),
                 Catch::Matchers::WithinAbs(foo_length, tolerance));

    for (auto it = copied_lattice.get_elements().begin();
         it != copied_lattice.get_elements().end();
         ++it) {
        CHECK(it->has_lattice());
        CHECK(&(it->get_lattice()) == &copied_lattice);
    }
}

TEST_CASE("lattice_reference_particle")
{
    Lattice lattice(name);
    Reference_particle reference_particle(1, 1.0, 1.25);
    lattice.set_reference_particle(reference_particle);
    CHECK(lattice.get_reference_particle().get_mass() == 1.0);
    CHECK(lattice.get_reference_particle().get_gamma() == 1.25);
    CHECK(lattice.get_reference_particle().get_total_energy() == 1.25);
    // set new energy, beta=13/84, gamma=85/84
    lattice.get_reference_particle().set_total_energy(85.0 / 84.0);
    REQUIRE_THAT(lattice.get_reference_particle().get_gamma(),
                 Catch::Matchers::WithinAbs(85.0 / 84.0, tolerance));
    REQUIRE_THAT(lattice.get_reference_particle().get_beta(),
                 Catch::Matchers::WithinAbs(13.0 / 85.0, tolerance));
}

TEST_CASE("lattice_energy")
{
    Lattice lattice(name);
    Reference_particle reference_particle(1, 1.0, 1.25);
    lattice.set_reference_particle(reference_particle);
    CHECK(lattice.get_lattice_energy() == 1.25);
    lattice.set_lattice_energy(7.0);
    CHECK(lattice.get_lattice_energy() == 7.0);
}

TEST_CASE("test_lsexpr")
{
#if 0
    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element d("quadrupole", "d");
    d.set_double_attribute("l", quad_length);

    Lattice lattice(name);
    lattice.append(f);
    lattice.append(o);
    lattice.append(d);
    lattice.append(o);
    Reference_particle reference_particle(charge, mass, total_energy);
    lattice.set_reference_particle(reference_particle);

    Lsexpr lsexpr(lattice.as_lsexpr());
    Lattice loaded(lsexpr);

    auto it = (loaded.get_elements()).begin();
    CHECK(it != loaded.get_elements().end());
    CHECK(it->get_name() == "f");
    CHECK(it->get_type_name() == "quadrupole");
    CHECK(it->get_double_attribute("l") == Approx(quad_length).margin(tolerance));
    ++it;

    CHECK(it->get_name() == "o");
    CHECK(it->get_type_name() == "drift");
    CHECK(it->get_double_attribute("l") == Approx(drift_length).margin(tolerance));
    ++it;

    CHECK(it->get_name() == "d");
    CHECK(it->get_type_name() == "quadrupole");
    CHECK(it->get_double_attribute("l") == Approx(quad_length).margin(tolerance));
    ++it;

    CHECK(it->get_name() == "o");
    CHECK(it->get_type_name() == "drift");
    CHECK(it->get_double_attribute("l") == Approx(drift_length).margin(tolerance));
#endif
}

TEST_CASE("test_serialize1")
{
    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element d("quadrupole", "d");
    d.set_double_attribute("l", quad_length);

    Lattice lattice(name);
    lattice.append(f);
    lattice.append(o);
    lattice.append(d);
    lattice.append(o);
    json_save<Lattice>(lattice, "lattice1.json");

    Lattice loaded;
    json_load<Lattice>(loaded, "lattice1.json");

    auto it = loaded.get_elements().begin();
    CHECK(it != loaded.get_elements().end());
    CHECK(it->get_name() == "f");
    CHECK(it->get_type_name() == "quadrupole");
    REQUIRE_THAT(it->get_double_attribute("l"),
                 Catch::Matchers::WithinAbs(quad_length, tolerance));
    ++it;

    CHECK(it->get_name() == "o");
    CHECK(it->get_type_name() == "drift");
    REQUIRE_THAT(it->get_double_attribute("l"),
                 Catch::Matchers::WithinAbs(drift_length, tolerance));
    ++it;

    CHECK(it->get_name() == "d");
    CHECK(it->get_type_name() == "quadrupole");
    REQUIRE_THAT(it->get_double_attribute("l"),
                 Catch::Matchers::WithinAbs(quad_length, tolerance));
    ++it;

    CHECK(it->get_name() == "o");
    CHECK(it->get_type_name() == "drift");
    REQUIRE_THAT(it->get_double_attribute("l"),
                 Catch::Matchers::WithinAbs(drift_length, tolerance));
}

TEST_CASE("test_serialize2")
{
    Lattice lattice(name);
    Reference_particle reference_particle(charge, mass, total_energy);
    lattice.set_reference_particle(reference_particle);
    json_save<Lattice>(lattice, "lattice2.json");

    Lattice loaded;
    json_load<Lattice>(loaded, "lattice2.json");

    CHECK(loaded.get_reference_particle().equal(reference_particle, tolerance));
}
