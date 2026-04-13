
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "synergia/lattice/madx_reader.h"

#include "synergia/utils/cereal_files.h"

TEST_CASE("print")
{
#if 0
    Dynamic_lattice dl;
    dl.set_variable("x", "1");
    dl.set_variable("a", "x+3");
    //dl.print();
#endif
}

TEST_CASE("dynamic lattice")
{
    std::string str = R"(

        start = 0.35;

        x = 1.0;
        y = {1, 2, 3, 4};
        z = "abc";

        o: drift, l=0.2;

        a: quadrupole, l=0.0, k1=x+1.0;
        b: quadrupole, l=o->l;
        c: quadrupole, l=0.4;
        d: quadrupole, l=o->l*4, k1=o->l*5;

        seq: sequence, l=3.0;
        a, at=0;
        b, at=0;
        c, at=start+0.1;
        d, at=0.6, from=c;
        endsequence;
    )";

    MadX_reader reader;
    reader.parse(str);

    Logger screen(0, LoggerV::DEBUG);

    auto lattice = reader.get_dynamic_lattice("seq");
    auto& elms = lattice.get_elements();

    // original a->k1
    REQUIRE_THAT(lattice.get_elements().front().get_double_attribute("k1"),
                 Catch::Matchers::WithinAbs(2.0, 1e-12));

    // set x
    lattice.get_lattice_tree().set_variable("x", 3.0);

    // a->k1 after setting a new x
    REQUIRE_THAT(lattice.get_elements().front().get_double_attribute("k1"),
                 Catch::Matchers::WithinAbs(4.0, 1e-12));

    // find element d
    auto it = elms.end();
    --it;
    --it;

    // original values
    REQUIRE_THAT(it->get_double_attribute("l"),
                 Catch::Matchers::WithinAbsMatcher(0.8, 1e-12));
    REQUIRE_THAT(it->get_double_attribute("k1"),
                 Catch::Matchers::WithinAbsMatcher(1.0, 1e-12));

    // set o->l to 0.3
    lattice.get_lattice_tree().set_element_attribute("o", "l", 0.3);
    lattice.get_lattice_tree().print();

    // updated values
    REQUIRE_THAT(it->get_double_attribute("l"),
                 Catch::Matchers::WithinAbsMatcher(1.2, 1e-12));
    REQUIRE_THAT(it->get_double_attribute("k1"),
                 Catch::Matchers::WithinAbsMatcher(1.5, 1e-12));

    // element c
    --it;

    // catch parsing error
    REQUIRE_THROWS(it->set_double_attribute("k1", "o->l*3+"));

    // set k1
    REQUIRE_NOTHROW(it->set_double_attribute("k1", "o->l*3"));

    // check value 0.3*3 = 0.9
    REQUIRE_THAT(it->get_double_attribute("k1"),
                 Catch::Matchers::WithinAbsMatcher(0.9, 1e-12));
}

TEST_CASE("serialization")
{
    {
        std::string str = R"(

            start = 0.35;

            x = 1.0;
            y = {1, 2, 3, 4};
            z = "abc";

            a: quadrupole, l=0.0, k1=x+1.0;
            b: quadrupole, l=0.2;
            c: quadrupole, l=0.4;
            d: quadrupole, l=0.8;

            seq: sequence, l=3.0;
            a, at=0;
            b, at=0;
            c, at=start+0.1;
            d, at=0.6, from=c;
            endsequence;
        )";

        MadX_reader reader;
        reader.parse(str);

        auto lattice = reader.get_dynamic_lattice("seq");
        json_save(lattice, "dyn_lattice.json");

        std::cout << lattice.get_lattice_tree().mx.to_madx() << "\n";
        std::cout << lattice.as_string() << "\n";
    }

    {
        std::cout << "\n\nload from file:\n";

        Lattice lattice;
        json_load(lattice, "dyn_lattice.json");

        lattice.get_lattice_tree().set_variable("x", 5.0);
        REQUIRE_THAT(lattice.get_elements().front().get_double_attribute("k1"),
                     Catch::Matchers::WithinAbs(6.0, 1e-12));

        std::cout << lattice.get_lattice_tree().mx.to_madx() << "\n";
        std::cout << lattice.as_string() << "\n";
    }
}
