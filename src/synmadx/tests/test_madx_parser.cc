
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <boost/math/constants/constants.hpp>

#include "synergia/foundation/physical_constants.h"
#include "synergia/lattice/mx_parse.h"

using namespace synergia;
const double tolerance = 1.0e-12;

TEST_CASE("line_expansion")
{
    using namespace synergia;

    std::string str = R"(
        a: rbend; 
        b: rbend; 
        c: rbend; 
        d: rbend; 
        e: rbend; 
        f: rbend; 
        g: rbend; 
        h: rbend;
	    r: line=(g,h); 
        s: line=(c,r,d); 
        t: line=(,,, 2*s 2*(e,f),,-s,,,, ,-(a,b),,);
    )";

    MadX mx;

    REQUIRE_NOTHROW(parse_madx(str, mx));
    REQUIRE(mx.line_count() == 3);

    MadX_line line = mx.line("t");

    CHECK(line.element_count() == 18);
    CHECK(line.element_name(0) == "c");
    CHECK(line.element_name(1) == "g");
    CHECK(line.element_name(2) == "h");
    CHECK(line.element_name(3) == "d");
    CHECK(line.element_name(4) == "c");
    CHECK(line.element_name(5) == "g");
    CHECK(line.element_name(6) == "h");
    CHECK(line.element_name(7) == "d");
    CHECK(line.element_name(8) == "e");
    CHECK(line.element_name(9) == "f");
    CHECK(line.element_name(10) == "e");
    CHECK(line.element_name(11) == "f");
    CHECK(line.element_name(12) == "d");
    CHECK(line.element_name(13) == "h");
    CHECK(line.element_name(14) == "g");
    CHECK(line.element_name(15) == "c");
    CHECK(line.element_name(16) == "b");
    CHECK(line.element_name(17) == "a");
}

TEST_CASE("line_expansion2")
{
    using namespace synergia;

    std::string str = R"(
        bpm: drift, l=0.01;

        fq1: quadrupole, l=1.0, k1=0.02;
        d1: drift, l=1.0;
        dq2: quadrupole, l=1.0, k1=-0.02;
        d2: drift, l=1.0;
        d3: drift, l=1.0;
        bpm2: bpm;

        line1 : line=(fq1, bpm1, d1);
        bpm1: bpm;

        line2: line=(dq2, d2, bpm1, d3);

        trouble: line=(line1, line2);
    )";

    MadX mx;

    REQUIRE_NOTHROW(parse_madx(str, mx));

    MadX_line line = mx.line("trouble");
    // line.print();

    CHECK(line.element_count() == 7);
    CHECK(line.element_name(0) == "fq1");
    CHECK(line.element_name(1) == "bpm1");
    CHECK(line.element_name(2) == "d1");
    CHECK(line.element_name(3) == "dq2");
    CHECK(line.element_name(4) == "d2");
    CHECK(line.element_name(5) == "bpm1");
    CHECK(line.element_name(6) == "d3");
}

TEST_CASE("sequence")
{
    using namespace synergia;

    std::string str = R"(

        start = 0.35;

        a: quadrupole, l=0.0;
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

    MadX mx;

    REQUIRE_NOTHROW(parse_madx(str, mx));

    MadX_sequence seq = mx.sequence("seq");
    seq.print();

    CHECK(seq.label() == "seq");
    CHECK(seq.length() == 3.0);
    CHECK(seq.element_count() == 4);
    CHECK(seq.refer() == SEQ_REF_CENTRE);
    CHECK(seq.refpos() == "");

    CHECK(seq.element(0).label() == "a");
    CHECK(seq.element(1).label() == "b");
    CHECK(seq.element(2).label() == "c");
    CHECK(seq.element(3).label() == "d");

    CHECK(seq.element(0).name() == "quadrupole");
    CHECK(seq.element(1).name() == "quadrupole");
    CHECK(seq.element(2).name() == "quadrupole");
    CHECK(seq.element(3).name() == "quadrupole");

    CHECK(seq.element(0).attribute_as_number("l") == 0.0);
    CHECK(seq.element(1).attribute_as_number("l") == 0.2);
    CHECK(seq.element(2).attribute_as_number("l") == 0.4);
    CHECK(seq.element(3).attribute_as_number("l") == 0.8);

    CHECK(seq.element_at(0) == 0.0);
    CHECK(seq.element_at(1) == 0.0);
    REQUIRE_THAT(seq.element_at(2),
                 Catch::Matchers::WithinAbs(0.45, tolerance));
    REQUIRE_THAT(seq.element_at(3), Catch::Matchers::WithinAbs(0.6, tolerance));

    CHECK(seq.element_from(0) == 0.0);
    CHECK(seq.element_from(1) == 0.0);
    CHECK(seq.element_from(2) == 0.0);
    REQUIRE_THAT(seq.element_from(3),
                 Catch::Matchers::WithinAbs(0.45, tolerance));
}

TEST_CASE("foil element")
{
    using namespace synergia;

    std::string str = R"(
        f: foil, xmin=-0.05, xmax=0.05, ymin=-0.03, ymax=0.03, thick=600;
    )";

    MadX mx;

    REQUIRE_NOTHROW(parse_madx(str, mx));

    auto foil = mx.command("f");

    CHECK(foil.is_element());
    CHECK(!foil.is_reference());
    CHECK(!foil.is_command());

    REQUIRE_THAT(foil.attribute_as_number("xmin"),
                 Catch::Matchers::WithinAbs(-0.05, tolerance));
    REQUIRE_THAT(foil.attribute_as_number("xmax"),
                 Catch::Matchers::WithinAbs(0.05, tolerance));
    REQUIRE_THAT(foil.attribute_as_number("ymin"),
                 Catch::Matchers::WithinAbs(-0.03, tolerance));
    REQUIRE_THAT(foil.attribute_as_number("ymax"),
                 Catch::Matchers::WithinAbs(0.03, tolerance));
    REQUIRE_THAT(foil.attribute_as_number("thick"),
                 Catch::Matchers::WithinAbs(600, tolerance));
}

TEST_CASE("construct")
{
    MadX mx;
}

TEST_CASE("blank_line")
{
    std::string str = "";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
}

TEST_CASE("garbage")
{
    std::string str = "your momma uses portions->of madx syntax";
    MadX mx;

    CHECK_THROWS_AS(parse_madx(str, mx), std::runtime_error);
}

TEST_CASE("comment1")
{
    std::string str = "//your momma uses portions->of madx syntax";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
}

TEST_CASE("comment2")
{
    std::string str = "!your momma uses portions->of madx syntax";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
}

TEST_CASE("variable_assignment")
{
    std::string str = "x=1;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.variable_as_number("x") == 1);
}

TEST_CASE("mod_variable_assignment")
{
    std::string str = "x:=1;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.variable_as_number("x") == 1);
}

TEST_CASE("keyword_leading_values")
{
    // n.b. "pine" starts with "pi"
    std::string str = "xpine = 3; v = xpine; endmark: marker;";
    str += "fodo: sequence, refer=entry, l=12.0; endmark, at=5.0; endsequence;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));

    MadX_sequence seq = mx.sequence("fodo");
    CHECK(seq.element_count() == 1);
    REQUIRE_THAT(seq.element_at(0), Catch::Matchers::WithinAbs(5.0, tolerance));

    MadX_command cmd = seq.element(0);
    CHECK(cmd.name() == "marker");
}

TEST_CASE("mad_constants")
{
    std::string str =
        "a = pi; b = twopi; c = degrad; d = raddeg; ee = e; "
        "f = emass; g = pmass; h = mumass; i = clight; j = qelect; ";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));

    REQUIRE_THAT(mx.variable_as_number("a"),
                 Catch::Matchers::WithinAbs(
                     boost::math::constants::pi<double>(), tolerance));
    REQUIRE_THAT(mx.variable_as_number("b"),
                 Catch::Matchers::WithinAbs(
                     boost::math::constants::two_pi<double>(), tolerance));
    REQUIRE_THAT(mx.variable_as_number("c"),
                 Catch::Matchers::WithinAbs(
                     180.0 / boost::math::constants::pi<double>(), tolerance));
    REQUIRE_THAT(mx.variable_as_number("d"),
                 Catch::Matchers::WithinAbs(
                     boost::math::constants::pi<double>() / 180.0, tolerance));
    REQUIRE_THAT(mx.variable_as_number("ee"),
                 Catch::Matchers::WithinAbs(boost::math::constants::e<double>(),
                                            tolerance));

    REQUIRE_THAT(mx.variable_as_number("f"),
                 Catch::Matchers::WithinAbs(pconstants::me, tolerance));
    REQUIRE_THAT(mx.variable_as_number("g"),
                 Catch::Matchers::WithinAbs(pconstants::mp, tolerance));
    REQUIRE_THAT(mx.variable_as_number("h"),
                 Catch::Matchers::WithinAbs(pconstants::mmu, tolerance));
    REQUIRE_THAT(mx.variable_as_number("i"),
                 Catch::Matchers::WithinAbs(pconstants::c, tolerance));
    REQUIRE_THAT(mx.variable_as_number("j"),
                 Catch::Matchers::WithinAbs(pconstants::e, tolerance));
}

TEST_CASE("newline_separation")
{
    std::string str = "x=1;\n y=2;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.variable_as_number("x") == 1);
    CHECK(mx.variable_as_number("y") == 2);
}

TEST_CASE("semicolon_separation")
{
    std::string str = "x=1;y=2;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.variable_as_number("x") == 1);
    CHECK(mx.variable_as_number("y") == 2);
}

TEST_CASE("floating_point")
{
    // 7/32 has an exact floating point representation
    std::string str =
        "x=1.234;y=0.21875; ze2=.21875e2; wep2=.21875E2; wep02=.21875e+02;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    REQUIRE_THAT(mx.variable_as_number("x"),
                 Catch::Matchers::WithinAbs(1.234, tolerance));
    CHECK(mx.variable_as_number("y") == 0.21875);
    CHECK(mx.variable_as_number("ze2") == 21.875);
    CHECK(mx.variable_as_number("wep2") == 21.875);
    CHECK(mx.variable_as_number("wep02") == 21.875);
}

TEST_CASE("variable_assignment_expression")
{
    std::string str = "foo.bar=pi*sin(1.2e-4)^0.69;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    REQUIRE_THAT(
        mx.variable_as_number("foo.bar"),
        Catch::Matchers::WithinAbs(boost::math::constants::pi<double>() *
                                       pow(sin(1.2e-4), 0.69),
                                   tolerance));
}

TEST_CASE("caps_variable_assignment")
{
    std::string str = "X=1;Y=X;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.variable_as_number("y") == 1);
}

// ==========================================================================
// commands

TEST_CASE("command")
{
    std::string str = "beam;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.label_count() == 1);

    MadX_command cmd = mx.command("beam");
    CHECK(cmd.name() == "beam");
    CHECK(cmd.attribute_count() == 0);
}

#if 0
// This is defined away because upper-case commands do not work.
TEST_CASE("upper_command")
{
  string str = "BEAM;";
  MadX   mx;

  CHECK_NOTHROW( parse_madx( str, mx ) );
  CHECK( mx.command_count() == 1 );

  MadX_command cmd = mx.command(0);
  CHECK( cmd.name() == "beam" );
  CHECK( cmd.attribute_count() == 0 );
}
#endif

TEST_CASE("command_attrs")
{
    std::string str = "beam,a=1,B=3*(4+5);";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.label_count() == 1);

    MadX_command cmd = mx.command("beam");
    CHECK(cmd.name() == "beam");
    CHECK(cmd.attribute_count() == 2);
    CHECK(cmd.attribute_as_number("a") == 1);
    CHECK(cmd.attribute_as_number("b") == 3 * (4 + 5));
}

TEST_CASE("command_particle_attrs")
{
    std::string str = "beAM, particle=proton;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.label_count() == 1);

    MadX_command cmd = mx.command("beam");
    CHECK(cmd.name() == "beam");
    CHECK(cmd.attribute_count() == 5);
    CHECK(cmd.attribute_as_string("particle") == "proton");
}

TEST_CASE("command_special_attrs1")
{
    std::string str = "mp: multipole, knl:={0, 1, 1}, type=octpn;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.label_count() == 1);

    MadX_command cmd = mx.command("mp");
    CHECK(cmd.name() == "multipole");
    CHECK(cmd.attribute_count() == 2);
    CHECK(cmd.attribute_as_string("type") == "octpn");
}

TEST_CASE("command_special_attrs2")
{
    std::string str = "mp: multipole, knl:={0, 1, 1}, TYPE=wgl;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.label_count() == 1);

    MadX_command cmd = mx.command("mp");
    CHECK(cmd.name() == "multipole");
    CHECK(cmd.attribute_count() == 2);
    CHECK(cmd.attribute_as_string("type") == "wgl");
}

TEST_CASE("command_special_attrs3")
{
    std::string str = "mp: multipole, knl:={0, 1, 1}, type=\"special\";";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.label_count() == 1);

    MadX_command cmd = mx.command("mp");
    CHECK(cmd.name() == "multipole");
    CHECK(cmd.attribute_count() == 2);
    CHECK(cmd.attribute_as_string("type") == "special");
}

TEST_CASE("command_omitted_comma")
{
    std::string str = "call file = './foo.dbx';";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.command_count() == 0);

    CHECK(mx.variable_as_number("a") == 3);
    CHECK(mx.variable_as_number("b") == 2);
}

TEST_CASE("command_beam_particle")
{
    std::string str = " BEAM, PARTICLE=Proton, MASS=0.93827, CHARGE=1., "
                      "ENERGY=0.93827 + 0.160;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));

    MadX_command cmd = mx.command("beam");
    CHECK(cmd.name() == "beam");
    CHECK(cmd.attribute_count() == 6);
    REQUIRE_THAT(cmd.attribute_as_number("mass"),
                 Catch::Matchers::WithinRel(0.93827));
    CHECK(cmd.attribute_as_number("charge") == 1);
    REQUIRE_THAT(cmd.attribute_as_number("energy"),
                 Catch::Matchers::WithinRel(0.93827 + 0.160));
}

TEST_CASE("command_beam_particle_abbreviate")
{
    std::string str = " BEAM, PARTICLE=Prot, MASS=0.93827, CHARGE=1., "
                      "ENERGY=0.93827 + 0.160;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));

    MadX_command cmd = mx.command("beam");
    CHECK(cmd.name() == "beam");
    CHECK(cmd.attribute_count() == 6);
    REQUIRE_THAT(cmd.attribute_as_number("mass"),
                 Catch::Matchers::WithinRel(0.93827));
    CHECK(cmd.attribute_as_number("charge") == 1);
    REQUIRE_THAT(cmd.attribute_as_number("energy"),
                 Catch::Matchers::WithinRel(0.93827 + 0.160));
}

TEST_CASE("command_assign")
{
    std::string str = "q1: quadrupole,l=3.14;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.label_count() == 1);

    MadX_command cmd = mx.command("q1");
    CHECK(cmd.name() == "quadrupole");
    CHECK(cmd.attribute_count() == 1);
    REQUIRE_THAT(cmd.attribute_as_number("l"),
                 Catch::Matchers::WithinAbs(3.14, tolerance));
}

TEST_CASE("subscripted_ident")
{
    std::string str = "foo: quadrupole, a=1; x=foo->a;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.variable_as_number("x") == 1);
}

TEST_CASE("vector_attribute")
{
    std::string str = "mpfi1: multipole,knl:={ 0, 1.1, 2.2, 3.3, 4.4, 5.5 };";
    MadX mx;
    CHECK_NOTHROW(parse_madx(str, mx));
    MadX_command cmd = mx.command("mpfi1");
    std::vector<double> knl(cmd.attribute_as_number_seq("knl"));

    CHECK(knl.size() == 6);
    const double tolerance = 1.0e-12;
    REQUIRE_THAT(knl.at(0), Catch::Matchers::WithinAbs(0.0, tolerance));
    REQUIRE_THAT(knl.at(1), Catch::Matchers::WithinAbs(1.1, tolerance));
    REQUIRE_THAT(knl.at(2), Catch::Matchers::WithinAbs(2.2, tolerance));
    REQUIRE_THAT(knl.at(3), Catch::Matchers::WithinAbs(3.3, tolerance));
    REQUIRE_THAT(knl.at(4), Catch::Matchers::WithinAbs(4.4, tolerance));
    REQUIRE_THAT(knl.at(5), Catch::Matchers::WithinAbs(5.5, tolerance));
}

TEST_CASE("continuation")
{
    std::string str = "q1: quadrupole,l=\n3.14,k1=0.2;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.label_count() == 1);

    MadX_command cmd = mx.command("q1");
    CHECK(cmd.name() == "quadrupole");
    CHECK(cmd.attribute_count() == 2);
    REQUIRE_THAT(cmd.attribute_as_number("l"),
                 Catch::Matchers::WithinAbs(3.14, tolerance));

    REQUIRE_THAT(cmd.attribute_as_number("k1"),
                 Catch::Matchers::WithinAbs(0.2, tolerance));
}

TEST_CASE("continuation2")
{
    std::string str = "q2: quadrupole,l=! )junk)\n3.14,k1=0.2;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.label_count() == 1);

    MadX_command cmd = mx.command("q2");
    CHECK(cmd.name() == "quadrupole");
    CHECK(cmd.attribute_count() == 2);
    REQUIRE_THAT(cmd.attribute_as_number("l"),
                 Catch::Matchers::WithinAbs(3.14, tolerance));
    REQUIRE_THAT(cmd.attribute_as_number("k1"),
                 Catch::Matchers::WithinAbs(0.2, tolerance));
}

TEST_CASE("continuation3")
{
    std::string str = "q3: quadrupole,l=3.14,\nk1=0.2;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.label_count() == 1);

    MadX_command cmd = mx.command("q3");
    CHECK(cmd.name() == "quadrupole");
    CHECK(cmd.attribute_count() == 2);
    REQUIRE_THAT(cmd.attribute_as_number("l"),
                 Catch::Matchers::WithinAbs(3.14, tolerance));
    REQUIRE_THAT(cmd.attribute_as_number("k1"),
                 Catch::Matchers::WithinAbs(0.2, tolerance));
}

TEST_CASE("tkicker")
{
    std::string str = "tk: tkicker, hkick=0.01;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.label_count() == 1);

    MadX_command cmd = mx.command("tk");
    CHECK(cmd.name() == "tkicker");
    CHECK(cmd.attribute_count() == 1);
    REQUIRE_THAT(cmd.attribute_as_number("hkick"),
                 Catch::Matchers::WithinAbs(0.01, tolerance));
}

TEST_CASE("constfoc")
{
    std::string str =
        "cf: constfoc, betaH=0.01, betaV=0.02, betaL=0.03, nuL=0.04;";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.label_count() == 1);

    MadX_command cmd = mx.command("cf");
    CHECK(cmd.name() == "constfoc");
    CHECK(cmd.attribute_count() == 4);

    REQUIRE_THAT(cmd.attribute_as_number("betaH"),
                 Catch::Matchers::WithinAbs(0.01, tolerance));
    REQUIRE_THAT(cmd.attribute_as_number("betaV"),
                 Catch::Matchers::WithinAbs(0.02, tolerance));
    REQUIRE_THAT(cmd.attribute_as_number("betaL"),
                 Catch::Matchers::WithinAbs(0.03, tolerance));
    REQUIRE_THAT(cmd.attribute_as_number("nuL"),
                 Catch::Matchers::WithinAbs(0.04, tolerance));
}

TEST_CASE("comment_at_tail")
{
    std::string str = "tk: tkicker, hkick=0.01;!afdsak";
    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
}

TEST_CASE("parse_file")
{
    MadX mx;

    // parse_madx_file( "./PS.madx", mx );
    // CHECK_NOTHROW( parse_madx_file( "./pstest/PS.madx", mx ) );
}

TEST_CASE("matrix")
{
    std::string str =
        "m1: matrix, type=abc, L=1.2, kick1=0.001, kick2=0.001, kick3=0.002, "
        "kick4=0.002, kick5=0.003, kick6=0.003, rm11=1.1, rm12=1.2, rm32=3.2, "
        "rm54=5.4, tm111=1.11, tm321=3.21;";

    MadX mx;

    CHECK_NOTHROW(parse_madx(str, mx));
    CHECK(mx.label_count() == 1);

    MadX_command cmd = mx.command("m1");
    CHECK(cmd.name() == "matrix");
    CHECK(cmd.attribute_count() == 14);
    CHECK(cmd.attribute_as_string("type") == "abc");
    REQUIRE_THAT(cmd.attribute_as_number("l"),
                 Catch::Matchers::WithinAbs(1.2, tolerance));
    REQUIRE_THAT(cmd.attribute_as_number("kick1"),
                 Catch::Matchers::WithinAbs(0.001, tolerance));
    REQUIRE_THAT(cmd.attribute_as_number("kick2"),
                 Catch::Matchers::WithinAbs(0.001, tolerance));
    REQUIRE_THAT(cmd.attribute_as_number("kick3"),
                 Catch::Matchers::WithinAbs(0.002, tolerance));
    REQUIRE_THAT(cmd.attribute_as_number("kick4"),
                 Catch::Matchers::WithinAbs(0.002, tolerance));
    REQUIRE_THAT(cmd.attribute_as_number("kick5"),
                 Catch::Matchers::WithinAbs(0.003, tolerance));
    REQUIRE_THAT(cmd.attribute_as_number("kick6"),
                 Catch::Matchers::WithinAbs(0.003, tolerance));
    REQUIRE_THAT(cmd.attribute_as_number("rm11"),
                 Catch::Matchers::WithinAbs(1.1, tolerance));
    REQUIRE_THAT(cmd.attribute_as_number("rm12"),
                 Catch::Matchers::WithinAbs(1.2, tolerance));
    REQUIRE_THAT(cmd.attribute_as_number("rm32"),
                 Catch::Matchers::WithinAbs(3.2, tolerance));
    REQUIRE_THAT(cmd.attribute_as_number("rm54"),
                 Catch::Matchers::WithinAbs(5.4, tolerance));
    REQUIRE_THAT(cmd.attribute_as_number("tm111"),
                 Catch::Matchers::WithinAbs(1.11, tolerance));
    REQUIRE_THAT(cmd.attribute_as_number("tm321"),
                 Catch::Matchers::WithinAbs(3.21, tolerance));
}
