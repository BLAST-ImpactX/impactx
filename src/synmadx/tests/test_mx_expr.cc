
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <iostream>

#include "synergia/lattice/mx_expr.h"
#include "synergia/lattice/mx_parse.h"

using namespace synergia;

const double tolerance = 1.0e-12;

TEST_CASE("mx_expr")
{
    {
        mx_expr expr;
        CHECK(parse_expression("3+1", expr));
        REQUIRE_NOTHROW(mx_eval(expr));
        REQUIRE_THAT(mx_eval(expr), Catch::Matchers::WithinAbs(4.0, tolerance));
        CHECK(mx_expr_is_number(expr));
    }

    {
        mx_expr expr;
        CHECK(parse_expression("3+a", expr));
        REQUIRE_THROWS(mx_eval(expr));
        CHECK(not mx_expr_is_number(expr));
    }
}

TEST_CASE("mx_expr_writer")
{
    {
        mx_expr expr;
        CHECK(parse_expression("3", expr));
        CHECK(mx_expr_str(expr) == "3");
    }

    {
        // "-3" is parsed as a real number, not a uop(-) + 3
        mx_expr expr;
        CHECK(parse_expression("-3", expr));
        CHECK(mx_expr_str(expr) == "-3");
    }

    {
        // "+3" is parsed as a real number, not a uop(+) + 3
        mx_expr expr;
        CHECK(parse_expression("+3", expr));
        CHECK(mx_expr_str(expr) == "3");
    }

    {
        // precision
        double d = 3.1234567890123456;

        mx_expr expr;
        CHECK(parse_expression("3.1234567890123456", expr));

        std::string s = mx_expr_str(expr);
        REQUIRE_THAT(std::stod(s), Catch::Matchers::WithinAbs(d, 1e-14));
    }

    {
        mx_expr expr;
        CHECK(parse_expression("(3)", expr));
        CHECK(mx_expr_str(expr) == "(3)");
    }

    {
        mx_expr expr;
        CHECK(parse_expression("((3))", expr));
        CHECK(mx_expr_str(expr) == "((3))");
    }

    {
        mx_expr expr;
        CHECK(parse_expression("3+1", expr));
        CHECK(mx_expr_str(expr) == "3+1");
    }

    {
        mx_expr expr;
        CHECK(parse_expression("3+1*2", expr));
        CHECK(mx_expr_str(expr) == "3+1*2");
    }

    {
        mx_expr expr;
        CHECK(parse_expression("(3+1)*2", expr));
        CHECK(mx_expr_str(expr) == "(3+1)*2");
    }

    {
        mx_expr expr;
        CHECK(parse_expression("x", expr));
        CHECK(mx_expr_str(expr) == "x");
    }

    {
        mx_expr expr;
        CHECK(parse_expression("+x", expr));
        CHECK(mx_expr_str(expr) == "+x");
    }

    {
        mx_expr expr;
        CHECK(parse_expression("-x", expr));
        CHECK(mx_expr_str(expr) == "-x");
    }

    {
        mx_expr expr;
        CHECK(parse_expression("x->y", expr));
        CHECK(mx_expr_str(expr) == "x->y");
    }

    {
        mx_expr expr;
        CHECK(parse_expression("sin(x)", expr));
        CHECK(mx_expr_str(expr) == "sin(x)");
    }

    {
        mx_expr expr;
        CHECK(parse_expression("sin(-x)", expr));
        CHECK(mx_expr_str(expr) == "sin(-x)");
    }

    {
        mx_expr expr;
        CHECK(
            parse_expression("(-3+1.1234556790123456)*2+sin(a->c) + pi", expr));

        std::cout << mx_expr_str(expr) << "\n";
    }
}
