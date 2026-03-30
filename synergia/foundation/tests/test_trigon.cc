
#include <catch2/catch_test_macros.hpp>

#include "synergia/foundation/trigon.h"

#include <complex>

template <unsigned int P>
using Trig = Trigon<std::complex<double>, P, 6>;
static const std::complex<double> complex_zero = {0.0, 0.0};

TEST_CASE("factorials")
{
    static_assert(factorial(0) == 1);
    static_assert(factorial(1) == 1);
    static_assert(factorial(2) == 2);
    static_assert(factorial(3) == 6);
    static_assert(factorial(4) == 24);
    static_assert(factorial(5) == 120);
    static_assert(factorial(6) == 720);
    static_assert(factorial(7) == 5040);
    static_assert(factorial(8) == 40320);
    static_assert(factorial(9) == 362880);
    static_assert(factorial(10) == 3628800);
    static_assert(factorial(11) == 39916800);
    static_assert(factorial(12) == 479001600);
    // factorial(13) is too large to represent in a 32-bit integer, and so
    // it returns 0.
    static_assert(factorial(13) == 0);
}

TEST_CASE("arr_t")
{
    // Default-constructed arr_t objects should be zero-initialized.
    arr_t<double, 6> x1;
    for (double v : x1) {
        CHECK(v == 0.0);
    }

    arr_t<std::complex<double>, 6> x2;
    for (std::complex<double> v : x2) {
        CHECK(v.real() == 0.0);
        CHECK(v.imag() == 0.0);
    }
}

TEST_CASE("Power 0")
{
    using trig_t = Trig<0>;
    static_assert(std::is_same_v<trig_t::data_type, std::complex<double>>);
    static_assert(trig_t::dim == 6);
    static_assert(trig_t::count == 1);
    trig_t x;
    CHECK(x.value() == complex_zero);
    trig_t y(std::complex<double>(2.5, 3.5));
    CHECK(y.value() == std::complex<double>(2.5, 3.5));
    y.set(complex_zero);
    CHECK(y.value() == complex_zero);
}

TEST_CASE("Power 1")
{
    using trig_t = Trig<1>;
    static_assert(std::is_same_v<trig_t::data_type, std::complex<double>>);
    static_assert(trig_t::dim == 6);
    static_assert(trig_t::count == 6);
    trig_t x;
    CHECK(x.value() == complex_zero);
}

TEST_CASE("Power 2")
{
    using trig_t = Trig<2>;
    static_assert(std::is_same_v<trig_t::data_type, std::complex<double>>);
    static_assert(trig_t::dim == 6);
    static_assert(trig_t::count == 21);
    trig_t x;
    CHECK(x.value() == complex_zero);
}
