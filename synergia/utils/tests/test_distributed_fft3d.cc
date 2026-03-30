#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "synergia/utils/distributed_fft3d.h"

// set DBGPRINT to 1 to print values for tolerance failures
#define DBGPRINT 0

#if DBGPRINT
#include <iostream>
#endif

// define FAILME to force failures
#define FAILME 0

// n.b. We use 0,1,2 here instead of x,y,z because
// we may use z,y,x ordering of arrays.
const int shape0 = 32;
const int shape1 = 16;
const int shape2 = 8;

TEST_CASE("construct")
{
    REQUIRE_NOTHROW(Distributed_fft3d());
}

TEST_CASE("methods")
{
    // construct the fft3d object
    Distributed_fft3d fft;

    // construct the communicator so that every rank does
    // the full calculation, by dividing the world into
    // subgroups of size 1
    auto comm_world = std::make_shared<Commxx>(Commxx::World);
    auto comm = comm_world->divide(1);

    fft.construct({shape0, shape1, shape2}, comm);

    SECTION("roundtrip normalization")
    {
        double normalization = fft.get_roundtrip_normalization();
        CHECK(normalization == 1.0 / (shape0 * shape1 * shape2));
    }

    SECTION("get shape")
    {
        auto s = fft.get_shape();

        CHECK(s.size() == 3);
        CHECK(s[0] == shape0);
        CHECK(s[1] == shape1);
        CHECK(s[2] == shape2);
    }

    SECTION("get comm")
    {
        CHECK(fft.get_comm() == comm);
    }

    SECTION("get lower")
    {
        CHECK(fft.get_lower() == 0);
        CHECK(fft.get_upper() == shape2);
    }
}
