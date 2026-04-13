#include <catch2/catch_test_macros.hpp>

#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/bunch_simulator_impl.h"

TEST_CASE("divide bunches train(1, 0)", "[Bunch_simulator]")
{
    size_t nb_pt = 1;
    size_t nb_st = 0;

    int mpi_size = 1;

    std::vector<int> p_ranks;
    std::vector<int> s_ranks;

    {
        mpi_size = 1;
        REQUIRE_NOTHROW(
            impl::divide_bunches(mpi_size, nb_pt, nb_st, p_ranks, s_ranks));

        CHECK(p_ranks == std::vector<int>{0});
        CHECK(s_ranks == std::vector<int>{});
    }

    {
        mpi_size = 2;
        REQUIRE_NOTHROW(
            impl::divide_bunches(mpi_size, nb_pt, nb_st, p_ranks, s_ranks));

        CHECK(p_ranks == std::vector<int>{0, 1});
        CHECK(s_ranks == std::vector<int>{});
    }

    {
        mpi_size = 3;
        REQUIRE_NOTHROW(
            impl::divide_bunches(mpi_size, nb_pt, nb_st, p_ranks, s_ranks));

        CHECK(p_ranks == std::vector<int>{0, 1, 2});
        CHECK(s_ranks == std::vector<int>{});
    }
}

TEST_CASE("divide bunches train(3, 0)", "[Bunch_simulator]")
{
    size_t nb_pt = 3;
    size_t nb_st = 0;

    int mpi_size = 1;

    std::vector<int> p_ranks;
    std::vector<int> s_ranks;

    {
        mpi_size = 1;
        REQUIRE_NOTHROW(
            impl::divide_bunches(mpi_size, nb_pt, nb_st, p_ranks, s_ranks));

        CHECK(p_ranks == std::vector<int>{0});
        CHECK(s_ranks == std::vector<int>{});
    }

    {
        mpi_size = 2;
        CHECK_THROWS(
            impl::divide_bunches(mpi_size, nb_pt, nb_st, p_ranks, s_ranks));
    }

    {
        mpi_size = 3;
        REQUIRE_NOTHROW(
            impl::divide_bunches(mpi_size, nb_pt, nb_st, p_ranks, s_ranks));

        CHECK(p_ranks == std::vector<int>{0, 1, 2});
        CHECK(s_ranks == std::vector<int>{});
    }

    {
        mpi_size = 4;
        CHECK_THROWS(
            impl::divide_bunches(mpi_size, nb_pt, nb_st, p_ranks, s_ranks));
    }

    {
        mpi_size = 5;
        CHECK_THROWS(
            impl::divide_bunches(mpi_size, nb_pt, nb_st, p_ranks, s_ranks));
    }

    {
        mpi_size = 6;
        REQUIRE_NOTHROW(
            impl::divide_bunches(mpi_size, nb_pt, nb_st, p_ranks, s_ranks));

        CHECK(p_ranks == std::vector<int>{0, 1, 2, 3, 4, 5});
        CHECK(s_ranks == std::vector<int>{});
    }
}

TEST_CASE("divide bunches train(4, 2)", "[Bunch_simulator]")
{
    size_t nb_pt = 4;
    size_t nb_st = 2;

    int mpi_size = 1;

    std::vector<int> p_ranks;
    std::vector<int> s_ranks;

    {
        mpi_size = 1;
        REQUIRE_NOTHROW(
            impl::divide_bunches(mpi_size, nb_pt, nb_st, p_ranks, s_ranks));

        CHECK(p_ranks == std::vector<int>{0});
        CHECK(s_ranks == std::vector<int>{0});
    }

    {
        mpi_size = 2;
        CHECK_THROWS(
            impl::divide_bunches(mpi_size, nb_pt, nb_st, p_ranks, s_ranks));
    }

    {
        mpi_size = 3;
        REQUIRE_NOTHROW(
            impl::divide_bunches(mpi_size, nb_pt, nb_st, p_ranks, s_ranks));

        CHECK(p_ranks == std::vector<int>{0, 1});
        CHECK(s_ranks == std::vector<int>{2});
    }

    {
        mpi_size = 4;
        CHECK_THROWS(
            impl::divide_bunches(mpi_size, nb_pt, nb_st, p_ranks, s_ranks));
    }

    {
        mpi_size = 6;
        REQUIRE_NOTHROW(
            impl::divide_bunches(mpi_size, nb_pt, nb_st, p_ranks, s_ranks));

        CHECK(p_ranks == std::vector<int>{0, 1, 2, 3});
        CHECK(s_ranks == std::vector<int>{4, 5});
    }

    {
        mpi_size = 12;
        REQUIRE_NOTHROW(
            impl::divide_bunches(mpi_size, nb_pt, nb_st, p_ranks, s_ranks));

        CHECK(p_ranks == std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7});
        CHECK(s_ranks == std::vector<int>{8, 9, 10, 11});
    }
}
