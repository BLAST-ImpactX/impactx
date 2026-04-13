

#include <catch2/catch_test_macros.hpp>

#include "synergia/utils/cereal_files.h"
#include "synergia/utils/commxx.h"

int
global_rank()
{
    int rank = 0;
    int const status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (status != MPI_SUCCESS) return -1;
    return rank;
}

int
global_size()
{
    int size = 0;
    int const status = MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (status != MPI_SUCCESS) return -1;
    return size;
}

TEST_CASE("static functions", "[commxx]")
{
    CHECK(Commxx::world_rank() == global_rank());
    CHECK(Commxx::world_size() == global_size());
}

TEST_CASE("default construct", "[commxx]")
{
    REQUIRE_NOTHROW(Commxx());
    Commxx defaulted;
    CHECK(defaulted.get_type() == comm_type::world);
    CHECK(defaulted.get_rank() == global_rank());
    CHECK(defaulted.rank() == global_rank());
    CHECK(defaulted.get_size() == global_size());
    CHECK(defaulted.size() == global_size());
    CHECK(defaulted.has_this_rank());
    CHECK(defaulted.is_root());

    REQUIRE_THROWS_AS(defaulted.dup(), std::runtime_error);
    REQUIRE_THROWS_AS(defaulted.split(1), std::runtime_error);
    REQUIRE_THROWS_AS(defaulted.split(0, 1), std::runtime_error);
}

TEST_CASE("construct with type", "[commxx]")
{
    REQUIRE_NOTHROW(Commxx(comm_type::world));
    REQUIRE_NOTHROW(Commxx(comm_type::null));

    // creating a comm_type::regular is not allowed
    REQUIRE_THROWS(Commxx(comm_type::regular));
}

TEST_CASE("static construct", "[commxx]")
{
    SECTION("null")
    {
        auto null = Commxx::Null;
        Commxx comm(comm_type::null);

        REQUIRE(null == comm);
    }

    SECTION("world")
    {
        auto world = Commxx::World;
        Commxx comm(comm_type::world);

        REQUIRE(world == comm);
    }
}

TEST_CASE("construct world", "[commxx]")
{
    SECTION("world")
    {
        Commxx comm(comm_type::world);
        CHECK(comm == MPI_COMM_WORLD);
    }

    SECTION("default world")
    {
        Commxx comm;
        CHECK(comm == MPI_COMM_WORLD);
    }
}

TEST_CASE("construct null", "[commxx]")
{
    Commxx comm(comm_type::null);
    CHECK(comm == MPI_COMM_NULL);
}

TEST_CASE("get_type")
{
    SECTION("world")
    {
        Commxx comm(comm_type::world);
        CHECK(comm.get_type() == comm_type::world);
    }

    SECTION("null")
    {
        Commxx comm(comm_type::null);
        CHECK(comm.get_type() == comm_type::null);
    }
}

TEST_CASE("compare")
{
    SECTION("null to null")
    {
        auto comm1 = Commxx::Null;
        auto comm2 = Commxx::Null;
        REQUIRE(comm1 == comm2);
    }

    SECTION("null to world")
    {
        auto comm1 = Commxx::Null;
        auto comm2 = Commxx::World;
        REQUIRE(comm1 != comm2);
    }

    SECTION("world to world")
    {
        auto comm1 = Commxx::World;
        auto comm2 = Commxx::World;
        REQUIRE(comm1 == comm2);
    }
}

TEST_CASE("get_size")
{
    SECTION("null")
    {
        Commxx comm(comm_type::null);
        REQUIRE_THROWS(comm.get_size());
        REQUIRE_THROWS(comm.size());
    }

    SECTION("world")
    {
        // get size from MPI_COMM_WORLD
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        Commxx comm(comm_type::world);
        REQUIRE_NOTHROW(comm.get_size());
        REQUIRE_NOTHROW(comm.size());

        CHECK(world_size == comm.get_size());
        CHECK(world_size == comm.size());
    }
}

TEST_CASE("get_rank")
{
    SECTION("null")
    {
        Commxx comm(comm_type::null);
        REQUIRE_THROWS(comm.get_rank());
        REQUIRE_THROWS(comm.rank());
    }

    SECTION("world")
    {
        // get size from MPI_COMM_WORLD
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

        Commxx comm(comm_type::world);
        REQUIRE_NOTHROW(comm.get_rank());
        REQUIRE_NOTHROW(comm.rank());

        REQUIRE(world_rank == comm.get_rank());
        REQUIRE(world_rank == comm.rank());
    }
}

TEST_CASE("is_null")
{
    SECTION("null")
    {
        Commxx comm(comm_type::null);
        REQUIRE(comm.is_null());
    }

    SECTION("world")
    {
        Commxx comm(comm_type::world);
        REQUIRE(!comm.is_null());
    }
}

TEST_CASE("is_root")
{
    SECTION("null")
    {
        Commxx comm(comm_type::null);
        REQUIRE(comm.is_root());
    }

    SECTION("world")
    {
        Commxx comm(comm_type::world);
        REQUIRE(comm.is_root());
    }
}

TEST_CASE("has_this_rank")
{
    SECTION("null")
    {
        Commxx comm(comm_type::null);
        REQUIRE(!comm.has_this_rank());
    }

    SECTION("world")
    {
        Commxx comm(comm_type::world);
        REQUIRE(comm.has_this_rank());
    }
}

TEST_CASE("dup")
{
    SECTION("null")
    {
        // cannot dup on a null communicator
        auto comm = Commxx::Null;
        REQUIRE_THROWS(comm.dup());
    }

    SECTION("world non-shared_ptr")
    {
        // dup on a non shared_ptr is not allowed
        auto comm = Commxx::World;
        REQUIRE_THROWS(comm.dup());
    }

    SECTION("world")
    {
        auto comm1_ptr = std::make_shared<Commxx>(comm_type::world);
        REQUIRE_NOTHROW(comm1_ptr->dup());

        auto comm2 = comm1_ptr->dup();
        auto const& comm1 = *comm1_ptr;

        // duplicated communicator is not the same as the original one
        CHECK(comm1 != comm2);

        // but the size and rank of each processor should be the same
        CHECK(comm1.size() == comm2.size());
        CHECK(comm1.rank() == comm2.rank());

        // compare the groups of two communicator
        MPI_Group g1;
        MPI_Group g2;

        MPI_Comm_group(comm1, &g1);
        MPI_Comm_group(comm2, &g2);

        int result;
        int res = MPI_Group_compare(g1, g2, &result);

        // the groups should be MPI_IDENT, meaning the member and
        // ordering are exactly the same
        REQUIRE(res == MPI_SUCCESS);
        CHECK(result == MPI_IDENT);

        MPI_Group_free(&g1);
        MPI_Group_free(&g2);
    }
}

TEST_CASE("split with color only")
{
    SECTION("null")
    {
        // cannot split a null comm
        auto comm = Commxx::Null;
        REQUIRE_THROWS(comm.split(0));
    }

    SECTION("world non-shared_ptr")
    {
        // split on a non shared_ptr is not allowed
        auto comm = Commxx::World;
        REQUIRE_THROWS(comm.split(0));
    }

    SECTION("world")
    {
        auto comm1_ptr = std::make_shared<Commxx>(comm_type::world);
        auto const& comm1 = *comm1_ptr;

        auto world_size = Commxx::world_size();
        auto world_rank = Commxx::world_rank();

        SECTION("no throw")
        {
            REQUIRE_NOTHROW(comm1_ptr->split(0));
        }

        SECTION("same color")
        {
            // this is effective just a dup
            REQUIRE_NOTHROW(comm1_ptr->split(0));
            auto comm2 = comm1_ptr->split(0);

            CHECK(comm2.size() == world_size);
            CHECK(comm2.rank() == world_rank);
        }

        SECTION("one color per rank")
        {
            int color = world_rank;

            REQUIRE_NOTHROW(comm1_ptr->split(color));
            auto comm2 = comm1_ptr->split(color);

            CHECK(comm2.size() == 1);
            CHECK(comm2.rank() == 0);
        }

        if (world_size == 3) {
            SECTION("split color [0 1 1]")
            {
                // workd_rank:       0  1  2
                // color:            0  1  1
                // size after split: 1  2  2
                // rank after split: 0  0  1
                int color[3] = {0, 1, 1};
                int size[3] = {1, 2, 2};
                int rank[3] = {0, 0, 1};

                REQUIRE_NOTHROW(comm1_ptr->split(color[world_rank]));
                auto comm2 = comm1_ptr->split(color[world_rank]);

                CHECK(comm2.size() == size[world_rank]);
                CHECK(comm2.rank() == rank[world_rank]);
            }
        } else if (world_size == 4) {
            SECTION("split color [0 0 1 1]")
            {
                // workd_rank:       0  1  2  3
                // color:            0  0  1  1
                // size after split: 2  2  2  2
                // rank after split: 0  1  0  1
                int color[4] = {0, 0, 1, 1};
                int size[4] = {2, 2, 2, 2};
                int rank[4] = {0, 1, 0, 1};

                REQUIRE_NOTHROW(comm1_ptr->split(color[world_rank]));
                auto comm2 = comm1_ptr->split(color[world_rank]);

                CHECK(comm2.size() == size[world_rank]);
                CHECK(comm2.rank() == rank[world_rank]);
            }

            SECTION("split color [0 1 1 1]")
            {
                // workd_rank:       0  1  2  3
                // color:            0  1  1  1
                // size after split: 1  3  3  3
                // rank after split: 0  0  1  2
                int color[4] = {0, 1, 1, 1};
                int size[4] = {1, 3, 3, 3};
                int rank[4] = {0, 0, 1, 2};

                REQUIRE_NOTHROW(comm1_ptr->split(color[world_rank]));
                auto comm2 = comm1_ptr->split(color[world_rank]);

                CHECK(comm2.size() == size[world_rank]);
                CHECK(comm2.rank() == rank[world_rank]);
            }
        }
    }
}

TEST_CASE("split with color and key")
{
    SECTION("null")
    {
        // cannot split a null comm
        auto comm = Commxx::Null;
        REQUIRE_THROWS(comm.split(0, 0));
    }

    SECTION("world non-shared_ptr")
    {
        // split on a non shared_ptr is not allowed
        auto comm = Commxx::World;
        REQUIRE_THROWS(comm.split(0, 0));
    }

    SECTION("world")
    {
        auto comm1_ptr = std::make_shared<Commxx>(comm_type::world);
        auto const& comm1 = *comm1_ptr;

        auto world_size = Commxx::world_size();
        auto world_rank = Commxx::world_rank();

        SECTION("no throw")
        {
            REQUIRE_NOTHROW(comm1_ptr->split(0, world_rank));
        }

        SECTION("same color, ascending key")
        {
            // this is effective just a dup
            REQUIRE_NOTHROW(comm1_ptr->split(0, world_rank));
            auto comm2 = comm1_ptr->split(0, world_rank);

            CHECK(comm2.size() == world_size);
            CHECK(comm2.rank() == world_rank);
        }

        SECTION("same color, descending key")
        {
            // this is effective just a dup
            REQUIRE_NOTHROW(comm1_ptr->split(0, world_size - world_rank - 1));
            auto comm2 = comm1_ptr->split(0, world_size - world_rank - 1);

            CHECK(comm2.size() == world_size);
            CHECK(comm2.rank() == world_size - world_rank - 1);
        }

        SECTION("one color per rank")
        {
            int color = world_rank;

            REQUIRE_NOTHROW(comm1_ptr->split(color, world_rank));
            auto comm2 = comm1_ptr->split(color, world_rank);

            CHECK(comm2.size() == 1);
            CHECK(comm2.rank() == 0);
        }

        if (world_size == 3) {
            SECTION("split color/key [0/0 1/2 1/1]")
            {
                // workd_rank:       0  1  2
                // color:            0  1  1
                // key:              0  2  1
                // size after split: 1  2  2
                // rank after split: 0  1  0
                int color[3] = {0, 1, 1};
                int key[3] = {0, 2, 1};
                int size[3] = {1, 2, 2};
                int rank[3] = {0, 1, 0};

                REQUIRE_NOTHROW(
                    comm1_ptr->split(color[world_rank], key[world_rank]));
                auto comm2 =
                    comm1_ptr->split(color[world_rank], key[world_rank]);

                CHECK(comm2.size() == size[world_rank]);
                CHECK(comm2.rank() == rank[world_rank]);
            }
        } else if (world_size == 4) {
            SECTION("split color/key [0/0 0/1 1/2 1/1]")
            {
                // workd_rank:       0  1  2  3
                // color:            0  0  1  1
                // key:              0  1  2  1
                // size after split: 2  2  2  2
                // rank after split: 0  1  1  0
                int color[4] = {0, 0, 1, 1};
                int key[4] = {0, 1, 2, 1};
                int size[4] = {2, 2, 2, 2};
                int rank[4] = {0, 1, 1, 0};

                REQUIRE_NOTHROW(
                    comm1_ptr->split(color[world_rank], key[world_rank]));
                auto comm2 =
                    comm1_ptr->split(color[world_rank], key[world_rank]);

                CHECK(comm2.size() == size[world_rank]);
                CHECK(comm2.rank() == rank[world_rank]);
            }

            SECTION("split color/key [0/0 1/3 1/2 1/1]")
            {
                // workd_rank:       0  1  2  3
                // color:            0  1  1  1
                // key:              0  3  2  1
                // size after split: 1  3  3  3
                // rank after split: 0  2  1  0
                int color[4] = {0, 1, 1, 1};
                int key[4] = {0, 3, 2, 1};
                int size[4] = {1, 3, 3, 3};
                int rank[4] = {0, 2, 1, 0};

                REQUIRE_NOTHROW(
                    comm1_ptr->split(color[world_rank], key[world_rank]));
                auto comm2 =
                    comm1_ptr->split(color[world_rank], key[world_rank]);

                CHECK(comm2.size() == size[world_rank]);
                CHECK(comm2.rank() == rank[world_rank]);
            }
        }
    }
}

TEST_CASE("divide")
{
    SECTION("null")
    {
        // cannot divide a null comm
        auto comm = Commxx::Null;
        REQUIRE_THROWS(comm.divide(1));
    }

    SECTION("world non-shared_ptr")
    {
        // split on a non shared_ptr is not allowed
        auto comm = Commxx::World;
        REQUIRE_THROWS(comm.divide(1));
    }

    SECTION("world")
    {
        auto comm1_ptr = std::make_shared<Commxx>(comm_type::world);
        auto const& comm1 = *comm1_ptr;

        auto world_size = Commxx::world_size();
        auto world_rank = Commxx::world_rank();

        SECTION("subgroup size > world_size")
        {
            // effectively a dup
            REQUIRE_NOTHROW(comm1_ptr->divide(world_size + 1));
            auto comm2 = comm1_ptr->divide(world_size + 1);

            CHECK(comm2.size() == world_size);
            CHECK(comm2.rank() == world_rank);
        }

        SECTION("subgroup size = world_size")
        {
            // also effectively a dup
            REQUIRE_NOTHROW(comm1_ptr->divide(world_size));
            auto comm2 = comm1_ptr->divide(world_size);

            CHECK(comm2.size() == world_size);
            CHECK(comm2.rank() == world_rank);
        }

        if (world_size == 2) {
            SECTION("divide into 2 groups")
            {
                int subgroup_size = world_size / 2;

                REQUIRE_NOTHROW(comm1_ptr->divide(subgroup_size));
                auto comm2 = comm1_ptr->divide(subgroup_size);

                CHECK(comm2.size() == subgroup_size);
                CHECK(comm2.rank() == 0);
            }
        } else if (world_size == 3) {
            SECTION("uneven divide")
            {
                REQUIRE_THROWS(comm1_ptr->divide(2));
            }
        } else if (world_size == 4) {
            SECTION("divide into 2 groups")
            {
                int subgroup_size = world_size / 2;

                REQUIRE_NOTHROW(comm1_ptr->divide(subgroup_size));
                auto comm2 = comm1_ptr->divide(subgroup_size);

                int ranks[4] = {0, 1, 0, 1};
                CHECK(comm2.size() == subgroup_size);
                CHECK(comm2.rank() == ranks[world_rank]);
            }
        }
    }
}

TEST_CASE("group")
{
    SECTION("null")
    {
        // cannot group a null comm
        auto comm = Commxx::Null;
        REQUIRE_THROWS(comm.group({}));
    }

    SECTION("world non-shared_ptr")
    {
        // group on a non shared_ptr is not allowed
        auto comm = Commxx::World;
        REQUIRE_THROWS(comm.group({}));
    }

    SECTION("world")
    {
        auto comm1_ptr = std::make_shared<Commxx>(comm_type::world);
        auto const& comm1 = *comm1_ptr;

        auto world_size = Commxx::world_size();
        auto world_rank = Commxx::world_rank();

        SECTION("empty group")
        {
            std::vector<int> group = {};

            REQUIRE_NOTHROW(comm1_ptr->group(group));
            auto comm2 = comm1_ptr->group(group);

            CHECK(comm2.is_null());
        }

        SECTION("group with ranks out of range")
        {
            std::vector<int> group = {100, 102};

            REQUIRE_NOTHROW(comm1_ptr->group(group));
            auto comm2 = comm1_ptr->group(group);

            CHECK(comm2.is_null());
        }

        SECTION("group of 1")
        {
            std::vector<int> group = {2};

            REQUIRE_NOTHROW(comm1_ptr->group(group));
            auto comm2 = comm1_ptr->group(group);

            if (world_rank == 2) {
                CHECK(comm2.size() == 1);
                CHECK(comm2.rank() == 0);
            } else {
                CHECK(comm2.is_null());
            }
        }

        SECTION("group of 2")
        {
            std::vector<int> ranks = {2, 1};

            REQUIRE_NOTHROW(comm1_ptr->group(ranks));
            auto comm2 = comm1_ptr->group(ranks);

            int size[5] = {0, 0, 1, 2, 2};

            if (world_rank == 1) {
                CHECK(comm2.size() == size[world_size]);
                CHECK(comm2.rank() == 0);
            } else if (world_rank == 2) {
                CHECK(comm2.size() == size[world_size]);
                CHECK(comm2.rank() == 1);
            } else {
                CHECK(comm2.is_null());
            }
        }
    }
}

TEST_CASE("combined")
{
    auto world_size = Commxx::world_size();
    auto world_rank = Commxx::world_rank();

    if (world_size == 4) {
        auto comm1_ptr = std::make_shared<Commxx>(comm_type::world);
        auto const& comm1 = *comm1_ptr;

        int color = world_rank < 2 ? 0 : 1;

        auto comm2_ptr = std::make_shared<Commxx>(comm1_ptr->split(color));
        auto comm3_ptr = std::make_shared<Commxx>(comm2_ptr->divide(1));

        CHECK(comm3_ptr->size() == 1);
        CHECK(comm3_ptr->rank() == 0);
    }
}
