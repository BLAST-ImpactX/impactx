#include <catch2/catch_test_macros.hpp>

#include "synergia/utils/hdf5_file.h"
#include "synergia/utils/kokkos_views.h"

TEST_CASE("hdf5_append", "[Hdf5_file_append]")
{
    CHECK(true);
}

TEST_CASE("hdf5_append_scalar", "[Hdf5_file_append]")
{
    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    std::stringstream ss;
    ss << "hdf5_append_scalar_" << mpi_size << ".h5";
    std::string fname = ss.str();

    {
        Hdf5_file file(fname, Hdf5_file::Flag::truncate, Commxx());

        // v1: [3, 4]
        file.append("v1", 3, false);
        file.append("v1", 4, false);

        // v2: [3.2, 4.2]
        file.append("v2", 3.2, false);
        file.append("v2", 4.2, false);

        // v3: [[3, 3, 3, 3], [4, 4, 4, 4]]
        file.append("v3", 3.0, true);
        file.append("v3", 4.0, true);

        // v4: [[10, 11, 12, 13], [20, 21, 22, 23]]
        file.append("v4", 10.0 + mpi_rank, true);
        file.append("v4", 20.0 + mpi_rank, true);
    }

    {
        Hdf5_file file(fname, Hdf5_file::Flag::read_only, Commxx());

        // v1
        auto v1 = file.read<karray1i_row>("v1");
        REQUIRE(v1.extent(0) == 2);
        CHECK(v1(0) == 3);
        CHECK(v1(1) == 4);

        // v2
        auto v2 = file.read<karray1d_row>("v2");
        REQUIRE(v2.extent(0) == 2);
        CHECK(v2(0) == 3.2);
        CHECK(v2(1) == 4.2);

        // v3
        auto v3 = file.read<karray2d_row>("v3");
        REQUIRE(v3.extent(0) == 2);
        REQUIRE(v3.extent(1) == mpi_size);
        for (int i = 0; i < mpi_size; ++i) {
            CHECK(v3(0, i) == 3.0);
            CHECK(v3(1, i) == 4.0);
        }

        // v4
        auto v4 = file.read<karray2d_row>("v4");
        REQUIRE(v4.extent(0) == 2);
        REQUIRE(v4.extent(1) == mpi_size);
        for (int i = 0; i < mpi_size; ++i) {
            CHECK(v4(0, i) == 10.0 + i);
            CHECK(v4(1, i) == 20.0 + i);
        }
    }
}

TEST_CASE("hdf5_append_kv", "[Hdf5_file_append]")
{
    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    std::stringstream ss;
    ss << "hdf5_append_kv_" << mpi_size << ".h5";
    std::string fname = ss.str();

    {
        Hdf5_file file(fname, Hdf5_file::Flag::truncate, Commxx());

        // v1_mem : [1, 2]
        //          [3, 4]
        // v1_file: [ [1, 2], [3, 4] ]
        karray1d_row v1("v", 2);
        v1(0) = 1.0;
        v1(1) = 2.0;
        file.append("v1", v1, false);

        v1(0) = 3.0;
        v1(1) = 4.0;
        file.append("v1", v1, false);

        // v2_mem : [1, 2]
        //          [3, 4]
        // v2_file: [ [1, 2, 1, 2, 1, 2, ... ],
        //            [3, 4, 3, 4, 3, 4, ... ] ]
        karray1d_row v2("v", 2);
        v2(0) = 1.0;
        v2(1) = 2.0;
        file.append("v2", v2, true);

        v2(0) = 3.0;
        v2(1) = 4.0;
        file.append("v2", v2, true);

        // v3_mem:  [[1.1, 1.2], [2.1, 2.2]]
        //          [[3.1, 3.2], [4.1, 4.2]]
        // v3_file: [ [[1.1, 1.2], [2.1, 2.2]],
        //            [[3.1, 3.2], [4.1, 4.2]] ]
        karray2d_row v3("v", 2, 2);
        v3(0, 0) = 1.1;
        v3(0, 1) = 1.2;
        v3(1, 0) = 2.1;
        v3(1, 1) = 2.2;
        file.append("v3", v3, false);

        v3(0, 0) = 3.1;
        v3(0, 1) = 3.2;
        v3(1, 0) = 4.1;
        v3(1, 1) = 4.2;
        file.append("v3", v3, false);

        // v3_mem:  [[1.1, 1.2], [2.1, 2.2]]
        //          [[3.1, 3.2], [4.1, 4.2]]
        // v3_file: [ [[1.1, 1.2], [2.1, 2.2], [1.1, 1.2], [2.1, 2.2], ... ],
        //            [[3.1, 3.2], [4.1, 4.2], [3.1, 3.2], [4.1, 4.2], ... ] ]
        karray2d_row v4("v", 2, 2);
        v4(0, 0) = 1.1;
        v4(0, 1) = 1.2;
        v4(1, 0) = 2.1;
        v4(1, 1) = 2.2;
        file.append("v4", v4, true);

        v4(0, 0) = 3.1;
        v4(0, 1) = 3.2;
        v4(1, 0) = 4.1;
        v4(1, 1) = 4.2;
        file.append("v4", v4, true);
    }

    {
        Hdf5_file file(fname, Hdf5_file::Flag::read_only, Commxx());

        // v1
        auto v1 = file.read<karray2d_row>("v1");
        REQUIRE(v1.extent(0) == 2);
        REQUIRE(v1.extent(1) == 2);
        CHECK(v1(0, 0) == 1.0);
        CHECK(v1(0, 1) == 2.0);
        CHECK(v1(1, 0) == 3.0);
        CHECK(v1(1, 1) == 4.0);

        // v1
        auto v2 = file.read<karray2d_row>("v2");
        REQUIRE(v2.extent(0) == 2);
        REQUIRE(v2.extent(1) == mpi_size * 2);
        for (int i = 0; i < mpi_size; ++i) {
            CHECK(v2(0, i * 2 + 0) == 1.0);
            CHECK(v2(0, i * 2 + 1) == 2.0);
            CHECK(v2(1, i * 2 + 0) == 3.0);
            CHECK(v2(1, i * 2 + 1) == 4.0);
        }

        // v3
        auto v3 = file.read<karray3d_row>("v3");
        REQUIRE(v3.extent(0) == 2);
        REQUIRE(v3.extent(1) == 2);
        REQUIRE(v3.extent(2) == 2);
        CHECK(v3(0, 0, 0) == 1.1);
        CHECK(v3(0, 0, 1) == 1.2);
        CHECK(v3(0, 1, 0) == 2.1);
        CHECK(v3(0, 1, 1) == 2.2);
        CHECK(v3(1, 0, 0) == 3.1);
        CHECK(v3(1, 0, 1) == 3.2);
        CHECK(v3(1, 1, 0) == 4.1);
        CHECK(v3(1, 1, 1) == 4.2);

        // v4
        auto v4 = file.read<karray3d_row>("v4");
        REQUIRE(v4.extent(0) == 2);
        REQUIRE(v4.extent(1) == mpi_size * 2);
        REQUIRE(v4.extent(2) == 2);
        for (int i = 0; i < mpi_size; ++i) {
            CHECK(v4(0, i * 2 + 0, 0) == 1.1);
            CHECK(v4(0, i * 2 + 0, 1) == 1.2);
            CHECK(v4(0, i * 2 + 1, 0) == 2.1);
            CHECK(v4(0, i * 2 + 1, 1) == 2.2);
            CHECK(v4(1, i * 2 + 0, 0) == 3.1);
            CHECK(v4(1, i * 2 + 0, 1) == 3.2);
            CHECK(v4(1, i * 2 + 1, 0) == 4.1);
            CHECK(v4(1, i * 2 + 1, 1) == 4.2);
        }
    }
}

TEST_CASE("hdf5_append_resume", "[Hdf5_file_append]")
{
    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    std::stringstream ss;
    ss << "hdf5_append_resume_" << mpi_size << ".h5";
    std::string fname = ss.str();

    {
        // create a new file
        Hdf5_file file(fname, Hdf5_file::Flag::truncate, Commxx());

        // v1: [3, 4]
        file.append("v1", 3, false);
        file.append("v1", 4, false);
    }

    {
        // close and re-open to append more
        Hdf5_file file(fname, Hdf5_file::Flag::read_write, Commxx());

        // v1: [3, 4, 5, 6]
        file.append("v1", 5, false);
        file.append("v1", 6, false);

        // v2: [5.0, 6.0, 7.0]
        file.append("v2", 5.0, false);
        file.append("v2", 6.0, false);
        file.append("v2", 7.0, false);
    }

    {
        Hdf5_file file(fname, Hdf5_file::Flag::read_only, Commxx());

        // v1
        auto v1 = file.read<karray1i_row>("v1");
        REQUIRE(v1.extent(0) == 4);
        CHECK(v1(0) == 3);
        CHECK(v1(1) == 4);
        CHECK(v1(2) == 5);
        CHECK(v1(3) == 6);

        // v2
        auto v2 = file.read<karray1d_row>("v2");
        REQUIRE(v2.extent(0) == 3);
        CHECK(v2(0) == 5.0);
        CHECK(v2(1) == 6.0);
        CHECK(v2(2) == 7.0);
    }
}
