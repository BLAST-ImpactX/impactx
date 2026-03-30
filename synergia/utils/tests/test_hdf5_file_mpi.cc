
#include <catch2/catch_test_macros.hpp>

#include "synergia/utils/hdf5_file.h"
#include "synergia/utils/kokkos_views.h"

TEST_CASE("hdf5_file_scalar_read", "[Hdf5_file]")
{
    {
        Hdf5_file file(
            "hdf5_file_scalar_read.h5", Hdf5_file::Flag::truncate, Commxx());

        int i = 3;
        CHECK_NOTHROW(file.write("int", i, false));
    }

    {
        Hdf5_file file(
            "hdf5_file_scalar_read.h5", Hdf5_file::Flag::read_only, Commxx());

        auto i = file.read<int>("int");

        CHECK(i == 3);
    }
}

TEST_CASE("hdf5_file_kv_read", "[Hdf5_file]")
{
    {
        Hdf5_file file(
            "hdf5_file_kv_read.h5", Hdf5_file::Flag::truncate, Commxx());

        karray2d_row arr1("arr1", 4, 3);
        arr1(0, 0) = 10.0;
        arr1(1, 0) = 11.0;
        arr1(2, 0) = 12.0;
        arr1(3, 0) = 13.0;

        arr1(1, 1) = 1.0;
        arr1(2, 1) = 2.0;
        arr1(3, 2) = 3.0;

        CHECK_NOTHROW(file.write("arr1", arr1, false));
    }

    {
        Hdf5_file file(
            "hdf5_file_kv_read.h5", Hdf5_file::Flag::read_only, Commxx());

        auto arr1 = file.read<karray2d_row>("arr1");

        REQUIRE(arr1.extent(0) == 4);
        REQUIRE(arr1.extent(1) == 3);

        CHECK(arr1(1, 1) == 1.0);
        CHECK(arr1(2, 1) == 2.0);
        CHECK(arr1(3, 2) == 3.0);
    }

    {
        int mpi_size = Commxx::world_size();
        int mpi_rank = Commxx::world_rank();

        if (mpi_size <= 4) {
            Hdf5_file file(
                "hdf5_file_kv_read.h5", Hdf5_file::Flag::read_only, Commxx());

            std::vector<int> rows;

            switch (mpi_size) {
                case 1:
                    rows = {4};
                    break;
                case 2:
                    rows = {1, 3};
                    break;
                case 3:
                    rows = {1, 1, 2};
                    break;
                case 4:
                    rows = {1, 1, 0, 2};
                    break;
                default:
                    break;
            }

            auto arr1 = file.read<karray2d_row>("arr1", rows[mpi_rank]);

            REQUIRE(arr1.extent(0) == rows[mpi_rank]);
            REQUIRE(arr1.extent(1) == 3);

            for (int i = 0; i < rows[mpi_rank]; ++i) {
                int off = 0;
                for (int r = 0; r < mpi_rank; ++r)
                    off += rows[r];

                CHECK(arr1(i, 0) == (10.0 + off + i));
            }
        }
    }
}

TEST_CASE("hdf5_file_append", "[Hdf5_file_seq_writer]")
{
    {
        Hdf5_file file("hdf5_file_seq.h5", Hdf5_file::Flag::truncate, Commxx());

        auto mpi_rank = Commxx::world_rank();
        auto mpi_size = Commxx::world_size();

        file.append_single("int", 5);
        file.append_single("int", 6);
        file.append_single("int", 7);

        file.append_collective("int2", 5);
        file.append_collective("int2", 6);
        file.append_collective("int2", 7);

        file.append_collective("int3", mpi_rank * 10 + 1);
        file.append_collective("int3", mpi_rank * 10 + 2);
        file.append_collective("int3", mpi_rank * 10 + 3);

        karray1d arr("arr", 2);
        arr[0] = mpi_rank * 10 + 0;
        arr[1] = mpi_rank * 10 + 1;

        file.append_collective("arr1", arr);

        arr[0] = mpi_rank * 10 + 5;
        arr[1] = mpi_rank * 10 + 6;

        file.append_collective("arr1", arr);

        arr[0] = mpi_rank * 10 + 8;
        arr[1] = mpi_rank * 10 + 9;

        file.append_collective("arr1", arr);

        karray2d arr2("arr2", 2, 3);
        file.append_collective("arr2", arr2);

        karray2d arr3("arr3", 3, 2);
        file.append_single("arr3", arr3);
    }
}

TEST_CASE("hdf5_file", "[Hdf5_file]")
{
    {
        Hdf5_file file(
            "hdf5_file_test.h5", Hdf5_file::Flag::truncate, Commxx());

        file.write_collective("int", 5);
        file.write_single("int2", 6);

        int mpi_rank = Commxx::world_rank();
        std::vector<int> vi(6);
        for (int i = 0; i < vi.size(); ++i)
            vi[i] = i + mpi_rank * 100;

        file.write("vi", vi.data(), vi.size(), true);
        file.write("vi2", vi.data(), vi.size(), false);

        karray2d_row arr1("arr1", 4, 3);
        arr1(2, 1) = 2.0;
        arr1(3, 2) = 3.0;

        file.write("arr1", arr1, true);
        file.write("arr2", arr1, false);
    }

    {
        Hdf5_file file(
            "hdf5_file_test.h5", Hdf5_file::Flag::read_only, Commxx());

        int mpi_rank = Commxx::world_rank();
        int mpi_size = Commxx::world_size();

        int sz = mpi_rank == 0 ? mpi_size * 6 : 0;

        std::vector<int> vi(sz, 3);
        file.read("vi", vi.data(), vi.size());
    }
}
