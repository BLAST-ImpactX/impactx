#include <numeric>

#include <catch2/catch_test_macros.hpp>

#include "synergia/utils/hdf5_file.h"
#include "synergia/utils/kokkos_views.h"

TEST_CASE("hdf5_file_dim_check_array_nothrow", "[Hdf5_file_write]")
{
    std::stringstream ss;
    ss << "hdf5_file_test_dim_array_nothrow_" << Commxx::world_size() << ".h5";

    Hdf5_file file(ss.str(), Hdf5_file::Flag::truncate, Commxx());

    std::vector<int> vi(6);
    std::iota(vi.begin(), vi.end(), Commxx::world_rank() * 100);

    CHECK_NOTHROW(file.write("vi", vi.data(), vi.size(), true));
    CHECK_NOTHROW(file.write("vi2", vi.data(), vi.size(), false));
}

TEST_CASE("hdf5_file_dim_check_array", "[Hdf5_file_write]")
{
    std::stringstream ss;
    ss << "hdf5_file_test_dim_array_" << Commxx::world_size() << ".h5";

    Hdf5_file file(ss.str(), Hdf5_file::Flag::truncate, Commxx());

    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    std::vector<int> vi(6 + mpi_rank);
    std::iota(vi.begin(), vi.end(), mpi_rank * 100);

    CHECK_NOTHROW(file.write("vi", vi.data(), vi.size(), true));

    if (mpi_size > 1)
        CHECK_THROWS(file.write("vi2", vi.data(), vi.size(), false));
    else
        CHECK_NOTHROW(file.write("vi2", vi.data(), vi.size(), false));
}

TEST_CASE("hdf5_file_dim_check_zero_sized", "[Hdf5_file_write]")
{
    std::stringstream ss;
    ss << "hdf5_file_test_dim_array_zero_sized_" << Commxx::world_size()
       << ".h5";

    Hdf5_file file(ss.str(), Hdf5_file::Flag::truncate, Commxx());

    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    std::vector<int> vi(mpi_rank);
    std::iota(vi.begin(), vi.end(), mpi_rank * 100);

    CHECK_NOTHROW(file.write("vi", vi.data(), vi.size(), true));

    if (mpi_size > 1)
        CHECK_THROWS(file.write("vi2", vi.data(), vi.size(), false));
    else
        CHECK_NOTHROW(file.write("vi2", vi.data(), vi.size(), false));
}

TEST_CASE("hdf5_file_dim_check_kv_nothrow", "[Hdf5_file_write]")
{
    std::stringstream ss;
    ss << "hdf5_file_test_dim_kv_nothrow_" << Commxx::world_size() << ".h5";

    Hdf5_file file(ss.str(), Hdf5_file::Flag::truncate, Commxx());

    karray2d_row arr1("arr1", 4, 3);
    arr1(2, 1) = 2.0;
    arr1(3, 2) = 3.0;

    CHECK_NOTHROW(file.write("arr1", arr1, true));
    CHECK_NOTHROW(file.write("arr2", arr1, false));
}

TEST_CASE("hdf5_file_dim_check_kv_1", "[Hdf5_file_write]")
{
    std::stringstream ss;
    ss << "hdf5_file_test_dim_kv_1_" << Commxx::world_size() << ".h5";

    Hdf5_file file(ss.str(), Hdf5_file::Flag::truncate, Commxx());

    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    karray2d_row arr1("arr1", 4 + mpi_rank, 3);
    arr1(2, 1) = 2.0;
    arr1(3, 2) = 3.0;

    CHECK_NOTHROW(file.write("arr1", arr1, true));

    if (mpi_size > 1)
        CHECK_THROWS(file.write("arr2", arr1, false));
    else
        CHECK_NOTHROW(file.write("arr2", arr1, false));
}

TEST_CASE("hdf5_file_dim_check_kv_2", "[Hdf5_file_write]")
{
    std::stringstream ss;
    ss << "hdf5_file_test_dim_kv_2_" << Commxx::world_size() << ".h5";

    Hdf5_file file(ss.str(), Hdf5_file::Flag::truncate, Commxx());

    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    karray2d_row arr1("arr1", 4 + mpi_rank, 3 + mpi_rank);
    arr1(2, 1) = 2.0;
    arr1(3, 2) = 3.0;

    if (mpi_size > 1)
        CHECK_THROWS(file.write("arr1", arr1, true));
    else
        CHECK_NOTHROW(file.write("arr1", arr1, true));

    if (mpi_size > 1)
        CHECK_THROWS(file.write("arr2", arr1, false));
    else
        CHECK_NOTHROW(file.write("arr2", arr1, false));
}

TEST_CASE("hdf5_file_dim_check_kv_zero_dim", "[Hdf5_file_write]")
{
    std::stringstream ss;
    ss << "hdf5_file_test_dim_kv_zero_dim_" << Commxx::world_size() << ".h5";

    Hdf5_file file(ss.str(), Hdf5_file::Flag::truncate, Commxx());

    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    karray2d_row arr1("arr1", mpi_rank, 3);
    // arr1(2, 1) = 2.0;
    // arr1(3, 2) = 3.0;

    CHECK_NOTHROW(file.write("arr1", arr1, true));

    if (mpi_size > 1)
        CHECK_THROWS(file.write("arr2", arr1, false));
    else
        CHECK_NOTHROW(file.write("arr2", arr1, false));
}

TEST_CASE("hdf5_file_dim_check_kv_diff_rank", "[Hdf5_file_write]")
{
    std::stringstream ss;
    ss << "hdf5_file_test_dim_kv_diff_rank_" << Commxx::world_size() << ".h5";

    Hdf5_file file(ss.str(), Hdf5_file::Flag::truncate, Commxx());

    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    if (mpi_rank == 0) {
        karray2d_row arr1("arr1", 4, 3);

        if (mpi_size == 1) {
            CHECK_NOTHROW(file.write("arr1", arr1, true));
            CHECK_NOTHROW(file.write("arr2", arr1, false));
        } else {
            CHECK_THROWS(file.write("arr1", arr1, true));
            CHECK_THROWS(file.write("arr2", arr1, false));
        }
    } else {
        karray3d_row arr1("arr1", 4, 3, 2);

        CHECK_THROWS(file.write("arr1", arr1, true));
        CHECK_THROWS(file.write("arr2", arr1, false));
    }
}

TEST_CASE("hdf5_file_scalar_write", "[Hdf5_file_write]")
{
    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    {
        std::stringstream ss;
        ss << "hdf5_file_scalar_write_" << mpi_size << ".h5";
        Hdf5_file file(ss.str(), Hdf5_file::Flag::truncate, Commxx());

        int i = 3;
        // i1: 3, i2: [3,3,3,3]
        CHECK_NOTHROW(file.write("i1", i, false));
        CHECK_NOTHROW(file.write("i2", i, true));

        double d = 4.2;
        // d1: 4.2, d2: [4.2, 4.2, 4.2, 4.2]
        CHECK_NOTHROW(file.write("d1", d, false));
        CHECK_NOTHROW(file.write("d2", d, true));

        double d2 = 5.2 + mpi_rank;
        // d3: 8.2 (root_rank = 3 in mpi_4)
        // d4: [5.2, 6.2, 7.2, 8.2]
        CHECK_NOTHROW(file.write("d3", d2, false));
        CHECK_NOTHROW(file.write("d4", d2, true));
    }

    {
        std::stringstream ss;
        ss << "hdf5_file_scalar_write_" << mpi_size << ".h5";
        Hdf5_file file(ss.str(), Hdf5_file::Flag::read_only, Commxx());

        auto i1 = file.read<int>("i1");
        CHECK(i1 == 3);

        auto d1 = file.read<double>("d1");
        CHECK(d1 == 4.2);

        auto i2 = file.read<karray1i_row>("i2");
        REQUIRE(i2.extent(0) == mpi_size);
        for (int i = 0; i < mpi_size; ++i)
            CHECK(i2(i) == 3);

        auto d2 = file.read<karray1d_row>("d2");
        REQUIRE(d2.extent(0) == mpi_size);
        for (int i = 0; i < mpi_size; ++i)
            CHECK(d2(i) == 4.2);

        auto d3 = file.read<double>("d3");
        // if(mpi_rank == file.get_root_rank()) CHECK(d3 == 5.2 + mpi_rank);

        auto d4 = file.read<karray1d_row>("d4");
        REQUIRE(d4.extent(0) == mpi_size);
        for (int i = 0; i < mpi_size; ++i)
            CHECK(d4(i) == 5.2 + i);
    }
}

TEST_CASE("hdf5_file_kv1d_write", "[Hdf5_file_write]")
{
    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    {
        std::stringstream ss;
        ss << "hdf5_file_kv1d_write_" << mpi_size << ".h5";
        Hdf5_file file(ss.str(), Hdf5_file::Flag::truncate, Commxx());

        // r0: [0, 1, 2]
        // r1: [0, 1, 2]
        karray1d_row ka1("ka", 3);
        ka1(0) = 0.0;
        ka1(1) = 1.0;
        ka1(2) = 2.0;

        // k1: [0, 1, 2]
        CHECK_NOTHROW(file.write("k1", ka1, false));

        // k2: [0, 1, 2, 0, 1, 2]
        CHECK_NOTHROW(file.write("k2", ka1, true));

        // r0: [ 0,  1,  2]
        // r1: [10, 11, 12, 0]
        // r2: [20, 21, 22, 0, 0]
        karray1d_row ka2("ka", 3 + mpi_rank);
        ka2(0) = 0.0 + mpi_rank * 10;
        ka2(1) = 1.0 + mpi_rank * 10;
        ka2(2) = 2.0 + mpi_rank * 10;

        // k3: [0, 1, 2] if mpi_size=1
        if (mpi_size > 1)
            CHECK_THROWS(file.write("k3", ka2, false));
        else
            CHECK_NOTHROW(file.write("k3", ka2, false));

        // k4: [0, 1, 2, 10, 11, 12, 0, ... ]
        CHECK_NOTHROW(file.write("k4", ka2, true));

        // r0: []
        // r1: []
        karray1d_row ka3("ka", 0);

        // k5: []
        CHECK_NOTHROW(file.write("k5", ka3, false));

        // k6: []
        CHECK_NOTHROW(file.write("k6", ka3, true));

        // r0: []
        // r1: [0]
        karray1d_row ka4("ka", mpi_rank);

        // k7: []
        if (mpi_size > 1)
            CHECK_THROWS(file.write("k7", ka4, false));
        else
            CHECK_NOTHROW(file.write("k7", ka4, false));

        // k8: [0]
        CHECK_NOTHROW(file.write("k8", ka4, true));
    }

    {
        std::stringstream ss;
        ss << "hdf5_file_kv1d_write_" << mpi_size << ".h5";
        Hdf5_file file(ss.str(), Hdf5_file::Flag::read_only, Commxx());

        // k1
        auto k1 = file.read<karray1d_row>("k1");
        REQUIRE(k1.extent(0) == 3);
        CHECK(k1(0) == 0.0);
        CHECK(k1(1) == 1.0);
        CHECK(k1(2) == 2.0);

        // k2
        auto k2 = file.read<karray1d_row>("k2");
        REQUIRE(k2.extent(0) == 3 * mpi_size);
        for (int r = 0; r < mpi_size; ++r) {
            CHECK(k2(r * 3 + 0) == 0.0);
            CHECK(k2(r * 3 + 1) == 1.0);
            CHECK(k2(r * 3 + 2) == 2.0);
        }

        // k3
        if (mpi_size == 1) {
            auto k3 = file.read<karray1d_row>("k3");
            REQUIRE(k3.extent(0) == 3);
            CHECK(k3(0) == 0.0);
            CHECK(k3(1) == 1.0);
            CHECK(k3(2) == 2.0);
        } else {
            CHECK_THROWS(file.read<karray1d_row>("k3"));
        }

        // k4
        auto k4 = file.read<karray1d_row>("k4");
        REQUIRE(k4.extent(0) == 3 * mpi_size + (mpi_size - 1) * mpi_size / 2);
        for (int r = 0; r < mpi_size; ++r) {
            int off = 0;
            for (int i = 0; i < r; ++i)
                off += (3 + i);

            CHECK(k4(off + 0) == 0.0 + r * 10);
            CHECK(k4(off + 1) == 1.0 + r * 10);
            CHECK(k4(off + 2) == 2.0 + r * 10);
        }

        // k5
        auto k5 = file.read<karray1d_row>("k5");
        REQUIRE(k5.extent(0) == 0);

        // k6
        auto k6 = file.read<karray1d_row>("k6");
        REQUIRE(k6.extent(0) == 0);

        // k7
        if (mpi_size == 1) {
            auto k7 = file.read<karray1d_row>("k7");
            REQUIRE(k7.extent(0) == 0);
        } else {
            CHECK_THROWS(file.read<karray1d_row>("k7"));
        }

        // k8
        auto k8 = file.read<karray1d_row>("k8");
        int k8_size = (mpi_size - 1) * mpi_size / 2;
        REQUIRE(k8.extent(0) == k8_size);
        for (int i = 0; i < k8_size; ++i)
            CHECK(k8(i) == 0);
    }
}

TEST_CASE("hdf5_file_kv2d_write", "[Hdf5_file_write]")
{
    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    {
        std::stringstream ss;
        ss << "hdf5_file_kv2d_write_" << mpi_size << ".h5";
        Hdf5_file file(ss.str(), Hdf5_file::Flag::truncate, Commxx());

        // r0: [0, 1, 2]
        // r1: [0, 1, 2]
        karray2d_row ka1("ka", 3, 2);
        ka1(0, 1) = 0.0;
        ka1(1, 1) = 1.0;
        ka1(2, 1) = 2.0;

        // k1: [0, 1, 2]
        CHECK_NOTHROW(file.write("k1", ka1, false));

        // k2: [0, 1, 2, 0, 1, 2]
        CHECK_NOTHROW(file.write("k2", ka1, true));

        // r0: [ 0,  1,  2]
        // r1: [10, 11, 12, 0]
        // r2: [20, 21, 22, 0, 0]
        karray2d_row ka2("ka", 3 + mpi_rank, 2);
        ka2(0, 1) = 0.0 + mpi_rank * 10;
        ka2(1, 1) = 1.0 + mpi_rank * 10;
        ka2(2, 1) = 2.0 + mpi_rank * 10;

        // k3: [0, 1, 2] if mpi_size=1
        if (mpi_size > 1)
            CHECK_THROWS(file.write("k3", ka2, false));
        else
            CHECK_NOTHROW(file.write("k3", ka2, false));

        // k4: [0, 1, 2, 10, 11, 12, 0, ... ]
        CHECK_NOTHROW(file.write("k4", ka2, true));

        // r0: []
        // r1: []
        karray2d_row ka3("ka", 0, 2);

        // k5: []
        CHECK_NOTHROW(file.write("k5", ka3, false));

        // k6: []
        CHECK_NOTHROW(file.write("k6", ka3, true));

        // r0: []
        // r1: [0]
        karray2d_row ka4("ka", mpi_rank, 2);

        // k7: []
        if (mpi_size > 1)
            CHECK_THROWS(file.write("k7", ka4, false));
        else
            CHECK_NOTHROW(file.write("k7", ka4, false));

        // k8: [0]
        CHECK_NOTHROW(file.write("k8", ka4, true));
    }

    {
        std::stringstream ss;
        ss << "hdf5_file_kv2d_write_" << mpi_size << ".h5";
        Hdf5_file file(ss.str(), Hdf5_file::Flag::read_only, Commxx());

        // k1
        auto k1 = file.read<karray2d_row>("k1");
        REQUIRE(k1.extent(0) == 3);
        REQUIRE(k1.extent(1) == 2);
        CHECK(k1(0, 1) == 0.0);
        CHECK(k1(1, 1) == 1.0);
        CHECK(k1(2, 1) == 2.0);

        // k2
        auto k2 = file.read<karray2d_row>("k2");
        REQUIRE(k2.extent(0) == 3 * mpi_size);
        REQUIRE(k2.extent(1) == 2);
        for (int r = 0; r < mpi_size; ++r) {
            CHECK(k2(r * 3 + 0, 1) == 0.0);
            CHECK(k2(r * 3 + 1, 1) == 1.0);
            CHECK(k2(r * 3 + 2, 1) == 2.0);
        }

        // k3
        if (mpi_size == 1) {
            auto k3 = file.read<karray2d_row>("k3");
            REQUIRE(k3.extent(0) == 3);
            REQUIRE(k3.extent(1) == 2);
            CHECK(k3(0, 1) == 0.0);
            CHECK(k3(1, 1) == 1.0);
            CHECK(k3(2, 1) == 2.0);
        } else {
            CHECK_THROWS(file.read<karray2d_row>("k3"));
        }

        // k4
        auto k4 = file.read<karray2d_row>("k4");
        REQUIRE(k4.extent(0) == 3 * mpi_size + (mpi_size - 1) * mpi_size / 2);
        REQUIRE(k4.extent(1) == 2);
        for (int r = 0; r < mpi_size; ++r) {
            int off = 0;
            for (int i = 0; i < r; ++i)
                off += (3 + i);

            CHECK(k4(off + 0, 1) == 0.0 + r * 10);
            CHECK(k4(off + 1, 1) == 1.0 + r * 10);
            CHECK(k4(off + 2, 1) == 2.0 + r * 10);
        }

        // k5
        auto k5 = file.read<karray2d_row>("k5");
        REQUIRE(k5.extent(0) == 0);
        REQUIRE(k5.extent(1) == 2);

        // k6
        auto k6 = file.read<karray2d_row>("k6");
        REQUIRE(k6.extent(0) == 0);
        REQUIRE(k6.extent(1) == 2);

        // k7
        if (mpi_size == 1) {
            auto k7 = file.read<karray2d_row>("k7");
            REQUIRE(k7.extent(0) == 0);
            REQUIRE(k7.extent(1) == 2);
        } else {
            CHECK_THROWS(file.read<karray2d_row>("k7"));
        }

        // k8
        auto k8 = file.read<karray2d_row>("k8");
        int k8_size = (mpi_size - 1) * mpi_size / 2;
        REQUIRE(k8.extent(0) == k8_size);
        REQUIRE(k8.extent(1) == 2);
        for (int i = 0; i < k8_size; ++i)
            CHECK(k8(i, 1) == 0);
    }
}
