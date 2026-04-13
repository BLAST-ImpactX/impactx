#include <catch2/catch_test_macros.hpp>

#include "synergia/utils/hdf5_file.h"
#include "synergia/utils/kokkos_views.h"

std::vector<int>
get_rows(int mpi_size)
{
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

    return rows;
}

std::vector<int>
get_offsets(std::vector<int> const& rows)
{
    std::vector<int> offs(rows.size(), 0);
    for (int i = 0; i < rows.size() - 1; ++i)
        offs[i + 1] = offs[i] + rows[i];
    return offs;
}

void
prepare_file(std::string const& fname)
{
    Hdf5_file file(fname, Hdf5_file::Flag::truncate, Commxx());

    // v1: 3
    file.write("v1", 3, false);

    // v2: 2.0
    file.write("v2", 2.0, false);

    // v3: []
    karray1d_row v3("v", 0);
    file.write("v3", v3, false);

    // v4: [2, 3, 4, 5]
    karray1d_row v4("v", 4);
    for (int i = 0; i < 4; ++i)
        v4(i) = i + 2.0;
    file.write("v4", v4, false);

    // v5: []
    karray2d_row v5("v", 0, 2);
    file.write("v5", v5, false);

    // v6: [[0, 2], [0, 3], [0, 4], [0, 5]]
    karray2d_row v6("v", 4, 2);
    for (int i = 0; i < 4; ++i)
        v6(i, 1) = i + 2.0;
    file.write("v6", v6, false);
}

TEST_CASE("hdf5_read_get_dims", "[Hdf5_file_read]")
{
    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    std::stringstream ss;
    ss << "hdf5_read_" << mpi_size << ".h5";

    prepare_file(ss.str());

    Hdf5_file file(ss.str(), Hdf5_file::Flag::read_only, Commxx());

    // names
    auto names = file.get_dataset_names();
    REQUIRE(names.size() == 6);
    CHECK(names[0] == "v1");
    CHECK(names[1] == "v2");
    CHECK(names[2] == "v3");
    CHECK(names[3] == "v4");
    CHECK(names[4] == "v5");
    CHECK(names[5] == "v6");

    // dataset exists
    CHECK(file.has_dataset("v1"));
    CHECK(file.has_dataset("v2"));
    CHECK(file.has_dataset("v3"));
    CHECK(file.has_dataset("v4"));
    CHECK(file.has_dataset("v5"));
    CHECK(file.has_dataset("v6"));
    CHECK_FALSE(file.has_dataset("vx"));

    // v1
    auto d1 = file.get_dims("v1");
    CHECK(d1.size() == 0);

    // v2
    auto d2 = file.get_dims("v2");
    CHECK(d2.size() == 0);

    // v3
    auto d3 = file.get_dims("v3");
    REQUIRE(d3.size() == 1);
    CHECK(d3[0] == 0);

    // v4
    auto d4 = file.get_dims("v4");
    REQUIRE(d4.size() == 1);
    CHECK(d4[0] == 4);

    // v5
    auto d5 = file.get_dims("v5");
    REQUIRE(d5.size() == 2);
    CHECK(d5[0] == 0);
    CHECK(d5[1] == 2);

    // v6
    auto d6 = file.get_dims("v6");
    REQUIRE(d6.size() == 2);
    CHECK(d6[0] == 4);
    CHECK(d6[1] == 2);

    // vx
    CHECK_THROWS(file.get_dims("vx"));
}

TEST_CASE("hdf5_read_single", "[Hdf5_file_read]")
{
    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    std::stringstream ss;
    ss << "hdf5_read_" << mpi_size << ".h5";

    prepare_file(ss.str());

    Hdf5_file file(ss.str(), Hdf5_file::Flag::read_only, Commxx());

    // v1 as int
    auto v1 = file.read<int>("v1");
    CHECK(v1 == 3);

    // v1 as double
    CHECK_THROWS(file.read<double>("v1"));

    // v2 as double
    auto v2 = file.read<double>("v2");
    CHECK(v2 == 2.0);

    // v2 as int
    CHECK_THROWS(file.read<int>("v2"));

    // v3
    auto v3 = file.read<karray1d_row>("v3");
    CHECK(v3.extent(0) == 0);

    // v4
    auto v4 = file.read<karray1d_row>("v4");
    REQUIRE(v4.extent(0) == 4);
    CHECK(v4(0) == 2);
    CHECK(v4(1) == 3);
    CHECK(v4(2) == 4);
    CHECK(v4(3) == 5);

    // v5
    auto v5 = file.read<karray2d_row>("v5");
    REQUIRE(v5.extent(0) == 0);
    REQUIRE(v5.extent(1) == 2);

    // v6
    auto v6 = file.read<karray2d_row>("v6");
    REQUIRE(v6.extent(0) == 4);
    REQUIRE(v6.extent(1) == 2);
    CHECK(v6(0, 1) == 2);
    CHECK(v6(1, 1) == 3);
    CHECK(v6(2, 1) == 4);
    CHECK(v6(3, 1) == 5);
}

TEST_CASE("hdf5_read_collective", "[Hdf5_file_read]")
{
    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    std::stringstream ss;
    ss << "hdf5_read_" << mpi_size << ".h5";

    auto rows = get_rows(mpi_size);
    auto offs = get_offsets(rows);

    prepare_file(ss.str());

    Hdf5_file file(ss.str(), Hdf5_file::Flag::read_only, Commxx());

    // v1 as int
    CHECK_THROWS(file.read<int>("v1", 1));

    // v2 as double
    CHECK_THROWS(file.read<double>("v2", 1));

    // v3
    auto v3 = file.read<karray1d_row>("v3", 0);
    REQUIRE(v3.extent(0) == 0);

    // v4 read with karray
    auto v4 = file.read<karray1d_row>("v4", rows[mpi_rank]);
    REQUIRE(v4.extent(0) == rows[mpi_rank]);
    for (int i = 0; i < rows[mpi_rank]; ++i)
        CHECK(v4(i) == offs[mpi_rank] + i + 2.0);

    // v4 read with raw array
    std::vector<double> rv4(rows[mpi_rank]);
    file.read("v4", rv4.data(), rows[mpi_rank]);
    for (int i = 0; i < rows[mpi_rank]; ++i)
        CHECK(rv4[i] == offs[mpi_rank] + i + 2.0);

    // v5
    auto v5 = file.read<karray2d_row>("v5", 0);
    REQUIRE(v5.extent(0) == 0);
    REQUIRE(v5.extent(1) == 2);

    // v6
    auto v6 = file.read<karray2d_row>("v6", rows[mpi_rank]);
    REQUIRE(v6.extent(0) == rows[mpi_rank]);
    REQUIRE(v6.extent(1) == 2);
    for (int i = 0; i < rows[mpi_rank]; ++i)
        CHECK(v6(i, 1) == offs[mpi_rank] + i + 2.0);
}
