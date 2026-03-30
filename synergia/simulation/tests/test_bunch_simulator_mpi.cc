#include <catch2/catch_test_macros.hpp>

#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/bunch_simulator_impl.h"

const int charge = 1;
const double mass = 100.0;
const double total_energy = 125.0;

const Reference_particle ref(charge, mass, total_energy);

const int total_num = 1024;
const double real_num = 2.0e12;

void
check_simulator(Bunch_simulator const& sim,
                int nb_pt,
                int nb_st,
                int size,
                int rank)
{
    CHECK(sim[0].get_num_bunches() == nb_pt);
    CHECK(sim[1].get_num_bunches() == nb_st);
}

TEST_CASE("create single bunch", "[Bunch_simulator]")
{
    int mpi_size = Commxx::world_size();
    int mpi_rank = Commxx::world_rank();

    const size_t nb_pt = 1;
    const size_t nb_st = 0;

    if (mpi_size == 1) {
        auto sim = Bunch_simulator::create_single_bunch_simulator(
            ref, total_num, real_num, Commxx());

        CHECK(sim[0].get_num_bunches() == nb_pt);
        CHECK(sim[1].get_num_bunches() == nb_st);

        REQUIRE(sim[0].get_bunch_array_size() == 1);
        REQUIRE(sim[1].get_bunch_array_size() == 0);

        REQUIRE(sim[0][0].get_array_index() == 0);
        REQUIRE(sim[0][0].get_bunch_index() == 0);
        REQUIRE(sim[0][0].get_bucket_index() == 0);

        REQUIRE(sim.has_local_bunch(0, 0));

        CHECK(sim.get_bunch(0, 0).get_array_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_bunch_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_bucket_index() == 0);

        int local_num = total_num / 1;
        CHECK(sim.get_bunch(0, 0).get_local_num() == local_num);

        CHECK_THROWS(sim.get_bunch(0, 1));
    } else if (mpi_size == 2) {
        auto sim = Bunch_simulator::create_single_bunch_simulator(
            ref, total_num, real_num, Commxx());

        CHECK(sim[0].get_num_bunches() == nb_pt);
        CHECK(sim[1].get_num_bunches() == nb_st);

        REQUIRE(sim[0].get_bunch_array_size() == 1);
        REQUIRE(sim[1].get_bunch_array_size() == 0);

        REQUIRE(sim[0][0].get_array_index() == 0);
        REQUIRE(sim[0][0].get_bunch_index() == 0);
        REQUIRE(sim[0][0].get_bucket_index() == 0);

        REQUIRE(sim.has_local_bunch(0, 0));

        CHECK(sim.get_bunch(0, 0).get_array_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_bunch_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_bucket_index() == 0);

        int local_num = total_num / 2;
        CHECK(sim.get_bunch(0, 0).get_local_num() == local_num);

        CHECK_THROWS(sim.get_bunch(0, 1));
    } else if (mpi_size == 3) {
        auto sim = Bunch_simulator::create_single_bunch_simulator(
            ref, total_num, real_num, Commxx());

        CHECK(sim[0].get_num_bunches() == nb_pt);
        CHECK(sim[1].get_num_bunches() == nb_st);

        REQUIRE(sim[0].get_bunch_array_size() == 1);
        REQUIRE(sim[1].get_bunch_array_size() == 0);

        REQUIRE(sim[0][0].get_array_index() == 0);
        REQUIRE(sim[0][0].get_bunch_index() == 0);
        REQUIRE(sim[0][0].get_bucket_index() == 0);

        REQUIRE(sim.has_local_bunch(0, 0));

        CHECK(sim.get_bunch(0, 0).get_array_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_bunch_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_bucket_index() == 0);

        int local_num = 341 + (mpi_rank < 2 ? 0 : 1);
        CHECK(sim.get_bunch(0, 0).get_local_num() == local_num);

        CHECK_THROWS(sim.get_bunch(0, 1));
    } else if (mpi_size == 4) {
        auto sim = Bunch_simulator::create_single_bunch_simulator(
            ref, total_num, real_num, Commxx());

        CHECK(sim[0].get_num_bunches() == nb_pt);
        CHECK(sim[1].get_num_bunches() == nb_st);

        REQUIRE(sim[0].get_bunch_array_size() == 1);
        REQUIRE(sim[1].get_bunch_array_size() == 0);

        REQUIRE(sim[0][0].get_array_index() == 0);
        REQUIRE(sim[0][0].get_bunch_index() == 0);
        REQUIRE(sim[0][0].get_bucket_index() == 0);

        REQUIRE(sim.has_local_bunch(0, 0));

        CHECK(sim.get_bunch(0, 0).get_array_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_bunch_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_bucket_index() == 0);

        int local_num = total_num / 4;
        CHECK(sim.get_bunch(0, 0).get_local_num() == local_num);

        CHECK_THROWS(sim.get_bunch(0, 1));
    } else if (mpi_size == 5) {
        auto sim = Bunch_simulator::create_single_bunch_simulator(
            ref, total_num, real_num, Commxx());

        CHECK(sim[0].get_num_bunches() == nb_pt);
        CHECK(sim[1].get_num_bunches() == nb_st);

        REQUIRE(sim[0].get_bunch_array_size() == 1);
        REQUIRE(sim[1].get_bunch_array_size() == 0);

        REQUIRE(sim[0][0].get_array_index() == 0);
        REQUIRE(sim[0][0].get_bunch_index() == 0);
        REQUIRE(sim[0][0].get_bucket_index() == 0);

        REQUIRE(sim.has_local_bunch(0, 0));

        CHECK(sim.get_bunch(0, 0).get_array_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_bunch_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_bucket_index() == 0);

        int local_num = 204 + (mpi_rank < 1 ? 0 : 1);
        CHECK(sim.get_bunch(0, 0).get_local_num() == local_num);

        CHECK_THROWS(sim.get_bunch(0, 1));
    } else if (mpi_size == 6) {
        auto sim = Bunch_simulator::create_single_bunch_simulator(
            ref, total_num, real_num, Commxx());

        CHECK(sim[0].get_num_bunches() == nb_pt);
        CHECK(sim[1].get_num_bunches() == nb_st);

        REQUIRE(sim[0].get_bunch_array_size() == 1);
        REQUIRE(sim[1].get_bunch_array_size() == 0);

        REQUIRE(sim[0][0].get_array_index() == 0);
        REQUIRE(sim[0][0].get_bunch_index() == 0);
        REQUIRE(sim[0][0].get_bucket_index() == 0);

        REQUIRE(sim.has_local_bunch(0, 0));

        CHECK(sim.get_bunch(0, 0).get_array_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_bunch_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_bucket_index() == 0);

        int local_num = 170 + (mpi_rank < 2 ? 0 : 1);
        CHECK(sim.get_bunch(0, 0).get_local_num() == local_num);

        CHECK_THROWS(sim.get_bunch(0, 1));
    } else {
        CHECK(false);
    }
}

TEST_CASE("create bunch train (3 bunches)", "[Bunch_simulator]")
{
    int mpi_size = Commxx::world_size();
    int mpi_rank = Commxx::world_rank();

    const size_t nb_pt = 3;
    const size_t nb_st = 0;

    const double spacing = 1.0;

    if (mpi_size == 1) {
        auto sim = Bunch_simulator::create_bunch_train_simulator(
            ref, total_num, real_num, nb_pt, spacing, Commxx());

        CHECK(sim[0].get_num_bunches() == 3);
        CHECK(sim[1].get_num_bunches() == 0);

        REQUIRE(sim[0].get_bunch_array_size() == 3);
        REQUIRE(sim[1].get_bunch_array_size() == 0);

        REQUIRE(sim[0][0].get_array_index() == 0);
        REQUIRE(sim[0][0].get_bunch_index() == 0);
        REQUIRE(sim[0][0].get_bucket_index() == 0);

        REQUIRE(sim[0][1].get_array_index() == 1);
        REQUIRE(sim[0][1].get_bunch_index() == 1);
        REQUIRE(sim[0][1].get_bucket_index() == 1);

        REQUIRE(sim[0][2].get_array_index() == 2);
        REQUIRE(sim[0][2].get_bunch_index() == 2);
        REQUIRE(sim[0][2].get_bucket_index() == 2);

        REQUIRE(sim.has_local_bunch(0, 0));
        REQUIRE(sim.has_local_bunch(0, 1));
        REQUIRE(sim.has_local_bunch(0, 2));

        CHECK_THROWS(sim.get_bunch(0, 3));

        CHECK(sim.get_bunch(0, 0).get_array_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_bunch_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_bucket_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_local_num() == total_num);

        CHECK(sim.get_bunch(0, 1).get_array_index() == 1);
        CHECK(sim.get_bunch(0, 1).get_bunch_index() == 1);
        CHECK(sim.get_bunch(0, 1).get_bucket_index() == 1);
        CHECK(sim.get_bunch(0, 1).get_local_num() == total_num);

        CHECK(sim.get_bunch(0, 2).get_array_index() == 2);
        CHECK(sim.get_bunch(0, 2).get_bunch_index() == 2);
        CHECK(sim.get_bunch(0, 2).get_bucket_index() == 2);
        CHECK(sim.get_bunch(0, 2).get_local_num() == total_num);
    } else if (mpi_size == 2) {
        CHECK_THROWS(Bunch_simulator::create_bunch_train_simulator(
            ref, total_num, real_num, nb_pt, spacing, Commxx()));
    } else if (mpi_size == 3) {
        auto sim = Bunch_simulator::create_bunch_train_simulator(
            ref, total_num, real_num, nb_pt, spacing, Commxx());

        CHECK(sim[0].get_num_bunches() == 3);
        CHECK(sim[1].get_num_bunches() == 0);

        REQUIRE(sim[0].get_bunch_array_size() == 1);
        REQUIRE(sim[1].get_bunch_array_size() == 0);

        REQUIRE(sim[0][0].get_array_index() == 0);
        REQUIRE(sim[0][0].get_bunch_index() == mpi_rank);
        REQUIRE(sim[0][0].get_bunch_index() == mpi_rank);

        REQUIRE(sim.has_local_bunch(0, 0) == (mpi_rank == 0));
        REQUIRE(sim.has_local_bunch(0, 1) == (mpi_rank == 1));
        REQUIRE(sim.has_local_bunch(0, 2) == (mpi_rank == 2));

        CHECK_THROWS(sim.get_bunch(0, 3));

        CHECK(sim.get_bunch(0, mpi_rank).get_array_index() == 0);
        CHECK(sim.get_bunch(0, mpi_rank).get_bunch_index() == mpi_rank);
        CHECK(sim.get_bunch(0, mpi_rank).get_bucket_index() == mpi_rank);

        int local_num = total_num / 1;
        CHECK(sim.get_bunch(0, mpi_rank).get_local_num() == local_num);
    } else if (mpi_size == 4) {
        CHECK_THROWS(Bunch_simulator::create_bunch_train_simulator(
            ref, total_num, real_num, nb_pt, spacing, Commxx()));
    } else if (mpi_size == 5) {
        CHECK_THROWS(Bunch_simulator::create_bunch_train_simulator(
            ref, total_num, real_num, nb_pt, spacing, Commxx()));
    } else if (mpi_size == 6) {
        auto sim = Bunch_simulator::create_bunch_train_simulator(
            ref, total_num, real_num, nb_pt, spacing, Commxx());

        CHECK(sim[0].get_num_bunches() == 3);
        CHECK(sim[1].get_num_bunches() == 0);

        REQUIRE(sim[0].get_bunch_array_size() == 1);
        REQUIRE(sim[1].get_bunch_array_size() == 0);

        int bunch_idx = mpi_rank / 2;

        REQUIRE(sim[0][0].get_array_index() == 0);
        REQUIRE(sim[0][0].get_bunch_index() == bunch_idx);
        REQUIRE(sim[0][0].get_bunch_index() == bunch_idx);

        bool has_bunch_0 = ((int)(mpi_rank / 2) == 0); // true if rank is 0 or 1
        bool has_bunch_1 = ((int)(mpi_rank / 2) == 1); // true if rank is 2 or 3
        bool has_bunch_2 = ((int)(mpi_rank / 2) == 2); // true if rank is 4 or 5

        REQUIRE(sim.has_local_bunch(0, 0) == has_bunch_0);
        REQUIRE(sim.has_local_bunch(0, 1) == has_bunch_1);
        REQUIRE(sim.has_local_bunch(0, 2) == has_bunch_2);

        CHECK_THROWS(sim.get_bunch(0, 3));

        CHECK(sim.get_bunch(0, bunch_idx).get_array_index() == 0);
        CHECK(sim.get_bunch(0, bunch_idx).get_bunch_index() == bunch_idx);
        CHECK(sim.get_bunch(0, bunch_idx).get_bucket_index() == bunch_idx);

        int local_num = total_num / 2;
        CHECK(sim.get_bunch(0, bunch_idx).get_local_num() == local_num);
    } else {
        CHECK(false);
    }
}

TEST_CASE("create two bunch trains (4 and 2 bunches)", "[Bunch_simulator]")
{
    int mpi_size = Commxx::world_size();
    int mpi_rank = Commxx::world_rank();

    const size_t nb_pt = 4;
    const size_t nb_st = 2;

    const double spacing = 1.0;

    if (mpi_size == 1) {
        auto sim = Bunch_simulator::create_two_trains_simulator(ref,
                                                                ref,
                                                                total_num,
                                                                real_num,
                                                                nb_pt,
                                                                nb_st,
                                                                spacing,
                                                                spacing,
                                                                Commxx());

        CHECK(sim[0].get_num_bunches() == 4);
        CHECK(sim[1].get_num_bunches() == 2);

        REQUIRE(sim[0].get_bunch_array_size() == 4);
        REQUIRE(sim[1].get_bunch_array_size() == 2);

        REQUIRE(sim[0][0].get_array_index() == 0);
        REQUIRE(sim[0][0].get_bunch_index() == 0);
        REQUIRE(sim[0][0].get_bucket_index() == 0);

        REQUIRE(sim[0][1].get_array_index() == 1);
        REQUIRE(sim[0][1].get_bunch_index() == 1);
        REQUIRE(sim[0][1].get_bucket_index() == 1);

        REQUIRE(sim[0][2].get_array_index() == 2);
        REQUIRE(sim[0][2].get_bunch_index() == 2);
        REQUIRE(sim[0][2].get_bucket_index() == 2);

        REQUIRE(sim[0][3].get_array_index() == 3);
        REQUIRE(sim[0][3].get_bunch_index() == 3);
        REQUIRE(sim[0][3].get_bucket_index() == 3);

        REQUIRE(sim[1][0].get_array_index() == 0);
        REQUIRE(sim[1][0].get_bunch_index() == 0);
        REQUIRE(sim[1][0].get_bucket_index() == 0);

        REQUIRE(sim[1][1].get_array_index() == 1);
        REQUIRE(sim[1][1].get_bunch_index() == 1);
        REQUIRE(sim[1][1].get_bucket_index() == 1);

        REQUIRE(sim.has_local_bunch(0, 0));
        REQUIRE(sim.has_local_bunch(0, 1));
        REQUIRE(sim.has_local_bunch(0, 2));
        REQUIRE(sim.has_local_bunch(0, 3));

        REQUIRE(sim.has_local_bunch(1, 0));
        REQUIRE(sim.has_local_bunch(1, 1));

        CHECK_THROWS(sim.get_bunch(0, 4));
        CHECK_THROWS(sim.get_bunch(1, 2));

        CHECK(sim.get_bunch(0, 0).get_array_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_bunch_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_bucket_index() == 0);
        CHECK(sim.get_bunch(0, 0).get_local_num() == total_num);

        CHECK(sim.get_bunch(0, 1).get_array_index() == 1);
        CHECK(sim.get_bunch(0, 1).get_bunch_index() == 1);
        CHECK(sim.get_bunch(0, 1).get_bucket_index() == 1);
        CHECK(sim.get_bunch(0, 1).get_local_num() == total_num);

        CHECK(sim.get_bunch(0, 2).get_array_index() == 2);
        CHECK(sim.get_bunch(0, 2).get_bunch_index() == 2);
        CHECK(sim.get_bunch(0, 2).get_bucket_index() == 2);
        CHECK(sim.get_bunch(0, 2).get_local_num() == total_num);
    } else if (mpi_size == 2) {
        CHECK_THROWS(Bunch_simulator::create_two_trains_simulator(ref,
                                                                  ref,
                                                                  total_num,
                                                                  real_num,
                                                                  nb_pt,
                                                                  nb_st,
                                                                  spacing,
                                                                  spacing,
                                                                  Commxx()));
    } else if (mpi_size == 3) {
        auto sim = Bunch_simulator::create_two_trains_simulator(ref,
                                                                ref,
                                                                total_num,
                                                                real_num,
                                                                nb_pt,
                                                                nb_st,
                                                                spacing,
                                                                spacing,
                                                                Commxx());

        CHECK(sim[0].get_num_bunches() == 4);
        CHECK(sim[1].get_num_bunches() == 2);

        if (mpi_rank < 2) {
            REQUIRE(sim[0].get_bunch_array_size() == 2);
            REQUIRE(sim[1].get_bunch_array_size() == 0);

            REQUIRE(sim[0][0].get_array_index() == 0);
            REQUIRE(sim[0][0].get_bunch_index() == mpi_rank * 2);
            REQUIRE(sim[0][0].get_bucket_index() == mpi_rank * 2);
            REQUIRE(sim[0][0].get_local_num() == total_num);

            REQUIRE(sim[0][1].get_array_index() == 1);
            REQUIRE(sim[0][1].get_bunch_index() == mpi_rank * 2 + 1);
            REQUIRE(sim[0][1].get_bucket_index() == mpi_rank * 2 + 1);
            REQUIRE(sim[0][1].get_local_num() == total_num);

            REQUIRE(sim.has_local_bunch(0, 0) == (mpi_rank == 0));
            REQUIRE(sim.has_local_bunch(0, 1) == (mpi_rank == 0));
            REQUIRE(sim.has_local_bunch(0, 2) == (mpi_rank == 1));
            REQUIRE(sim.has_local_bunch(0, 3) == (mpi_rank == 1));

            REQUIRE(sim.has_local_bunch(1, 0) == false);
            REQUIRE(sim.has_local_bunch(1, 1) == false);

            REQUIRE(sim[0].get_array_idx_of_bunch(0) ==
                    (mpi_rank == 0 ? 0 : -1));
            REQUIRE(sim[0].get_array_idx_of_bunch(1) ==
                    (mpi_rank == 0 ? 1 : -1));
            REQUIRE(sim[0].get_array_idx_of_bunch(2) ==
                    (mpi_rank == 1 ? 0 : -1));
            REQUIRE(sim[0].get_array_idx_of_bunch(3) ==
                    (mpi_rank == 1 ? 1 : -1));
        } else {
            REQUIRE(sim[0].get_bunch_array_size() == 0);
            REQUIRE(sim[1].get_bunch_array_size() == 2);

            REQUIRE(sim[1][0].get_array_index() == 0);
            REQUIRE(sim[1][0].get_bunch_index() == (mpi_rank - 2) * 2);
            REQUIRE(sim[1][0].get_bunch_index() == (mpi_rank - 2) * 2);
            REQUIRE(sim[1][0].get_local_num() == total_num);

            REQUIRE(sim[1][1].get_array_index() == 1);
            REQUIRE(sim[1][1].get_bunch_index() == (mpi_rank - 2) * 2 + 1);
            REQUIRE(sim[1][1].get_bunch_index() == (mpi_rank - 2) * 2 + 1);
            REQUIRE(sim[1][0].get_local_num() == total_num);

            REQUIRE(sim.has_local_bunch(0, 0) == false);
            REQUIRE(sim.has_local_bunch(0, 1) == false);
            REQUIRE(sim.has_local_bunch(0, 2) == false);
            REQUIRE(sim.has_local_bunch(0, 3) == false);

            REQUIRE(sim.has_local_bunch(1, 0) == true);
            REQUIRE(sim.has_local_bunch(1, 1) == true);

            REQUIRE(sim[1].get_array_idx_of_bunch(0) == 0);
            REQUIRE(sim[1].get_array_idx_of_bunch(1) == 1);
        }
    } else if (mpi_size == 4) {
        CHECK_THROWS(Bunch_simulator::create_two_trains_simulator(ref,
                                                                  ref,
                                                                  total_num,
                                                                  real_num,
                                                                  nb_pt,
                                                                  nb_st,
                                                                  spacing,
                                                                  spacing,
                                                                  Commxx()));
    } else if (mpi_size == 5) {
        CHECK_THROWS(Bunch_simulator::create_two_trains_simulator(ref,
                                                                  ref,
                                                                  total_num,
                                                                  real_num,
                                                                  nb_pt,
                                                                  nb_st,
                                                                  spacing,
                                                                  spacing,
                                                                  Commxx()));
    } else if (mpi_size == 6) {
        auto sim = Bunch_simulator::create_two_trains_simulator(ref,
                                                                ref,
                                                                total_num,
                                                                real_num,
                                                                nb_pt,
                                                                nb_st,
                                                                spacing,
                                                                spacing,
                                                                Commxx());

        CHECK(sim[0].get_num_bunches() == 4);
        CHECK(sim[1].get_num_bunches() == 2);

        if (mpi_rank < 4) {
            REQUIRE(sim[0].get_bunch_array_size() == 1);
            REQUIRE(sim[1].get_bunch_array_size() == 0);

            REQUIRE(sim[0][0].get_array_index() == 0);
            REQUIRE(sim[0][0].get_bunch_index() == mpi_rank);
            REQUIRE(sim[0][0].get_bucket_index() == mpi_rank);
            REQUIRE(sim[0][0].get_local_num() == total_num);

            REQUIRE(sim.has_local_bunch(0, 0) == (mpi_rank == 0));
            REQUIRE(sim.has_local_bunch(0, 1) == (mpi_rank == 1));
            REQUIRE(sim.has_local_bunch(0, 2) == (mpi_rank == 2));
            REQUIRE(sim.has_local_bunch(0, 3) == (mpi_rank == 3));

            REQUIRE(sim.has_local_bunch(1, 0) == false);
            REQUIRE(sim.has_local_bunch(1, 1) == false);

            REQUIRE(sim[0].get_array_idx_of_bunch(0) ==
                    (mpi_rank == 0 ? 0 : -1));
            REQUIRE(sim[0].get_array_idx_of_bunch(1) ==
                    (mpi_rank == 1 ? 0 : -1));
            REQUIRE(sim[0].get_array_idx_of_bunch(2) ==
                    (mpi_rank == 2 ? 0 : -1));
            REQUIRE(sim[0].get_array_idx_of_bunch(3) ==
                    (mpi_rank == 3 ? 0 : -1));
        } else {
            REQUIRE(sim[0].get_bunch_array_size() == 0);
            REQUIRE(sim[1].get_bunch_array_size() == 1);

            REQUIRE(sim[1][0].get_array_index() == 0);
            REQUIRE(sim[1][0].get_bunch_index() == mpi_rank - 4);
            REQUIRE(sim[1][0].get_bunch_index() == mpi_rank - 4);
            REQUIRE(sim[1][0].get_local_num() == total_num);

            REQUIRE(sim.has_local_bunch(0, 0) == false);
            REQUIRE(sim.has_local_bunch(0, 1) == false);
            REQUIRE(sim.has_local_bunch(0, 2) == false);
            REQUIRE(sim.has_local_bunch(0, 3) == false);

            REQUIRE(sim.has_local_bunch(1, 0) == (mpi_rank == 4));
            REQUIRE(sim.has_local_bunch(1, 1) == (mpi_rank == 5));

            REQUIRE(sim[1].get_array_idx_of_bunch(0) ==
                    (mpi_rank == 4 ? 0 : -1));
            REQUIRE(sim[1].get_array_idx_of_bunch(1) ==
                    (mpi_rank == 5 ? 0 : -1));
        }
    } else {
        CHECK(false);
    }
}

TEST_CASE("serialization of two bunch trains (4 and 2 bunches)",
          "[Bunch_simulator]")
{
    int mpi_size = Commxx::world_size();
    int mpi_rank = Commxx::world_rank();

    const size_t nb_pt = 4;
    const size_t nb_st = 2;

    const double spacing = 1.0;

    std::string data;

    if (mpi_size == 1 || mpi_size == 3 || mpi_size == 6) {
        {
            // save
            // auto comm = Commxx().split(0);
            auto comm = Commxx();
            auto sim = Bunch_simulator::create_two_trains_simulator(ref,
                                                                    ref,
                                                                    total_num,
                                                                    real_num,
                                                                    nb_pt,
                                                                    nb_st,
                                                                    spacing,
                                                                    spacing,
                                                                    comm);

            data = sim.dump();
        }

        {
            // load
            auto sim = Bunch_simulator::load_from_string(data);
        }
    }
}
