
#include <catch2/catch_test_macros.hpp>

#include "synergia/bunch/bunch_particles.h"

// number of particles
constexpr int np = 11;

// idx of lost particles
constexpr auto losts = {3, 5};

void
init_particle_values(BunchParticles& bp)
{
    bp.checkout_particles();

    for (int i = 0; i < np; ++i)
        for (int j = 0; j < 6; ++j)
            bp.hparts(i, j) = i + j * 0.1;

    for (auto l : losts)
        bp.hmasks(l) = 0;

    bp.checkin_particles();
    bp.update_valid_num();
}

void
check_particle_values(BunchParticles& bp)
{
    bp.checkout_particles();

    for (int i = 0; i < np; ++i) {
        for (int j = 0; j < 6; ++j)
            CHECK(bp.hparts(i, j) == i + j * 0.1);

        auto it = std::find(losts.begin(), losts.end(), i);
        if (it == losts.end())
            CHECK(bp.hmasks(i) == 1);
        else
            CHECK(bp.hmasks(i) == 0);
    }
}

TEST_CASE("BunchParticles", "[BunchParticles]")
{
    // create object
    BunchParticles bp(ParticleGroup::regular, np, -1, Commxx::World);

    // init particle data
    init_particle_values(bp);

    // check
    REQUIRE(bp.size() == np);
    REQUIRE(bp.capacity() >= np);
    REQUIRE(bp.num_valid() == np - losts.size());

    SECTION("verify values")
    {
        check_particle_values(bp);
    }

    SECTION("resizing bigger")
    {
        bp.reserve(np + 6, Commxx::World);

        REQUIRE(bp.size() == np);
        REQUIRE(bp.capacity() >= np + 6);
        REQUIRE(bp.num_valid() == np - losts.size());

        check_particle_values(bp);
    }

    SECTION("resizing smaller")
    {
        bp.reserve(np - 6, Commxx::World);

        REQUIRE(bp.size() == np);
        REQUIRE(bp.capacity() >= np);
        REQUIRE(bp.num_valid() == np - losts.size());

        check_particle_values(bp);
    }

#if not defined SYNERGIA_HAVE_OPENPMD
    SECTION("write/read hdf5 file")
    {
        {
            Hdf5_file file(
                "bp_test.h5", Hdf5_file::Flag::truncate, Commxx::World);
            bp.write_file(file, -1, 0, Commxx::World);
        }

        SECTION("read into 0 sized bunch particle")
        {
            BunchParticles bp2(ParticleGroup::regular, 0, 0, Commxx::World);

            REQUIRE(bp2.size() == 0);
            REQUIRE(bp2.capacity() >= 0);
            REQUIRE(bp2.num_valid() == 0);

            Hdf5_file file(
                "bp_test.h5", Hdf5_file::Flag::read_only, Commxx::World);
            bp2.read_file(file, Commxx::World);

            REQUIRE(bp2.size() == np);
            REQUIRE(bp2.capacity() >= np);
            REQUIRE(bp2.num_valid() == np - losts.size());

            check_particle_values(bp2);
        }

        SECTION("read into 0 size but reserved bunch particle")
        {
            BunchParticles bp2(
                ParticleGroup::regular, 0, np + 6, Commxx::World);

            REQUIRE(bp2.size() == 0);
            REQUIRE(bp2.capacity() >= np + 6);
            REQUIRE(bp2.num_valid() == 0);

            Hdf5_file file(
                "bp_test.h5", Hdf5_file::Flag::read_only, Commxx::World);
            bp2.read_file(file, Commxx::World);

            REQUIRE(bp2.size() == np);
            REQUIRE(bp2.capacity() >= np + 6);
            REQUIRE(bp2.num_valid() == np - losts.size());

            check_particle_values(bp2);
        }

        SECTION("read into non-zero sized bunch particle")
        {
            BunchParticles bp2(
                ParticleGroup::regular, np + 6, -1, Commxx::World);

            REQUIRE(bp2.size() == np + 6);
            REQUIRE(bp2.capacity() >= np + 6);
            REQUIRE(bp2.num_valid() == np + 6);

            Hdf5_file file(
                "bp_test.h5", Hdf5_file::Flag::read_only, Commxx::World);
            bp2.read_file(file, Commxx::World);

            REQUIRE(bp2.size() == np);
            REQUIRE(bp2.capacity() >= np + 6);
            REQUIRE(bp2.num_valid() == np - losts.size());

            check_particle_values(bp2);
        }
    }
#endif
}
