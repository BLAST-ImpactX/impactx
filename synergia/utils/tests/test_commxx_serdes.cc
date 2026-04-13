

#include <catch2/catch_test_macros.hpp>

#include "synergia/utils/cereal_files.h"
#include "synergia/utils/commxx.h"

TEST_CASE("null commxx")
{
    std::stringstream ss;
    ss << "commxx_serdes_null_" << Commxx::world_size() << "_"
       << Commxx::world_rank() << ".json";

    // create and save
    {
        auto comm = Commxx::Null;
        REQUIRE_NOTHROW(json_save(comm, ss.str()));
    }

    // reload
    {
        Commxx comm;
        REQUIRE_NOTHROW(json_load(comm, ss.str()));

        CHECK(comm.is_null());
        CHECK(comm.get_type() == comm_type::null);
    }
}

TEST_CASE("world commxx")
{
    std::stringstream ss;
    ss << "commxx_serdes_world_" << Commxx::world_size() << "_"
       << Commxx::world_rank() << ".json";

    // create and save
    {
        auto comm = Commxx::World;
        REQUIRE_NOTHROW(json_save(comm, ss.str()));
    }

    // reload
    {
        Commxx comm;
        REQUIRE_NOTHROW(json_load(comm, ss.str()));

        CHECK(!comm.is_null());
        CHECK(comm.get_type() == comm_type::world);

        CHECK(comm.size() == Commxx::world_size());
        CHECK(comm.rank() == Commxx::world_rank());
    }
}

TEST_CASE("world commxx shared_ptr")
{
    std::stringstream ss;
    ss << "commxx_serdes_world_sp_" << Commxx::world_size() << "_"
       << Commxx::world_rank() << ".json";

    // create and save
    {
        auto comm = std::make_shared<Commxx>(comm_type::world);
        REQUIRE_NOTHROW(json_save(comm, ss.str()));
    }

    // reload
    {
        auto comm = std::make_shared<Commxx>();
        REQUIRE_NOTHROW(json_load(comm, ss.str()));

        CHECK(!comm->is_null());
        CHECK(comm->get_type() == comm_type::world);

        CHECK(comm->size() == Commxx::world_size());
        CHECK(comm->rank() == Commxx::world_rank());
    }
}

TEST_CASE("combo1")
{
    auto world_size = Commxx::world_size();
    auto world_rank = Commxx::world_rank();

    if (world_size != 4) return;

    std::stringstream ss;
    ss << "commxx_serdes_combo1_" << Commxx::world_size() << "_"
       << Commxx::world_rank() << ".json";

    // create and save
    {
        auto comm1_ptr = std::make_shared<Commxx>(comm_type::world);
        auto const& comm1 = *comm1_ptr;

        auto comm2_ptr = std::make_shared<Commxx>(comm1_ptr->divide(2));
        auto comm3_ptr = std::make_shared<Commxx>(comm2_ptr->divide(1));

        REQUIRE_NOTHROW(json_save(comm3_ptr, ss.str()));
    }

    // reload
    {
        auto comm = std::make_shared<Commxx>();
        REQUIRE_NOTHROW(json_load(comm, ss.str()));

        CHECK(comm->size() == 1);
        CHECK(comm->rank() == 0);
    }
}
