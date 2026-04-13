#include <catch2/catch_test_macros.hpp>

#include "Kokkos_Core.hpp"

TEST_CASE("kokkos_padding", "[BunchParticles]")
{
    auto alloc = Kokkos::view_alloc("test", Kokkos::AllowPadding);
    auto darr = Kokkos::View<double**, Kokkos::LayoutLeft>(alloc, 11, 7);
    auto harr = Kokkos::create_mirror_view(darr);

    INFO("This is a test for Kokkos::AllowPadding issue on the "
         "device memory. If it fails, this means the issue "
         " has been fixed in Kokkos.");

    REQUIRE(darr.stride(0) == harr.stride(0));
    REQUIRE(darr.stride(1) != harr.stride(1));

    CHECK_THROWS(Kokkos::deep_copy(harr, darr));
}
