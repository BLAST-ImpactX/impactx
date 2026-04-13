#include <Kokkos_Core.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

// For debugging, set to 1.
#define IDXPRINT 0

KOKKOS_INLINE_FUNCTION
void
get_leftmost_indices_offset(double pos,
                            double left,
                            double inv_cell_size,
                            int& idx,
                            double& off)
{
    double scaled_location = (pos - left) * inv_cell_size - 0.5;
    // The static cast is likely spurious, but it might be better to
    // have it be optimized away rather than not having it here.
    idx = static_cast<int>(Kokkos::floor(scaled_location));
    off = scaled_location - idx;

#if IDXPRINT == 1
    std::cout << "particle with pos : " << pos << ", left : " << left
              << ", inv_cell_size : " << inv_cell_size << ", idx : " << idx
              << ", off : " << off << std::endl;
#endif

    return;
}
