#ifndef KOKKOS_TOOLS_H
#define KOKKOS_TOOLS_H

#include "synergia/utils/kokkos_types.h"
#include "synergia/utils/kokkos_views.h"
#include "synergia/utils/logger.h"

namespace ku {
  struct alg_zeroer {
    karray1d_dev arr;

    KOKKOS_INLINE_FUNCTION
    void
    operator()(const int i) const
    {
      arr(i) = 0.0;
    }
  };

  inline void
  zero_karray(karray1d_dev const& arr)
  {
    alg_zeroer alg{arr};
    Kokkos::parallel_for(arr.extent(0), alg);
  }

  inline void
  print_arr_sum(Logger& logger,
                karray1d_hst const& arr,
                int offset = 0,
                int length = -1)
  {
    double sum = 0;
    for (int i = offset; i < offset + length; ++i) sum += arr(i);

    logger << std::setprecision(12);
    logger << "      " << arr.label() << " = " << sum << "\n";
  }

  inline void
  print_arr_sum(Logger& logger,
                karray1d_dev const& arr,
                int offset = 0,
                int length = -1)
  {
    if (length < 0) length = arr.extent(0);

    karray1d_hst harr = Kokkos::create_mirror_view(arr);
    Kokkos::deep_copy(harr, arr);

    print_arr_sum(logger, harr, offset, length);
  }

}

#endif
