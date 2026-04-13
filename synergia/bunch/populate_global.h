#ifndef POPULATE_GLOBAL_H
#define POPULATE_GLOBAL_H

#include "synergia/bunch/populate.h"

void populate_global_6d(uint64_t seed,
                        Bunch& bunch,
                        const_karray1d means,
                        const_karray2d_row covariances);

void populate_global_6d_truncated(uint64_t seed,
                                  Bunch& bunch,
                                  const_karray1d means,
                                  const_karray2d_row covariances,
                                  const_karray1d limits);

#endif
