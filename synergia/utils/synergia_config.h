#ifndef SYNERGIA_CONFIG_H_
#define SYNERGIA_CONFIG_H_

#include <climits>
#include <cstdint>

#if defined SYNERGIA_HAVE_OPENPMD
#include <openPMD/openPMD.hpp>
#else
#include "synergia/utils/hdf5_file.h"
#endif

#ifdef SYNERGIA_HAVE_OPENPMD
using io_device = openPMD::Series;
#else
using io_device = Hdf5_file;
#endif

#if SIZE_MAX == UCHAR_MAX
#define MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
#define MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
#define MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
#define MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
#define MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
#error "Unable to define MPI_SIZE_T!"
#endif

#endif
