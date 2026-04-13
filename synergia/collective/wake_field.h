#ifndef WAKE_FIELD_H_
#define WAKE_FIELD_H_

#include "synergia/utils/kokkos_views.h"

#include <mpi.h>

#include <string>
#include <vector>

struct Wake_field {
  const std::string wake_file;
  const std::string wake_type;

  /// assume the  wake functions are stored using a quadratic grid
  ///  z[i]= (i-istart)^2*delta_z+zstart for i>istart
  ///  z[i]= -(i-istart)^2*delta_z+zstart  for i< istart
  int istart;
  double zstart;
  double delta_z;

  // number of terms
  int size_wake;

   // z_coord, z_wake, xw_lead, xw_trail, yw_lead, yw_trail
  // all in a single buffer.
  //
  // total size = size_wake * 6. Fortran ordering.
  //
  // x/yw_lead: wake term proportional with the displacement of
  //     the leading (source) particle
  //
  // x/yw_trail: wake term proportional with the displacement of
  //     the trailing (affected) particle
  //
  karray1d_dev terms;
  karray1d_hst h_terms;

#if 0
    // wake terms
    karray1d_dev z_coord;
    karray1d_dev z_wake;

    // wake term proportional with the displacement of the leading (source) particle
    karray1d_dev xw_lead;

    // wake term proportional with the displacement of the trailing (affected) particle
    karray1d_dev xw_trail;

    // wake term proportional with the displacement of the leading (source) particle
    karray1d_dev yw_lead;

    // wake term proportional with the displacement of the trail particle
    karray1d_dev yw_trail;

    // host mirrors
    karray1d_hst h_z_coord;
    karray1d_hst h_z_wake;
    karray1d_hst h_xw_lead;
    karray1d_hst h_xw_trail;
    karray1d_hst h_yw_lead;
    karray1d_hst h_yw_trail;
#endif

  Wake_field(std::string const& wake_file, std::string const& wake_type,
    const double& mwf_xlead,
    const double& mwf_xtrail,
    const double& mwf_ylead,
    const double& mwf_ytrail,
    const double& mwf_zwake);

  void multiply_xw_lead(double mltp);
  void multiply_xw_trail(double mltp);
  void multiply_yw_lead(double mltp);
  void multiply_yw_trail(double mltp);
  void multiply_z_wake(double mltp);
};

#endif /* WAKE_FIELD_H_ */
