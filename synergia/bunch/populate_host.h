#ifndef BUNCH_POPULATE_HOST_H
#define BUNCH_POPULATE_HOST_H

#include <array>

void adjust_moments_host(double const* means,
                         double const* covariances,
                         double const* bunch_mean,
                         double const* bunch_mom2,
                         int num_particles,
                         int num_particles_slots,
                         double* particles);

void get_correlation_matrix_host(double* correlation_matrix,
                                 double const* one_turn_map,
                                 double arms,
                                 double brms,
                                 double crms,
                                 double beta,
                                 std::array<int, 3> const& rms_index);

#endif
