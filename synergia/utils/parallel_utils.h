#ifndef PARALLEL_UTILS_H_
#define PARALLEL_UTILS_H_

#include <vector>

#include "commxx.h"

/// See decompose_1d_raw. The number of processors is extracted from the
/// Commxx object.
void decompose_1d(Commxx const& comm,
                  int length,
                  std::vector<int>& offsets,
                  std::vector<int>& counts);

/// See decompose_1d_raw. In this simpler version, the number of processors is
/// extracted from the Commxx object and the number of counts assigned to the
/// current processor is returned. The lengths on other processors and offsets
/// are discarded.
int decompose_1d_local(Commxx const& comm, int length);

/// Distribute elements across multiple processors. Where decompose assumes that
/// each element will be assigned to a single processor, distribute assumes that
/// elements may be shared among multiple processors.  The return value is a
/// vector of length elements, each member of which is a vector containing the
/// ranks for each element.
/// @param processors is the number of processors.
/// @param elements is the number of elements.
/// Example1: elements=5, processors=3... vector[0]=0, vector[1]=1, vector[2]=1,
///                                       vector[3]=2, vector[4]=2
/// Example2: elements=5, processors=12... vector[0]=(0,1) vector[1]=(2,3),
/// vector[2]=(4,5)
///                                       vector[3]=(6,7,8) vector[4]=(9,10,11)
std::vector<std::vector<int>> distribute_1d_raw(int processors, int elements);

/// See distribute_1d_raw. The number of processors is extracted from the
/// Commxx object.
std::vector<std::vector<int>> distribute_1d(Commxx const& comm, int elements);

void counts_and_offsets_for_impedance_raw(unsigned int processors,
                                          int length,
                                          std::vector<int>& offsets,
                                          std::vector<int>& counts);

void counts_and_offsets_for_impedance(Commxx const& comm,
                                      int length,
                                      std::vector<int>& offsets,
                                      std::vector<int>& counts);

#endif /* PARALLEL_UTILS_H_ */
