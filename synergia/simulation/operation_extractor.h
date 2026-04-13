#ifndef OPERATION_EXTRACTOR_H_
#define OPERATION_EXTRACTOR_H_

#include <memory>
#include <string>
#include <vector>

class Independent_operation;
class Lattice;
class Lattice_element_slice;

void extract_independent_operations(
  std::string const& extractor_type,
  Lattice const& lattice,
  std::vector<Lattice_element_slice> const& slices,
  std::vector<std::unique_ptr<Independent_operation>>& operations);

std::unique_ptr<Independent_operation> extract_aperture_operation(
  std::string const& aperture_type,
  Lattice_element_slice const& slice);

#endif /* OPERATION_EXTRACTOR_H_ */
