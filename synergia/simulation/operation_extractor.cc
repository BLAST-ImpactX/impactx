
#include "synergia/simulation/operation_extractor.h"

#include "synergia/simulation/aperture_operation.h"
#include "synergia/simulation/independent_operation.h"

namespace {
  void
  chef_map_operation_extract(
    Lattice const& lattice,
    std::vector<Lattice_element_slice> const& slices,
    std::vector<std::unique_ptr<Independent_operation>>& operations)
  {}

  void
  chef_propagator_operation_extract(
    Lattice const& lattice,
    std::vector<Lattice_element_slice> const& slices,
    std::vector<std::unique_ptr<Independent_operation>>& operations)
  {}

  void
  libff_operation_extract(
    Lattice const& lattice,
    std::vector<Lattice_element_slice> const& slices,
    std::vector<std::unique_ptr<Independent_operation>>& operations)
  {
    operations.push_back(std::make_unique<LibFF_operation>(slices));
  }

} // namespace

void
extract_independent_operations(
  std::string const& extractor_type,
  Lattice const& lattice,
  std::vector<Lattice_element_slice> const& slices,
  std::vector<std::unique_ptr<Independent_operation>>& operations)
{
  if (extractor_type == "chef_map") {
    chef_map_operation_extract(lattice, slices, operations);
  } else if (extractor_type == "chef_propagator") {
    chef_propagator_operation_extract(lattice, slices, operations);
  }
#if 0
    else if (extractor_type == "chef_mixed")
    {
        return chef_mixed_operation_extract(lattice, slices);
    }
#endif
  else if (extractor_type == "libff" || extractor_type == "default") {
    libff_operation_extract(lattice, slices, operations);
  } else {
    throw std::runtime_error("unknown extractor_type: " + extractor_type);
  }
}

std::unique_ptr<Independent_operation>
extract_aperture_operation(std::string const& aperture_type,
                           Lattice_element_slice const& slice)
{
  if (aperture_type == Finite_aperture::type) {
    return std::make_unique<Aperture_operation<Finite_aperture>>(slice);
  } else if (aperture_type == Circular_aperture::type) {
    return std::make_unique<Aperture_operation<Circular_aperture>>(slice);
  } else if (aperture_type == Elliptical_aperture::type) {
    return std::make_unique<Aperture_operation<Elliptical_aperture>>(slice);
  } else if (aperture_type == Rectangular_aperture::type) {
    return std::make_unique<Aperture_operation<Rectangular_aperture>>(slice);
  } else if (aperture_type == Polygon_aperture::type) {
    return std::make_unique<Aperture_operation<Polygon_aperture>>(slice);
  } else {
    // return std::make_unique<Aperture_operation<Dummy_aperture>>(slice);
    throw std::runtime_error("unknown aperture_type " + aperture_type);
  }
}
