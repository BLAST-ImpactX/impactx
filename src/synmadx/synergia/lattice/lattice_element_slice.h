#ifndef LATTICE_ELEMENT_SLICE_H_
#define LATTICE_ELEMENT_SLICE_H_

#include "synergia/lattice/lattice_element.h"

class Lattice_element_slice {
private:
  Lattice_element const* element;

  bool whole;
  bool left_edge;
  bool right_edge;

  double left;
  double right;
  double ref_ct;

public:
  // lattice functions
  latt_func_t lf;

public:
  explicit Lattice_element_slice(Lattice_element const& element);

  Lattice_element_slice(Lattice_element const& element,
                        double left,
                        double right);

  bool
  is_whole() const
  {
    return whole;
  }
  bool
  has_left_edge() const
  {
    return left_edge;
  }
  bool
  has_right_edge() const
  {
    return right_edge;
  }

  double
  get_left() const
  {
    return left;
  }
  double
  get_right() const
  {
    return right;
  }

  double
  get_length() const
  {
    return right - left;
  }

  void
  set_reference_ct(double ct)
  {
    ref_ct = ct;
  }
  double
  get_reference_ct() const
  {
    return ref_ct;
  }

  Lattice_element const&
  get_lattice_element() const
  {
    return *element;
  }

  std::string as_string() const;

  void print() const;
};

#endif /* LATTICE_ELEMENT_SLICE_H_ */
