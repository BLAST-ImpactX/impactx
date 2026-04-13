#ifndef SYNERIGA_LATTICE_LATTICE_TREE_H
#define SYNERIGA_LATTICE_LATTICE_TREE_H

#include "synergia/lattice/madx.h"
#include "synergia/lattice/mx_parse.h"

#include "synergia/utils/cereal.h"

class Lattice_tree {
public:
  // default ctor for serialization
  Lattice_tree() : mx() {}

  explicit Lattice_tree(synergia::MadX const& madx) : mx(madx) {}

  // set the value of a variable
  void set_variable(std::string const& name, double val);
  void set_variable(std::string const& name, std::string const& val);

  // set the attribute value of an element
  void set_element_attribute(std::string const& label,
                             std::string const& attr,
                             double val);

  void set_element_attribute(std::string const& label,
                             std::string const& attr,
                             std::string const& val);

  void print() const;

public:
  synergia::MadX mx;

private:
  friend class cereal::access;

  template <class Archive>
  void
  save(Archive& ar) const
  {
    std::string madx = mx.to_madx();
    ar(madx);
  }

  template <class Archive>
  void
  load(Archive& ar)
  {
    std::string madx;
    ar(madx);

    parse_madx(madx, mx);
  }
};

#endif
