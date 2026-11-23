#ifndef SYNERGIA_LATTICE_DYNAMIC_LATTICE_H
#define SYNERGIA_LATTICE_DYNAMIC_LATTICE_H

#include "synergia/lattice/madx.h"

class Lattice;

class Dynamic_lattice {

public:
  Dynamic_lattice(/*possible DL options*/) {}

  void read_madx_file(std::string const& filename);
  void read_madx_string(std::string const& str);

  Lattice construct_lattice(std::string const& line_name);

  void set_variable(std::string const& name, double val);
  void set_variable(std::string const& name, std::string const& val);

  void print() const;

private:
  synergia::MadX mx;
};

#endif
