#include "synergia/lattice/dynamic_lattice.h"
#include "synergia/lattice/madx_reader.h"
#include "synergia/lattice/mx_parse.h"

void
Dynamic_lattice::read_madx_file(std::string const& filename)
{
  synergia::parse_madx_file(filename, mx);
}

void
Dynamic_lattice::read_madx_string(std::string const& str)
{
  synergia::parse_madx(str, mx);
}

Lattice
Dynamic_lattice::construct_lattice(std::string const& line_name)
{
  return MadX_reader::get_lattice(line_name, mx);
}

void
Dynamic_lattice::set_variable(std::string const& name, double val)
{
  synergia::mx_expr expr;
  synergia::parse_expression(std::to_string(val), expr);

  mx.insert_variable(name, expr);
}

void
Dynamic_lattice::set_variable(std::string const& name, std::string const& val)
{
  synergia::mx_expr expr;
  synergia::parse_expression(val, expr);

  mx.insert_variable(name, expr);
}

void
Dynamic_lattice::print() const
{
  mx.print();
}
