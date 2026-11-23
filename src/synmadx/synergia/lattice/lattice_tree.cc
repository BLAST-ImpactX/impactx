#include "synergia/lattice/lattice_tree.h"
#include "synergia/lattice/mx_parse.h"

using namespace synergia;

void
Lattice_tree::set_variable(std::string const& name, double val)
{
  mx.insert_variable(name, mx_expr(val));
}

void
Lattice_tree::set_variable(std::string const& name, std::string const& val)
{
  synergia::mx_expr expr;
  synergia::parse_expression(val, expr);

  mx.insert_variable(name, expr);
}

void
Lattice_tree::set_element_attribute(std::string const& label,
                                    std::string const& attr,
                                    double val)
{
  mx.command_ref(label).insert_attribute(attr, mx_expr(val));
}

void
Lattice_tree::set_element_attribute(std::string const& label,
                                    std::string const& attr,
                                    std::string const& val)
{
  synergia::mx_expr expr;
  synergia::parse_expression(val, expr);

  mx.command_ref(label).insert_attribute(attr, expr);
}

void
Lattice_tree::print() const
{
  mx.print();
}
