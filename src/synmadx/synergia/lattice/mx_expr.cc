#include "mx_expr.h"
#include "madx.h"
#include "mx_parse.h"

#include <cmath>
#include <limits>
#include <stdexcept>

#include <iomanip>
#include <sstream>

using namespace synergia;
using namespace boost;

// mx_calculator
double mx_calculator::nan = std::numeric_limits<double>::quiet_NaN();

double
mx_calculator::operator()(double val) const
{
  return val;
}

double
mx_calculator::operator()(std::string const& ref) const
{
  if (mx == NULL) {
    if (std::isnan(def))
      throw std::runtime_error("Unable to locate reference " + ref);

    return def;
  }

  return mx->variable_as_number(ref, def);
}

double
mx_calculator::operator()(string_pair_t const& ref) const
{
  if (mx == NULL) {
    if (std::isnan(def))
      throw std::runtime_error("Unable to locate reference " + ref.first +
                               "->" + ref.second);

    return def;
  }

  return mx->command(ref.first).attribute_as_number(ref.second, def);
}

double
mx_calculator::operator()(nop_t const& n) const
{
  return boost::apply_visitor(*this, n.expr);
}

double
mx_calculator::operator()(uop_t const& u) const
{
  return u.func.op(boost::apply_visitor(*this, u.param));
}

double
mx_calculator::operator()(bop_t const& b) const
{
  return b.func.op(boost::apply_visitor(*this, b.lhs),
                   boost::apply_visitor(*this, b.rhs));
}

// mx_ref_checker
bool
mx_expr_ref_checker::operator()(double val) const
{
  return false;
}

bool
mx_expr_ref_checker::operator()(std::string const& ref) const
{
  return true;
}

bool
mx_expr_ref_checker::operator()(string_pair_t const& ref) const
{
  return true;
}

bool
mx_expr_ref_checker::operator()(nop_t const& n) const
{
  return boost::apply_visitor(*this, n.expr);
}

bool
mx_expr_ref_checker::operator()(uop_t const& u) const
{
  return boost::apply_visitor(*this, u.param);
}

bool
mx_expr_ref_checker::operator()(bop_t const& b) const
{
  return boost::apply_visitor(*this, b.lhs) ||
         boost::apply_visitor(*this, b.rhs);
}

// mx_expr_writer
std::string
mx_expr_writer::operator()(double val) const
{
  std::stringstream ss;
  ss << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
     << val;
  return ss.str();
}

std::string
mx_expr_writer::operator()(std::string const& ref) const
{
  return ref;
}

std::string
mx_expr_writer::operator()(string_pair_t const& ref) const
{
  return ref.first + "->" + ref.second;
}

std::string
mx_expr_writer::operator()(nop_t const& n) const
{
  if (n.primary) {
    return "(" + boost::apply_visitor(*this, n.expr) + ")";
  } else {
    return boost::apply_visitor(*this, n.expr);
  }
}

std::string
mx_expr_writer::operator()(uop_t const& u) const
{
  switch (u.func.tag) {
    case op_tag::pos: return "+" + apply_visitor(*this, u.param);

    case op_tag::neg: return "-" + apply_visitor(*this, u.param);

    case op_tag::abs: return "abs(" + apply_visitor(*this, u.param) + ")";

    case op_tag::acos: return "acos(" + apply_visitor(*this, u.param) + ")";

    case op_tag::asin: return "asin(" + apply_visitor(*this, u.param) + ")";

    case op_tag::atan: return "atan(" + apply_visitor(*this, u.param) + ")";

    case op_tag::ceil: return "ceil(" + apply_visitor(*this, u.param) + ")";

    case op_tag::cos: return "cos(" + apply_visitor(*this, u.param) + ")";

    case op_tag::cosh: return "cosh(" + apply_visitor(*this, u.param) + ")";

    case op_tag::exp: return "exp(" + apply_visitor(*this, u.param) + ")";

    case op_tag::floor: return "floor(" + apply_visitor(*this, u.param) + ")";

    case op_tag::log: return "log(" + apply_visitor(*this, u.param) + ")";

    case op_tag::log10: return "log10(" + apply_visitor(*this, u.param) + ")";

    case op_tag::sin: return "sin(" + apply_visitor(*this, u.param) + ")";

    case op_tag::sinh: return "sinh(" + apply_visitor(*this, u.param) + ")";

    case op_tag::sqrt: return "sqrt(" + apply_visitor(*this, u.param) + ")";

    case op_tag::tan: return "tan(" + apply_visitor(*this, u.param) + ")";

    case op_tag::tanh: return "tanh(" + apply_visitor(*this, u.param) + ")";

    default:
      throw std::runtime_error("mx_expr_writer(): invalid uop_t operator");
  }
}

std::string
mx_expr_writer::operator()(bop_t const& b) const
{
  switch (b.func.tag) {
    case op_tag::add:
      return apply_visitor(*this, b.lhs) + "+" + apply_visitor(*this, b.rhs);

    case op_tag::sub:
      return apply_visitor(*this, b.lhs) + "-" + apply_visitor(*this, b.rhs);

    case op_tag::mul:
      return apply_visitor(*this, b.lhs) + "*" + apply_visitor(*this, b.rhs);

    case op_tag::div:
      return apply_visitor(*this, b.lhs) + "/" + apply_visitor(*this, b.rhs);

    case op_tag::pow:
      return "pow(" + apply_visitor(*this, b.lhs) + "," +
             apply_visitor(*this, b.rhs) + ")";

    case op_tag::atan2:
      return "atan2(" + apply_visitor(*this, b.lhs) + "," +
             apply_visitor(*this, b.rhs) + ")";

    default:
      throw std::runtime_error("mx_expr_writer(): invalid bop_t operator");
  }
}

// util functions
double
synergia::mx_eval(mx_expr const& expr, double def)
{
  return boost::apply_visitor(mx_calculator(def), expr);
}

double
synergia::mx_eval(mx_expr const& expr, MadX const& mx, double def)
{
  return boost::apply_visitor(mx_calculator(mx, def), expr);
}

std::string
synergia::mx_expr_refstr(mx_expr const& expr)
{
  mx_expr ex = get<nop_t>(get<nop_t>(get<nop_t>(expr).expr).expr).expr;

  if (ex.which() != 1) // string
    throw std::runtime_error("unable to get ref string from mx_expr");

  return get<std::string>(ex);
}

bool
synergia::mx_expr_is_number(mx_expr const& expr)
{
  return not boost::apply_visitor(mx_expr_ref_checker(), expr);
}

std::string
synergia::mx_expr_str(mx_expr const& expr)
{
  return boost::apply_visitor(mx_expr_writer(), expr);
}

bool
synergia::parse_expr(std::string const& s, mx_expr& expr)
{
  return synergia::parse_expression(s, expr);
}
