#include "lattice_element.h"
#include "lattice.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

// synergia::mx_expr
using namespace synergia;

namespace {
  using type_map_t = std::map<std::string, element_type>;

  const type_map_t type_map = {
    {element_type_name::generic, element_type::generic},
    {element_type_name::drift, element_type::drift},
    {element_type_name::rbend, element_type::rbend},
    {element_type_name::sbend, element_type::sbend},
    {element_type_name::quadrupole, element_type::quadrupole},
    {element_type_name::multipole, element_type::multipole},
    {element_type_name::rfcavity, element_type::rfcavity},
    {element_type_name::hkicker, element_type::hkicker},
    {element_type_name::vkicker, element_type::vkicker},
    {element_type_name::kicker, element_type::kicker},
    {element_type_name::monitor, element_type::monitor},
    {element_type_name::hmonitor, element_type::hmonitor},
    {element_type_name::vmonitor, element_type::vmonitor},
    {element_type_name::sextupole, element_type::sextupole},
    {element_type_name::octupole, element_type::octupole},
    {element_type_name::marker, element_type::marker},
    {element_type_name::instrument, element_type::instrument},
    {element_type_name::rcollimator, element_type::rcollimator},
    {element_type_name::nllens, element_type::nllens},
    {element_type_name::solenoid, element_type::solenoid},
    {element_type_name::elens, element_type::elens},
    {element_type_name::foil, element_type::foil},
    {element_type_name::dipedge, element_type::dipedge},
    {element_type_name::matrix, element_type::matrix},
  };

  element_type
  find_type(std::string const& stype)
  {
    auto r = type_map.find(stype);
    if (r == type_map.end())
      throw std::runtime_error("invalid element type " + stype);
    return r->second;
  }
}

Lattice_element::Lattice_element()
  : name("")
  , format(element_format::madx)
  , stype("generic")
  , type(find_type(stype))
  , ancestors()
  , string_attributes()
  , lazy_double_attributes()
  , lazy_vector_attributes()
  , length_attribute_name("l")
  , bend_angle_attribute_name("angle")
  , revision(0)
  , lattice_ptr(nullptr)
  , markers{}
{}

Lattice_element::Lattice_element(std::string const& type,
                                 std::string const& name,
                                 element_format format)
  : name(name)
  , format(format)
  , stype(type)
  , type(find_type(stype))
  , ancestors()
  , string_attributes()
  , lazy_double_attributes()
  , lazy_vector_attributes()
  , length_attribute_name("l")
  , bend_angle_attribute_name("angle")
  , revision(0)
  , lattice_ptr(nullptr)
  , markers{}
{}

Lattice_element::Lattice_element(Lsexpr const& lsexpr)
  : name("")
  , format(element_format::madx)
  , stype("generic")
  , type(find_type(stype))
  , ancestors()
  , string_attributes()
  , lazy_double_attributes()
  , lazy_vector_attributes()
  , length_attribute_name("l")
  , bend_angle_attribute_name("angle")
  , revision(0)
  , lattice_ptr(nullptr)
  , markers{}
{
  using namespace synergia;

  for (auto const& lse : lsexpr) {
    if (lse.is_labeled()) {
      if (lse.get_label() == "type") {
        stype = lse.get_string();
        type = find_type(stype);
      } else if (lse.get_label() == "name") {
        name = lse.get_string();
      } else if (lse.get_label() == "ancestors") {
        auto ancestors_vector = lse.get_string_vector();
        std::copy(ancestors_vector.begin(),
                  ancestors_vector.end(),
                  std::back_inserter(ancestors));
      } else if (lse.get_label() == "double_attributes") {
        for (auto const& attr : lse)
          lazy_double_attributes[attr.get_label()] = mx_expr(attr.get_double());
      } else if (lse.get_label() == "string_attributes") {
        for (auto const& attr : lse)
          string_attributes[attr.get_label()] = attr.get_string();
      } else if (lse.get_label() == "vector_attributes") {
        for (auto const& attr : lse) {
          std::vector<double> vd = attr.get_double_vector();

          std::vector<mx_expr> ve;
          for (double d : vd) ve.push_back(mx_expr(d));

          lazy_vector_attributes[attr.get_label()] = ve;
        }
      }
    } else {
      if (!lse.is_atomic()) {
        for (auto const& attr : lse)
          lazy_double_attributes[attr.get_label()] = mx_expr(attr.get_double());
      }
    }
  }
}

Lsexpr
Lattice_element::as_lsexpr() const
{

    Lsexpr retval;
#if 0 // as_lsexpr() might come back some day
    retval.push_back(Lsexpr(type, "type"));
    retval.push_back(Lsexpr(name, "name"));
    if (double_attributes.size() > 0) {
        Lsexpr attrs;
        for (std::map<std::string, double>::const_iterator it =
                 double_attributes.begin();
             it != double_attributes.end(); ++it) {
            Lsexpr attr(it->second);
            attr.set_label(it->first);
            attrs.push_back(attr);
        }
        if (!((string_attributes.size() == 0) &&
              (vector_attributes.size() == 0))) {
            attrs.set_label("double_attributes");
        }
    }
    if (double_attributes.size() > 0) {
        Lsexpr attrs;
        for (std::map<std::string, double>::const_iterator it =
                 double_attributes.begin();
             it != double_attributes.end(); ++it) {
            Lsexpr attr(it->second);
            attr.set_label(it->first);
            attrs.push_back(attr);
        }
        if (!((string_attributes.size() == 0) &&
              (vector_attributes.size() == 0))) {
            attrs.set_label("double_attributes");
        }
        retval.push_back(attrs);
    }
    if (string_attributes.size() > 0) {
        Lsexpr attrs;
        for (std::map<std::string, std::string>::const_iterator it =
                 string_attributes.begin();
             it != string_attributes.end(); ++it) {
            Lsexpr attr(it->second);
            attr.set_label(it->first);
            attrs.push_back(attr);
        }
        attrs.set_label("string_attributes");
        retval.push_back(attrs);
    }
    if (vector_attributes.size() > 0) {
        Lsexpr attrs;
        for (std::map<std::string, std::vector<double> >::const_iterator it =
                 vector_attributes.begin();
             it != vector_attributes.end(); ++it) {
            Lsexpr attr(it->second);
            attr.set_label(it->first);
            attrs.push_back(attr);
        }
        attrs.set_label("vector_attributes");
        retval.push_back(attrs);
    }
    if (ancestors.size() > 0) {
        std::vector<std::string> ancestors_vector;
        std::copy(ancestors.begin(), ancestors.end(),
                  std::back_inserter(ancestors_vector));
        retval.push_back(Lsexpr(ancestors_vector, "ancestors"));
    }
#endif
  return retval;
}

std::vector<std::string>
Lattice_element::get_all_type_names()
{
  static std::vector<std::string> names;

  if (names.empty())
    for (auto const& type : type_map) names.push_back(type.first);

  return names;
}

std::string const&
Lattice_element::get_type_name() const
{
  return stype;
}

element_type
Lattice_element::get_type() const
{
  return type;
}

std::string const&
Lattice_element::get_name() const
{
  return name;
}

element_format
Lattice_element::get_format() const
{
  return format;
}

void
Lattice_element::add_ancestor(std::string const& ancestor)
{
  ancestors.push_back(ancestor);
}

std::list<std::string> const&
Lattice_element::get_ancestors() const
{
  return ancestors;
}

void
Lattice_element::duplicate_attribute(std::string const& name,
                                     std::string const& new_name,
                                     bool overwrite)
{
  auto duplicator = [](auto& attrs,
                       std::string const& name,
                       std::string const& new_name,
                       bool overwrite) {
    if (attrs.find(new_name) != attrs.end() && !overwrite)
      throw std::runtime_error("Lattice_element::duplicate_attribute(): "
                               "the target attribute " +
                               new_name + " already exists");

    if (attrs.find(name) != attrs.end()) attrs[new_name] = attrs[name];
  };

  duplicator(lazy_double_attributes, name, new_name, overwrite);
  duplicator(lazy_vector_attributes, name, new_name, overwrite);
  duplicator(string_attributes, name, new_name, overwrite);
}

void
Lattice_element::delete_attribute(std::string const& name)
{
  lazy_double_attributes.erase(name);
  lazy_vector_attributes.erase(name);
  string_attributes.erase(name);
}

void
Lattice_element::rename_attribute(std::string const& name,
                                  std::string const& new_name)
{
  duplicate_attribute(name, new_name);
  delete_attribute(name);
}

void
Lattice_element::copy_attributes_from(Lattice_element const& o)
{
  lazy_double_attributes = o.lazy_double_attributes;
  lazy_vector_attributes = o.lazy_vector_attributes;
  string_attributes = o.string_attributes;
}

void
Lattice_element::set_double_attribute(std::string const& name,
                                      double value,
                                      bool increment_revision)
{
  lazy_double_attributes[name] = mx_expr(value);
  if (increment_revision) ++revision;
}

void
Lattice_element::set_double_attribute(std::string const& name,
                                      mx_expr const& value,
                                      bool increment_revision)
{
  lazy_double_attributes[name] = value;
  if (increment_revision) ++revision;
}

void
Lattice_element::set_double_attribute(std::string const& name,
                                      std::string const& value,
                                      bool increment_revision)
{
  mx_expr expr;

  if (!parse_expr(value, expr)) {
    throw std::runtime_error("Lattice_element::set_double_attribute() "
                             "cannot parse the value " +
                             value +
                             " into "
                             "a valid expression");
  }

  lazy_double_attributes[name] = expr;
  if (increment_revision) ++revision;
}

void
Lattice_element::set_double_attribute(std::string const& name,
                                      const char* value,
                                      bool increment_revision)
{
  set_double_attribute(name, std::string(value), increment_revision);
}

void
Lattice_element::set_default_double_attribute(std::string const& name,
                                              double value,
                                              bool increment_revision)
{
  if (!has_double_attribute(name)) {
    lazy_double_attributes[name] = mx_expr(value);
    if (increment_revision) ++revision;
  }
}

void
Lattice_element::remove_double_attribute(std::string const& name)
{
  if (has_double_attribute(name)) {
    lazy_double_attributes.erase(name);
  }
}

bool
Lattice_element::has_double_attribute(std::string const& name) const
{
  return lazy_double_attributes.count(name) > 0;
}

double
Lattice_element::get_double_attribute(std::string const& name) const
{
  auto lr = lazy_double_attributes.find(name);

  // no such attribute
  if (lr == lazy_double_attributes.end()) {
    throw std::runtime_error("Lattice_element::get_double_attribute: element " +
                             this->name + " of type " + stype +
                             " has no double attribute '" + name + "'");
  }

  // found the attribute, now evaluate either with or without
  // dynamic lattice
  //
  // default the references to 0.0 when evaluating the lazy value.
  //
  // this will set undefined variables to 0.0 instead of throwing
  // an excpetion for "undefined reference".
  //
  if (lattice_ptr && lattice_ptr->is_dynamic_lattice()) {
    auto const& mx = lattice_ptr->get_lattice_tree().mx;
    return mx_eval(lr->second, mx, 0.0);
  } else {
    return mx_eval(lr->second, 0.0);
  }
}

double
Lattice_element::get_double_attribute(std::string const& name, double val) const
{
  auto lr = lazy_double_attributes.find(name);

  // no such attribute, returns the def value
  if (lr == lazy_double_attributes.end()) return val;

  // default the references to 0.0 when evaluating the lazy value.
  if (lattice_ptr && lattice_ptr->is_dynamic_lattice()) {
    auto const& mx = lattice_ptr->get_lattice_tree().mx;
    return mx_eval(lr->second, mx, 0.0);
  } else {
    return mx_eval(lr->second, 0.0);
  }
}

std::vector<std::string>
Lattice_element::get_double_attribute_names() const
{
  std::vector<std::string> names;

  for (auto const& attr : lazy_double_attributes) names.push_back(attr.first);

  return names;
}

std::map<std::string, double>
Lattice_element::get_double_attributes() const
{
  std::cout << "WARNING: it is not recommended to keep using "
               "Lattice_element::get_double_attributes() method. For the "
               "purposes of copying double attributes from one element "
               "to another, please use Lattice_element::copy_attributes_from() "
               "instead. For exploring individual double attributes, please "
               "call the Lattice_element::get_double_attribute(name). For "
               "iterating through the double attribute keys, please use "
               "Lattice_element::get_double_attribute_names().";

  std::map<std::string, double> attrs;

  for (auto const& attr : lazy_double_attributes) {
    if (lattice_ptr && lattice_ptr->is_dynamic_lattice()) {
      auto const& mx = lattice_ptr->get_lattice_tree().mx;
      attrs[attr.first] = mx_eval(attr.second, mx, 0.0);
    } else {
      attrs[attr.first] = mx_eval(attr.second, 0.0);
    }
  }

  return attrs;
}

void
Lattice_element::set_string_attribute(std::string const& name,
                                      std::string const& value,
                                      bool increment_revision)
{
  string_attributes[name] = value;
  if (increment_revision) ++revision;
}

void
Lattice_element::remove_string_attribute(std::string const& name)
{
  if (has_string_attribute(name)) {
    string_attributes.erase(name);
  }
}

void
Lattice_element::set_default_string_attribute(std::string const& name,
                                              std::string const& value,
                                              bool increment_revision)
{
  if (!has_string_attribute(name)) {
    string_attributes[name] = value;
    if (increment_revision) ++revision;
  }
}

bool
Lattice_element::has_string_attribute(std::string const& name) const
{
  bool retval = (string_attributes.count(name) > 0);
  return retval;
}

std::string const&
Lattice_element::get_string_attribute(std::string const& name) const
{
  auto r = string_attributes.find(name);
  if (r == string_attributes.end()) {
    throw std::runtime_error("Lattice_element::get_string_attribute: element " +
                             this->name + " of type " + stype +
                             " has no string attribute '" + name + "'");
  } else {
    return r->second;
  }
}

std::string const&
Lattice_element::get_string_attribute(std::string const& name,
                                      std::string const& val) const
{
  auto r = string_attributes.find(name);
  if (r == string_attributes.end())
    return val;
  else
    return r->second;
}

void
Lattice_element::set_vector_attribute(std::string const& name,
                                      std::vector<double> const& value,
                                      bool increment_revision)
{
  std::vector<mx_expr> ve;
  for (double d : value) ve.emplace_back(d);

  lazy_vector_attributes[name] = ve;
  if (increment_revision) ++revision;
}

bool
Lattice_element::has_vector_attribute(std::string const& name) const
{
  bool retval = (lazy_vector_attributes.count(name) > 0);
  return retval;
}

void
Lattice_element::remove_vector_attribute(std::string const& name)
{
  if (has_vector_attribute(name)) {
    lazy_vector_attributes.erase(name);
  }
}

std::vector<double>
Lattice_element::get_vector_attribute(std::string const& name) const
{
  auto r = lazy_vector_attributes.find(name);

  // no such attribute
  if (r == lazy_vector_attributes.end()) {
    throw std::runtime_error("Lattice_element::get_vector_attribute: element " +
                             this->name + " of type " + stype +
                             " has no vector attribute '" + name + "'");
  }

  // the returned array
  std::vector<double> vd;

  // default the references to 0.0 when evaluating the lazy value.
  if (lattice_ptr && lattice_ptr->is_dynamic_lattice()) {
    auto const& mx = lattice_ptr->get_lattice_tree().mx;
    for (auto const& expr : r->second) vd.push_back(mx_eval(expr, mx, 0.0));
  } else {
    for (auto const& expr : r->second) vd.push_back(mx_eval(expr, 0.0));
  }

  return vd;
}

std::vector<double>
Lattice_element::get_vector_attribute(std::string const& name,
                                      std::vector<double> const& val) const
{
  auto r = lazy_vector_attributes.find(name);

  // no such attribute, returns the default val
  if (r == lazy_vector_attributes.end()) return val;

  // the double array to be returned
  std::vector<double> vd;

  // default the references to 0.0 when evaluating the lazy values.
  if (lattice_ptr && lattice_ptr->is_dynamic_lattice()) {
    auto const& mx = lattice_ptr->get_lattice_tree().mx;
    for (auto const& expr : r->second) vd.push_back(mx_eval(expr, mx, 0.0));
  } else {
    for (auto const& expr : r->second) vd.push_back(mx_eval(expr, 0.0));
  }

  return vd;
}

void
Lattice_element::set_length_attribute_name(std::string const& attribute_name)
{
  length_attribute_name = attribute_name;
}

void
Lattice_element::set_bend_angle_attribute_name(
  std::string const& attribute_name)
{
  bend_angle_attribute_name = attribute_name;
}

double
Lattice_element::get_length() const
{
  return get_double_attribute(length_attribute_name, 0.0);
}

double
Lattice_element::get_bend_angle() const
{
  return get_double_attribute(bend_angle_attribute_name, 0.0);
}

long int
Lattice_element::get_revision() const
{
  return revision;
}

bool
Lattice_element::has_lattice() const
{
  return (lattice_ptr != 0);
}

void
Lattice_element::set_lattice(Lattice& lattice)
{
  lattice_ptr = &lattice;
}

Lattice const&
Lattice_element::get_lattice() const
{
  if (!has_lattice()) {
    throw std::runtime_error(
      "Lattice_element::get_lattice: element not part of any lattice");
  }
  return *lattice_ptr;
}

namespace {
  void
  validate_tunes_corrector(Lattice_element const& e)
  {
    // TODO: quads, CFbends, multipoles
    // for now, quads without tilt or skew
    // throw if not a valid corrector
  }

  void
  validate_chrom_corrector(Lattice_element const& e)
  {
    // TODO: sexts, CFbends, multipoles
    // for now, sextupole and thin-sextupoles without tilt or skew
    // throw if not a valid corrector
  }
}

void
Lattice_element::set_marker(marker_type t)
{
  switch (t) {
    case marker_type::h_tunes_corrector:

      if (has_marker(marker_type::v_tunes_corrector))
        throw std::runtime_error(
          "Lattice_element::set_marker(): "
          "v_tunes_corrector has been set for the element " +
          name + ", unable to set the h_tunes_corrector.");

      validate_tunes_corrector(*this);

      break;

    case marker_type::v_tunes_corrector:

      if (has_marker(marker_type::h_tunes_corrector))
        throw std::runtime_error(
          "Lattice_element::set_marker(): "
          "h_tunes_corrector has been set for the element " +
          name + ", unable to set the v_tunes_corrector.");

      validate_tunes_corrector(*this);

      break;

    case marker_type::h_chrom_corrector:

      if (has_marker(marker_type::v_chrom_corrector))
        throw std::runtime_error(
          "Lattice_element::set_marker(): "
          "v_chrom_corrector has been set for the element " +
          name + ", unable to set the h_chrom_corrector.");

      validate_chrom_corrector(*this);

      break;

    case marker_type::v_chrom_corrector:

      if (has_marker(marker_type::h_chrom_corrector))
        throw std::runtime_error(
          "Lattice_element::set_marker(): "
          "h_chrom_corrector has been set for the element " +
          name + ", unable to set the v_chrom_corrector.");

      validate_chrom_corrector(*this);

      break;

    default: break;
  }

  // set the marker
  markers[(int)t] = true;
}

std::string
Lattice_element::as_string() const
{
  std::stringstream sstream;
  sstream.precision(15);

  for (auto const& ancestor : ancestors) sstream << ancestor << ":";

  sstream << " " << stype << " ";
  sstream << name << ": ";
  bool first_attr = true;

  for (auto const& attr : lazy_double_attributes) {
    if (first_attr)
      first_attr = false;
    else
      sstream << ", ";

    sstream << attr.first << "=" << mx_expr_str(attr.second);
  }

  for (auto const& attr : string_attributes) {
    if (first_attr)
      first_attr = false;
    else
      sstream << ", ";

    sstream << attr.first << "=" << attr.second;
  }

  for (auto const& attr : lazy_vector_attributes) {
    if (first_attr)
      first_attr = false;
    else
      sstream << ", ";

    sstream << attr.first << "={";
    for (size_t i = 0; i < attr.second.size(); ++i) {
      if (i) sstream << ", ";
      sstream << mx_expr_str((attr.second)[i]);
    }
    sstream << "}";
  }

  return sstream.str();
}

void
Lattice_element::print() const
{
  std::cout << as_string() << std::endl;
}

std::unordered_set<std::string> const non_madx_double_attributes = {
  "entry_edge_kick", "exit_edge_kick", "kicks", "yoshida_order", "kl",
  "a1", "a2", "a3", "a4", "a5", "a6", "a7",
  "b1", "b2", "b3", "b4", "b5", "b6", "b7"
};

std::unordered_set<std::string> const non_madx_string_attributes = {
  "extractor_type", "propagator_type"
};

std::string
Lattice_element::as_madx(bool sanitize) const
{
  std::stringstream ss;
  ss << name << ": " << stype;

  for (auto const& attr : lazy_double_attributes) {
    if (sanitize && non_madx_double_attributes.count(attr.first)) {
      continue;  // skip this attribute if it is not madx kosher
    }
    ss << ", " << attr.first << "=" << std::setprecision(16) << mx_expr_str(attr.second);
  }

  for (auto const& attr : string_attributes) {
    if (sanitize && non_madx_string_attributes.count(attr.first)) {
      continue; // skip this attriute if it is not madx kosher
    }
    ss << ", " << attr.first << "=" << attr.second;
  }

  for (auto const& attr : lazy_vector_attributes) {
    ss << ", " << attr.first << "={";

    for (int i = 0; i < attr.second.size(); ++i) {
      if (i) ss << ", ";
      ss << std::setprecision(16) << mx_expr_str(attr.second[i]);
    }

    ss << "}";
  }

  ss << ";";

  return ss.str();
}
