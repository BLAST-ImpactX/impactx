#ifndef LATTICE_ELEMENT_H_
#define LATTICE_ELEMENT_H_

#include <array>
#include <list>
#include <map>
#include <string>
#include <vector>

#include "synergia/lattice/mx_expr.h"
#include "synergia/utils/cereal.h"
#include "synergia/utils/lsexpr.h"

enum class element_format {
  mad8,
  madx,
};

enum class element_type {
  generic,

  drift,
  rbend,
  sbend,
  quadrupole,
  multipole,
  rfcavity,

  hkicker,
  vkicker,
  kicker,

  sextupole,
  octupole,
  monitor,
  hmonitor,
  vmonitor,
  marker,
  instrument,
  rcollimator,

  nllens,
  solenoid,
  elens,

  foil,
  dipedge,
  matrix
};

namespace element_type_name {
  constexpr char const* generic = "generic";
  constexpr char const* drift = "drift";
  constexpr char const* rbend = "rbend";
  constexpr char const* sbend = "sbend";
  constexpr char const* quadrupole = "quadrupole";
  constexpr char const* multipole = "multipole";
  constexpr char const* rfcavity = "rfcavity";
  constexpr char const* hkicker = "hkicker";
  constexpr char const* vkicker = "vkicker";
  constexpr char const* kicker = "kicker";
  constexpr char const* monitor = "monitor";
  constexpr char const* hmonitor = "hmonitor";
  constexpr char const* vmonitor = "vmonitor";
  constexpr char const* sextupole = "sextupole";
  constexpr char const* octupole = "octupole";
  constexpr char const* marker = "marker";
  constexpr char const* instrument = "instrument";
  constexpr char const* rcollimator = "rcollimator";
  constexpr char const* nllens = "nllens";
  constexpr char const* solenoid = "solenoid";
  constexpr char const* elens = "elens";
  constexpr char const* foil = "foil";
  constexpr char const* dipedge = "dipedge";
  constexpr char const* matrix = "matrix";
}

enum class marker_type {
  h_tunes_corrector,
  v_tunes_corrector,

  h_chrom_corrector,
  v_chrom_corrector,

  end
};

// Lattice Functions
struct latt_func_t {
  struct lf_val_t {
    double hor;
    double ver;
  };

  double arcLength;

  lf_val_t dispersion;
  lf_val_t dPrime;
  lf_val_t beta;
  lf_val_t alpha;
  lf_val_t psi;
};

// forward declaration
class Lattice;

/// The Lattice_element class contains the description of a single
/// lattice element. Each element has a name, a (string) type and
/// dictionaries of named double and string attributes.
/// Lattice structure is described by a list of ancestors stored in
/// an element.
class Lattice_element {

private:
  std::string name;
  element_format format;

  std::string stype;
  element_type type;

  std::list<std::string> ancestors;

  // std::map<std::string, double> double_attributes;
  std::map<std::string, std::string> string_attributes;
  // std::map<std::string, std::vector<double>> vector_attributes;

  // lazy attributes
  std::map<std::string, synergia::mx_expr> lazy_double_attributes;
  std::map<std::string, std::vector<synergia::mx_expr>> lazy_vector_attributes;

  std::string length_attribute_name;
  std::string bend_angle_attribute_name;

  long int revision;

  Lattice* lattice_ptr;

  // marked as mutable because this attribute is not a lattice
  // intrinsic attribute, but an attribute serves as the result of
  // bunch propagation through the lattice element (aperture
  // operation to be specific).
  // the constness of a lattice element is that the intrinsic
  // attributes of the element (e.g., the strength and length)
  // is unaffected by the bunch propagation
  mutable double deposited_charge = 0.0;

  // markers
  std::array<bool, (int)marker_type::end> markers;

public:
  // lattice functions
  latt_func_t lf;

public:
  // get all element type names
  static std::vector<std::string> get_all_type_names();

public:
  /// Construct a Lattice_element with an empty name and type.
  /// for serialization only
  Lattice_element();

  /// Construct a Lattice_element.
  /// @param name name
  /// @param type type
  Lattice_element(std::string const& type,
                  std::string const& name,
                  element_format format = element_format::madx);

  /// Construct a Lattice_element from the Lsexpr representation
  /// @param lsexpr representation
  explicit Lattice_element(Lsexpr const& lsexpr);

  /// Extract an Lsexpr representation of the Lattice_element
  Lsexpr as_lsexpr() const;

  /// Get the type
  std::string const& get_type_name() const;

  element_type get_type() const;

  /// Get the name
  std::string const& get_name() const;

  /// Get the version info
  element_format get_format() const;

  /// Add an ancestor to the list of ancestors
  /// @param ancestor ancestor name
  void add_ancestor(std::string const& ancestor);

  /// Get the list of ancestors
  std::list<std::string> const& get_ancestors() const;

  /// Duplicate the named attribute with the new_name.
  /// It looks for the named attribute in all categories: double,
  /// lazy_double, vector, lasy_vector, and strings
  /// Do nothing if the named attribute is not found.
  /// Throws an exception if the new name is already an attribute
  /// and the overwrite flag is false
  /// @param name name of the attribute to be duplicated
  /// @param new_name attribute name after the duplication
  void duplicate_attribute(std::string const& name,
                           std::string const& new_name,
                           bool overwrite = false);

  /// Delete the named attribute.
  /// It looks for the named attribute in all categories: double,
  /// lazy_double, vector, lazy_vector, and strings
  /// Do nothing if the named attribute is not found.
  /// @param name name of the attribute to be duplicated
  void delete_attribute(std::string const& name);

  /// Rename the named attribute to new_name
  /// It is equivalent of calling the duplicate_attribute(name, new_name)
  /// then delete_attribute(name).
  /// @param name name of the attribute to be renamed
  /// @param new_name attribute name after the renaming
  void rename_attribute(std::string const& name, std::string const& new_name);

  /// Copy all attributes from another element
  void copy_attributes_from(Lattice_element const& elem);

  /// Set the value of the named double attribute
  /// @param name attribute name
  /// @param value attribute value
  /// @param increment_revision can be set to false for
  /// attributes that do not affect dynamics
  void set_double_attribute(std::string const& name,
                            double value,
                            bool increment_revision = true);

  void set_double_attribute(std::string const& name,
                            synergia::mx_expr const& value,
                            bool increment_revision = true);

  void set_double_attribute(std::string const& name,
                            const char* value,
                            bool increment_revision = true);

  void set_double_attribute(std::string const& name,
                            std::string const& value,
                            bool increment_revision = true);

  void set_default_double_attribute(std::string const& name,
                                    double value,
                                    bool increment_revision = true);

  /// Remove the named double attribute from the element
  /// @param name attribute name

  void remove_double_attribute(std::string const& name);

  /// Check for the existence of the named double attribute
  /// @param name attribute name
  bool has_double_attribute(std::string const& name) const;

  /// Get the value of the named double attribute
  /// @param name attribute name
  double get_double_attribute(std::string const& name) const;

  /// Get the value of the named double attribute
  /// @param name attribute name
  /// @param val default value if the specified attribute doesnt exist
  double get_double_attribute(std::string const& name, double val) const;

  /// Set the value of the named string attribute
  /// @param name attribute name
  /// @param value attribute value
  /// @param increment_revision can be set to false for
  /// attributes that do not affect dynamics
  void set_string_attribute(std::string const& name,
                            std::string const& value,
                            bool increment_revision = true);

  void set_default_string_attribute(std::string const& name,
                                    std::string const& value,
                                    bool incremnt_revision = true);

  /// Remove the named string attribute
  /// @param name attribute name
  void remove_string_attribute(std::string const&name);

  /// Check for the existence of the named string attribute
  /// @param name attribute name
  bool has_string_attribute(std::string const& name) const;

  /// Get the value of the named string attribute
  /// @param name attribute name
  std::string const& get_string_attribute(std::string const& name) const;

  /// Get the value of the named string attribute
  /// @param name attribute name
  /// @param val default value if the specified attribute doesnt exist
  std::string const& get_string_attribute(std::string const& name,
                                          std::string const& val) const;

  /// Set the value of the named vector attribute
  /// @param name attribute name
  /// @param value attribute value
  /// @param increment_revision can be set to false for attributes that do not
  /// affect dynamics
  void set_vector_attribute(std::string const& name,
                            std::vector<double> const& value,
                            bool increment_revision = true);

  /// Check for the existence of the named vector attribute
  /// @param name attribute name
  bool has_vector_attribute(std::string const& name) const;

  /// remove the named vector attribute
  /// @param name attribute name
  void remove_vector_attribute(std::string const& name);

  /// Get the value of the named vector attribute
  /// @param name attribute name
  std::vector<double> get_vector_attribute(std::string const& name) const;

  /// Get the value of the named vector attribute
  /// @param name attribute name
  /// @param val default value if the specified attribute doesnt exist
  std::vector<double> get_vector_attribute(
    std::string const& name,
    std::vector<double> const& val) const;

  /// Set the attribute name to be used to determine the length
  /// of the Lattice_element
  /// @param attribute_name attribute name
  void set_length_attribute_name(std::string const& attribute_name);

  /// Set the attribute name to be used to determine the bend_angle
  /// of the Lattice_element
  /// @param attribute_name attribute name
  void set_bend_angle_attribute_name(std::string const& attribute_name);

  /// Get the Lattice_element's length
  double get_length() const;

  /// Get the Lattice_element's bend angle
  double get_bend_angle() const;

  /// (re)set a marker for the element. examples of supported
  /// markers are,
  ///   horizontal/vertical tunes corrector,
  ///   horizontal/vertical chromaticity crrector
  void set_marker(marker_type t);

  void
  reset_marker(marker_type t)
  {
    markers[(int)t] = false;
  }

  void
  reset_markers()
  {
    for (auto& m : markers) m = false;
  }

  bool
  has_marker(marker_type t) const
  {
    return markers[(int)t];
  }

  /// deposited charge
  double
  get_deposited_charge(int bunch = 0, int train = 0) const
  {
    return deposited_charge;
  }

  void
  set_deposited_charge(double charge, int bunch = 0, int train = 0) const
  {
    deposited_charge = charge;
  }

  void
  deposit_charge(double charge, int bunch = 0, int train = 0) const
  {
    deposited_charge += charge;
  }

  /// Get the Lattice_element's revision number
  long int get_revision() const;

  /// Check whether the element has a reference to a parent lattice
  bool has_lattice() const;

  /// Set the reference to the parent lattice
  void set_lattice(Lattice& lattice);

  /// Get a reference to the parent lattice
  Lattice const& get_lattice() const;

  /// Return a human-readable description of the Lattice_element
  std::string as_string() const;

  /// Return a madx string
  std::string as_madx(bool sanitize=false) const;

  /// Print a human-readable description of the Lattice_element
  /// The Python version of the function is named "print_".
  void print() const;

  // get the names of all double attributes
  std::vector<std::string> get_double_attribute_names() const;

  /// deprecated
  std::map<std::string, double> get_double_attributes() const;

  /// deprecated
  std::map<std::string, std::string> const&
  get_string_attributes() const
  {
    return string_attributes;
  }

private:
  friend class Lattice;
  friend class cereal::access;

  template <class Archive>
  void
  serialize(Archive& ar)
  {
    ar(CEREAL_NVP(name));
    ar(CEREAL_NVP(format));
    ar(CEREAL_NVP(stype));
    ar(CEREAL_NVP(type));
    ar(CEREAL_NVP(ancestors));

    ar(CEREAL_NVP(string_attributes));
    ar(CEREAL_NVP(lazy_double_attributes));
    ar(CEREAL_NVP(lazy_vector_attributes));

    ar(CEREAL_NVP(length_attribute_name));
    ar(CEREAL_NVP(bend_angle_attribute_name));
    ar(CEREAL_NVP(revision));
    ar(CEREAL_NVP(markers));
  }
};

#endif /* LATTICE_ELEMENT_H_ */
