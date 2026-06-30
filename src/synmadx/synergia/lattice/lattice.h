#ifndef LATTICE_H_
#define LATTICE_H_

#include <list>
#include <string>
#include <optional>
#include <stdexcept>

#include "synergia/foundation/reference_particle.h"
#include "synergia/lattice/lattice_element.h"
#include "synergia/lattice/lattice_tree.h"

#include <iostream>

/// The Lattice class contains an abstract representation of an ordered
/// set of objects of type Lattice_element.
/// Each element of the Lattice is unique.
class Lattice {
public:
  struct update_flags_t {
    bool ref;       // lattice reference particle updated
    bool structure; // add or remove any element
    bool element;   // changes to the element attributes
  };

private:
  std::string name;

  std::optional<Reference_particle> reference_particle;
  std::list<Lattice_element> elements;

  update_flags_t updated;

  // Lattice tree object for evaluating variables
  // in the lattice element attributes
  Lattice_tree tree;

public:
  /// Construct a Lattice object without a name.
  /// Defaults to interpreting elements as Mad8 elements
  Lattice();

  /// Copy, move, and assignment of Lattices contain copies of elements
  Lattice(Lattice const& lattice);
  Lattice(Lattice&& lattice) noexcept;
  Lattice& operator=(Lattice const& lattice);

  /// Construct a Lattice object with a name
  /// Defaults to interpreting elements as Mad8 elements
  /// @param name an arbitrary name
  explicit Lattice(std::string const& name);

  /// Construct a Lattice object with a name and a reference particle
  /// Defaults to interpreting elements as Mad8 elements
  /// @param name an arbitrary name
  Lattice(std::string const& name, Reference_particle const& ref);

  /// Construct a dynamic lattice object
  Lattice(std::string const& name, Lattice_tree const& tree);

  /// Construct a Lattice from the Lsexpr representation
  /// @param lsexpr representation
  explicit Lattice(Lsexpr const& lsexpr);

#if 0 // as_lsexpr() might come back some day
    /// Extract an Lsexpr representation of the Lattice
    Lsexpr
    as_lsexpr() const;
#endif

  /// Get the Lattice name
  std::string const&
  get_name() const
  {
    return name;
  }

  /// Always a dynamic lattice
  bool
  is_dynamic_lattice() const
  {
    return true;
  }

  /// Set the Lattice reference particle
  /// @param reference_particle a Reference_particle
  void
  set_reference_particle(Reference_particle const& ref)
  {
    reference_particle = ref;
    updated.ref = true;
  }


    // check whether the reference particle is valid, throw if not
    static void check_reference_particle_value(std::optional<Reference_particle> reference_particle)
    {
       if (!reference_particle.has_value()) {
            throw std::runtime_error("reference particle not set- did you forget a BEAM? statement?");
       }
    }

    /// Get the Lattice reference particle (const)
    Reference_particle const& get_reference_particle() const
    {
        check_reference_particle_value(reference_particle);
        // only returns if valid
        return reference_particle.value();
    }

    Reference_particle& get_reference_particle()
    {
        check_reference_particle_value(reference_particle);
        // only returns if valid
        return reference_particle.value();
    }

    double
    get_lattice_energy() const
    {
       check_reference_particle_value(reference_particle);
        // only returns if valid
      return (reference_particle.value()).get_total_energy();
    }

    void
    set_lattice_energy(double energy)
    {
       check_reference_particle_value(reference_particle);
        // only returns if valid
        // std::cout << "EGS: set_lattice_energy: " << energy << std::endl;
      (reference_particle.value()).set_total_energy(energy);
    }

    update_flags_t update();
    update_flags_t is_updated() const { return updated; }

  /// Append a copy of a Lattice_element.
  /// @param element a Lattice_element
  void append(Lattice_element const& element);

    /// Set the value of the named double attribute on all elements
    /// @param name attribute name
    /// @param value attribute value
    /// @param increment_revision can be set to false for attributes that do not affect dynamics
    void
    set_all_double_attribute(
            std::string const& name, double value,
            bool increment_revision = true);

  /// Set the value of the named string attribute on all elements
  /// @param name attribute name
  /// @param value attribute value
  /// @param increment_revision can be set to false for attributes that do not
  /// affect dynamics
  void set_all_string_attribute(std::string const& name,
                                std::string const& value,
                                bool increment_revision = true);

  /// Get the list of elements in the Lattice
  std::list<Lattice_element> const& get_elements() const;

  std::list<Lattice_element>& get_elements();

  /// Clear the h/v tunes and chromaticity markers for all lattice elements
  void reset_all_markers();

  /// Get the combined length of all the elements in the Lattice
  double get_length() const;

  /// Get the total angle in radians subtended by all the elements in the
  /// Lattice
  double get_total_angle() const;

  /// Return a human-readable summary of the elements in the Lattice.
  std::string as_string() const;

  /// Print a human-readable summary of the elements in the Lattice.
  /// The Python version of this function is named "print_".
  void print(std::ostream& os = std::cerr) const;

  /// Get the Lattice_tree object
  Lattice_tree& get_lattice_tree();

  Lattice_tree const& get_lattice_tree() const;

  /// Set the Lattice_tree object
  void set_lattice_tree(Lattice_tree const& tree);

  /// Set a variable in the Lattice_tree object
  void set_variable(std::string const& name, double val);

  void set_variable(std::string const& name, std::string const& val);

public:
  // export madx file
  void export_madx_file(std::string const& filename, bool const sanitize=false) const;

  // read from madx file
  static Lattice import_madx_file(std::string const& filename,
                                  std::string const& line);
};

#endif /* LATTICE_H_ */
