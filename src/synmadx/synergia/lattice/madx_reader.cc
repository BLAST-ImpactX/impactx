#include <algorithm>
#include <iomanip>
#include <stdexcept>

#include "synergia/foundation/four_momentum.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/foundation/reference_particle.h"
#include "synergia/utils/digits.h"

#include "madx_reader.h"
#include "mx_parse.h"

using namespace synergia;

namespace {

  void
  insert_line(Lattice& lattice,
              MadX const& mx,
              std::string const& line_name,
              bool dynamic = false)
  {
    MadX_line const& line = mx.line(line_name);

    for (int i = 0; i < line.element_count(); ++i) {
      MadX_command cmd = line.element(i, true);
      std::string name = line.element_name(i);
      std::string type = cmd.name();

      Lattice_element element(type, name);
      std::vector<std::string> attr_names(cmd.attribute_names());

      for (std::vector<std::string>::const_iterator it = attr_names.begin();
           it != attr_names.end();
           ++it)
        for (auto const& name : attr_names) {
          MadX_value_type vt = cmd.attribute_type(name);

          switch (vt) {
            case NONE:

            case STRING:
              element.set_string_attribute(name,
                                           cmd.attribute_as_string(name, ""));
              break;

            case NUMBER:
              if (dynamic)
                element.set_double_attribute(name, cmd.attribute_as_expr(name));
              else
                element.set_double_attribute(
                  name, cmd.attribute_as_number(name, 0.0));
              break;

            case ARRAY:
              element.set_vector_attribute(
                name, cmd.attribute_as_number_seq(name, 0.0));
              break;

            default:
              throw std::runtime_error("unable to process attribute " + name +
                                       " of element " + name);
          }
        }

      lattice.append(element);
    }
  }

  double
  insert_sequence(Lattice& lattice,
                  MadX const& mx,
                  std::string const& line_name,
                  bool dynamic = false)
  {
    MadX_sequence const& sequence = mx.sequence(line_name);

    double r = 0.5;
    double total_length = sequence.length();
    double current_pos = 0.0;
    const double min_drift_length = 1.0e-6;
    int drift_count = 0;
    int drift_digits = digits(sequence.element_count());

    MadX_sequence_refer ref = sequence.refer();
    if (ref == SEQ_REF_ENTRY)
      r = 1.0;
    else if (ref == SEQ_REF_CENTRE)
      r = 0.5;
    else
      r = 0.0;

    for (int i = 0; i < sequence.element_count(); ++i) {

      double at = sequence.element_at(i);
      double from = sequence.element_from(i);
      std::string label = sequence.element_label(i);

      if (sequence.element_type(i) == ENTRY_SEQUENCE) {
        double l = insert_sequence(lattice, mx, label, dynamic);
        current_pos =
          at + from + l; // sub-sequence always refer to the entry point

        continue;
      }

      MadX_command cmd = sequence.element(i);
      std::string type(cmd.name());
      std::string name(cmd.label());

      Lattice_element element(type, name);
      element.set_lattice(lattice);

      std::vector<string_t> attribute_names(cmd.attribute_names());

      for (std::vector<string_t>::iterator it = attribute_names.begin();
           it != attribute_names.end();
           ++it) {

        MadX_value_type vt = cmd.attribute_type(*it);

        switch (vt) {
          case NONE:
          case STRING:
            element.set_string_attribute(*it, cmd.attribute_as_string(*it, ""));
            break;
          case NUMBER:
            if (dynamic)
              element.set_double_attribute(*it, cmd.attribute_as_expr(*it));
            else
              element.set_double_attribute(*it,
                                           cmd.attribute_as_number(*it, 0.0));

            break;
          case ARRAY:
            element.set_vector_attribute(*it,
                                         cmd.attribute_as_number_seq(*it, 0.0));
            break;
          default:
            throw std::runtime_error("unable to process attribute " + *it +
                                     " of element " + name);
        }
      }

      double drift_length =
        at + from - current_pos - element.get_length() * (1.0 - r);
      if (drift_length > min_drift_length) {
        std::stringstream name_stream;
        name_stream << "auto_drift";
        name_stream << "_";
        name_stream << line_name;
        name_stream << "_";
        name_stream << std::setw(drift_digits);
        name_stream << std::setfill('0');
        name_stream << drift_count;
        Lattice_element drift("drift", name_stream.str());
        drift.set_double_attribute("l", drift_length, false);
        lattice.append(drift);
        ++drift_count;
      }
      lattice.append(element);
      current_pos = at + from + element.get_length() * r;
    }

    double final_drift_length = sequence.length() - current_pos;
    if (final_drift_length > min_drift_length) {
      std::stringstream name_stream;
      name_stream << "auto_drift";
      name_stream << "_";
      name_stream << line_name;
      name_stream << "_";
      name_stream << std::setw(drift_digits);
      name_stream << std::setfill('0');
      name_stream << drift_count;
      Lattice_element drift("drift", name_stream.str());
      drift.set_double_attribute("l", final_drift_length, false);
      lattice.append(drift);
      ++drift_count;
    }

    return sequence.length();
  }

  void
  extract_reference_particle(Lattice& lattice, MadX const& mx)
  {
    try {
      synergia::MadX_command command(mx.command("beam"));
      std::vector<std::string> attributes(command.attribute_names());
      double mass = 0, charge = 0, energy = 0, pc = 0, gamma = 0;
      for (std::vector<std::string>::const_iterator it = attributes.begin();
           it != attributes.end();
           ++it) {
        if (*it == "particle") {
          std::string particle(command.attribute_as_string(*it));
          if (particle == "proton") {
            mass = pconstants::mp;
            charge = pconstants::proton_charge;
          } else if (particle == "antiproton") {
            mass = pconstants::mp;
            charge = pconstants::antiproton_charge;
          } else if (particle == "electron") {
            mass = pconstants::me;
            charge = pconstants::electron_charge;
          } else if (particle == "positron") {
            mass = pconstants::me;
            charge = pconstants::positron_charge;
          } else if (particle == "negmuon") {
            mass = pconstants::mmu;
            charge = pconstants::muon_charge;
          } else if (particle == "posmuon") {
            mass = pconstants::mmu;
            charge = pconstants::antimuon_charge;
          }
        } else if (*it == "mass") {
          mass = command.attribute_as_number(*it);
        } else if (*it == "charge") {
          charge = command.attribute_as_number(*it);
        } else if (*it == "energy") {
          energy = command.attribute_as_number(*it);
        } else if (*it == "pc") {
          pc = command.attribute_as_number(*it);
        } else if (*it == "gamma") {
          gamma = command.attribute_as_number(*it);
        }
      }
      Four_momentum four_momentum(mass);
      if (energy > 0) { four_momentum.set_total_energy(energy); }
      if (pc > 0) { four_momentum.set_momentum(pc); }
      if (gamma > 0) { four_momentum.set_gamma(gamma); }
      Reference_particle reference_particle((int)charge, four_momentum);
      lattice.set_reference_particle(reference_particle);
    }
    catch (...) {
    }
  }
}

void
MadX_reader::parse(std::string const& string)
{
  synergia::parse_madx(string, madx);
}

void
MadX_reader::parse_file(std::string const& filename)
{
  synergia::parse_madx_file(filename, madx);
}

std::vector<std::string>
MadX_reader::get_line_names() const
{
  return madx.line_labels();
}

std::vector<std::string>
MadX_reader::get_sequence_names() const
{
  return madx.sequence_labels();
}

std::vector<std::string>
MadX_reader::get_all_names() const
{
  std::vector<std::string> retval(madx.line_labels());
  std::vector<std::string> sequence_labels(madx.sequence_labels());
  for (std::vector<std::string>::const_iterator it = sequence_labels.begin();
       it != sequence_labels.end();
       ++it) {
    retval.push_back(*it);
  }
  return retval;
}

double
MadX_reader::get_double_variable(std::string const& name) const
{
  return madx.variable_as_number(name);
}

std::string
MadX_reader::get_string_variable(std::string const& name) const
{
  return madx.variable_as_string(name);
}

// construct the regular lattice
Lattice
MadX_reader::get_lattice(std::string const& line_name)
{
  return MadX_reader::get_lattice(line_name, madx);
}

Lattice
MadX_reader::get_lattice(std::string const& line_name,
                         std::string const& filename)
{
  parse_file(filename);
  return get_lattice(line_name);
}

Lattice
MadX_reader::get_lattice(std::string const& line_name, synergia::MadX const& mx)
{
  Lattice lattice(line_name);

  std::vector<std::string> lines(mx.line_labels());
  if (std::find(lines.begin(), lines.end(), line_name) != lines.end()) {
    insert_line(lattice, mx, line_name);
  } else {
    std::vector<std::string> seq(mx.sequence_labels());
    if (std::find(seq.begin(), seq.end(), line_name) != seq.end()) {
      insert_sequence(lattice, mx, line_name);
    } else {
      throw std::runtime_error("MadX_reader::get_lattice: "
                               "cannot find line or sequence with label " +
                               line_name);
    }
  }

  extract_reference_particle(lattice, mx);
  return lattice;
}

// construct a dynamic lattice
Lattice
MadX_reader::get_dynamic_lattice(std::string const& line_name)
{
  return MadX_reader::get_dynamic_lattice(line_name, madx);
}

Lattice
MadX_reader::get_dynamic_lattice(std::string const& line_name,
                                 std::string const& filename)
{
  parse_file(filename);
  return get_dynamic_lattice(line_name);
}

Lattice
MadX_reader::get_dynamic_lattice(std::string const& line_name,
                                 synergia::MadX const& mx)
{
  Lattice lattice(line_name);
  lattice.set_lattice_tree(Lattice_tree(mx));

  std::vector<std::string> lines(mx.line_labels());
  if (std::find(lines.begin(), lines.end(), line_name) != lines.end()) {
    insert_line(lattice, mx, line_name, true);
  } else {
    std::vector<std::string> seq(mx.sequence_labels());
    if (std::find(seq.begin(), seq.end(), line_name) != seq.end()) {
      insert_sequence(lattice, mx, line_name, true);
    } else {
      throw std::runtime_error("MadX_reader::get_lattice: "
                               "cannot find line or sequence with label " +
                               line_name);
    }
  }

  extract_reference_particle(lattice, mx);
  return lattice;
}
