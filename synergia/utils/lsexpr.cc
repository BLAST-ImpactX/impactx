#include <cstring>
#include <iostream>
#include <limits>

#include "lsexpr.h"

Lsexpr::Lsexpr() : label(""), is_atom(false), atom(""), sequence() {}

Lsexpr::Lsexpr(std::string const& atom, std::string const& label)
  : label(label), is_atom(true), atom(atom), sequence()
{}

std::string
Lsexpr::to_string(int i) const
{
  std::stringstream ss;
  ss << i;
  return ss.str();
}

int
Lsexpr::stoi(std::string const& s) const
{
  int retval;
  std::istringstream ss(s);
  ss >> retval;
  return retval;
}

std::string
Lsexpr::to_string(double d) const
{
  std::stringstream ss;
  const int precision = std::numeric_limits<double>::max_digits10;
  ss.precision(precision);
  ss << d;
  return (ss.str());
}

double
Lsexpr::stod(std::string const& s) const
{
  double retval;
  std::istringstream ss(s);
  ss >> retval;
  return retval;
}

Lsexpr::Lsexpr(int atom, std::string const& label)
  : label(label), is_atom(true), atom(to_string(atom)), sequence()
{}

Lsexpr::Lsexpr(double atom, std::string const& label)
  : label(label), is_atom(true), atom(to_string(atom)), sequence()
{}

Lsexpr::Lsexpr(std::vector<std::string> const& sequence,
               std::string const& label)
  : label(label), is_atom(false), atom(), sequence()
{
  for (std::vector<std::string>::const_iterator it = sequence.begin();
       it != sequence.end();
       ++it) {
    this->sequence.push_back(Lsexpr(*it));
  }
}

Lsexpr::Lsexpr(std::vector<int> const& sequence, std::string const& label)
  : label(label), is_atom(false), atom(), sequence()
{
  for (std::vector<int>::const_iterator it = sequence.begin();
       it != sequence.end();
       ++it) {
    this->sequence.push_back(Lsexpr(*it));
  }
}

Lsexpr::Lsexpr(std::vector<double> const& sequence, std::string const& label)
  : label(label), is_atom(false), atom(), sequence()
{
  for (std::vector<double>::const_iterator it = sequence.begin();
       it != sequence.end();
       ++it) {
    this->sequence.push_back(Lsexpr(*it));
  }
}

Lsexpr::Lsexpr(std::istream& stream)
  : label(""), is_atom(false), atom(""), sequence()
{
  std::vector<Lsexpr> stack;
  Lsexpr root;
  stack.push_back(root);
  bool in_str(false);
  std::string current_word("");
  std::string current_label("");
  char c(stream.get());
  while (!stream.eof()) {
    if (c == '{' && !in_str) {
      Lsexpr lsexpr;
      if (!current_label.empty()) {
        lsexpr.set_label(current_label);
        current_label = "";
      }
      stack.push_back(lsexpr);
    } else if (c == '}' && !in_str) {
      if (!current_word.empty()) {
        Lsexpr lsexpr(current_word);
        if (!current_label.empty()) { lsexpr.set_label(current_label); }
        stack.back().push_back(lsexpr);
      }
      current_word = "";
      current_label = "";
      Lsexpr tmp(stack.back());
      stack.pop_back();
      stack.back().push_back(tmp);
    } else if ((c == ',') && !in_str) {
      if (!current_word.empty()) {
        Lsexpr lsexpr(current_word);
        if (!current_label.empty()) { lsexpr.set_label(current_label); }
        stack.back().push_back(lsexpr);
      }
      current_word = "";
      current_label = "";
    } else if (c == '"') {
      in_str = !in_str;
    } else if (((c == ' ') || (c == '\n') || (c == '\t')) && !in_str) {
      // pass
    } else if ((c == ':' && !in_str)) {
      current_label = current_word;
      current_word = "";
    } else {
      current_word += c;
    }
    c = stream.get();
  }
  if (!current_word.empty()) {
    Lsexpr lsexpr(current_word);
    if (!current_label.empty()) { lsexpr.set_label(current_label); }
    stack.back().push_back(lsexpr);
  }
  *this = stack.front().sequence.front();
}

void
Lsexpr::set_label(std::string const& label)
{
  this->label = label;
}

bool
Lsexpr::is_labeled() const
{
  return !label.empty();
}

std::string const&
Lsexpr::get_label() const
{
  return label;
}

bool
Lsexpr::is_atomic() const
{
  return is_atom;
}

std::string
Lsexpr::get_string() const
{
  if (!is_atom) {
    throw std::runtime_error("Lsexpr::get_string: non-atomic Lsexpr");
  }
  return atom;
}

int
Lsexpr::get_int() const
{
  if (!is_atom) {
    throw std::runtime_error("Lsexpr::get_int: non-atomic Lsexpr");
  }
  return stoi(atom);
}

double
Lsexpr::get_double() const
{
  if (!is_atom) {
    throw std::runtime_error("Lsexpr::get_double: non-atomic Lsexpr " + label);
  }
  return stod(atom);
}

void
Lsexpr::push_back(Lsexpr const& lsexpr)
{
  sequence.push_back(lsexpr);
}

Lsexpr::const_iterator_t
Lsexpr::begin() const
{
  return sequence.begin();
}

Lsexpr::const_iterator_t
Lsexpr::end() const
{
  return sequence.end();
}

Lsexpr::iterator_t
Lsexpr::begin()
{
  return sequence.begin();
}

Lsexpr::iterator_t
Lsexpr::end()
{
  return sequence.end();
}

std::vector<std::string>
Lsexpr::get_string_vector() const
{
  if (is_atom) {
    throw std::runtime_error("Lsexpr::get_string_vector: atomic Lsexpr");
  }
  std::vector<std::string> retval;
  for (const_iterator_t it = begin(); it != end(); ++it) {
    retval.push_back(it->get_string());
  }
  return retval;
}

std::vector<int>
Lsexpr::get_int_vector() const
{
  if (is_atom) {
    throw std::runtime_error("Lsexpr::get_int_vector: atomic Lsexpr");
  }
  std::vector<int> retval;
  for (const_iterator_t it = begin(); it != end(); ++it) {
    retval.push_back(stoi(it->get_string()));
  }
  return retval;
}

std::vector<double>
Lsexpr::get_double_vector() const
{
  if (is_atom) {
    throw std::runtime_error("Lsexpr::get_double_vector: atomic Lsexpr");
  }
  std::vector<double> retval;
  for (const_iterator_t it = begin(); it != end(); ++it) {
    retval.push_back(stod(it->get_string()));
  }
  return retval;
}

namespace {
  const char* specials = " \t\n:{},\"";

  void
  write_atom(std::string const& atom, std::ostream& stream)
  {
    if (std::strcspn(atom.c_str(), specials) != atom.length()) {
      stream << '"' << atom << '"';
    } else {
      stream << atom;
    }
  }

  void
  write_label(std::string const& label, std::ostream& stream)
  {
    write_atom(label, stream);
    stream << ": ";
  }

  void
  write_lsexpr_internal(Lsexpr const& lsexpr,
                        std::ostream& stream,
                        int indent,
                        bool first)
  {
    if (lsexpr.is_labeled()) {
      if (!lsexpr.is_atomic() && !first) {
        stream << "\n";
        if (indent > -1) { stream << std::string(indent + 1, ' '); }
      }
      write_label(lsexpr.get_label(), stream);
    }
    if (lsexpr.is_atomic()) {
      write_atom(lsexpr.get_string(), stream);
    } else {
      if (!first || lsexpr.is_labeled()) {
        if (!lsexpr.is_atomic()) {
          if (lsexpr.is_labeled()) {
            indent += 4;
          } else {
            indent += 1;
          }
        }
        stream << "\n" << std::string(indent + 1, ' ');
      }
      stream << "{";
      bool first_in_sequence(true);
      for (Lsexpr::const_iterator_t it = lsexpr.begin(); it != lsexpr.end();
           ++it) {
        write_lsexpr_internal(*it, stream, indent, first_in_sequence);
        first_in_sequence = false;
        if ((it + 1) != lsexpr.end()) { stream << ", "; }
      }
      stream << "}";
      if (!first || lsexpr.is_labeled()) {
        if (lsexpr.is_labeled()) {
          indent -= 4;
        } else {
          indent -= 1;
        }
      }
    }
  }
}

void
Lsexpr::write(std::ostream& stream) const
{
  write_lsexpr_internal(*this, stream, 0, true);
}

Lsexpr
read_lsexpr_file(std::string const& filename)
{
  std::ifstream f(filename.c_str());
  if (!f.is_open()) {
    throw std::runtime_error("read_lsexpr_file:: unable to open '" + filename +
                             "'");
  }
  Lsexpr retval(f);
  return retval;
}

void
write_lsexpr_file(Lsexpr const& lsexpr, std::string const& filename)
{
  std::ofstream f(filename.c_str());
  if (!f.is_open()) {
    throw std::runtime_error("write_lsexpr_file:: unable to open '" + filename +
                             "'");
  }
  lsexpr.write(f);
}
