#ifndef LSEXPR_H
#define LSEXPR_H

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

class Lsexpr;

class Lsexpr {
public:
  typedef std::vector<Lsexpr> sequence_t;
  typedef sequence_t::const_iterator const_iterator_t;
  typedef sequence_t::iterator iterator_t;

  Lsexpr();

  Lsexpr(std::string const& atom, std::string const& label = "");

  Lsexpr(int atom, std::string const& label = "");

  Lsexpr(double atom, std::string const& label = "");

  Lsexpr(std::vector<std::string> const& sequence,
         std::string const& label = "");

  Lsexpr(std::vector<int> const& sequence, std::string const& label = "");

  Lsexpr(std::vector<double> const& sequence, std::string const& label = "");

  Lsexpr(std::istream& stream);

  void set_label(std::string const& label);

  bool is_labeled() const;

  std::string const& get_label() const;

  bool is_atomic() const;

  std::string get_string() const;

  int get_int() const;

  double get_double() const;

  void push_back(Lsexpr const& lsexpr);

  const_iterator_t begin() const;

  const_iterator_t end() const;

  iterator_t begin();

  iterator_t end();

  std::vector<std::string> get_string_vector() const;

  std::vector<int> get_int_vector() const;

  std::vector<double> get_double_vector() const;

  void write(std::ostream& stream) const;

private:
  std::string to_string(int i) const;
  int stoi(std::string const& s) const;
  std::string to_string(double d) const;
  double stod(std::string const& s) const;

  std::string label;
  bool is_atom;
  std::string atom;
  sequence_t sequence;
};

Lsexpr read_lsexpr_file(std::string const& filename);

void write_lsexpr_file(Lsexpr const& lsexpr, std::string const& filename);

#endif // LSEXPR_H
