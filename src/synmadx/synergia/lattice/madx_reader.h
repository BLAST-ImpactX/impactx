#ifndef MADX_READER_CC_
#define MADX_READER_CC_

#include "lattice.h"
#include "madx.h"

class MadX_reader {
private:
  synergia::MadX madx;

public:
  void parse(std::string const& string);
  void parse_file(std::string const& filename);

  std::vector<std::string> get_line_names() const;
  std::vector<std::string> get_sequence_names() const;
  std::vector<std::string> get_all_names() const;

  double get_double_variable(std::string const& name) const;
  std::string get_string_variable(std::string const& name) const;

  Lattice get_lattice(std::string const& line_name);
  Lattice get_lattice(std::string const& line_name,
                      std::string const& filename);

  Lattice get_dynamic_lattice(std::string const& line_name);
  Lattice get_dynamic_lattice(std::string const& line_name,
                              std::string const& filename);

  static Lattice get_lattice(std::string const& line_name,
                             synergia::MadX const& mx);

  static Lattice get_dynamic_lattice(std::string const& line_name,
                                     synergia::MadX const& mx);
};
#endif /* MADX_READER_CC_ */
