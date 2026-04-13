#ifndef SYNERGIA_LATTICE_MAD8_H
#define SYNERGIA_LATTICE_MAD8_H

#include <any>
#include <string>
#include <vector>
#include <map>



namespace synergia
{
  struct Mad8_value;
  class  Mad8_command;
  class  Mad8_line;
  class  Mad8;

  enum Mad8_value_type
  { NONE, INT, DOUBLE, NUMBER, STRING };

  enum Mad8_line_element_type
  { LABEL, MULTIPLIER, MINUS, LINE };
}

struct synergia::Mad8_value
{
  std::any      value;
  Mad8_value_type type;
};

typedef std::string                              string_t;
typedef std::map<string_t, synergia::Mad8_value> value_map_t;

class synergia::Mad8_command
{
public:
  Mad8_command() : name_(), attributes_() { }

  // accessor
  string_t name() const;
  size_t   attribute_count() const;
  string_t attribute_as_string(string_t const & name) const;  
  double   attribute_as_number(string_t const & name) const;  

  // modifier
  void set_name(string_t const & name);
  void insert_attribute(string_t const & name, double val);
  void insert_attribute(string_t const & name, string_t const & val);

private:
  string_t     name_;
  value_map_t  attributes_;
};

typedef std::vector<synergia::Mad8_command>        commands_v_t;
typedef std::map<string_t, synergia::Mad8_command> commands_m_t;

class synergia::Mad8_line
{
  typedef Mad8_line_element_type ele_t;

public:
  Mad8_line() : elements_(), lines_(), index_() { }

  // accessor
  size_t    element_count() const;
  ele_t     element_type(size_t idx) const;
  string_t  element_as_string(size_t idx) const;
  Mad8_line element_as_line  (size_t idx) const;

  // modifier
  void insert_operator(string_t const & op);
  void insert_element(string_t const & ele);
  void insert_subline(Mad8_line const & line);

private:

  std::vector<string_t>  elements_;
  std::vector<Mad8_line> lines_;
  std::vector<std::pair<ele_t, size_t> > index_;  // index_[i].first = element type
                                                  // index_[i].second = index in respective container
};

typedef std::map<string_t, synergia::Mad8_line>    lines_m_t;

class synergia::Mad8
{
public:
  Mad8() : variables_(), cmd_seq_(), cmd_map_(), lines_() {}

  // accessor
  string_t variable_as_string(string_t const & name) const;
  double   variable_as_number(string_t const & name) const;

  size_t command_count() const;  // un-labeled commands
  Mad8_command const & command(size_t idx) const;

  size_t label_count() const;    // labeled commands
  Mad8_command const & label(string_t const & l) const;

  size_t line_count() const;     // labeled lines
  Mad8_line const & line(string_t const & l) const;

  // modifier
  void insert_variable(string_t const & name, double value);
  void insert_variable(string_t const & name, string_t const & value);
  void insert_label   (string_t const & name, Mad8_command const & cmd);
  void insert_line    (string_t const & name, Mad8_line const & line);
  void insert_command (Mad8_command const & cmd);

private:
  value_map_t  variables_;
  commands_v_t cmd_seq_;
  commands_m_t cmd_map_;
  lines_m_t    lines_;
};



#endif
