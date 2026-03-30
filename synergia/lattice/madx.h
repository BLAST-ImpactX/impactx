#ifndef SYNERGIA_LATTICE_MAD8_H
#define SYNERGIA_LATTICE_MAD8_H

#include <any>
#include <string>
#include <vector>
#include <map>

#include "mx_expr.h"

namespace synergia
{
  struct MadX_value;
  class  MadX_command;
  class  MadX_line;
  class  MadX_sequence;
  class  MadX;

  enum MadX_value_type
  { NONE
  //, BOOLEAN
  //, INT
  //, DOUBLE
  , NUMBER
  , STRING
  , DEFFERED_NUMBER
  , ARRAY
  , DEFFERED_ARRAY 
  };

  enum MadX_command_type
  { NULL_COMMAND
  , ELEMENT
  , ELEMENT_REF
  , EXECUTABLE
  , SEQUENCE };

  enum MadX_sequence_refer
  { SEQ_REF_ENTRY
  , SEQ_REF_CENTRE
  , SEQ_REF_EXIT };

  enum MadX_entry_type
  { ENTRY_NULL
  , ENTRY_VARIABLE
  , ENTRY_COMMAND
  , ENTRY_LINE
  , ENTRY_SEQUENCE };


  std::string 
  to_string(MadX_value const& val);
}

struct synergia::MadX_value
{
  std::any      value;
  MadX_value_type type;
};

typedef std::string                              string_t;
typedef std::map<string_t, synergia::MadX_value> value_map_t;

class synergia::MadX_command
{
public:
  MadX_command(MadX * parent = NULL)
    : mx(parent)
    , name_()
    , label_()
    , type_(NULL_COMMAND)
    , attributes_()
  { }

  // accessor
  string_t name() const;
  string_t label() const;

  size_t                attribute_count() const;
  std::vector<string_t> attribute_names() const;
  MadX_value_type       attribute_type(string_t const & name) const;

  string_t attribute_as_string(string_t const & name) const;
  string_t attribute_as_string(string_t const & name, string_t const & def) const;
  double   attribute_as_number(string_t const & name) const;
  double   attribute_as_number(string_t const & name, double def) const;
  bool     attribute_as_boolean(string_t const & name) const;

  mx_expr  attribute_as_expr(string_t const& name) const;
  std::vector<mx_expr> attribute_as_expr_seq(string_t const & name) const;

  std::vector<double> attribute_as_number_seq(string_t const & name) const;
  std::vector<double> attribute_as_number_seq(string_t const & name, double def) const;

  MadX_command_type type() const;
  bool is_element() const;
  bool is_reference() const;
  bool is_command() const;

  // modifier
  void set_parent(MadX const & parent);
  void set_name(string_t const & name, MadX_command_type cmd_type);
  void set_label(string_t const & label);
  void insert_attribute(string_t const & name, string_t const & val);
  void insert_attribute(string_t const & name, mx_expr  const & val);
  void insert_attribute(string_t const & name, mx_exprs const & val);
  void merge_with_overwrite(MadX_command const & other);
  void merge(MadX_command const & other);

  std::string to_string() const;

private:
  const MadX *      mx;
  string_t          name_;
  string_t          label_;
  MadX_command_type type_;
  value_map_t       attributes_;
};

typedef std::pair<string_t, synergia::MadX_command> labeled_cmd_t;
typedef std::vector<labeled_cmd_t>                  labeled_cmd_v_t;
typedef std::vector<synergia::MadX_command>         commands_v_t;
typedef std::map<string_t, synergia::MadX_command>  commands_m_t;

class synergia::MadX_line
{
public:
  MadX_line(MadX const & parent) 
    : parent(parent), elements_() { }

  // accessor
  size_t       element_count() const;
  string_t     element_name(size_t idx) const;
  MadX_command element(size_t idx, bool resolve=true) const;

  // modifier
  void insert_element(string_t const & ele);

  // print
  void print() const;

private:

  MadX const & parent;
  std::vector<string_t>  elements_;
};

typedef std::map<string_t, synergia::MadX_line>    lines_m_t;

class synergia::MadX_sequence
{
private:
  struct seq_element
  {
    seq_element(std::string const & label_, 
            mx_expr const& at_, 
            std::string const & from_)
      : label(label_)
      , from_str(from_)
      , at(at_)
      , from(0.0) 
    { }

    std::string label;
    std::string from_str;
    mx_expr     at;
    double      from;
  };

  typedef std::vector<seq_element> seq_ele_v_t;

public:
  MadX_sequence(MadX const & parent)
    : parent(parent), lbl(), l(0.0), r(SEQ_REF_CENTRE), rp(), seq_() { }

  // accessor
  string_t label() const;
  double   length() const;
  size_t   element_count() const;
  MadX_command      element(size_t idx) const;
  double            element_at(size_t idx) const;
  double            element_from(size_t idx) const;
  string_t          element_label(size_t idx) const;
  MadX_entry_type   element_type(size_t idx) const;
  MadX_sequence_refer refer() const;
  string_t            refpos() const;

  // modifier
  void set_label(string_t const & label);
  void set_length(double length);
  void set_refer(MadX_sequence_refer ref);
  void set_refpos(string_t const & refpos);
  void add_element(string_t const & label, mx_expr const& at, string_t const & from);
  void finalize();
  void reset();

  // print
  void print() const;

private:
  MadX const & parent;
  string_t lbl;
  double l;
  MadX_sequence_refer r;
  string_t rp;
  seq_ele_v_t   seq_;
};

typedef std::map<string_t, synergia::MadX_sequence> sequences_m_t;

class synergia::MadX
{
public:

  // ctor
  MadX();

  // both copy and assignment ignores the lines and sequences
  // for the reason that the copy/assignment only happens when
  // mx is used for dynamic lattice, and no references can be
  // made into the inside of lines or sequences.
  MadX(MadX const& mx);
  MadX& operator=(MadX const& mx);

  // accessor
  string_t variable_as_string (string_t const & name) const;
  string_t variable_as_string (string_t const & name, string_t const & def) const;
  double   variable_as_number (string_t const & name) const;
  double   variable_as_number (string_t const & name, double def) const;
  bool     variable_as_boolean(string_t const & name) const;
  std::vector<double> variable_as_number_seq(string_t const & name) const;
  std::vector<double> variable_as_number_seq(string_t const & name, double def) const;

  size_t command_count() const;  // un-labeled commands
  std::vector<string_t > commands() const;
  MadX_command command(size_t idx, bool resolve = true) const;

  size_t label_count() const;    // labeled commands
  std::vector<string_t > command_labels() const;
  MadX_command command(string_t const & l, bool resolve = true) const;
  MadX_command& command_ref(string_t const& l);

  size_t line_count() const;     // labeled lines
  std::vector<string_t > line_labels() const;
  MadX_line const & line(string_t const & l) const;

  size_t sequence_count() const; // labeled sequences
  std::vector<string_t> sequence_labels() const;
  MadX_sequence const & sequence(string_t const & s) const;
  MadX_sequence const & current_sequence( ) const;
  MadX_sequence & current_sequence( );

  bool building_sequence() const;

  MadX_entry_type entry_type(string_t const & entry) const;

  void print() const;

  // modifier
  void insert_variable(string_t const & name, string_t const & value);
  void insert_variable(string_t const & name, mx_expr  const & value);
  void insert_variable(string_t const & name, mx_exprs const & value);
  void insert_label   (string_t const & name, MadX_command const & cmd);
  void insert_line    (string_t const & name, MadX_line const & line);
  void insert_command (MadX_command const & cmd);
  void fuse_command   (string_t const & name, MadX_command const & cmd);
                          // fuse all attr in cmd into label with the provided name

  void start_sequence( string_t const & label, double length
                     , string_t const & refer, string_t const & refpos );
  void append_sequence_element( string_t const & label, mx_expr const& at
                              , string_t const & from );
  void end_sequence();

  // alias
  void insert_attribute(string_t const & name, string_t const & value)
    { insert_variable(name, value); }
  void insert_attribute(string_t const & name, mx_expr  const & value)
    { insert_variable(name, value); }
  void insert_attribute(string_t const & name, mx_exprs const & value)
    { insert_variable(name, value); }

  // export
  std::string to_madx() const;

private:

  std::string export_variables() const;
  std::string export_unnamed_cmds() const;
  std::string export_labled_cmds() const;

public:
  static double   nan;  // not a number
  static string_t nst;  // not a string

private:
  value_map_t   variables_;
  commands_v_t  cmd_seq_;
  commands_m_t  cmd_map_;
  lines_m_t     lines_;
  sequences_m_t seqs_;
  MadX_sequence cur_seq_;       // sequence thats being built currently
  bool          building_seq_;  // currently building sequence?
};


#endif
