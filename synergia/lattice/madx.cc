
#include "madx.h"
#include "mx_expr.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>

using namespace synergia;
using namespace std;


//===========================================================================
// Static initializer

double MadX::nan = std::numeric_limits<double>::quiet_NaN();
string MadX::nst = string("!!!!*&^%&*IMANULLSTRING(*&*&^^^");


//===========================================================================
// Helper functions

namespace
{
  double madx_nan = MadX::nan;
  string madx_nst = MadX::nst;

  mx_expr
    retrieve_expr_from_map( value_map_t const & m
                            , string_t const & k )
  {
    string_t key(k);
    std::transform(key.begin(), key.end(), key.begin(), ::tolower);

    value_map_t::const_iterator it = m.find(key);

    if( it == m.end() )
      throw std::runtime_error( "retrieve number: cannot find attribute with name " + key);

   if( it->second.type!=NUMBER && it->second.type!=DEFFERED_NUMBER )
      throw std::runtime_error( "the requested key '" + k + "' cannot be retrieved as an expr" );

    return std::any_cast<mx_expr>(it->second.value);
  }

  std::vector<mx_expr>
    retrieve_expr_seq_from_map( value_map_t const & m
                            , string_t const & k )
  {
    string_t key(k);
    std::transform(key.begin(), key.end(), key.begin(), ::tolower);

    value_map_t::const_iterator it = m.find(key);

    if( it == m.end() )
      throw std::runtime_error( "retrieve number: cannot find attribute with name " + key);

   if( it->second.type!=ARRAY && it->second.type!=DEFFERED_ARRAY)
      throw std::runtime_error( "the requested key '" + k + "' cannot be retrieved as an expr_seq" );

    return std::any_cast<std::vector<mx_expr>>(it->second.value);
  }


  string_t
    retrieve_string_from_map( value_map_t const & m
                            , string_t const & k 
                            , string_t const & def )
  {
    string_t key(k);
    std::transform(key.begin(), key.end(), key.begin(), ::tolower);

    value_map_t::const_iterator it = m.find(key);
    if( it!=m.end() )
    {
      return std::any_cast<string_t>(it->second.value);
    }
    else
    {
      if( def==madx_nst )
        throw std::runtime_error( "retrieve string: cannot find attribute with name " + key);
      else
        return def;
    }
  }

  double
    retrieve_number_from_map( value_map_t const & m
                            , string_t const & k
                            , MadX const & global
                            , double def )
  {
    string_t key(k);
    std::transform(key.begin(), key.end(), key.begin(), ::tolower);

    value_map_t::const_iterator it = m.find(key);
    if( it!=m.end() )
    {
      if( it->second.type == NUMBER )
      {
        mx_expr e = std::any_cast<mx_expr>(it->second.value);
        return boost::apply_visitor(mx_calculator(global, def), e);
      }
      else
      {
        throw std::runtime_error( "the requested key '" + k + "' cannot be retrieved as a number" );
      }
    }
    else
    {
      if( std::isnan(def) )
        throw std::runtime_error( "retrieve number: cannot find attribute with name " + key);
      else
        return def;
    }
  }

  std::vector<double>
    retrieve_number_seq_from_map( value_map_t const & m
                                , string_t const & k
                                , MadX const & global
                                , double def )
  {
    string_t key(k);
    std::transform(key.begin(), key.end(), key.begin(), ::tolower);

    value_map_t::const_iterator it = m.find(key);
    if( it!=m.end() )
    {
     if( it->second.type == ARRAY )
      {
        mx_exprs es = std::any_cast<mx_exprs>(it->second.value);
        std::vector<double> vd;

        for( mx_exprs::const_iterator it = es.begin()
           ; it != es.end(); ++it )
        {
          vd.push_back( boost::apply_visitor(mx_calculator(global, def), *it) );
        }
        return vd;
      }
      else
      {
        throw std::runtime_error( "the requested key '" + k + "' cannot be retrieved as a sequence of number" );
      }
    }
    else
    {
      if( std::isnan(def) )
        throw std::runtime_error( "retrieve number sequence: cannot find attribute with name " + key);
      else
      {
        vector<double> r(1); 
        r[0] = def;
        return r;
      }
    }
  }

  MadX_command
    resolve_command(MadX_command const & cmd, MadX const & mx, bool resolve)
  {
    if( !resolve || !cmd.is_reference() )  return cmd;

    MadX_command result = mx.command(cmd.name());
    result.merge_with_overwrite(cmd);
    return result;
  }

}


//===========================================================================
// free functions
//
std::string
synergia::to_string(MadX_value const& val)
{
    std::string str;

    switch(val.type)
    {
    case STRING:
        str += "\"" + std::any_cast<std::string>(val.value) + "\"";
        break;

    case NUMBER:
    case DEFFERED_NUMBER:
        str += mx_expr_str(std::any_cast<mx_expr>(val.value));
        break;

    case ARRAY:
    case DEFFERED_ARRAY:
        str += "{";
        for(auto const& expr : std::any_cast<mx_exprs>(val.value))
            str += mx_expr_str(expr) + ",";
        if(str.back() == ',') str.pop_back();
        str += "}";

    default:
        break;
    }

    return str;
}


//===========================================================================
// MadX_command

string_t
  MadX_command::name() const
{
  return name_;
}

string_t
  MadX_command::label() const
{
  return label_;
}

size_t
  MadX_command::attribute_count() const
{
  return attributes_.size();
}

std::vector<string_t>
  MadX_command::attribute_names() const
{
  std::vector<string_t> names;
  for( value_map_t::const_iterator it = attributes_.begin()
     ; it!= attributes_.end(); ++it )
    names.push_back(it->first);
  return names;
}

MadX_value_type
  MadX_command::attribute_type( string_t const & name ) const
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  value_map_t::const_iterator it = attributes_.find(key);
  if( it!=attributes_.end() )
  {
    return it->second.type;
  }
  else
  {
    throw std::runtime_error( "MadX_command::attribute_type:"
        " cannot find attribute with name " + key);
  }
}

string_t
  MadX_command::attribute_as_string( string_t const & name ) const
{
  return retrieve_string_from_map( attributes_, name, madx_nst );
}

string_t
  MadX_command::attribute_as_string( string_t const & name, string_t const & def ) const
{
  return retrieve_string_from_map( attributes_, name, def );
}

double
  MadX_command::attribute_as_number( string_t const & name ) const
{
  return retrieve_number_from_map( attributes_, name, *mx, madx_nan );
}

double
  MadX_command::attribute_as_number( string_t const & name, double def ) const
{
  return retrieve_number_from_map(attributes_, name, *mx, def);
}

bool
  MadX_command::attribute_as_boolean( string_t const & name ) const
{
  return std::abs(retrieve_number_from_map(attributes_, name, *mx, madx_nan)) > 1e-10;
}

mx_expr
  MadX_command::attribute_as_expr( string_t const& name ) const
{
  return retrieve_expr_from_map(attributes_, name);
}

std::vector<mx_expr>
  MadX_command::attribute_as_expr_seq( string_t const& name ) const
{
  return retrieve_expr_seq_from_map(attributes_, name);
}

std::vector<double>
  MadX_command::attribute_as_number_seq( string_t const & name ) const
{
  return retrieve_number_seq_from_map(attributes_, name, *mx, madx_nan);
}

std::vector<double>
  MadX_command::attribute_as_number_seq( string_t const & name, double def ) const
{
  return retrieve_number_seq_from_map(attributes_, name, *mx, def);
}

void
  MadX_command::set_parent( MadX const & parent )
{
  mx = &parent;
}

void
  MadX_command::set_name( string_t const & name, MadX_command_type type )
{
  name_ = name;
  std::transform(name_.begin(), name_.end(), name_.begin(), ::tolower);
  type_ = type;
}

void
  MadX_command::set_label( string_t const & label)
{
  label_ = label;
}

MadX_command_type
  MadX_command::type() const
{
  return type_;
}

bool
  MadX_command::is_element() const
{
  return type_ == ELEMENT;
}

bool
  MadX_command::is_reference() const
{
  return type_ == ELEMENT_REF;
}

bool
  MadX_command::is_command() const
{
  return type_ == EXECUTABLE;
}

void
  MadX_command::insert_attribute( string_t const & name, string_t const & value )
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  MadX_value v;
  v.value = std::any(value);
  v.type  = value.empty() ? NONE : STRING;

  attributes_.insert_or_assign(key, v);
}

void
  MadX_command::insert_attribute( string_t const & name, mx_expr const & value )
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  MadX_value v;
  v.value = std::any(value);
  v.type  = NUMBER;

  attributes_.insert_or_assign(key, v);
}

void
  MadX_command::insert_attribute( string_t const & name, mx_exprs const & value )
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  MadX_value v;
  v.value = std::any(value);
  v.type  = ARRAY;

  attributes_.insert_or_assign(key, v);
}

void
  MadX_command::merge_with_overwrite(MadX_command const & other)
{
  value_map_t map(other.attributes_);
  map.insert(attributes_.begin(), attributes_.end());
  std::swap(attributes_, map);
}

void
  MadX_command::merge(MadX_command const & other)
{
  attributes_.insert(other.attributes_.begin(), other.attributes_.end());
}

std::string
MadX_command::to_string() const
{
    std::string cmd;
    cmd += name() + ",";

    for(auto const& var : attributes_)
        cmd += var.first + "=" + ::to_string(var.second) + ",";

    if (cmd.back() == ',') cmd.pop_back();

    return cmd;
}


//===========================================================================
// MadX_line
size_t
  MadX_line::element_count() const
{
  return elements_.size();
}

string_t
  MadX_line::element_name(size_t idx) const
{
  return elements_[idx];
}

MadX_command
  MadX_line::element(size_t idx, bool resolve) const
{
  return parent.command(elements_[idx], resolve);
}

void
  MadX_line::insert_element(string_t const & ele)
{
  string_t e(ele);
  std::transform(e.begin(), e.end(), e.begin(), ::tolower);

  elements_.push_back(e);
}

void
  MadX_line::print() const
{
  cout << "( ";
  for( vector<string>::const_iterator it=elements_.begin()
     ; it!=elements_.end(); ++it )
  {
    cout << *it;
    if( it+1 != elements_.end() ) cout << ", ";
  }
  cout << " )\n";
}

//===========================================================================
// MadX_sequence

string_t
  MadX_sequence::label() const
{
  return lbl;
}

double
  MadX_sequence::length() const
{
  return l;
}

size_t
  MadX_sequence::element_count() const
{
  return seq_.size();
}

MadX_command
  MadX_sequence::element(size_t idx) const
{
  return parent.command(seq_[idx].label);
}

double
  MadX_sequence::element_at(size_t idx) const
{
  return boost::apply_visitor(mx_calculator(parent), seq_[idx].at);
}

double
  MadX_sequence::element_from(size_t idx) const
{
  return seq_[idx].from;
}

string_t
  MadX_sequence::element_label(size_t idx) const
{
  return seq_[idx].label;
}

MadX_entry_type
  MadX_sequence::element_type(size_t idx) const
{
//  MadX_command const & cmd = element(idx);
//  std::string key = cmd.label();
//  if( key.empty() ) key = cmd.name();
  std::string key = seq_[idx].label;
  return parent.entry_type(key);
}

MadX_sequence_refer
  MadX_sequence::refer() const
{
  return r;
}

string
  MadX_sequence::refpos() const
{
  return rp;
}

void
  MadX_sequence::set_label(string_t const & label)
{
  lbl = label;
}

void
  MadX_sequence::set_length(double length)
{
  l = length;
}

void
  MadX_sequence::set_refer(MadX_sequence_refer refer)
{
  r = refer;
}

void
  MadX_sequence::set_refpos(string const & refpos)
{
  rp = refpos;
}

void
  MadX_sequence::add_element(string_t const & label, mx_expr const& at, string_t const & from)
{
  seq_.push_back( seq_element(label, at, from) );
}

void
  MadX_sequence::finalize()
{
  // resolve all the 'from' references
  for (seq_ele_v_t::iterator it = seq_.begin(); it != seq_.end(); ++it)
  {
    if (!it->from_str.empty())
    {
      bool found = false;

      for (seq_ele_v_t::const_iterator cit = seq_.begin(); cit != seq_.end(); ++cit)
      {
        if (cit == it) continue;

        if (cit->label == it->from_str)
        {
          double at = boost::apply_visitor(mx_calculator(parent), cit->at);
          it->from = at;
          found = true;
          break;
        }
      }

      if (!found)
        throw runtime_error("fatal: 'from' reference to unknown element: " + it->from_str);
    }
  }
}

void
  MadX_sequence::reset()
{
  lbl = string_t();
  l = 0.0;
  r = SEQ_REF_CENTRE;
  rp = string_t();
  seq_.clear();
}

void
  MadX_sequence::print() const
{
  cout << "\n";
  cout << "  l = " << length() << ", elements = "
       << element_count() << ", refer = ";

  switch(refer())
  {
  case SEQ_REF_ENTRY:  cout << "entry"; break;
  case SEQ_REF_CENTRE: cout << "centre"; break;
  case SEQ_REF_EXIT:   cout << "exit"; break;
  default: break;
  }

  cout << ", refpos = " << refpos();
  cout << ",\n";

  for(auto const& ele : seq_)
  {
      cout << "  " << ele.label 
           << ", at = " << boost::apply_visitor(mx_calculator(parent), ele.at)
           << ", from = " << ele.from
           << ", from_ref = " << ele.from_str
           << "\n";
  }

  cout << "\n";
}


//===========================================================================
// MadX

MadX::MadX()
  : variables_()
  , cmd_seq_()
  , cmd_map_()
  , lines_()
  , seqs_()
  , cur_seq_(*this)
  , building_seq_(false)
{ }

// lines and sequences are not copied for the moment
// for the reason that the copied MadX objects are only used in the
// Lattice_tree for referencing, and one cannot reference an element
// defined in a line or a sequence.
//
// It is, however, more appropriate to have a proper data structure
// that only contains the references (variables and commands) for
// the Lattice_tree, and leaves the MadX object as an intermediate
// data strucuture for madx parsing, and remains non-copyable
MadX::MadX(MadX const& o) 
    : variables_(o.variables_)
    , cmd_seq_(o.cmd_seq_)
    , cmd_map_(o.cmd_map_)
    , lines_()
    , seqs_()
    , cur_seq_(*this)
    , building_seq_(false)
{
    for(auto& cmd : cmd_seq_) cmd.set_parent(*this);
    for(auto& cmd : cmd_map_) cmd.second.set_parent(*this);
}

MadX&
MadX::operator=(MadX const& o) 
{
    variables_ = o.variables_;
    cmd_seq_ = o.cmd_seq_;
    cmd_map_ = o.cmd_map_;

    lines_.clear();
    seqs_.clear();

    cur_seq_.reset();
    building_seq_ = false;

    for(auto& cmd : cmd_seq_) cmd.set_parent(*this);
    for(auto& cmd : cmd_map_) cmd.second.set_parent(*this);

    return *this;
}


string_t
  MadX::variable_as_string( string_t const & name ) const
{
  return retrieve_string_from_map( variables_, name, madx_nst );
}

string_t
  MadX::variable_as_string( string_t const & name, string_t const & def ) const
{
  return retrieve_string_from_map( variables_, name, def );
}

double
  MadX::variable_as_number( string_t const & name ) const
{
  return retrieve_number_from_map( variables_, name, *this, madx_nan );
}

double
  MadX::variable_as_number( string_t const & name, double def ) const
{
  return retrieve_number_from_map( variables_, name, *this, def );
}

bool
  MadX::variable_as_boolean( string_t const & name ) const
{
  return std::abs(retrieve_number_from_map( variables_, name, *this, madx_nan )) > 1e-10;
}

std::vector<double>
  MadX::variable_as_number_seq( string_t const & name ) const
{
  return retrieve_number_seq_from_map( variables_, name, *this, madx_nan );
}

std::vector<double>
  MadX::variable_as_number_seq( string_t const & name, double def ) const
{
  return retrieve_number_seq_from_map( variables_, name, *this, def );
}

size_t
  MadX::command_count() const
{
  return cmd_seq_.size();
}

std::vector<string_t >
  MadX::commands() const
{
    std::vector<string_t > commands;
    for(commands_v_t::const_iterator it = cmd_seq_.begin();
            it != cmd_seq_.end(); ++it)
    {
        commands.push_back(it->name());
    }
    return commands;
}

MadX_command
  MadX::command( size_t idx, bool resolve ) const
{
  return resolve_command(cmd_seq_[idx], *this, resolve);
}

size_t
  MadX::label_count() const
{
  return cmd_map_.size();
}

std::vector<string_t>
  MadX::command_labels() const
{
  std::vector<string_t> labels;
  for( commands_m_t::const_iterator it = cmd_map_.begin()
     ; it!=cmd_map_.end(); ++it )
    labels.push_back(it->first);
  return labels;
}

MadX_command
  MadX::command( string_t const & label, bool resolve ) const
{
  string_t key(label);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  commands_m_t::const_iterator it = cmd_map_.find(key);
  if( it!=cmd_map_.end() )
  {
    return resolve_command(it->second, *this, resolve);
  }
  else
  {
    throw std::runtime_error( "cannot find command with label " + key );
  }
}

MadX_command&
MadX::command_ref(string_t const& label)
{
  string_t key(label);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  commands_m_t::iterator it = cmd_map_.find(key);
  if( it!=cmd_map_.end() )
  {
    return it->second;
  }
  else
  {
    throw std::runtime_error( "cannot find command with label " + key );
  }
}


size_t
  MadX::line_count() const
{
  return lines_.size();
}

std::vector<string_t>
  MadX::line_labels() const
{
  std::vector<string_t> labels;
  for( lines_m_t::const_iterator it = lines_.begin()
     ; it!=lines_.end(); ++it)
    labels.push_back(it->first);
  return labels;
}

MadX_line const &
  MadX::line( string_t const & line ) const
{
  string_t key(line);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  lines_m_t::const_iterator it = lines_.find(key);

  if( it!=lines_.end() )
    return it->second;
  else
    throw std::runtime_error( "cannot find line with label " + line );
}

size_t
  MadX::sequence_count() const
{
  return seqs_.size();
}

std::vector<string_t>
  MadX::sequence_labels() const
{
  std::vector<string_t> labels;
  for( sequences_m_t::const_iterator it = seqs_.begin()
     ; it!=seqs_.end(); ++it)
    labels.push_back(it->first);
  return labels;
}

MadX_sequence const &
  MadX::sequence( string_t const & seq ) const
{
  string_t key(seq);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  sequences_m_t::const_iterator it = seqs_.find(key);

  if( it!=seqs_.end() )
    return it->second;
  else
    throw std::runtime_error( "cannot find sequence with label " + seq );
}

MadX_sequence const &
  MadX::current_sequence( ) const
{
  return cur_seq_;
}

MadX_sequence &
  MadX::current_sequence( )
{
  return cur_seq_;
}

bool
  MadX::building_sequence() const
{
  return building_seq_;
}

MadX_entry_type
  MadX::entry_type(string_t const & entry) const
{
  string_t key(entry);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  if( variables_.find(key) != variables_.end() )
    return ENTRY_VARIABLE;

  if( seqs_.find(key) != seqs_.end() )
    return ENTRY_SEQUENCE;

  if( lines_.find(key) != lines_.end() )
    return ENTRY_LINE;

  if( cmd_map_.find(key) != cmd_map_.end() )
    return ENTRY_COMMAND;

  return ENTRY_NULL;
}

void
  MadX::insert_variable(string_t const & name, string_t const & value)
{
  string_t key(name);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  MadX_value v;
  v.value = std::any(value);
  v.type  = value.empty() ? NONE : STRING;

  variables_[key] = v;
}

void
  MadX::insert_variable(string_t const & name, mx_expr const & value)
{
  string_t key(name);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  MadX_value v;
  v.value = std::any(value);
  v.type  = NUMBER;

  variables_[key] = v;
}

void
  MadX::insert_variable(string_t const & name, mx_exprs const & value)
{
  string_t key(name);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  MadX_value v;
  v.value = std::any(value);
  v.type  = ARRAY;

  variables_[key] = v;
}

void
  MadX::insert_label(string_t const & name, MadX_command const & cmd)
{
  string_t key(name);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  cmd_map_[key] = cmd; // always overwrite
  cmd_map_[key].set_parent(*this);
}

void
  MadX::fuse_command(string_t const & name, MadX_command const & cmd)
{
  string_t key(name);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  commands_m_t::iterator it = cmd_map_.find( key );
  if( it!=cmd_map_.end() )
  {
    it->second.merge_with_overwrite( cmd );
  }
}

void
  MadX::insert_line(string_t const & name, MadX_line const & line)
{
  string_t key(name);
  std::transform( key.begin(), key.end(), key.begin(), ::tolower );

  lines_.insert( std::make_pair(key, line) );
}

void
  MadX::insert_command(MadX_command const & cmd)
{
  cmd_seq_.push_back(cmd);
  cmd_seq_.back().set_parent(*this);
}

void
  MadX::start_sequence( string_t const & label
                      , double length
                      , string_t const & refer
                      , string_t const & refpos )
{
  building_seq_ = true;
  cur_seq_.set_label(label);
  cur_seq_.set_length(length);
  cur_seq_.set_refpos(refpos);

  if( refer=="entry"  )      cur_seq_.set_refer(SEQ_REF_ENTRY);
  else if( refer=="center" ) cur_seq_.set_refer(SEQ_REF_CENTRE);
  else if( refer=="exit" )   cur_seq_.set_refer(SEQ_REF_EXIT);
  else                       cur_seq_.set_refer(SEQ_REF_CENTRE);
}

void
  MadX::append_sequence_element(string_t const & label, mx_expr const& at, string_t const & from)
{
  if( building_seq_ )
    cur_seq_.add_element(label, at, from);
}

void
  MadX::end_sequence()
{
  // done and insert the sequence to madx object
  building_seq_ = false;
  cur_seq_.finalize();
  seqs_.insert(std::make_pair(cur_seq_.label(), cur_seq_));
  cur_seq_.reset();
}

void
MadX::print() const
{
    for(auto const& var : variables_)
    { 
        std::cout << var.first << " = " 
            << to_string(var.second) << "\n";
    }

    for(auto const& cmd : cmd_map_)
    {
        std::cout << cmd.first << ": "
            << cmd.second.to_string() << "\n";
    }

    for(auto const& line : lines_)
    {
        std::cout << line.first << " : line = ";
        line.second.print();
    }

    for(auto const& seq : seqs_ )
    {
        cout << seq.first << " : sequence = ";
        seq.second.print();
    }
}

std::string
MadX::export_variables() const
{
    std::string vars;
    
    for(auto const& var : variables_)
        vars += var.first + "=" + to_string(var.second) + ";\n";

    return vars;
}

std::string 
MadX::export_unnamed_cmds() const
{
    std::string cmds;

    for(auto const& cmd : cmd_seq_)
    {
        if (cmd.type() == ELEMENT || cmd.type() == ELEMENT_REF)
            cmds += cmd.to_string() + ";\n";
    }

    return cmds;
}

std::string 
MadX::export_labled_cmds() const
{
    std::string cmds;

    for(auto const& cmd : cmd_map_)
    {
        if (cmd.second.type() == ELEMENT || cmd.second.type() == ELEMENT_REF)
            cmds += cmd.first + ":" + cmd.second.to_string() + ";\n";
    }

    return cmds;
}

std::string 
MadX::to_madx() const
{
    std::string madx;

    madx += export_variables();
    madx += export_unnamed_cmds();
    madx += export_labled_cmds();

    return madx;
}


