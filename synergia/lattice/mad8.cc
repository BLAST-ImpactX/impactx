#include <any>
#include <stdexcept>

#include "mad8.h"

using namespace synergia;

//===========================================================================
// Mad8_command

string_t
Mad8_command::name() const
{
  return name_;
}

size_t
Mad8_command::attribute_count() const
{
  return attributes_.size();
}

string_t
Mad8_command::attribute_as_string(string_t const& name) const
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  value_map_t::const_iterator it = attributes_.find(key);
  if (it != attributes_.end()) {
    return std::any_cast<string_t>(it->second.value);
  } else {
    throw std::runtime_error("cannot find attribute with name " + key);
  }
}

double
Mad8_command::attribute_as_number(string_t const& name) const
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  value_map_t::const_iterator it = attributes_.find(key);
  if (it != attributes_.end()) {
    return std::any_cast<double>(it->second.value);
  } else {
    throw std::runtime_error("cannot find attribute with name " + key);
  }
}

void
Mad8_command::set_name(string_t const& name)
{
  name_ = name;
  std::transform(name_.begin(), name_.end(), name_.begin(), ::tolower);
}

void
Mad8_command::insert_attribute(string_t const& name, double value)
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  Mad8_value v;
  v.value = std::any(value);
  v.type = NUMBER;

  attributes_.insert(std::make_pair(key, v));
}

void
Mad8_command::insert_attribute(string_t const& name, string_t const& value)
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  Mad8_value v;
  v.value = std::any(value);
  v.type = value.empty() ? NONE : STRING;

  attributes_.insert(std::make_pair(key, v));
}

//===========================================================================
// Mad8_line

size_t
Mad8_line::element_count() const
{
  return index_.size();
}

Mad8_line_element_type
Mad8_line::element_type(size_t idx) const
{
  return index_[idx].first;
}

string_t
Mad8_line::element_as_string(size_t idx) const
{
  return elements_[index_[idx].second];
}

Mad8_line
Mad8_line::element_as_line(size_t idx) const
{
  return lines_[index_[idx].second];
}

void
Mad8_line::insert_operator(string_t const& op)
{
  elements_.push_back(op);
  size_t idx = elements_.size() - 1;

  Mad8_line_element_type t;
  if (op.compare("-") == 0)
    t = MINUS;
  else if (op[op.size() - 1] == '*')
    t = MULTIPLIER;
  else
    throw std::runtime_error("Unrecognized Mad8 line operator: " + op);

  index_.push_back(std::make_pair(t, idx));
}

void
Mad8_line::insert_element(string_t const& ele)
{
  string_t e(ele);
  std::transform(e.begin(), e.end(), e.begin(), ::tolower);

  elements_.push_back(e);
  size_t idx = elements_.size() - 1;
  index_.push_back(std::make_pair(LABEL, idx));
}

void
Mad8_line::insert_subline(Mad8_line const& line)
{
  lines_.push_back(line);
  size_t idx = lines_.size() - 1;
  index_.push_back(std::make_pair(LINE, idx));
}

//===========================================================================
// Mad8

string_t
Mad8::variable_as_string(string_t const& name) const
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  value_map_t::const_iterator it = variables_.find(key);

  if (it != variables_.end())
    return std::any_cast<string_t>(it->second.value);
  else
    throw std::runtime_error("cannot find variable with name " + key);
}

double
Mad8::variable_as_number(string_t const& name) const
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  value_map_t::const_iterator it = variables_.find(key);

  if (it != variables_.end())
    return std::any_cast<double>(it->second.value);
  else
    throw std::runtime_error("cannot find variable with name " + key);
}

size_t
Mad8::command_count() const
{
  return cmd_seq_.size();
}

Mad8_command const&
Mad8::command(size_t idx) const
{
  return cmd_seq_[idx];
}

size_t
Mad8::label_count() const
{
  return cmd_map_.size();
}

Mad8_command const&
Mad8::label(string_t const& label) const
{
  string_t key(label);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  commands_m_t::const_iterator it = cmd_map_.find(key);

  if (it != cmd_map_.end())
    return it->second;
  else
    throw std::runtime_error("cannot find command with label " + key);
}

size_t
Mad8::line_count() const
{
  return lines_.size();
}

Mad8_line const&
Mad8::line(string_t const& line) const
{
  string_t key(line);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  lines_m_t::const_iterator it = lines_.find(line);

  if (it != lines_.end())
    return it->second;
  else
    throw std::runtime_error("cannot find line with label " + line);
}

void
Mad8::insert_variable(string_t const& name, double value)
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  Mad8_value v;
  v.value = std::any(value);
  v.type = NUMBER;

  variables_.insert(std::make_pair(key, v));
}

void
Mad8::insert_variable(string_t const& name, string_t const& value)
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  Mad8_value v;
  v.value = std::any(value);
  v.type = value.empty() ? NONE : STRING;

  variables_.insert(std::make_pair(key, v));
}

void
Mad8::insert_label(string_t const& name, Mad8_command const& cmd)
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  cmd_map_.insert(std::make_pair(key, cmd));
}

void
Mad8::insert_line(string_t const& name, Mad8_line const& line)
{
  string_t key(name);
  std::transform(key.begin(), key.end(), key.begin(), ::tolower);

  lines_.insert(std::make_pair(key, line));
}

void
Mad8::insert_command(Mad8_command const& cmd)
{
  cmd_seq_.push_back(cmd);
}
