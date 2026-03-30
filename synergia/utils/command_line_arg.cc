#include "command_line_arg.h"

Command_line_arg::Command_line_arg(const char* char_ptr) : str(char_ptr) {}

bool
Command_line_arg::is_equal_pair() const
{
  return str.find('=') != std::string::npos;
}

std::string
Command_line_arg::get_lhs() const
{
  return str.substr(0, str.find('='));
}

std::string
Command_line_arg::get_rhs() const
{
  return str.substr(str.find('=') + 1, std::string::npos);
}

template <>
bool
Command_line_arg::extract_value<bool>() const
{
  bool retval;
  if ((get_rhs() == "True") || (get_rhs() == "true") || (get_rhs() == "t") ||
      (get_rhs() == "1")) {
    retval = true;
  } else if ((get_rhs() == "False") || (get_rhs() == "false") ||
             (get_rhs() == "f") || (get_rhs() == "0") || (get_rhs() == "nil")) {
    retval = false;
  } else {
    throw std::runtime_error("Cannot convert " + get_rhs() + " to bool");
  }
  return retval;
}
