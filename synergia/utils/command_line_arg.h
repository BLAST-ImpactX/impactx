#ifndef COMMAND_LINE_ARG_H_
#define COMMAND_LINE_ARG_H_
#include <sstream>
#include <stdexcept>
#include <string>

class Command_line_arg {
private:
  std::string str;

public:
  Command_line_arg(const char* char_ptr);
  bool is_equal_pair() const;
  std::string get_lhs() const;
  std::string get_rhs() const;
  template <typename T>
  T
  extract_value() const
  {
    std::stringstream str_stream(get_rhs());
    T retval;
    str_stream >> retval;
    return retval;
  }
};

#endif /* COMMAND_LINE_ARG_ */
