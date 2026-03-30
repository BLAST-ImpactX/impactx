#ifndef STRING_UTILS_H_
#define STRING_UTILS_H_
#include <string>

inline bool
false_string(std::string const& val)
{
  return (val == "false") || (val == "False") || (val == "FALSE");
}

#endif /* STRING_UTILS_H_ */
