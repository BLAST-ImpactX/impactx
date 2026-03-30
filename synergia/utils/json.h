#ifndef SYNERGIA_UTILS_JSON_H
#define SYNERGIA_UTILS_JSON_H

#include "synergia/utils/json.hpp"

namespace syn {
  using json = nlohmann::json;

  struct dummy_json {
    int
    dump(int indent = -1)
    {
      return 0;
    }

    friend std::ostream&
    operator<<(std::ostream& os, dummy_json const&)
    {
      return os;
    }
  };
}

#endif
