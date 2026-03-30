
#include <fstream>

#include "synergia/utils/json.h"

namespace syn {
    void
    checkpoint_save_as_json(std::string const& prop_str,
                            std::string const& sims_str,
                            std::vector<int> const& displs,
                            std::vector<int> const& lens)
    {
        syn::json cp = syn::json::object();

        cp["propagator"] = syn::json::parse(prop_str);
        cp["simulator"] = syn::json::array();

        for (int i = 0; i < displs.size(); ++i) {
            auto begin = sims_str.begin() + displs[i];
            auto d = syn::json::parse(begin, begin + lens[i]);
            cp["simulator"].push_back(d);
        }

        // write states to file
        std::ofstream file("cp_state.json");
        if (!file.good())
            throw std::runtime_error("Error at creating checkpointing file");

        file << cp;
    }

    std::pair<std::string, std::string>
    checkpoint_load_json(std::vector<char> const& buf, int rank)
    {
        auto cp = syn::json::parse(buf.begin(), buf.end());
        return std::make_pair(cp["propagator"].dump(),
                              cp["simulator"][rank].dump());
    }
}
