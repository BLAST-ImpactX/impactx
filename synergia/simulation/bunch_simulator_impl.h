
#ifndef BUNCH_SIMULATOR_IMPL_H
#define BUNCH_SIMULATOR_IMPL_H

#include <vector>

namespace impl {

    void divide_bunches(int size,
                        size_t num_bunches_pri,
                        size_t num_bunches_sec,
                        std::vector<int>& p_ranks,
                        std::vector<int>& s_ranks);

}

#endif
