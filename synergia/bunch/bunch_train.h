#ifndef BUNCH_TRAIN_H_
#define BUNCH_TRAIN_H_

#include "synergia/bunch/bunch.h"

#include <cereal/types/array.hpp>
#include <cereal/types/vector.hpp>

class Bunch_train {

  private:
    std::shared_ptr<Commxx> comm;

    std::vector<Bunch> bunches;
    std::vector<double> spacings;

    int index;

    int num_bunches;
    int num_buckets;

    // bunch_idx_map[bunch] stores the index in the local bunch array.
    // the array has the dim of num_bunches, if the bunch is not present in
    // the current rank, the stored value is -1
    std::vector<int> bunch_idx_map;

    /// counts and offsets are needed for impedance,
    /// counts.size()=offsets.size()=num_procs they are meaningfull only on the
    /// local rank=0 of every bunch communicator for example:  3 bunches on 5
    /// processors: proc 2 and 4 has local rank different form zero
    ///              (count[0]=1, offset[0]=0), (count[1]=1, offset[3]=1),
    ///              (count[3]=1, offset[3]=2)
    /// for example:  5 bunches on 3 processors:
    ///                (count[0]=1, offset[0]=0), (count[1]=2, offset[3]=1),
    ///                (count[2]=2, offset[2]=3)
    std::vector<int> proc_counts_imped;
    std::vector<int> proc_offsets_imped;

    void set_bucket_indices();
    void find_parent_comm();
    void calculates_counts_and_offsets_for_impedance();

  public:
    // Bunch_train(Bunches const& bunches, double spacing);
    // Bunch_train(Bunches const& bunches, std::vector<double > const&
    // spacings);

    Bunch_train(Reference_particle const& ref,
                size_t num_bunches,
                size_t num_particles_per_bunch,
                double num_real_particles_per_bunch,
                double spacing,
                Commxx const& comm = Commxx(),
                size_t num_spectator_per_bunch = 0,
                int index = 0);

    // Get the communicator
    Commxx const&
    get_comm() const
    {
        return *comm;
    }

    // total number of bunches in the train. it is not always
    // the same as the size of the bunch array
    int
    get_num_bunches() const
    {
        return num_bunches;
    }

    // number of local bunches on this rank. alias for
    // get_bunch_array_size()
    int
    get_num_local_bunches() const
    {
        return get_bunch_array_size();
    }

    // number of buckets. it is not necessary for the bunches to
    // occupy consecutive buckets, so there might be some empty
    // buckets in-between bunches
    // for now it is always the same as the number of bunches
    int
    get_num_buckets() const
    {
        return num_buckets;
    }

    // index of the train (0 or 1)
    int
    get_index() const
    {
        return index;
    }

    // number of bunches on this rank (size of the bunch array)
    size_t
    get_bunch_array_size() const
    {
        return bunches.size();
    }

    // get the local array index from the bunch index in the train
    // returns -1 if the bunch is not available locally
    int
    get_array_idx_of_bunch(int bunch_idx) const
    {
        return bunch_idx_map[bunch_idx];
    }

    // get the bunch index in the train from the local array index
    int
    get_train_idx_of_bunch(int array_idx) const
    {
        return bunches[array_idx].get_bunch_index();
    }

    // access the local bunch
    Bunch&
    operator[](size_t array_idx)
    {
        return bunches[array_idx];
    }

    Bunch const&
    operator[](size_t array_idx) const
    {
        return bunches[array_idx];
    }

    // access the local bunch array
    std::vector<Bunch>&
    get_bunches()
    {
        return bunches;
    }

    std::vector<Bunch> const&
    get_bunches() const
    {
        return bunches;
    }

    std::vector<double>& get_spacings();

    // update the total particle number for all bunches in the bunch train
    // note that calling each bunch's update_total_num() wont do the actual
    // work if the caller's rank is not part of the bunch's communicator.
    // on the other hand, calling bunch_train's update_bunch_total_num()
    // gurantees that all bunches are updated for all ranks
    void update_bunch_total_num();

    void
    set_longitudinal_boundary(LongitudinalBoundary lb, double param = 0.0)
    {
        for (auto& b : bunches)
            b.set_longitudinal_boundary(lb, param);
    }

    std::vector<int>& get_proc_counts_for_impedance();
    std::vector<int>& get_proc_offsets_for_impedance();

  private:
    friend class cereal::access;

    template <class Archive>
    void
    serialize(Archive& ar)
    {
        ar(CEREAL_NVP(comm));
        ar(CEREAL_NVP(bunches));
        ar(CEREAL_NVP(spacings));
        ar(CEREAL_NVP(index));
        ar(CEREAL_NVP(num_bunches));
        ar(CEREAL_NVP(num_buckets));
        ar(CEREAL_NVP(bunch_idx_map));
        ar(CEREAL_NVP(proc_counts_imped));
        ar(CEREAL_NVP(proc_offsets_imped));
    }
};

#endif /* BUNCH_TRAIN_H_ */
