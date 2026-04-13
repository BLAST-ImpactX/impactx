#include "bunch_train.h"
#include "synergia/utils/parallel_utils.h"

#if 0
void 
Bunch_train::find_parent_comm()
{
  try{
    if (bunches.size()>0) {
        // check if all bunches has the same parent communicator
        MPI_Comm comm_test=bunches[0]->get_comm().get_parent_sptr()->get();
	for (size_t i = 1; i < bunches.size(); ++i) {
 	      int result;
 	      MPI_Comm_compare( comm_test, bunches[i]->get_comm().get_parent_sptr()->get(), &result);
 	      if (result != MPI_IDENT) {
             	throw std::runtime_error("Bunch_train, find_parent_comm_sptr: bunches have different parrent comunicator");
 	      } 
 	} 
 	parent_comm_sptr=bunches[0]->get_comm().get_parent_sptr();
    }
    else{
	 throw std::runtime_error(
	    "Bunch_train, find_parent_comm_sptr: number of bunches is zero, there is no commuicator");
    }   
    has_parent_comm=true;
  }
  catch (std::exception const& e) {
         std::cout<<e.what()<<std::endl;
         MPI_Abort(MPI_COMM_WORLD, 333);
  }  
}
#endif

void
Bunch_train::calculates_counts_and_offsets_for_impedance()
{
#if 0
   
  try{
     if (!has_parent_comm) find_parent_comm_sptr(); 
     int size_parent=parent_comm_sptr->get_size();
     proc_counts_imped.resize(size_parent);
     proc_offsets_imped.resize(size_parent); 
 	  counts_and_offsets_for_impedance(*parent_comm_sptr, bunches.size(), proc_offsets_imped, proc_counts_imped);
  }
  catch (std::exception const& e) {
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 333);
  }
#endif
}

#if 0
Commxx
Bunch_train::get_parent_comm()
{
  if (!has_parent_comm) find_parent_comm(); 
  return parent_comm;
}
#endif

void
Bunch_train::set_bucket_indices()
{
#if 0
    std::list<int > found_indices;
    for (size_t i = 0; i < bunches.size(); ++i) {
        if (!bunches[i]->is_bucket_index_assigned()) {
            bunches[i]->set_bucket_index(i);
        }
        for (std::list<int >::const_iterator it = found_indices.begin();
                it != found_indices.end(); ++it) {
            if (*it >= bunches[i]->get_bucket_index()) {
                throw std::runtime_error(
                        "Bunch_train: bunch bucket indices must be either in strictly increasing order or all zero; otherwise wake field does not work");
            }
        }
        found_indices.push_back(bunches[i]->get_bucket_index());
    }
#endif
}

#if 0
Bunch_train::Bunch_train(Bunches const& bunches, double spacing) 
    : bunches(bunches)
    , spacings(std::vector<double >(bunches.size() - 1, spacing))
    , has_parent_comm(false)
{
    set_bucket_indices();
    calculates_counts_and_offsets_for_impedance();
}

Bunch_train::Bunch_train(Bunches const& bunches,
        std::vector<double > const& spacings) :
                bunches(bunches),
                spacings(spacings),
                has_parent_comm(false)
{
    if (spacings.size() != bunches.size() - 1) {
        throw std::runtime_error(
                "Bunch_train:: spacings must have length (length(bunches)-1)");
    }
    set_bucket_indices();
    calculates_counts_and_offsets_for_impedance();
}
#endif

Bunch_train::Bunch_train(Reference_particle const& ref,
                         size_t num_bunches,
                         size_t num_particles_per_bunch,
                         double num_real_particles_per_bunch,
                         double spacing,
                         Commxx const& bt_comm,
                         size_t num_spectator_per_bunch,
                         int index)
    : comm(std::make_shared<Commxx>(bt_comm))
    , bunches()
    , spacings()
    , index(index)
    , num_bunches(num_bunches)
    , num_buckets(num_bunches)
    , bunch_idx_map(num_bunches, -1)
{
    for (auto i = 0; i < num_bunches; ++i)
        spacings.emplace_back(spacing);

    // empty train
    if (comm->is_null()) return;
    if (num_bunches == 0) return;

    // construct bunches
    int rank = comm->rank();
    int size = comm->size();

    int bunches_per_rank = 0;
    int color = 0;

    if (size < num_bunches) {
        if (num_bunches % size != 0) {
            throw std::runtime_error("Bunch_train::Bunch_train() number of "
                                     "bunches must be divisible "
                                     "by the number of ranks.");
        }

        bunches_per_rank = num_bunches / size;
        color = rank;
    } else {
        if (size % num_bunches != 0) {
            throw std::runtime_error(
                "Bunch_train::Bunch_train() number of ranks must be divisible "
                "by the number of bunches.");
        }

        bunches_per_rank = 1;
        color = rank / (size / num_bunches); // rank / ranks_per_bunch
    }

    for (int b = 0; b < bunches_per_rank; ++b) {
        // bunch_index and array_index
        int b_idx = b + color * bunches_per_rank;
        int a_idx = bunches.size();

        // construct the idx map
        bunch_idx_map[b_idx] = a_idx;

        // construct and push bunch object
        bunches.emplace_back(
            ref,
            num_particles_per_bunch,
            num_real_particles_per_bunch,
            comm->split(color),
            num_spectator_per_bunch,
            b_idx, // bunch index in the train
            b_idx, // bucket index set to the same of bunch index
            a_idx, // array index in the bunches array
            index  // train index
        );
    }
}

std::vector<double>&
Bunch_train::get_spacings()
{
    return spacings;
}

std::vector<int>&
Bunch_train::get_proc_counts_for_impedance()
{
    return proc_counts_imped;
}

std::vector<int>&
Bunch_train::get_proc_offsets_for_impedance()
{
    return proc_offsets_imped;
}

void
Bunch_train::update_bunch_total_num()
{
#if 0
    const size_t nb = get_size();
    if (nb == 0) return;

    std::vector<int> nums(nb);

    for (int i=0; i<nb; ++i)
    {
        nums[i] = bunches[i].get_local_num();
    }

    MPI_Allreduce(MPI_IN_PLACE, &nums[0], nb, MPI_INT, MPI_SUM, 
            get_parent_comm());

    for (int i=0; i<nb; ++i)
    {
        bunches[i].set_total_num(nums[i]);
    }
#endif
}

#if 0
template<class Archive>
    void
    Bunch_train::serialize(Archive & ar, const unsigned int version)
    {
        ar & CEREAL_NVP(bunches);
        ar & CEREAL_NVP(spacings);
        ar & CEREAL_NVP(has_parent_comm);
        ar & CEREAL_NVP(parent_comm_sptr);
        ar & CEREAL_NVP(proc_counts_imped);
        ar & CEREAL_NVP(proc_offsets_imped);
    }

template
void
Bunch_train::serialize<cereal::BinaryOutputArchive >(
        cereal::BinaryOutputArchive & ar, const unsigned int version);

template
void
Bunch_train::serialize<cereal::XMLOutputArchive >(
        cereal::XMLOutputArchive & ar, const unsigned int version);

template
void
Bunch_train::serialize<cereal::BinaryInputArchive >(
        cereal::BinaryInputArchive & ar, const unsigned int version);

template
void
Bunch_train::serialize<cereal::XMLInputArchive >(
        cereal::XMLInputArchive & ar, const unsigned int version);
#endif
