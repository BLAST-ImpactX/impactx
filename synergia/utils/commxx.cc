#include <stdexcept>

#include "commxx.h"

namespace {
  struct comm_free {
    void
    operator()(MPI_Comm* comm)
    {
      if (comm == nullptr) return;
      if (*comm == MPI_COMM_NULL) return;

      MPI_Comm_free(comm);
      delete comm;
    }
  };

}

const Commxx Commxx::World(comm_type::world);
const Commxx Commxx::Null(comm_type::null);

inline bool
valid_mpi_communicator(std::shared_ptr<const Commxx> const& pcomm)
{
  if (!pcomm) return false;   // there is no controlled Commxx
  return ! pcomm->is_null();  // a null MPI communicator is not valid
}

Commxx
Commxx::create_child(std::shared_ptr<const Commxx>&& parent, int color, int key)
{
  return Commxx(std::move(parent), color, key);
}


Commxx::Commxx(comm_type type)
  : mpi_comm(type == comm_type::null ? nullptr : new MPI_Comm(MPI_COMM_WORLD))
  , parent_comm()
  , type(type)
  , color(0)
  , key(0)
{
  if (type == comm_type::regular)
    throw std::runtime_error("only null or mpi_comm_world can be created");
}

void
Commxx::split_parent_and_set_mpi_comm()
{
  if (!valid_mpi_communicator(parent_comm))
    throw std::runtime_error("invalid parent communicator while in split_parent_and_set_mpi_comm");

  MPI_Comm newcomm;
  int status = MPI_Comm_split(*(parent_comm->mpi_comm), color, key, &newcomm);
  if (status != MPI_SUCCESS)
    throw std::runtime_error("MPI_Comm_split failed in split_parent_and_set_mpi_comm");

  // Do not split_parent_and_set_mpi_comm the null communicator.
  // Note that MPI_Commm_split sets 'newcomm' to MPI_COMM_NULL if 'color' is
  // set to MPI_UNDEFINED.
  if (newcomm != MPI_COMM_NULL) {
    // Replace mpi_comm with a shared_ptr that will free and MPI_Comm and then
    // delete the pointer.
    mpi_comm.reset(new MPI_Comm(newcomm), comm_free());
  }
}

Commxx::Commxx(std::shared_ptr<const Commxx>&& parent, int color, int key)
  : mpi_comm()
  , parent_comm(std::move(parent))
  , type(comm_type::regular)
  , color(color)
  , key(key)
{
  split_parent_and_set_mpi_comm();
}

int
Commxx::get_rank() const
{
  if (type == comm_type::null)
    throw std::runtime_error("Cannot get_rank() for a null commxx");

  int rank;
  int status = MPI_Comm_rank(*mpi_comm, &rank);
  if (status != MPI_SUCCESS) {
    throw std::runtime_error("MPI error in MPI_Comm_rank");
  }
  return rank;
}

int
Commxx::get_size() const
{
  if (type == comm_type::null)
    throw std::runtime_error("Cannot get_size() for a null commxx");

  int size;
  int status = MPI_Comm_size(*mpi_comm, &size);
  if (status != MPI_SUCCESS) {
    throw std::runtime_error("MPI error in MPI_Comm_size");
  }
  return size;
}

Commxx
Commxx::dup() const
{
  if (is_null()) throw std::runtime_error("dup from a null comm");
  try {
    return create_child(shared_from_this(), 0, rank());
  }
  catch (std::bad_weak_ptr const&) {
    throw std::runtime_error("Illegal attempt to dup a Commxx not owned by a shared_ptr");
  }
}

Commxx
Commxx::split(int color) const
{
  if (is_null()) throw std::runtime_error("split from a null comm");
  try {
    return create_child(shared_from_this(), color, rank());
  }
  catch (std::bad_weak_ptr const&) {
    throw std::runtime_error("Illegal attempt to split a Commxx not owned by a shared_ptr");
  }
}

Commxx
Commxx::split(int color, int key) const
{
  if (is_null()) throw std::runtime_error("split from a null comm");
  try {
    return create_child(shared_from_this(), color, key);
  }
  catch (std::bad_weak_ptr const&) {
    throw std::runtime_error("Illegal attempt to split a Commxx not owned by a shared_ptr");
  }
}

Commxx
Commxx::divide(int subgroup_size) const
{
  if (is_null()) throw std::runtime_error("divide from a null comm");

  if (size() < subgroup_size) return dup();
  if (size() % subgroup_size != 0) {
    throw std::runtime_error(
      "Commxx::divide(): size must be divisible by subgroup_size");
  }

  int color = rank() / subgroup_size;
  return split(color);
}

Commxx
Commxx::group(std::vector<int> const& ranks) const
{
  if (is_null()) throw std::runtime_error("group from a null comm");

  int r = rank();
  bool in_range = std::find(ranks.begin(), ranks.end(), r) != ranks.end();
  int color = (in_range ? 0 : MPI_UNDEFINED);
  try {
    return create_child(shared_from_this(), color, r);
  }
  catch (std::bad_weak_ptr const&) {
    throw std::runtime_error("Illegal attempt to group a Commxx not owned by a shared_ptr");
  }
}

bool
operator==(Commxx const& comm1, Commxx const& comm2)
{
  // both null
  if (comm1.is_null() && comm2.is_null()) return true;

  // one is null, the other is not
  if (comm1.is_null() || comm2.is_null()) return false;

  // both not null
  int result;
  int status = MPI_Comm_compare((MPI_Comm)comm1, (MPI_Comm)comm2, &result);
  if (status != MPI_SUCCESS) {
    throw std::runtime_error("MPI_Comm_compare failed while comparing Commxx objects");
  }
  return result == MPI_IDENT;
}
