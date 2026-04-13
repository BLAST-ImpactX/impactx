#include "cereal_files.h"
#include "commxx.h"
#include "digits.h"

#include <filesystem>

namespace fs = std::filesystem;

const std::string serialization_directory("serialization");

class Parallel_helper {
private:
  bool parallel;

public:
  explicit Parallel_helper(bool parallel) : parallel(parallel) {}

  void
  barrier()
  {
    if (parallel) { MPI_Barrier(Commxx()); }
  }

  bool
  operate_locally()
  {
    if (parallel) {
      return Commxx::world_rank() == 0;
    } else {
      return true;
    }
  }
};

void
copy_file_overwrite_if_exists(std::string const& source,
                              std::string const& dest)
{
  fs::path const src{source};
  fs::path const dst{dest};
  fs::rename(src, dst);
}

std::string
get_serialization_directory()
{
  return serialization_directory;
}

void
remove_directory(std::string const& name, bool parallel)
{
  Parallel_helper parallel_helper(parallel);
  parallel_helper.barrier();
  if (parallel_helper.operate_locally()) { fs::remove_all(name); }
  parallel_helper.barrier();
}

void
remove_serialization_directory(bool parallel)
{
  Parallel_helper parallel_helper(parallel);
  parallel_helper.barrier();
  if (parallel_helper.operate_locally()) {
    fs::path const path_to_delete{get_serialization_directory()};
    fs::remove_all(path_to_delete);
  }
  parallel_helper.barrier();
}

void
ensure_serialization_directory_exists(bool parallel)
{
  Parallel_helper parallel_helper(parallel);
  parallel_helper.barrier();
  if (parallel_helper.operate_locally()) {
    fs::path const path_to_create{get_serialization_directory()};
    fs::create_directories(path_to_create);
  }
  parallel_helper.barrier();
}

void
rename_serialization_directory(std::string const& new_name, bool parallel)
{
  Parallel_helper parallel_helper(parallel);
  parallel_helper.barrier();
  if (parallel_helper.operate_locally()) {
    fs::path const serialization_path{get_serialization_directory()};
    fs::path const new_path{new_name};
    fs::rename(serialization_path, new_path);
  }
  parallel_helper.barrier();
}

void
symlink_serialization_directory(std::string const& existing_dir, bool parallel)
{
  Parallel_helper parallel_helper(parallel);
  parallel_helper.barrier();
  if (parallel_helper.operate_locally()) {
    fs::path const serialization_path{get_serialization_directory()};
    fs::path const existing_path{existing_dir};
    fs::create_directory_symlink(existing_path, serialization_path);
  }
  parallel_helper.barrier();
}

void
unlink_serialization_directory(bool parallel)
{
  Parallel_helper parallel_helper(parallel);
  parallel_helper.barrier();
  if (parallel_helper.operate_locally()) {
    fs::path const serialization_path{get_serialization_directory()};
    fs::remove(serialization_path);
  }
  parallel_helper.barrier();
}

std::string
get_combined_path(std::string const& directory,
                  std::string const& base_name,
                  bool parallel)
{
  std::string full_name(base_name);
  if (parallel) {
    Commxx commxx;
    std::stringstream sstream;
    sstream << std::setw(digits(commxx.get_size()));
    sstream << std::setfill('0');
    sstream << commxx.get_rank();
    sstream << "_";
    sstream << base_name;
    full_name = sstream.str();
  }
  return directory + "/" + full_name;
}

std::string
get_serialization_path(std::string const& base_name, bool parallel)
{
  return std::string();
#if 0
#if (defined BOOST_FILESYSTEM_VERSION) && (BOOST_FILESYSTEM_VERSION > 2)
    return get_combined_path(serialization_directory,
                             path(base_name).filename().string(), parallel);
#else
    return get_combined_path(serialization_directory,
                             path(base_name).filename(), parallel);
#endif
#endif
}

void
copy_to_serialization_directory(std::string const& file_name)
{
  copy_file_overwrite_if_exists(file_name, get_serialization_path(file_name));
}

void
copy_from_serialization_directory(std::string const& file_name)
{
  copy_file_overwrite_if_exists(get_serialization_path(file_name), file_name);
}
