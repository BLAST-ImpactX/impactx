/* Copyright 2026 The ImpactX Community
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "FilePrefix.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <filesystem>
#include <string>


namespace impactx::diagnostics
{
    namespace
    {
        std::filesystem::path
        normalize_file_prefix (std::string const & file_prefix)
        {
            std::filesystem::path path = std::filesystem::path(file_prefix).lexically_normal();
            while (!path.empty() && path.filename().empty() && path.has_parent_path() &&
                   path != path.root_path())
            {
                path = path.parent_path();
            }
            return path;
        }

        bool
        is_path_prefix (
            std::filesystem::path const & prefix,
            std::filesystem::path const & path
        )
        {
            auto path_it = path.begin();
            for (auto prefix_it = prefix.begin(); prefix_it != prefix.end(); ++prefix_it, ++path_it)
            {
                if (path_it == path.end() || *path_it != *prefix_it) {
                    return false;
                }
            }
            return true;
        }

        void
        create_file_prefix_directory (std::filesystem::path const & file_prefix)
        {
            if (file_prefix.empty() || file_prefix == std::filesystem::path(".")) {
                return;
            }

            constexpr int permission_flag_rwxrxrx = 0755;
            std::string const file_prefix_string = file_prefix.generic_string();
            if (!amrex::UtilCreateDirectory(file_prefix_string, permission_flag_rwxrxrx)) {
                amrex::CreateDirectoryFailed(file_prefix_string);
            }
        }
    }

    std::string
    FilePrefix ()
    {
        std::string file_prefix = "diags";
        amrex::ParmParse("diag").queryAdd("file_prefix", file_prefix);
        return normalize_file_prefix(file_prefix).generic_string();
    }

    std::string
    FilePrefixPath (std::string const & name)
    {
        std::filesystem::path const file_prefix = normalize_file_prefix(FilePrefix());
        create_file_prefix_directory(file_prefix);

        if (name.empty()) {
            return file_prefix.generic_string();
        }

        bool const use_current_directory = file_prefix.empty() ||
            file_prefix == std::filesystem::path(".");
        std::filesystem::path const path = use_current_directory
            ? std::filesystem::path(name)
            : file_prefix / std::filesystem::path(name);

        return path.lexically_normal().generic_string();
    }

    bool
    FilePrefixCanClean ()
    {
        std::filesystem::path const file_prefix = normalize_file_prefix(FilePrefix());
        if (file_prefix.empty() || file_prefix == std::filesystem::path(".")) {
            return false;
        }

        std::filesystem::path const absolute_file_prefix =
            std::filesystem::absolute(file_prefix).lexically_normal();
        std::filesystem::path const current_path =
            std::filesystem::current_path().lexically_normal();

        if (absolute_file_prefix == current_path ||
            absolute_file_prefix == absolute_file_prefix.root_path())
        {
            return false;
        }
        return !is_path_prefix(absolute_file_prefix, current_path);
    }
} // namespace impactx::diagnostics
