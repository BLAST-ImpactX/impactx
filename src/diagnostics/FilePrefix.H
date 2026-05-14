/* Copyright 2026 The ImpactX Community
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DIAGNOSTICS_FILE_PREFIX_H
#define IMPACTX_DIAGNOSTICS_FILE_PREFIX_H

#include <string>


namespace impactx::diagnostics
{
    /** Root directory for diagnostic output.
     *
     * This queries ParmParse every time, so runtime changes to diag.file_prefix
     * are honored.
     * Empty string and "." mean the current working directory.
     */
    std::string FilePrefix ();

    /** Path to a diagnostic file or directory under diag.file_prefix.
     */
    std::string FilePrefixPath (std::string const & name);

    /** Whether an existing diag.file_prefix directory may be renamed aside at initialization.
     *
     * Current-directory, root-directory, and ancestor-directory paths are not cleaned.
     */
    bool FilePrefixCanClean ();
} // namespace impactx::diagnostics

#endif // IMPACTX_DIAGNOSTICS_FILE_PREFIX_H
