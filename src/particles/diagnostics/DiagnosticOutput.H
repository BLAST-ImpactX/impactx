/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DIAGNOSTIC_OUTPUT_H
#define IMPACTX_DIAGNOSTIC_OUTPUT_H

#include "particles/ImpactXParticleContainer.H"
#include "particles/CovarianceMatrix.H"

#include <string>


namespace impactx::diagnostics
{
    /** ASCII output diagnostics associated with the beam particles.
     *
     * This prints reduced beam characteristics derived from particles in ASCII.
     *
     * @param pc container of the particles use for diagnostics
     * @param file_name the file name to write to
     * @param step the global step
     * @param append open a new file with a fresh header (false) or append data to an existing file (true)
     */
    void DiagnosticOutput (
        ImpactXParticleContainer const & pc,
        std::string file_name,
        int step = 0,
        bool append = false
    );

    /** ASCII output diagnostics associated with the covariance matrix.
     *
     * This prints reduced beam characteristics derived from the envelope description in ASCII.
     *
     * @param cm covariance matrix
     * @param ref_part reference particle
     * @param file_name the file name to write to
     * @param step the global step
     * @param append open a new file with a fresh header (false) or append data to an existing file (true)
     */
    void DiagnosticOutput (
        Map6x6 const & cm,
        RefPart const & ref_part,
        std::string file_name,
        int step = 0,
        bool append = false
    );

    /** ASCII output diagnostics associated with the reference particle.
     *
     * This prints reference particle properties in ASCII.
     *
     * @param ref_part reference particle
     * @param file_name the file name to write to
     * @param step the global step
     * @param append open a new file with a fresh header (false) or append data to an existing file (true)
     */
    void DiagnosticOutput (
        RefPart const & ref_part,
        std::string file_name,
        int step = 0,
        bool append = false
    );

} // namespace impactx::diagnostics

#endif // IMPACTX_DIAGNOSTIC_OUTPUT_H
