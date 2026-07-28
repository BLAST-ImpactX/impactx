/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_BEAM_MOMENTS_SELECTION_H
#define IMPACTX_BEAM_MOMENTS_SELECTION_H

#include <string>
#include <vector>


namespace impactx::diagnostics
{
    /** Reduction profiles, ordered by increasing coverage; each is a superset of
     *  the previous one: Positions < Sizes < Twiss < Full.
     *
     * A profile is the set of weighted power sums a single-pass reduction has to
     * accumulate; lighter profiles read fewer SoA arrays and reduce fewer values.
     * The runtime resolver picks the smallest profile whose reduced quantities
     * cover a requested set of output names.
     *
     * - Positions: means and sigmas of x, y, t (and charge)
     * - Sizes:     + means and sigmas of px, py, pt
     * - Twiss:     + emittances, dispersion and (dispersion-corrected) Twiss
     * - Full:      + full 6x6 covariance for eigenemittances
     */
    enum class MomentsProfile : int { Positions = 0, Sizes, Twiss, Full };

    /** A resolved beam-moments selection.
     *
     * Describes which compile-time reduction profile to run, whether to also
     * reduce the spin and min/max quantities, whether to compute eigenemittances,
     * and which output keys to emit.
     *
     * @param profile the smallest reduction profile covering the request
     * @param spin    also reduce the spin first and second moments
     * @param minmax  also reduce the per-coordinate minima and maxima
     * @param eigen   also compute the eigenemittances (implies profile == Full)
     * @param keys    output keys to emit; if empty, every key the profile covers
     *                is emitted (used for the "all" and lighter-default selections)
     */
    struct MomentsSelection
    {
        MomentsProfile profile = MomentsProfile::Full;
        bool spin = true;
        bool minmax = true;
        bool eigen = false;
        std::vector<std::string> keys {};
    };

    /** Resolve a list of requested output names into the fastest selection that
     *  can deliver them.
     *
     * Deprecated aliases (e.g. ``x_mean``, ``sig_x``, ``pt_min``) are accepted and
     * map to their canonical requirement; the emitted key keeps the requested
     * spelling. The special token ``"all"`` requests the full set (empty ``keys``).
     * An unknown name throws ``std::runtime_error``.
     *
     * @param names        requested output names (or {"all"})
     * @param eigen_default eigenemittances flag used when ``"all"`` is requested
     */
    MomentsSelection
    resolve_beam_moments_selection (std::vector<std::string> const & names, bool eigen_default);

    /** The full selection (all outputs). Used by backward-compatible callers
     *  (space charge, ASCII diagnostics, the deprecated Python binding).
     */
    MomentsSelection
    all_beam_moments_selection (bool eigen);

    /** The default selection: every output except min/max, and except the spin
     *  moments when spin tracking is off. Eigenemittances follow the eigen flag.
     */
    MomentsSelection
    default_beam_moments_selection (bool spin_on, bool eigen);

} // namespace impactx::diagnostics

#endif // IMPACTX_BEAM_MOMENTS_SELECTION_H
