/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_TRACKING_PROGRESSBAR_H
#define IMPACTX_TRACKING_PROGRESSBAR_H

#include "elements/All.H"

#include <cstddef>
#include <iosfwd>
#include <list>
#include <memory>
#include <string>
#include <variant>


namespace impactx
{
    /** Total number of forward slice-steps of a tracking run.
     *
     * This equals ``num_periods * sum_over_elements(nslice)`` and is used as the
     * denominator for the progress bar. It mirrors the global step counter
     * ``TrackingState::m_step``, which increments once per slice on forward tracking.
     *
     * @param lattice the beamline elements to traverse
     * @param num_periods number of periods (turns/cycles) through the lattice
     * @return the total number of forward slice-steps (>= 0)
     */
    inline int
    count_total_steps (
        std::list<elements::KnownElements> const & lattice,
        int num_periods
    )
    {
        int per_period = 0;
        for (auto const & element_variant : lattice)
        {
            int nslice = 1;
            std::visit([&nslice](auto && element) { nslice = element.nslice(); },
                       element_variant);
            per_period += nslice;
        }
        return per_period * num_periods;
    }

    /** Total reference-trajectory path length ``s`` of a tracking run, in meters.
     *
     * This equals ``num_periods * sum_over_elements(ds)`` and is used as the
     * denominator for the progress bar's ``s`` read-out (and, when positive, drives
     * the bar fill so it advances with physical distance rather than slice count).
     *
     * @param lattice the beamline elements to traverse
     * @param num_periods number of periods (turns/cycles) through the lattice
     * @return the total path length in meters (>= 0)
     */
    inline double
    count_total_length (
        std::list<elements::KnownElements> const & lattice,
        int num_periods
    )
    {
        double per_period = 0.0;
        for (auto const & element_variant : lattice)
        {
            double ds = 0.0;
            std::visit([&ds](auto && element) { ds = static_cast<double>(element.ds()); },
                       element_variant);
            per_period += ds;
        }
        return per_period * static_cast<double>(num_periods);
    }

    /** Format a wall-clock duration in seconds as ``M:SS`` or ``H:MM:SS``. */
    std::string
    format_time (double seconds);

    /** Render a single progress-bar frame (pure, no I/O, no carriage return).
     *
     * Kept as a free, side-effect-free function so it is deterministically testable
     * without a terminal.
     *
     * @param frac fraction completed in [0, 1] (drives the bar fill)
     * @param s_cur current path length ``s`` in meters (shown if @p s_tot > 0)
     * @param s_tot total path length in meters (values <= 0 hide the ``s`` read-out)
     * @param bar_width width of the bar body in character cells
     * @param elapsed wall-clock seconds since the start of tracking
     * @param label short label for the current element (name or type)
     * @param ascii use ASCII glyphs instead of Unicode block characters
     * @return the rendered line (without a trailing newline)
     */
    std::string
    render_bar (
        double frac,
        double s_cur,
        double s_tot,
        int bar_width,
        double elapsed,
        std::string const & label,
        bool ascii
    );

    /** A terminal progress bar / status line for the tracking loops.
     *
     * The output is MPI-safe: it is only drawn on the I/O rank. The mode is chosen
     * once at construction from the verbosity (``impactx.verbose``), the
     * ``impactx.progress`` input and whether ``stdout`` is an interactive terminal:
     *  - ``Silent``: ``verbose == 0`` (or an empty lattice) -> no output.
     *  - ``Live``:   interactive terminal, default verbosity -> a single line pinned
     *                to the bottom of the terminal, redrawn in place.
     *  - ``Banner``: non-interactive (pipe/file/batch) or ``verbose >= 2`` -> one
     *                ``++++ Starting step=N of T (P%) label`` line per step.
     *
     * In ``Live`` mode the constructor installs a filtering stream buffer on
     * ``amrex::OutStream()`` (which ``amrex::Print`` writes to). Any other output
     * that arrives while the bar is up -- e.g. the MLMG space-charge solver's
     * per-slice residuals, or warnings -- transparently erases the bar, scrolls
     * above it, and the bar is redrawn on the (new) bottom line. The buffer is
     * restored in ``finish()`` / the destructor. This coexistence uses ANSI
     * "erase line" (``\x1b[2K``), which every modern terminal supports.
     */
    class ProgressBar
    {
    public:
        /** Output mode, selected once at construction. */
        enum class Mode
        {
            Silent, //!< no output
            Banner, //!< one "++++ Starting step=N of T ..." line per step (pipe/file)
            Live    //!< a single status line pinned to the bottom of the terminal
        };

        /** @param total_steps result of count_total_steps(); values <= 0 disable output
         *  @param total_length result of count_total_length() (meters); drives the ``s`` read-out
         *  @param verbose the impactx.verbose level */
        ProgressBar (int total_steps, double total_length, int verbose);

        /** Restore ``amrex::OutStream()``'s stream buffer if it is still installed. */
        ~ProgressBar ();

        ProgressBar (ProgressBar const &) = delete;
        ProgressBar & operator= (ProgressBar const &) = delete;

        /** Report progress for one forward slice-step.
         *
         * In ``Live`` mode this redraws the pinned status line (throttled); in
         * ``Banner`` mode it prints one progress line; in ``Silent`` mode nothing.
         *
         * @param step the global step counter (1 .. total)
         * @param slice_step the slice index within the current element (Banner detail)
         * @param s reference-trajectory path length in meters at this step
         * @param label short label of the current element (name or type)
         */
        void show (int step, int slice_step, double s, std::string const & label);

        /** Finalize the bar (``Live`` mode only): draw the 100% frame + newline and
         * restore the stream buffer. */
        void finish ();

        /** @return true if a live carriage-return bar is drawn */
        bool is_live () const { return m_mode == Mode::Live; }

        /** @return true if per-step banner lines are emitted */
        bool wants_banner () const { return m_mode == Mode::Banner; }

    private:
        /** Erase the pinned bar and restore ``amrex::OutStream()``'s stream buffer. */
        void uninstall ();

        /** Render a frame that fits into @p max_cols terminal columns.
         *
         * Never returns a line wider than ``max_cols`` (a wrapped line could not be
         * erased with a carriage return): the bar is shrunk first, then the element
         * label is dropped, then the line is hard-truncated.
         */
        std::string render_fit (
            double frac, double s_cur, double elapsed,
            std::string const & label, int max_cols
        ) const;

        Mode m_mode = Mode::Silent;
        int m_verbose = 1;
        int m_total = 0;         //!< total forward slice-steps (bar denominator fallback)
        double m_total_s = 0.0;  //!< total path length in meters (for the ``s`` read-out)
        int m_bar_width = 24;
        int m_last_bucket = -1;  //!< throttle: redraw only when this changes
        int m_last_cols = -1;    //!< terminal width at the last redraw (resize detection)
        bool m_ascii = false;
        bool m_io = false;       //!< this is the I/O rank (only it draws)
        double m_start_time = 0.0;             //!< amrex::second() at construction (Live)
        std::streambuf * m_saved_streambuf = nullptr; //!< amrex::OutStream() buffer to restore
        std::unique_ptr<std::streambuf> m_bar_buf;    //!< bottom-pinned filtering buffer (Live)
    };

} // namespace impactx

#endif // IMPACTX_TRACKING_PROGRESSBAR_H
