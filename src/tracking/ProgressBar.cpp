/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "tracking/ProgressBar.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#if defined(_WIN32)
#   include <io.h>
#   ifndef WIN32_LEAN_AND_MEAN
#       define WIN32_LEAN_AND_MEAN
#   endif
#   ifndef NOMINMAX
#       define NOMINMAX // windows.h: no min/max macros that break std::min/std::max
#   endif
#   include <windows.h>
#else
#   include <sys/ioctl.h>
#   include <unistd.h>
#endif

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <ostream>
#include <streambuf>
#include <string>


namespace impactx
{
namespace
{
    /** Return true if stdout is connected to an interactive terminal. */
    bool
    stdout_is_tty ()
    {
#if defined(_WIN32)
        return _isatty(_fileno(stdout)) != 0;
#else
        return isatty(fileno(stdout)) != 0;
#endif
    }

    /** Return true if the process locale indicates UTF-8 output. */
    bool
    locale_is_utf8 ()
    {
        for (char const * name : {"LC_ALL", "LC_CTYPE", "LANG"})
        {
            char const * value = std::getenv(name);
            if (value != nullptr && value[0] != '\0')
            {
                std::string v(value);
                std::transform(v.begin(), v.end(), v.begin(),
                               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                return v.find("utf-8") != std::string::npos ||
                       v.find("utf8") != std::string::npos;
            }
        }
        // conservative: assume a non-UTF-8 console when the locale is unknown
        return false;
    }

    /** Shorten a label to at most @p width characters (ASCII-safe ".." suffix). */
    std::string
    truncate_label (std::string label, std::size_t width)
    {
        if (width >= 2 && label.size() > width)
        {
            label = label.substr(0, width - 2) + "..";
        }
        return label;
    }

    /** Current terminal width in columns, or 0 if unknown (not a terminal). */
    int
    terminal_columns ()
    {
#if defined(_WIN32)
        CONSOLE_SCREEN_BUFFER_INFO info;
        if (GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &info) != 0)
        {
            return static_cast<int>(info.srWindow.Right - info.srWindow.Left + 1);
        }
        return 0;
#else
        struct winsize w;
        if (ioctl(fileno(stdout), TIOCGWINSZ, &w) == 0)
        {
            return static_cast<int>(w.ws_col);
        }
        return 0;
#endif
    }

    /** Number of terminal columns of a rendered bar line.
     *
     * Counts UTF-8 code points; all glyphs the bar uses occupy one column each.
     */
    std::size_t
    utf8_columns (std::string const & s)
    {
        std::size_t n = 0;
        for (unsigned char const c : s)
        {
            if ((c & 0xC0) != 0x80) { ++n; } // count code-point start bytes
        }
        return n;
    }

    /** The longest prefix of @p s that occupies at most @p ncols terminal columns. */
    std::string
    utf8_prefix (std::string const & s, std::size_t ncols)
    {
        std::size_t n = 0;
        for (std::size_t i = 0; i < s.size(); ++i)
        {
            if ((static_cast<unsigned char>(s[i]) & 0xC0) != 0x80)
            {
                if (n == ncols) { return s.substr(0, i); }
                ++n;
            }
        }
        return s;
    }

    /** A filtering stream buffer that keeps a status line pinned to the bottom.
     *
     * Installed on ``amrex::OutStream()`` while a live bar is up. Foreign output
     * (anything written through ``amrex::Print`` / the wrapped stream) transparently
     * erases the bar, is forwarded to the real terminal, and the bar is redrawn on
     * the new bottom line -- so solver residuals, warnings, etc. scroll above a
     * pinned bar instead of shredding it. The bar itself is drawn via ``set_bar()``,
     * which writes straight to the wrapped buffer and never re-enters this filter.
     */
    class BottomBar final : public std::streambuf
    {
    public:
        explicit BottomBar (std::streambuf * dest)
            : m_dest(dest)
        {}

        /** Redraw the pinned bar with new text (erases the previous one). */
        void set_bar (std::string const & bar)
        {
            erase();
            m_bar = bar;
            m_bar_cols = utf8_columns(bar);
            raw(m_bar);
            m_shown = !m_bar.empty();
            m_dest->pubsync();
        }

        /** Draw a final bar, end the line, and stop pinning. */
        void finish_bar (std::string const & bar)
        {
            erase();
            raw(bar);
            raw("\n");
            m_bar.clear();
            m_bar_cols = 0;
            m_shown = false;
            m_dest->pubsync();
        }

        /** Erase the pinned bar (if any) and stop pinning, leaving the cursor at col 0. */
        void remove_bar ()
        {
            erase();
            m_bar.clear();
            m_bar_cols = 0;
            m_shown = false;
            m_dest->pubsync();
        }

    protected:
        int_type overflow (int_type ch) override
        {
            if (traits_type::eq_int_type(ch, traits_type::eof())) { return ch; }
            char const c = traits_type::to_char_type(ch);
            before_foreign();
            int_type const r = m_dest->sputc(c);
            after_foreign(c == '\n');
            return r;
        }

        std::streamsize xsputn (char const * s, std::streamsize n) override
        {
            if (n <= 0) { return 0; }
            before_foreign();
            std::streamsize const r = m_dest->sputn(s, n);
            after_foreign(s[n - 1] == '\n');
            return r;
        }

        int sync () override { return m_dest->pubsync(); }

    private:
        void raw (std::string const & s)
        {
            m_dest->sputn(s.data(), static_cast<std::streamsize>(s.size()));
        }

        void erase ()
        {
            if (m_shown)
            {
                raw("\r\x1b[2K"); // carriage return + ANSI "erase entire line"

                // If the terminal was shrunk below the drawn bar's width, the bar has
                // re-wrapped onto multiple rows; the erase above only cleared the
                // bottom one. Clear the remaining wrapped rows above it as well.
                int const cols = terminal_columns();
                if (cols > 0 && m_bar_cols > static_cast<std::size_t>(cols))
                {
                    std::size_t const extra =
                        (m_bar_cols - 1) / static_cast<std::size_t>(cols);
                    for (std::size_t i = 0; i < extra; ++i)
                    {
                        raw("\x1b[1A\x1b[2K"); // cursor up + erase that row
                    }
                    // return to the bar's bottom row: drawing at the top row instead
                    // would migrate the bar one row up on every shrink, until it
                    // reaches (and glitches at) the top of the terminal
                    raw("\x1b[" + std::to_string(extra) + "B"); // cursor down
                }
                m_shown = false;
            }
        }

        /** Called before forwarding foreign output: clear the bar off the line. */
        void before_foreign ()
        {
            erase();
        }

        /** Called after forwarding foreign output: repin the bar once a line ended. */
        void after_foreign (bool ended_line)
        {
            if (ended_line && !m_bar.empty())
            {
                raw(m_bar);
                m_shown = true;
            }
        }

        std::streambuf * m_dest;      //!< the real terminal buffer we forward to
        std::string m_bar;            //!< current pinned bar text (no trailing newline)
        std::size_t m_bar_cols = 0;   //!< terminal columns the pinned bar occupies
        bool m_shown = false;         //!< is the bar currently on the terminal's last line?
    };
} // namespace


    std::string
    format_time (double seconds)
    {
        if (seconds < 0.0) { seconds = 0.0; }
        auto const total = static_cast<long>(seconds + 0.5);
        long const h = total / 3600;
        long const m = (total % 3600) / 60;
        long const s = total % 60;

        char buffer[32];
        if (h > 0)
        {
            std::snprintf(buffer, sizeof(buffer), "%ld:%02ld:%02ld", h, m, s);
        }
        else
        {
            std::snprintf(buffer, sizeof(buffer), "%ld:%02ld", m, s);
        }
        return std::string(buffer);
    }


    std::string
    render_bar (
        double frac,
        double s_cur,
        double s_tot,
        int bar_width,
        double elapsed,
        std::string const & label,
        bool ascii
    )
    {
        frac = std::min(std::max(frac, 0.0), 1.0);
        if (bar_width < 1) { bar_width = 1; }

        // filled length, measured in eighths of a character cell
        int const eighths = static_cast<int>(std::lround(frac * 8.0 * bar_width));
        int const full = eighths / 8;
        int const rem = eighths % 8;

        std::string bar;
        std::string left_cap;
        std::string right_cap;
        if (ascii)
        {
            left_cap = "[";
            right_cap = "]";
            int cells = 0;
            bar.append(static_cast<std::size_t>(full), '#');
            cells += full;
            if (full < bar_width && rem > 0)
            {
                bar.push_back('+');
                cells += 1;
            }
            bar.append(static_cast<std::size_t>(bar_width - cells), '-');
        }
        else
        {
            left_cap = "▕";  // right one-eighth block (thin left cap)
            right_cap = "▏"; // left one-eighth block (thin right cap)
            static constexpr std::array<char const *, 8> partial = {
                " ", "▏", "▎", "▍",
                "▌", "▋", "▊", "▉"
            };
            int cells = 0;
            for (int i = 0; i < full; ++i) { bar += "█"; } // full block
            cells += full;
            if (full < bar_width && rem > 0)
            {
                bar += partial[static_cast<std::size_t>(rem)];
                cells += 1;
            }
            for (int i = cells; i < bar_width; ++i) { bar += " "; }
        }

        int const pct = static_cast<int>(std::lround(frac * 100.0));

        // progress read-out: path length s (if known), percent, and ETA / total time
        std::string suffix;
        if (s_tot > 0.0)
        {
            // pad the running s to the width of the total, so the field width (and
            // with it the whole line layout) stays constant while the digits grow
            char tot[32];
            int const tot_width = std::snprintf(tot, sizeof(tot), "%.2f", s_tot);
            char sb[96];
            std::snprintf(sb, sizeof(sb), "s=%*.2f/%s m  ", tot_width, s_cur, tot);
            suffix += sb;
        }
        {
            char pb[16];
            std::snprintf(pb, sizeof(pb), "%3d%%  ", pct);
            suffix += pb;
        }
        if (frac < 1.0)
        {
            double const eta = (frac > 0.0) ? elapsed * (1.0 - frac) / frac : 0.0;
            suffix += "ETA " + format_time(eta) + "  ";
        }
        else
        {
            suffix += "time " + format_time(elapsed) + "  ";
        }
        suffix += label;

        return "Tracking " + left_cap + bar + right_cap + "  " + suffix;
    }


    ProgressBar::ProgressBar (int total_steps, double total_length, int verbose)
        : m_verbose(verbose), m_total(total_steps), m_total_s(total_length)
    {
        m_io = amrex::ParallelDescriptor::IOProcessor();

        std::string progress = "auto";
        amrex::ParmParse("impactx").queryAdd("progress", progress);

        if (verbose <= 0 || m_total <= 0)
        {
            m_mode = Mode::Silent;
        }
        else if (progress == "on")
        {
            m_mode = Mode::Live; // explicit user override, even into a pipe
        }
        else if (progress == "off")
        {
            m_mode = Mode::Banner;
        }
        else // "auto"
        {
            if (verbose >= 2)
            {
                m_mode = Mode::Banner;
            }
            else
            {
                m_mode = stdout_is_tty() ? Mode::Live : Mode::Banner;
            }
        }

        bool ascii = !locale_is_utf8();
        amrex::ParmParse("impactx").queryAdd("progress_ascii", ascii);
        m_ascii = ascii;

        if (m_mode == Mode::Live && m_io)
        {
            m_start_time = amrex::second();

            // pin the bar to the bottom line: install a filter on amrex::Print()'s
            // stream so any foreign output (solver residuals, warnings, ...) scrolls
            // above the bar instead of shredding it (see BottomBar).
            std::ostream & out = amrex::OutStream();
            m_saved_streambuf = out.rdbuf();
            m_bar_buf = std::make_unique<BottomBar>(m_saved_streambuf);
            out.rdbuf(m_bar_buf.get());
        }
    }


    ProgressBar::~ProgressBar ()
    {
        uninstall();
    }


    void
    ProgressBar::uninstall ()
    {
        if (!m_bar_buf) { return; }

        // erase a still-visible bar (e.g. on an exception before finish())
        static_cast<BottomBar *>(m_bar_buf.get())->remove_bar();

        // restore amrex::OutStream()'s buffer BEFORE the filter buffer is freed
        if (m_saved_streambuf != nullptr)
        {
            amrex::OutStream().rdbuf(m_saved_streambuf);
            m_saved_streambuf = nullptr;
        }
        m_bar_buf.reset();
    }


    void
    ProgressBar::show (int step, int slice_step, double s, std::string const & label)
    {
        if (m_mode == Mode::Silent) { return; }

        // fraction completed: by path length s when known, else by slice-step count
        double frac = (m_total_s > 0.0)
            ? s / m_total_s
            : (m_total > 0 ? static_cast<double>(step) / m_total : 1.0);
        frac = std::min(std::max(frac, 0.0), 1.0);

        if (m_mode == Mode::Banner)
        {
            int const denom = (m_total > 0) ? m_total : 1;
            int const pct = (std::min(step, denom) * 100) / denom;
            amrex::Print() << "\n++++ Starting step=" << step
                           << " of " << m_total
                           << " (" << pct << "%) " << label;
            if (m_total_s > 0.0)
            {
                char sb[80];
                std::snprintf(sb, sizeof(sb), "  s=%.2f/%.2f m", s, m_total_s);
                amrex::Print() << sb;
            }
            if (m_verbose >= 2)
            {
                amrex::Print() << " slice_step=" << slice_step;
            }
            return;
        }

        // Mode::Live (I/O rank only): throttle redraws to when the rendered state
        // changes (clock-free), then repaint the bottom-pinned bar
        if (!m_io || !m_bar_buf) { return; }

        // re-query the terminal width every redraw: it tracks window resizes
        int const cols = terminal_columns();

        int const bucket = static_cast<int>(frac * 8.0 * m_bar_width);
        if (bucket == m_last_bucket && cols == m_last_cols && frac < 1.0) { return; }
        m_last_bucket = bucket;
        m_last_cols = cols;

        double const elapsed = amrex::second() - m_start_time;
        static_cast<BottomBar *>(m_bar_buf.get())->set_bar(
            render_fit(frac, s, elapsed, label, cols));
    }


    void
    ProgressBar::finish ()
    {
        if (m_mode != Mode::Live || !m_io || !m_bar_buf) { return; }

        double const elapsed = amrex::second() - m_start_time;
        static_cast<BottomBar *>(m_bar_buf.get())->finish_bar(
            render_fit(1.0, m_total_s, elapsed, "done", terminal_columns()));

        uninstall();
    }


    std::string
    ProgressBar::render_fit (
        double frac, double s_cur, double elapsed,
        std::string const & label, int max_cols
    ) const
    {
        // unknown width (e.g. a forced bar in a pipe): assume a classic 80-column line
        if (max_cols <= 0) { max_cols = 80; }
        // never draw into the last column: writing there triggers auto-wrap, and a
        // wrapped line cannot be erased with a carriage return anymore
        int const usable = std::max(max_cols - 1, 8);

        // The layout is computed from content-independent field reserves, NOT from the
        // rendered text: element labels, the running s digits and the ETA all change
        // from step to step, and a layout derived from them would make the bar body
        // change its width from frame to frame (a "fluctuating" bar).
        int const overhead = 9 + 2 + 2; // "Tracking " + bar caps + "  " separator
        int const pct_cols = 6;         // "100%  "
        int const time_cols = 12;       // "time H:MM:SS" / "ETA ..." reserve
        int s_cols = 0;
        if (m_total_s > 0.0)
        {
            char tot[32];
            // "s=" + total + "/" + total + " m" + "  "
            s_cols = 2 * std::snprintf(tot, sizeof(tot), "%.2f", m_total_s) + 7;
        }
        int const fixed = overhead + pct_cols + time_cols;

        // prefer the full-width bar; the element label gets the leftover columns
        double s_tot = m_total_s;
        int const label_cols = usable - fixed - s_cols - m_bar_width;
        std::string lab;
        if (label_cols >= 4)
        {
            lab = truncate_label(label, static_cast<std::size_t>(std::min(label_cols, 24)));
        }

        // narrow: drop the label, then the s read-out, before shrinking the bar
        int bar_width = std::min(m_bar_width, usable - fixed - s_cols);
        if (bar_width < 8)
        {
            s_tot = 0.0;
            bar_width = std::min(m_bar_width, usable - fixed);
        }
        bar_width = std::max(bar_width, 8);

        std::string line = render_bar(frac, s_cur, s_tot, bar_width, elapsed, lab, m_ascii);

        // last resort (very narrow terminal): hard-truncate at the column limit
        if (utf8_columns(line) > static_cast<std::size_t>(usable))
        {
            line = utf8_prefix(line, static_cast<std::size_t>(usable));
        }
        return line;
    }

} // namespace impactx
