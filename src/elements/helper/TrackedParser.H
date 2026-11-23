/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_HELPER_TRACKEDPARSER_H
#define IMPACTX_ELEMENTS_HELPER_TRACKEDPARSER_H

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
#include <AMReX_Vector.H>

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>


namespace impactx::elements
{
    /** A user-supplied formula, compiled for host and device on demand.
     *
     * The formula is host data that outlives AMReX sessions, like the arrays in
     * @c amrex::Gpu::TrackedVector. Its compiled form belongs to one AMReX session:
     * it is built on first use and dropped by @see release or, at the latest, by
     * @c amrex::Finalize, so that a later session compiles it again.
     *
     * Names other than the variables are resolved as constants from the ParmParse
     * table when the formula is compiled.
     *
     * @tparam N number of variables
     */
    template<int N>
    class TrackedParser
    {
    public:
        /** A formula
         *
         * @param expression the formula, e.g., ``"2.0*x*y"``
         * @param variables the names of the N variables, in the order of evaluation
         */
        TrackedParser (std::string expression, amrex::Vector<std::string> variables)
          : m_variables(std::move(variables))
        {
            set_expression(std::move(expression));
        }

        /** The formula as given */
        [[nodiscard]] std::string const & expression () const { return m_expression; }

        /** Replace the formula
         *
         * Checks the syntax first, so that a rejected formula leaves the old one in place.
         *
         * @param expression the new formula
         */
        void set_expression (std::string expression)
        {
            // an empty formula parses, but compiles to nothing
            if (expression.find_first_not_of(" \t\n") == std::string::npos) {
                throw std::invalid_argument("TrackedParser: the formula must not be empty");
            }
            try {
                // parsing needs no AMReX session
                [[maybe_unused]] amrex::Parser const syntax_check(expression);
            } catch (std::runtime_error const & e) {
                throw std::invalid_argument(
                    "TrackedParser: cannot parse '" + expression + "': " + e.what());
            }
            m_expression = std::move(expression);
            release();
        }

        /** The compiled formula, callable on host and device
         *
         * Compiles on first use in an AMReX session. The executor does not own the
         * compiled form: it is valid until @see release, @see set_expression or
         * @c amrex::Finalize.
         */
        [[nodiscard]] amrex::ParserExecutor<N> executor () const
        {
            if (!m_parser || !*m_parser) {
                m_parser = std::make_shared<amrex::Parser>(
                    amrex::ParmParse{}.makeParser(m_expression, m_variables));

                // drop the compiled form while AMReX is still up, also if no element
                // releases it, so that a later session compiles it again
                std::weak_ptr<amrex::Parser> weak_parser = m_parser;
                amrex::ExecOnFinalize([weak_parser]() {
                    if (auto parser = weak_parser.lock()) { *parser = amrex::Parser{}; }
                });
            }
            return m_parser->compile<N>();
        }

        /** Drop the compiled form; the formula is kept */
        void release () { m_parser.reset(); }

    private:
        std::string m_expression;                //!< the formula as given
        amrex::Vector<std::string> m_variables;  //!< names of the N variables
        //! compiled on demand, shared by copies until one of them changes its formula
        mutable std::shared_ptr<amrex::Parser> m_parser;
    };

} // namespace impactx::elements

#endif // IMPACTX_ELEMENTS_HELPER_TRACKEDPARSER_H
