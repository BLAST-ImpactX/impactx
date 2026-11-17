/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_DIAGNOSTICS_INSITU_H
#define IMPACTX_DIAGNOSTICS_INSITU_H

#include "particles/ImpactXParticleContainer.H"

#include <AMReX_Config.H>

#ifdef AMREX_USE_CONDUIT
#   include <conduit/conduit.hpp>
#endif

#include <memory>
#include <string>


namespace impactx::diagnostics::insitu
{
    /** Name of the Blueprint topology, coordinate set and field prefix of the beam */
    static constexpr auto beam_topology = "beam";

#ifdef AMREX_USE_CONDUIT
    /** Describe the particle beam as a Conduit Mesh Blueprint
     *
     * The particle data is not copied: the Blueprint arrays point to the particle
     * arrays, in host or device memory, wherever the particles are stored. Each
     * particle tile becomes one domain of a multi-domain mesh. The coordinate set
     * holds the positions x, y and t of the particles relative to the reference
     * particle, at fixed s. Every particle component is a field, named as in
     * ImpactXParticleContainer::GetRealSoANames, e.g., ``momentum_t``.
     *
     * The state of each domain holds the step as ``cycle``, the reference particle
     * position s in meters as ``time`` and ``s``, and the ``period``.
     *
     * The returned nodes are only valid until the particles change, e.g., until
     * the next push.
     *
     * This is MPI-collective.
     *
     * @param[in] pc particle beam, at fixed s
     * @param[in] step global step
     * @param[in] period current period (turn or cycle)
     * @param[out] mesh Blueprint mesh, empty on input
     */
    void
    beam_to_blueprint (
        ImpactXParticleContainer const & pc,
        int step,
        int period,
        conduit::Node & mesh
    );
#endif

    /** Ascent in situ visualization of the particle beam
     *
     * Keeps one Ascent instance open between calls to execute(), so the
     * Ascent pipelines are only rebuilt when the actions change.
     */
    class Ascent
    {
    public:
        /** Open Ascent
         *
         * The MPI communicator is set by ImpactX.
         * Unless set in options_json, Ascent forwards errors as exceptions, only
         * converts the fields that the actions use, executes exactly the
         * actions passed to execute() (``actions_file`` is empty), and does not
         * save its query results to a session file at exit (``default_dir`` and
         * ``session_name`` are empty: images without a directory in their name
         * are written to the current directory).
         *
         * @param[in] options_json Ascent options as JSON, e.g., {"default_dir": "path"}
         * @throws std::runtime_error if ImpactX is built without Ascent
         */
        explicit Ascent (std::string const & options_json = "{}");

        Ascent (Ascent const &) = delete;
        Ascent & operator= (Ascent const &) = delete;
        Ascent (Ascent &&) noexcept;
        Ascent & operator= (Ascent &&) noexcept;

        ~Ascent ();

        /** Publish the particle beam and execute Ascent actions
         *
         * This is MPI-collective. Rendered images are written by MPI rank 0.
         *
         * @param[in] pc particle beam, at fixed s, or nullptr to publish only the geometry
         * @param[in] actions_json Ascent actions as JSON (a list of actions)
         * @param[in] step global step
         * @param[in] period current period (turn or cycle)
         * @param[in] geometry_json optional: a domain of a Conduit Mesh Blueprint as JSON,
         *            e.g., of the lattice, that MPI rank 0 publishes next to the beam;
         *            its names must differ from the ones of the beam. Only the one of
         *            MPI rank 0 is used: other ranks can pass an empty string.
         * @return Ascent info as JSON, e.g., with the names of written images
         * @throws std::invalid_argument without a beam and without a geometry on MPI rank 0
         */
        std::string
        execute (
            ImpactXParticleContainer const * pc,
            std::string const & actions_json,
            int step,
            int period,
            std::string const & geometry_json = ""
        );

#ifdef AMREX_USE_CONDUIT
        /** Publish the particle beam and execute Ascent actions
         *
         * As execute() above, with the geometry as a Conduit node: large geometries
         * skip the conversion to and from JSON.
         *
         * @param[in] pc particle beam, at fixed s, or nullptr to publish only the geometry
         * @param[in] actions_json Ascent actions as JSON (a list of actions)
         * @param[in] step global step
         * @param[in] period current period (turn or cycle)
         * @param[in] geometry a domain of a Conduit Mesh Blueprint, e.g., of the lattice,
         *            that MPI rank 0 publishes (a copy of it) next to the beam; its names
         *            must differ from the ones of the beam. nullptr or an empty node for
         *            none. Only the one of MPI rank 0 is used.
         * @return Ascent info as JSON, e.g., with the names of written images
         * @throws std::invalid_argument without a beam and without a geometry on MPI rank 0
         */
        std::string
        execute (
            ImpactXParticleContainer const * pc,
            std::string const & actions_json,
            int step,
            int period,
            conduit::Node const * geometry
        );
#endif

        /** Close Ascent and release the published data
         *
         * This is MPI-collective. Calling it again does nothing.
         */
        void close ();

        /** Is Ascent open? */
        [[nodiscard]] bool is_open () const;

        struct Impl;

    private:
        std::unique_ptr<Impl> m_impl;
    };

    /** Close all open in situ visualization instances
     *
     * Called before AMReX and MPI are finalized. Instances that are still
     * referenced afterwards stay closed.
     *
     * This is MPI-collective.
     */
    void close_all ();

} // namespace impactx::diagnostics::insitu

#endif // IMPACTX_DIAGNOSTICS_INSITU_H
