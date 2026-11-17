/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#include "InSitu.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_String.H>
#include <AMReX_Vector.H>

#ifdef AMREX_USE_CONDUIT
#   include <AMReX_Conduit_Blueprint.H>
#endif

#ifdef AMREX_USE_ASCENT
#   include <ascent.hpp>
#endif

#include <exception>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


namespace impactx::diagnostics::insitu
{
#ifdef AMREX_USE_CONDUIT
    void
    beam_to_blueprint (
        ImpactXParticleContainer const & pc,
        int step,
        int period,
        conduit::Node & mesh
    )
    {
        BL_PROFILE("impactx::diagnostics::insitu::beam_to_blueprint");

        // the coordinate set describes particles at fixed s
        if (pc.GetCoordSystem() != CoordSystem::s)
        {
            throw std::runtime_error(
                "In situ visualization: the particles must be at fixed s. "
                "Call coordinate_transformation(beam, CoordSystem.s) or leave the "
                "at_fixed_t() block first.");
        }
        if (mesh.number_of_children() != 0)
        {
            throw std::invalid_argument(
                "In situ visualization: the Blueprint mesh must be empty.");
        }

        std::vector<std::string> const real_soa_names = pc.GetRealSoANames();
        std::vector<std::string> const int_soa_names = pc.GetIntSoANames();
        amrex::Vector<std::string> const real_names(real_soa_names.begin(), real_soa_names.end());
        amrex::Vector<std::string> const int_names(int_soa_names.begin(), int_soa_names.end());

        // zero-copy: x, y, t of the particles are the coordinates, all components are fields
        amrex::ParticleContainerToBlueprint(pc, real_names, int_names, mesh, beam_topology);

        RefPart const & ref_part = pc.GetRefParticle();
        for (conduit::index_t d = 0; d < mesh.number_of_children(); ++d)
        {
            conduit::Node & state = mesh.child(d)["state"];
            state["cycle"] = step;
            state["time"] = static_cast<double>(ref_part.s);
            state["s"] = static_cast<double>(ref_part.s);
            state["period"] = period;
        }

        // Ascent and Catalyst read the particle arrays outside of the AMReX GPU stream
        amrex::Gpu::streamSynchronize();
    }
#endif

    namespace
    {
        /** All instances that are open, to close them before AMReX and MPI are finalized */
        std::set<Ascent::Impl *> &
        open_ascent_instances ()
        {
            static std::set<Ascent::Impl *> instances;
            return instances;
        }
    } // namespace

#ifdef AMREX_USE_ASCENT
    struct Ascent::Impl
    {
        ::ascent::Ascent ascent;

        /** The last published mesh
         *
         * Ascent refers to the published data until the next publish, including
         * the point connectivity that this node owns.
         */
        conduit::Node mesh;

        bool is_open = false;

        void close ()
        {
            if (is_open)
            {
                // deregister first: if close throws, do not try again
                open_ascent_instances().erase(this);
                is_open = false;
                ascent.close();
                mesh.reset();
            }
        }
    };
#else
    struct Ascent::Impl
    {
        bool is_open = false;

        void close () { is_open = false; }
    };
#endif

    Ascent::Ascent (std::string const & options_json)
    {
#ifdef AMREX_USE_ASCENT
        conduit::Node options;
        options.parse(options_json, "json");

        // report errors to the caller instead of printing and continuing
        if (!options.has_child("exceptions"))
        {
            options["exceptions"] = "forward";
        }
        // convert only the fields that the actions use, e.g., not the 64bit particle ids
        if (!options.has_child("field_filtering"))
        {
            options["field_filtering"] = "true";
        }
        // execute the actions passed to execute(): an empty name stops Ascent from
        // replacing them with an ascent_actions.json/.yaml file in the current directory
        if (!options.has_child("actions_file"))
        {
            options["actions_file"] = "";
        }
        // Ascent keeps the results of queries in a cache for the whole process and saves
        // it at exit to <default_dir>/<session_name>.yaml, with the directory of the first
        // Ascent opened in the process. Only both empty disable that file. Rendered images
        // and extracts are written to paths with a directory or to the current directory.
        if (!options.has_child("default_dir") && !options.has_child("session_name"))
        {
            options["default_dir"] = "";
            options["session_name"] = "";
        }
#   ifdef AMREX_USE_MPI
        options["mpi_comm"] = MPI_Comm_c2f(amrex::ParallelDescriptor::Communicator());
#   endif

        m_impl = std::make_unique<Impl>();
        m_impl->ascent.open(options);
        m_impl->is_open = true;
        open_ascent_instances().insert(m_impl.get());
#else
        amrex::ignore_unused(options_json);
        throw std::runtime_error(
            "Ascent: ImpactX was built without Ascent support (ImpactX_ASCENT=OFF).");
#endif
    }

    Ascent::Ascent (Ascent &&) noexcept = default;
    Ascent & Ascent::operator= (Ascent && other) noexcept
    {
        if (this != &other)
        {
            // a move-assigned instance must not leave its own Ascent open
            try { close(); }
            catch (std::exception const & e)
            {
                amrex::Print() << "Ascent: could not close: " << e.what() << "\n";
            }
            m_impl = std::move(other.m_impl);
        }
        return *this;
    }

    Ascent::~Ascent ()
    {
        try { close(); }
        catch (std::exception const & e)
        {
            amrex::Print() << "Ascent: could not close: " << e.what() << "\n";
        }
        catch (...)
        {
            amrex::Print() << "Ascent: could not close.\n";
        }
    }

    std::string
    Ascent::execute (
        [[maybe_unused]] ImpactXParticleContainer const * pc,
        [[maybe_unused]] std::string const & actions_json,
        [[maybe_unused]] int step,
        [[maybe_unused]] int period,
        [[maybe_unused]] std::string const & geometry_json
    )
    {
#ifdef AMREX_USE_CONDUIT
        // parsed where it is given, e.g., on MPI rank 0
        conduit::Node geometry;
        if (!geometry_json.empty())
        {
            geometry.parse(geometry_json, "json");
        }
        return execute(pc, actions_json, step, period,
                       geometry_json.empty() ? nullptr : &geometry);
#else
        // without Conduit, there is no Ascent to open
        if (!is_open())
        {
            throw std::runtime_error("Ascent: execute() after close().");
        }
        return "{}";
#endif
    }

#ifdef AMREX_USE_CONDUIT
    std::string
    Ascent::execute (
        [[maybe_unused]] ImpactXParticleContainer const * pc,
        [[maybe_unused]] std::string const & actions_json,
        [[maybe_unused]] int step,
        [[maybe_unused]] int period,
        [[maybe_unused]] conduit::Node const * geometry
    )
    {
        if (!is_open())
        {
            throw std::runtime_error("Ascent: execute() after close().");
        }
        // the geometry of MPI rank 0 is published: all ranks need to know if there is one
        int has_geometry = (geometry == nullptr || geometry->dtype().is_empty()) ? 0 : 1;
        amrex::ParallelDescriptor::Bcast(&has_geometry, 1, 0);
        if (pc == nullptr && has_geometry == 0)
        {
            throw std::invalid_argument("Ascent: execute() needs a beam or a geometry.");
        }

#ifdef AMREX_USE_ASCENT
        BL_PROFILE("impactx::diagnostics::insitu::Ascent::execute");

        conduit::Node actions;
        actions.parse(actions_json, "json");

        conduit::Node mesh;
        if (pc != nullptr)
        {
            beam_to_blueprint(*pc, step, period, mesh);
        }

        // optional geometry, e.g., of the lattice: one more domain, on MPI rank 0
        if (has_geometry != 0)
        {
            // the domains of the beam are numbered 0, 1, ... over all MPI ranks
            amrex::Long num_beam_domains = mesh.number_of_children();
            amrex::ParallelDescriptor::ReduceLongSum(num_beam_domains);
            if (amrex::ParallelDescriptor::MyProc() == 0)
            {
                // a named child, as the domains of the beam: append() would turn the
                // mesh into a list and drop them
                conduit::Node & domain =
                    mesh[amrex::Concatenate("domain_", static_cast<int>(num_beam_domains), 6)];
                // a copy: Ascent refers to the published data until the next publish
                domain.set(*geometry);
                domain["state/domain_id"] = static_cast<int>(num_beam_domains);
                domain["state/cycle"] = step;
                domain["state/period"] = period;
                if (pc != nullptr)
                {
                    RefPart const & ref_part = pc->GetRefParticle();
                    domain["state/time"] = static_cast<double>(ref_part.s);
                    domain["state/s"] = static_cast<double>(ref_part.s);
                }
            }
        }

        // Ascent refers to the published data until the next publish: keep it
        // alive and release the previously published mesh only afterwards
        m_impl->ascent.publish(mesh);
        m_impl->mesh.swap(mesh);

        m_impl->ascent.execute(actions);

        conduit::Node info;
        m_impl->ascent.info(info);
        return info.to_json();
#else
        return "{}";
#endif
    }
#endif

    void
    Ascent::close ()
    {
        if (m_impl)
        {
            m_impl->close();
        }
    }

    bool
    Ascent::is_open () const
    {
        return m_impl && m_impl->is_open;
    }

    void
    close_all ()
    {
        // close() deregisters, so iterate over a copy
        std::set<Ascent::Impl *> const instances = open_ascent_instances();
        for (Ascent::Impl * instance : instances)
        {
            instance->close();
        }
    }

} // namespace impactx::diagnostics::insitu
