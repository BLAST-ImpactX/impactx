/* Copyright 2021-2023 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"

#include <particles/ImpactXParticleContainer.H>
#include <diagnostics/ReducedBeamCharacteristics.H>

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParticleContainer.H>

#include <algorithm>
#include <string>
#include <type_traits>
#include <vector>

namespace py = pybind11;
using namespace impactx;


namespace
{
    /** Register the ImpactX beam particle container and its iterators at one
     *  concrete precision T_PC = ImpactXParticleContainerT<float|double>.
     *
     * The Python class names are suffixed with @p suffix (empty for the build's
     * native/default precision so the historic names stay stable; e.g. "_sp" for
     * the single-precision storage container in a multi-precision binary). The
     * AMReX base classes are matched by C++ type, so they resolve to the right
     * (precision-distinct) pyAMReX registrations.
     */
    template <class T_PC>
    void register_particle_container (py::module& m, std::string const & suffix)
    {
        using iterator = typename T_PC::iterator;
        using const_iterator = typename T_PC::const_iterator;

        py::class_<
            iterator,
            typename iterator::base_t
        > py_pariter_soa(m, ("ImpactXParIter" + suffix).c_str());
        py_pariter_soa
            .def(py::init<typename iterator::ContainerType&, int>(),
                 py::arg("particle_container"), py::arg("level"))
            .def(py::init<typename iterator::ContainerType&, int, amrex::MFItInfo&>(),
                 py::arg("particle_container"), py::arg("level"), py::arg("info"))
        ;

        py::class_<
            const_iterator,
            typename const_iterator::base_t
        > py_parconstiter_soa(m, ("ImpactXParConstIter" + suffix).c_str());
        py_parconstiter_soa
            .def(py::init<typename const_iterator::ContainerType&, int>(),
                 py::arg("particle_container"), py::arg("level"))
            .def(py::init<typename const_iterator::ContainerType&, int, amrex::MFItInfo&>(),
                 py::arg("particle_container"), py::arg("level"), py::arg("info"))
        ;

        py::class_<
            T_PC,
            typename T_PC::base_t
        >(m, ("ImpactXParticleContainer" + suffix).c_str())
            //.def(py::init<>())

            .def_property_readonly("coord_system",
                &T_PC::GetCoordSystem,
                "Get the current coordinate system of particles in this container"
            )

            // simpler particle iterator loops: return types of this particle box
            // note: overwritten to return ImpactX instead of (py)AMReX iterators
            .def_property_readonly_static(
                "Iterator",
                [](py::object /* pc */){ return py::type::of<iterator>(); },
                "ImpactX iterator for particle boxes"
            )
            .def_property_readonly_static(
                "ConstIterator",
                [](py::object /* pc */){ return py::type::of<const_iterator>(); },
                "ImpactX constant iterator for particle boxes (read-only)"
            )

            .def("clear", &T_PC::clear,
                 py::arg("keep_mass")=false, py::arg("keep_charge")=false,
                 "Empty the container and reset the reference particle"
            )
            // internal: there is a high-level Python wrapper that accepts more types
            .def("_add_n_particles",
                 &T_PC::AddNParticles,
                 py::arg("x"), py::arg("y"), py::arg("t"),
                 py::arg("px"), py::arg("py"), py::arg("pt"),
                 py::arg("qm"), py::arg("bunch_charge")=py::none(), py::arg("w")=py::none(),
                 py::arg("sx")=py::none(), py::arg("sy")=py::none(), py::arg("sz")=py::none(),
                 "Add new particles to the container for fixed s.\n\n"
                 "Either the total charge (bunch_charge) or the weight of each\n"
                 "particle (w) must be provided.\n\n"
                 "Note: This can only be used *after* the initialization (grids) have\n"
                 "      been created, meaning after the call to ImpactX.init_grids\n"
                 "      has been made in the ImpactX class.\n\n"
                 ":param x: positions in x\n"
                 ":param y: positions in y\n"
                 ":param t: positions as time-of-flight in c*t\n"
                 ":param px: momentum in x\n"
                 ":param py: momentum in y\n"
                 ":param pt: momentum in t\n"
                 ":param qm: charge over mass in 1/eV\n"
                 ":param bunch_charge: total charge within a bunch in C"
                 ":param w: weight of each particle: how many real particles to represent"
                 ":param sx: spin component in x\n"
                 ":param sy: spin component in y\n"
                 ":param sz: spin component in z\n"
            )
            // Getter-only property is intentional: it returns the live mutable
            // RefPart by reference.
            // A writable property (with assignment) would be ambiguous:
            // alias (Pythonic) or copy-in (safe) for the simulation-owned RefPart.
            .def_property_readonly("ref",
                [](T_PC & pc) -> RefPart & {
                    return pc.GetRefParticle();
                },
                "Access the reference particle."
            )
            .def("ref_particle",
                [](T_PC & pc) -> RefPart & {
                    py::warnings::warn(
                        "ref_particle() is deprecated. Use beam.ref instead.",
                        PyExc_DeprecationWarning,
                        2
                    );
                    return pc.GetRefParticle();
                },
                py::return_value_policy::reference_internal,
                "Access the reference particle.\n\n"
                "Deprecated: use ``beam.ref``."
            )
            .def("set_ref_particle",
                 &T_PC::SetRefParticle,
                 py::arg("refpart"),
                 "Set reference particle attributes."
            )
            .def("set_bucket_length",
                 &T_PC::SetBucketLength,
                 py::arg("bucket_length"),
                 "Set bucket length for particle boundary condition."
            )
            .def("min_and_max_positions",
                 &T_PC::MinAndMaxPositions,
                 "Compute the min and max of the particle position in each dimension.\n\n"
                 ":return: x_min, y_min, z_min, x_max, y_max, z_max"
            )
            .def("mean_and_std_positions",
                 &T_PC::MeanAndStdPositions,
                 "Compute the mean and std of the particle position in each dimension.\n\n"
                 ":return: x_mean, x_std, y_mean, y_std, z_mean, z_std"
            )
            .def("reduced_beam_characteristics",
                 [](T_PC & pc) {
                     py::warnings::warn(
                        "WARNING: reduced_beam_characteristics() is deprecated. "
                        "Use beam_moments() instead.",
                        PyExc_DeprecationWarning,
                        2
                     );
                     return diagnostics::reduced_beam_characteristics(pc);
                 },
                 "Compute reduced beam characteristics like the position and momentum moments of the particle distribution, as well as emittance and Twiss parameters."
            )
            .def("beam_moments",
                 [](T_PC & pc) {
                     return pc.beam_moments();
                 },
                 "Calculate beam moments at current ``s`` like the position and momentum moments of the particle "
                 "distribution, as well as emittance and Twiss parameters."
            )
            .def("beam_moments_history_list",
                 [](T_PC & pc) {
                     return pc.beam_moments_history();
                 },
                 "Return the history of the beam moments on every step."
            )
            .def("record_beam_moments",
                 [](T_PC & pc) {
                     return pc.record_beam_moments();
                 },
                 "Calculate & record the beam moments at current s"
            )
            .def("reset_beam_moments_history",
                 [](T_PC & pc) {
                     return pc.reset_beam_moments_history();
                 },
                 "Reset the history of the beam moments."
            )
            .def_property("store_beam_moments",
                [](T_PC & pc){ return pc.store_beam_moments; },
                [](T_PC & pc, bool record){ pc.store_beam_moments = record; },
                "In situ calculate and store the beam moments for every simulation step."
            )
        ;

        py_pariter_soa.def("pc", &iterator::pc);
        py_parconstiter_soa.def("pc", &const_iterator::pc);
    }
} // namespace


void init_impactxparticlecontainer(py::module& m)
{
    py::enum_<CoordSystem>(m, "CoordSystem")
        .value("s", CoordSystem::s)
        .value("t", CoordSystem::t)
        .export_values();

    // Register the beam container for each compiled precision. The build's
    // native precision (== amrex::ParticleReal) keeps the historic, unsuffixed
    // class names; an additionally-compiled precision gets a suffix.
#ifdef IMPACTX_COMPILE_DOUBLE
    register_particle_container<ImpactXParticleContainerT<double>>(
        m, std::is_same_v<amrex::ParticleReal, double> ? "" : "_dp");
#endif
#ifdef IMPACTX_COMPILE_SINGLE
    register_particle_container<ImpactXParticleContainerT<float>>(
        m, std::is_same_v<amrex::ParticleReal, float> ? "" : "_sp");
#endif
}
