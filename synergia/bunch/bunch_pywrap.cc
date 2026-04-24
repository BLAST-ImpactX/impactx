
#include <cstddef>
#include <optional>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "synergia/bunch/bunch.h"
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/bunch/populate.h"

#include "synergia/bunch/diagnostics_bulk_track.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/bunch/diagnostics_loss.h"
#include "synergia/bunch/diagnostics_particles.h"
#include "synergia/bunch/diagnostics_worker.h"

#include "synergia/bunch/diagnostics_py.h"

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE(bunch, m)
{
    py::enum_<ParticleGroup>(m, "ParticleGroup")
        .value("regular", ParticleGroup::regular)
        .value("spectator", ParticleGroup::spectator);

    py::enum_<LongitudinalBoundary>(m, "LongitudinalBoundary")
        .value("open", LongitudinalBoundary::open)
        .value("periodic", LongitudinalBoundary::periodic)
        .value("aperture", LongitudinalBoundary::aperture)
        .value("bucket_barrier", LongitudinalBoundary::bucket_barrier);

    py::class_<HostParticles>(m, "Particles", py::buffer_protocol())
        .def_buffer([](HostParticles const& p) -> py::buffer_info {
            return py::buffer_info(
                p.data(),       // pointer to buffer
                sizeof(double), // size of one scalar
                py::format_descriptor<double>::format(),
                2,                          // num of dimensions
                {p.extent(0), p.extent(1)}, // dimensions
                {p.stride(0) * sizeof(double), p.stride(1) * sizeof(double)}
                // strides (in bytes)
            );
        });

    py::class_<HostParticleMasks>(m, "ParticleMasks", py::buffer_protocol())
        .def_buffer([](HostParticleMasks const& p) -> py::buffer_info {
            return py::buffer_info(
                p.data(),       // pointer to buffer
                sizeof(uint8_t), // size of one scalar
                py::format_descriptor<uint8_t>::format(),
                1,                          // num of dimensions
                {p.extent(0)}, // dimensions
                {sizeof(uint8_t)}
                // strides (in bytes)
            );
        });

    py::class_<karray1d>(m, "karray1d", py::buffer_protocol())
        .def_buffer([](karray1d const& p) -> py::buffer_info {
            return py::buffer_info(p.data(),
                                   sizeof(double),
                                   py::format_descriptor<double>::format(),
                                   1,
                                   {p.extent(0)},
                                   {p.stride(0) * sizeof(double)});
        });

    py::class_<karray2d_row>(m, "karray2d_row", py::buffer_protocol())
        .def_buffer([](karray2d_row const& p) -> py::buffer_info {
            return py::buffer_info(
                p.data(),
                sizeof(double),
                py::format_descriptor<double>::format(),
                2,
                {p.extent(0), p.extent(1)},
                {p.stride(0) * sizeof(double), p.stride(1) * sizeof(double)});
        });

    // Bunch
    py::class_<Bunch>(m, "Bunch")
        .def(py::init<Reference_particle const&,
                      int,
                      double,
                      Commxx,
                      int,
                      int,
                      int,
                      int>(),
             "Construct a Bunch object.",
             "reference_particle"_a,
             "total_num"_a,
             "real_num"_a,
             "comm"_a = Commxx(),
             "total_spectator_num"_a = 0,
             "bunch_index"_a = 0,
             "bucket_index"_a = 0,
             "array_index"_a = 0)

        .def("read_file_legacy",
             &Bunch::read_file_legacy,
             "Read particle data from a lagecy file (produced by Synergia2).",
             "filename"_a)

        .def("read_file",
             &Bunch::read_file,
             "Read particle data from a file.",
             "filename"_a)

        .def("write_file",
             py::overload_cast<std::string const&, int, int, int, int>(
                 &Bunch::write_file, py::const_),
             "filename"_a,
             "num_part"_a = -1,
             "offset"_a = 0,
             "num_part_spec"_a = 0,
             "offset_spec"_a = 0,
             "Write particle data to a file.")

#if defined SYNERGIA_HAVE_OPENPMD
        .def("read_openpmd_file",
             &Bunch::read_openpmd_file,
             "Read particle data from an OpenPMD file.",
             "filename"_a,
             "idx"_a = std::optional<std::size_t>(std::nullopt))

        .def("write_openpmd_file",
             &Bunch::write_openpmd_file,
             "filename"_a,
             "num_part"_a = -1,
             "offset"_a = 0,
             "num_part_spec"_a = 0,
             "offset_spec"_a = 0,
             "Write particle data to an OpenPMD file.")

#endif

        .def("get_comm", &Bunch::get_comm, "Get the communicator.")

        .def("get_reference_particle",
             (Reference_particle & (Bunch::*)()) &
                 Bunch::get_reference_particle,
             py::return_value_policy::reference,
             "Get the reference particle of the bunch. Returns reference to "
             "object")

        .def("get_design_reference_particle",
             (Reference_particle & (Bunch::*)()) &
                 Bunch::get_design_reference_particle,
             py::return_value_policy::reference,
             "Get the design reference particle of the bunch. Returns "
             "reference to object")

        .def("checkout_particles",
             &Bunch::checkout_particles,
             "Copy the particle array from device memory to host memory.",
             "particle_group"_a = ParticleGroup::regular)

        .def("checkin_particles",
             &Bunch::checkin_particles,
             "Copy the particle array from device memory to host memory.",
             "particle_group"_a = ParticleGroup::regular)

        .def("get_host_particles",
             (HostParticles(Bunch::*)(ParticleGroup)) &
                 Bunch::get_host_particles,
             // py::overload_cast<ParticleGroup>(&Bunch::get_host_particles),
             "Get host particles.",
             "particle_group"_a = ParticleGroup::regular)

        .def("get_host_particle_masks",
             (HostParticleMasks(Bunch::*)(ParticleGroup)) &
                 Bunch::get_host_particle_masks,
             // py::overload_cast<ParticleGroup>(&Bunch::get_host_particles),
             "Get host particle masks.",
             "particle_group"_a = ParticleGroup::regular)

        .def("size",
             &Bunch::size,
             "Get the number of particles (regular and lost particles)",
             "particle_group"_a = ParticleGroup::regular)

        .def(
            "capacity",
            &Bunch::capacity,
            "Get the number of all particle slots (regualr, lost, and reserved "
            "particles)",
            "particle_group"_a = ParticleGroup::regular)

        .def("reserve",
             &Bunch::reserve,
             "Reserve the particle slots to the new capcity",
             "capacity"_a,
             "particle_group"_a = ParticleGroup::regular)

        .def("get_local_num",
             &Bunch::get_local_num,
             "Get the number of valid particles in current rank",
             "particle_group"_a = ParticleGroup::regular)

        .def("get_total_num",
             &Bunch::get_total_num,
             "Get the number of particles in the bunch across all rank",
             "particle_group"_a = ParticleGroup::regular)

        .def("get_real_num",
             &Bunch::get_real_num,
             "Get the charge represented by particles in this bunch")

        .def("set_real_num",
            &Bunch::set_real_num,
            "Set the real number of charges represented by this bunch",
            "num"_a)

        .def("set_longitudinal_boundary",
             &Bunch::set_longitudinal_boundary,
             "Set the longitudinal boundary of the bunch",
             "boundary"_a,
             "param"_a = 0.0)

        .def("get_longitudinal_boundary",
             &Bunch::get_longitudinal_boundary,
             "Get the longitudinal boundary of the bunch")

        .def("adjust_bunch_particles_reference_energy",
            &Bunch::adjust_bunch_particles_reference_energy,
            "Set the reference energy of the bunch and adjust the normalization of the particles to match",
            "newE"_a)

        .def("offset_bunch_particles_cdt",
            &Bunch::offset_bunch_particles_cdt,
            "Add offset to the cdt of all bunch particles",
            "offset"_a)

        .def(
            "add_diagnostics",
            [](Bunch& self, std::shared_ptr<Diagnostics> const& diag) {
                // if the diagnostics is an inherited python type,
                // reg a ref of the python instance to keep it alive
                // so the __dict__ object which contains the actual
                // python methods is available for serializing and
                // calling from C++
                PyDiagnostics* p = dynamic_cast<PyDiagnostics*>(diag.get());

                if (p) p->reg_self();

                return self.add_diagnostics(diag);
            },
            "Add a diagnostics to the bunch object.",
            "diag"_a)

        .def("set_diag_loss_aperture",
             &Bunch::set_diag_loss_aperture,
             "Set the loss diagnostics for apertures.",
             "filename"_a)

        .def("set_diag_loss_zcut",
             &Bunch::set_diag_loss_zcut,
             "Set the loss diagnostics for z-cut.",
             "filename"_a)

        .def("diag_type",
             &Bunch::diag_type,
             "Get the type of the named diagnostics.",
             "name"_a)

        .def("diag_update",
             &Bunch::diag_update,
             "Performs the update on the named diagnostics.",
             "name"_a)

        .def("diag_update_and_write",
             &Bunch::diag_update_and_write,
             "Performs the update and write on the named diagnostics.",
             "name"_a)

        .def("set_bucket_index", &Bunch::set_bucket_index, "index"_a)

        .def("is_bucket_index_assigned", &Bunch::is_bucket_index_assigned)

        .def("inject", &Bunch::inject)

        .def("print_statistics", &Bunch::print_statistics)

        .def("print_particle",
             &Bunch::print_particle,
             "idx"_a,
             "logger"_a,
             "particle_group"_a = ParticleGroup::regular)

        .def("dump", &Bunch::dump, "Dump.")

        .def("load", &Bunch::load, "Load.");

    // Core_diagnostics
    py::class_<Core_diagnostics>(m, "Core_diagnostics")
        .def_static("calculate_mean_ka",
                    &Core_diagnostics::calculate_mean,
                    "Calculate the mean for the bunch.",
                    "bunch"_a)

        .def_static("calculate_abs_mean_ka",
                    &Core_diagnostics::calculate_abs_mean,
                    "Calculate the mean of absolute values for the bunch.",
                    "bunch"_a)

        .def_static(
            "calculate_std_ka",
            [](Bunch const& bunch, py::buffer b) {
                py::buffer_info info = b.request();
                Kokkos::View<double*,
                             Kokkos::HostSpace,
                             Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                    array((double*)info.ptr, info.shape[0]);
                return Core_diagnostics::calculate_std(bunch, array);
            },
            "Calculate the standard deviation for the bunch.",
            "bunch"_a,
            "mean"_a)

        .def_static(
            "calculate_mom2_ka",
            [](Bunch const& bunch, py::buffer b) {
                py::buffer_info info = b.request();
                Kokkos::View<double*,
                             Kokkos::HostSpace,
                             Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                    array((double*)info.ptr, info.shape[0]);
                return Core_diagnostics::calculate_mom2(bunch, array);
            },
            "Calculate the second order moments for the bunch.",
            "bunch"_a,
            "mean"_a);

    // Diagnostics_handler
    py::class_<Diagnostics_handler>(m, "Diagnostics_handler")
        .def("type",
             &Diagnostics_handler::type,
             "Get the type of the registered diagnostics.")

        .def("update", &Diagnostics_handler::type, "Update the diagnostics.")

        .def("write", &Diagnostics_handler::write, "Write the diagnostics.")

        .def("update_and_write",
             &Diagnostics_handler::update_and_write,
             "Update and write the diagnostics.");

    // Diagnostics_calculator base class
    py::class_<Diagnostics, PyDiagnostics, std::shared_ptr<Diagnostics>>(
        m, "Diagnostics")
        .def(py::init<std::string const&, std::string const&, bool>(),
             "Construct a Diagnostics object.",
             "type"_a = "py_diag",
             "filename"_a = "py_diag.h5",
             "serial"_a = true)

        .def("type", &Diagnostics::type, "Get the type of the Diagnostics.")

        .def("filename",
             &Diagnostics::filename,
             "Returns the filename of the Diagnostics.")

#if 0
        .def( "__reduce__",
                [](py::object const& self) {
                    std::cout << "diag reduce\n";
                    return py::tuple(py::cpp_function([](){
                        return PyDiagnostics("reduced pydiag");
                    }), py::tuple() );
                } )
#endif

        .def(py::pickle(
            [](py::object self) {
                return py::make_tuple(self.attr("__dict__"));
            },
            [](py::tuple const& t) {
                PyDiagnostics cpp_diag("restored pydiag", "py_diag.h5", true);
                auto py_diag = t[0].cast<py::dict>();
                return std::make_pair(cpp_diag, py_diag);
            }))

        ;

    // Diagnostics_dummy
    py::class_<Diagnostics_dummy,
               Diagnostics,
               std::shared_ptr<Diagnostics_dummy>>(m, "Diagnostics_dummy")
        .def(py::init<>());

    py::class_<Diagnostics_full2,
               Diagnostics,
               std::shared_ptr<Diagnostics_full2>>(m, "Diagnostics_full2")
        .def(py::init<std::string const&>(),
             "Construct a Diagnostics_full2 object.",
             "filename"_a = "diag_full2.h5");

    py::class_<Diagnostics_bulk_track,
               Diagnostics,
               std::shared_ptr<Diagnostics_bulk_track>>(
        m, "Diagnostics_bulk_track")
        .def(py::init<std::string const&, int, int, ParticleGroup>(),
             "Construct a Diagnostics_bulk_track object.",
             "filename"_a = "diag_bulk_track.h5",
             "num_tracks"_a = 0,
             "offset"_a = 0,
             "particlegroup"_a = ParticleGroup::regular);

    py::class_<Diagnostics_particles,
               Diagnostics,
               std::shared_ptr<Diagnostics_particles>>(m,
                                                       "Diagnostics_particles")
        .def(py::init<std::string const&, int, int, int, int>(),
             "Construct a Diagnostics_particles object.",
             "filename"_a = "diag_particles.h5",
             "num_part"_a = -1,
             "offset"_a = 0,
             "num_spec_part"_a = 0,
             "spec_offset"_a = 0);

    // populate
    // m.def( "populate_6d", populate_6d );

    m.def("populate_6d",
          [](Distribution& dist,
             Bunch& bunch,
             py::buffer p_means,
             py::buffer p_covars) {
              using ka1d_unmanaged =
                  Kokkos::View<double*,
                               Kokkos::HostSpace,
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

              using ka2d_unmanaged =
                  Kokkos::View<double**,
                               Kokkos::LayoutRight,
                               Kokkos::HostSpace,
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

              py::buffer_info pm_info = p_means.request();
              py::buffer_info pc_info = p_covars.request();

              if (pm_info.ndim != 1 || pm_info.shape[0] != 6) {
                  throw std::runtime_error(
                      "populate_6d: "
                      "mean must be an 1d array of 6 elements");
              }

              if (pc_info.ndim != 2 || pc_info.shape[0] != 6 ||
                  pc_info.shape[1] != 6) {
                  throw std::runtime_error(
                      "populate_6d: "
                      "covariances must be a 2d array of 6x6");
              }

              ka1d_unmanaged means((double*)pm_info.ptr, pm_info.shape[0]);

              ka2d_unmanaged covars(
                  (double*)pc_info.ptr, pc_info.shape[0], pc_info.shape[1]);

              populate_6d(dist, bunch, means, covars);
          });

    m.def("populate_6d_truncated",
          [](Distribution& dist,
             Bunch& bunch,
             py::buffer p_means,
             py::buffer p_covars,
             py::buffer p_limits) {
              using ka1d_unmanaged =
                  Kokkos::View<double*,
                               Kokkos::HostSpace,
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

              using ka2d_unmanaged =
                  Kokkos::View<double**,
                               Kokkos::LayoutRight,
                               Kokkos::HostSpace,
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

              py::buffer_info pm_info = p_means.request();
              py::buffer_info pl_info = p_limits.request();
              py::buffer_info pc_info = p_covars.request();

              if (pm_info.ndim != 1 || pm_info.shape[0] != 6) {
                  throw std::runtime_error(
                      "populate_6d_truncated: "
                      "mean must be an 1d array of 6 elements");
              }

              if (pl_info.ndim != 1 || pl_info.shape[0] != 6) {
                  throw std::runtime_error(
                      "populate_6d_truncated: "
                      "limits must be an 1d array of 6 elements");
              }

              if (pc_info.ndim != 2 || pc_info.shape[0] != 6 ||
                  pc_info.shape[1] != 6) {
                  throw std::runtime_error(
                      "populate_6d_truncated: "
                      "covariances must be a 2d array of 6x6");
              }

              ka1d_unmanaged means((double*)pm_info.ptr, pm_info.shape[0]);
              ka1d_unmanaged limits((double*)pl_info.ptr, pl_info.shape[0]);

              ka2d_unmanaged covars(
                  (double*)pc_info.ptr, pc_info.shape[0], pc_info.shape[1]);

              populate_6d_truncated(dist, bunch, means, covars, limits);
          });

    m.def("populate_transverse_gaussian",
          [](Distribution& dist,
             Bunch& bunch,
             py::buffer p_means,
             py::buffer p_covars,
             double cdt) {
              using ka1d_unmanaged =
                  Kokkos::View<double*,
                               Kokkos::HostSpace,
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

              using ka2d_unmanaged =
                  Kokkos::View<double**,
                               Kokkos::LayoutRight,
                               Kokkos::HostSpace,
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

              py::buffer_info pm_info = p_means.request();
              py::buffer_info pc_info = p_covars.request();

              if (pm_info.ndim != 1 || pm_info.shape[0] != 6) {
                  throw std::runtime_error(
                      "populate_transverse_gaussian: "
                      "mean must be an 1d array of 6 elements");
              }

              if (pc_info.ndim != 2 || pc_info.shape[0] != 6 ||
                  pc_info.shape[1] != 6) {
                  throw std::runtime_error(
                      "populate_transverse_gaussian: "
                      "covariances must be a 2d array of 6x6");
              }

              ka1d_unmanaged means((double*)pm_info.ptr, pm_info.shape[0]);

              ka2d_unmanaged covars(
                  (double*)pc_info.ptr, pc_info.shape[0], pc_info.shape[1]);

              populate_transverse_gaussian(dist, bunch, means, covars, cdt);
          });

    m.def("get_correlation_matrix",
          [](py::array_t<double> map,
             double arms,
             double brms,
             double crms,
             double beta,
             std::array<int, 3> const& rms_index = {0, 2, 4}) {
              if (map.ndim() != 2 || map.shape(0) != 6 || map.shape(1) != 6)
                  throw std::runtime_error(
                      "get_correlation_matrix(): "
                      "one_turn_map must be an array of (6, 6)");

              using ka2d_unmanaged =
                  Kokkos::View<double**,
                               Kokkos::LayoutRight,
                               Kokkos::HostSpace,
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

              ka2d_unmanaged ka_map(map.mutable_data(), 6, 6);

              auto matrix = get_correlation_matrix(
                  ka_map, arms, brms, crms, beta, rms_index);

              auto arr = py::array_t<double>({6, 6});

              for (int i = 0; i < 6; ++i)
                  for (int j = 0; j < 6; ++j)
                      arr.mutable_at(i, j) = matrix(i, j);

              return arr;
          });
}
