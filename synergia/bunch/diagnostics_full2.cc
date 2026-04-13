
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/bunch/diagnostics_io.h"
#include "synergia/utils/simple_timer.h"

void calculate_emittances(double* mom2,
                          double& emitx,
                          double& emity,
                          double& emitz,
                          double& emitxy,
                          double& emitxyz);

void
Diagnostics_full2::do_update(Bunch const& bunch)
{
    scoped_simple_timer timer("diag_full2_update");

    ref = bunch.get_reference_particle();

    num_particles = bunch.get_total_num();
    real_num_particles = bunch.get_real_num();

    min = Core_diagnostics::calculate_min(bunch);
    max = Core_diagnostics::calculate_max(bunch);
    mean = Core_diagnostics::calculate_mean(bunch);
    mom2 = Core_diagnostics::calculate_mom2(bunch, mean);

    for (int i = 0; i < 6; ++i) {
        std(i) = std::sqrt(mom2(i, i));

        for (int j = 0; j < 6; ++j)
            corr(i, j) = mom2(i, j) / std::sqrt(mom2(i, i) * mom2(j, j));
    }

    calculate_emittances(mom2.data(), emitx, emity, emitz, emitxy, emitxyz);
}

void
Diagnostics_full2::do_first_write(io_device& file)
{
#ifdef SYNERGIA_HAVE_OPENPMD
    file.setAttribute("charge", ref.get_charge());
    file.setAttribute("mass", ref.get_four_momentum().get_mass());
    file.setAttribute("abs_offset", ref.get_bunch_abs_offset());
    file.flush();
#else
    file.write("charge", ref.get_charge());
    file.write("mass", ref.get_four_momentum().get_mass());
    file.write("abs_offset", ref.get_bunch_abs_offset());
#endif
    return;
}

void
Diagnostics_full2::do_write(io_device& file, size_t iteration)
{
    scoped_simple_timer timer("diag_full2_write");
#ifdef SYNERGIA_HAVE_OPENPMD
    auto i = file.iterations[iteration];
    i.setAttribute("s", ref.get_s());
    i.setAttribute("s_n", ref.get_s_n());
    i.setAttribute("repetition", ref.get_repetition());
    i.setAttribute("abs_time", ref.get_bunch_abs_time());
    i.setAttribute("num_particles", num_particles);
    i.setAttribute("real_num_particles", real_num_particles);
    i.setAttribute("pz", ref.get_momentum());

    // Once we move to C++20, we could request a span based attribute
    // from OpenPMD and prevent making copies!
    i.setAttribute("mean", Core_diagnostics::kokkos_view_to_stl_vector(mean));
    i.setAttribute("std", Core_diagnostics::kokkos_view_to_stl_vector(std));
    i.setAttribute("min", Core_diagnostics::kokkos_view_to_stl_vector(min));
    i.setAttribute("max", Core_diagnostics::kokkos_view_to_stl_vector(max));
    i.setAttribute("mom2", Core_diagnostics::kokkos_view_to_stl_vector(mom2));
    i.setAttribute("corr", Core_diagnostics::kokkos_view_to_stl_vector(corr));

    i.setAttribute("emitx", emitx);
    i.setAttribute("emity", emity);
    i.setAttribute("emitz", emitz);
    i.setAttribute("emitxy", emitxy);
    i.setAttribute("emitxyz", emitxyz);

    file.flush();
#else
    // write serial
    file.append("s", ref.get_s());
    file.append("s_n", ref.get_s_n());
    file.append("repetition", ref.get_repetition());
    file.append("abs_time", ref.get_bunch_abs_time());
    file.append("num_particles", num_particles);
    file.append("real_num_particles", real_num_particles);
    file.append("pz", ref.get_momentum());

    file.append("mean", mean);
    file.append("std", std);
    file.append("min", min);
    file.append("max", max);
    file.append("mom2", mom2);
    file.append("corr", corr);

    file.append("emitx", emitx);
    file.append("emity", emity);
    file.append("emitz", emitz);
    file.append("emitxy", emitxy);
    file.append("emitxyz", emitxyz);
#endif
    return;
}
