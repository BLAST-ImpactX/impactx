#include "impedance.h"

#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/physical_constants.h"

#include "synergia/bunch/period.h"
#include "synergia/utils/kokkos_utils.h"

#include "synergia/utils/simple_timer.h"

#include <Kokkos_ScatterView.hpp>

typedef Kokkos::TeamPolicy<> team_policy;
typedef typename team_policy::member_type team_member;

using scatter_t =
    Kokkos::Experimental::ScatterView<double*, Kokkos::LayoutLeft>;

namespace {

    void
    zero_karray(karray1d_dev const& arr)
    {
        ku::alg_zeroer alg{arr};
        Kokkos::parallel_for(arr.extent(0), alg);
    }

    struct alg_write_bps {
        Bunch_props bps;
        karray1d_dev vbi_buf;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            int idx = bps.get_write_index(i);

            bps.xmean(idx) = vbi_buf(i * 5 + 0);
            bps.ymean(idx) = vbi_buf(i * 5 + 1);
            bps.zmean(idx) = vbi_buf(i * 5 + 2);
            bps.realnum(idx) = vbi_buf(i * 5 + 3);
            bps.bucket_index(idx) = (int)vbi_buf(i * 5 + 4);
        }
    };

    struct alg_z_binning {
        typedef double value_type[];

        const int value_count;
        const int z_grid;

        ConstParticles p;
        ConstParticleMasks masks;

        const double z_left;
        const double recip_h;

        alg_z_binning(ConstParticles const& p,
                      ConstParticleMasks const& masks,
                      int z_grid,
                      double z_left,
                      double h)
            : value_count(z_grid * 3)
            , z_grid(z_grid)
            , p(p)
            , masks(masks)
            , z_left(z_left)
            , recip_h(1.0 / h)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i, value_type sum) const
        {
            if (masks(i)) {
                int bin = (p(i, 4) - z_left) * recip_h;

                if (bin == z_grid) {
                    sum[z_grid * 0 + bin - 1] += 1;
                } else {
                    sum[z_grid * 0 + bin] += 1;       // zdensity
                    sum[z_grid * 1 + bin] += p(i, 0); // xmom
                    sum[z_grid * 2 + bin] += p(i, 2); // ymom
                }
            }
        }
    };

    struct alg_z_binning_sv {
        ConstParticles p;
        ConstParticleMasks masks;
        scatter_t scatter;

        const int z_grid;
        const double z_left;
        const double recip_h;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            if (masks(i)) {
                auto access = scatter.access();

                int bin = (p(i, 4) - z_left) * recip_h;
                if (bin == z_grid) {
                    access(z_grid * 0 + bin - 1) += 1;
                } else {
                    access(z_grid * 0 + bin) += 1;       // zdensity
                    access(z_grid * 1 + bin) += p(i, 0); // xmom
                    access(z_grid * 2 + bin) += p(i, 2); // ymom
                }
            }
        }
    };

    struct alg_z_normalize {
        karray1d_dev binning;
        const int z_grid;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            double zden = binning(z_grid * 0 + i);

            if (zden) {
                binning(z_grid * 1 + i) /= zden;
                binning(z_grid * 2 + i) /= zden;
            }
        }
    };

    KOKKOS_INLINE_FUNCTION
    int
    get_zindex_for_wake(double z, double dz, int istart, double zstart)
    {
        // if  (z< (-istart*istart*dz+zstart)) return -100;
        if (z >= zstart)
            return (static_cast<int>(floor(sqrt((z - zstart) / dz)))) + istart;
        else
            return (-static_cast<int>(floor(sqrt((zstart - z) / dz)))) + istart;
    }

    struct alg_z_wake_reduce {
        typedef kt::array_type<double, 5> value_type;

        const int i;
        const int z_grid;

        const double cell_size_z;
        const double N_factor;

        const int size_wake;
        const int istart;
        const double zstart;
        const double delta_z;

        karray1d_dev const& zbins;
        karray1d_dev const& wf;

        KOKKOS_INLINE_FUNCTION
        alg_z_wake_reduce(int i,
                          int z_grid,
                          Bunch_params const& bp,
                          Wake_field const& wf,
                          karray1d_dev const& zbins)
            : i(i)
            , z_grid(z_grid)
            , cell_size_z(bp.cell_size_z)
            , N_factor(bp.N_factor)
            , size_wake(wf.size_wake)
            , istart(wf.istart)
            , zstart(wf.zstart)
            , delta_z(wf.delta_z)
            , zbins(zbins)
            , wf(wf.terms)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int j, value_type& sum) const
        {
            double* zdensity = &zbins(z_grid * 0);
            double* xmom = &zbins(z_grid * 1);
            double* ymom = &zbins(z_grid * 2);

            double* z_coord = &wf(size_wake * 0);
            double* z_wake = &wf(size_wake * 1);
            double* xw_lead = &wf(size_wake * 2);
            double* xw_trail = &wf(size_wake * 3);
            double* yw_lead = &wf(size_wake * 4);
            double* yw_trail = &wf(size_wake * 5);

            double zji = (j - i) * cell_size_z;
            if (zji < z_coord[0]) return;

            // below it is assumed the wake function is stored using a quadratic
            // grid
            int iz = get_zindex_for_wake(zji, delta_z, istart, zstart);

            if (iz + 1 < size_wake) {
                double z1 = zji - z_coord[iz];
                double recip_z2 = 1.0 / (z_coord[iz + 1] - z_coord[iz]);

                double xwl = xw_lead[iz] +
                             z1 * (xw_lead[iz + 1] - xw_lead[iz]) * recip_z2;
                sum.data[0] += zdensity[j] * N_factor * xmom[j] * xwl;

                double xwt = xw_trail[iz] +
                             z1 * (xw_trail[iz + 1] - xw_trail[iz]) * recip_z2;
                sum.data[1] += zdensity[j] * N_factor * xwt;

                double ywl = yw_lead[iz] +
                             z1 * (yw_lead[iz + 1] - yw_lead[iz]) * recip_z2;
                sum.data[2] += zdensity[j] * N_factor * ymom[j] * ywl;

                double ywt = yw_trail[iz] +
                             z1 * (yw_trail[iz + 1] - yw_trail[iz]) * recip_z2;
                sum.data[3] += zdensity[j] * N_factor * ywt;

                double zw =
                    z_wake[iz] + z1 * (z_wake[iz + 1] - z_wake[iz]) * recip_z2;
                sum.data[4] += zdensity[j] * N_factor * zw;
            }
        }
    };

    struct alg_z_wake {
        Bunch_params bp;
        Wake_field wf;
        karray1d_dev zbins;
        karray1d_dev wakes;

        alg_z_wake(Bunch_params const& bp,
                   Wake_field const& wf,
                   karray1d_dev const& zbins,
                   karray1d_dev const& wakes)
            : bp(bp), wf(wf), zbins(zbins), wakes(wakes)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const team_member& member) const
        {
            // league size is the z_grid
            int z_grid = member.league_size();

            // li is the index i in zji
            int li = member.league_rank();
            int ti = member.team_rank();

            typedef kt::array_type<double, 5> value_t;
            typedef kt::SumArray<double, 5> array_sum_res_t;

            value_t sum;

            alg_z_wake_reduce z_wake_reduce{li, z_grid, bp, wf, zbins};

            Kokkos::parallel_reduce(Kokkos::TeamThreadRange(member, z_grid),
                                    z_wake_reduce,
                                    array_sum_res_t(sum));

            if (ti == 0) {
                wakes(z_grid * 0 + li) = sum.data[0];
                wakes(z_grid * 1 + li) = sum.data[1];
                wakes(z_grid * 2 + li) = sum.data[2];
                wakes(z_grid * 3 + li) = sum.data[3];
                wakes(z_grid * 4 + li) = sum.data[4];
            }
        }
    };

    KOKKOS_INLINE_FUNCTION
    void
    sum_over_bunch(double* sum,
                   int iz,
                   int size_wake,
                   double zji,
                   double xmean,
                   double ymean,
                   double realnum,
                   karray1d_dev const& wf)
    {
        double* z_coord = &wf(size_wake * 0);
        double* z_wake = &wf(size_wake * 1);
        double* xw_lead = &wf(size_wake * 2);
        double* xw_trail = &wf(size_wake * 3);
        double* yw_lead = &wf(size_wake * 4);
        double* yw_trail = &wf(size_wake * 5);

        if (iz + 1 < size_wake && iz > 0) {
            double z1 = zji - z_coord[iz];
            double recip_z2 = 1.0 / (z_coord[iz + 1] - z_coord[iz]);

            double xwl =
                xw_lead[iz] + z1 * (xw_lead[iz + 1] - xw_lead[iz]) * recip_z2;
            sum[0] += realnum * xmean * xwl;

            double xwt = xw_trail[iz] +
                         z1 * (xw_trail[iz + 1] - xw_trail[iz]) * recip_z2;
            sum[1] += realnum * xwt;

            double ywl =
                yw_lead[iz] + z1 * (yw_lead[iz + 1] - yw_lead[iz]) * recip_z2;
            sum[2] += realnum * ymean * ywl;

            double ywt = yw_trail[iz] +
                         z1 * (yw_trail[iz + 1] - yw_trail[iz]) * recip_z2;
            sum[3] += realnum * ywt;

            double zw =
                z_wake[iz] + z1 * (z_wake[iz + 1] - z_wake[iz]) * recip_z2;
            sum[4] += realnum * zw;
        }
    }

    struct alg_bunch_wake {
        Bunch_params bp;
        Wake_field wf;
        Bunch_props bps;
        karray1d_dev wakes;

        const int mean_bin;
        const double bunch_spacing;
        const double orbit_length;
        const bool full_machine;
        const int z_grid;

        alg_bunch_wake(Bunch_params const& bp,
                       Wake_field const& wf,
                       Bunch_props const& bps,
                       karray1d_dev const& wakes,
                       double bunch_spacing,
                       double orbit_length,
                       bool full_machine,
                       int z_grid)
            : bp(bp)
            , wf(wf)
            , bps(bps)
            , wakes(wakes)
            , mean_bin((bp.z_mean - bp.z_left) / bp.cell_size_z)
            , bunch_spacing(bunch_spacing)
            , orbit_length(orbit_length)
            , full_machine(full_machine)
            , z_grid(z_grid)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            double z_to_zmean = (mean_bin - i) * bp.cell_size_z;

            int num_bunches = bps.num_bunches;
            double sum[5] = {0, 0, 0, 0, 0};

            // current turn
            for (int j = 0; j < num_bunches; ++j) {
                // 0 is current turn
                int j_idx = bps.get_read_index(0, j);

                int j_bucket = bps.bucket_index[j_idx];
                if (j_bucket >= bp.bucket) continue;

                double zji = z_to_zmean +
                             bunch_spacing * (bp.bucket - j_bucket) +
                             (bps.zmean[j_idx] - bp.z_mean);

                // if (zji < z_coord[0]) continue;

                // below it is assumed the wake function is stored using a
                // quadratic grid
                int iz =
                    get_zindex_for_wake(zji, wf.delta_z, wf.istart, wf.zstart);

                // accumulate
                sum_over_bunch(sum,
                               iz,
                               wf.size_wake,
                               zji,
                               bps.xmean(j_idx),
                               bps.ymean(j_idx),
                               bps.realnum(j_idx),
                               wf.terms);
            }

            // full machine
            if (full_machine) {
                // TODO ...
            }

            // prev turn
            if (bps.registered_turns > 1) {
                for (int j = 0; j < num_bunches; ++j) {
                    // -1 is prev turn
                    int j_idx = bps.get_read_index(-1, j);

                    int j_bucket = bps.bucket_index[j_idx];
                    if (j_bucket < bp.bucket) continue;

                    double zji = z_to_zmean +
                                 bunch_spacing * (bp.bucket - j_bucket) +
                                 orbit_length + (bps.zmean(j_idx) - bp.z_mean);

                    // if (zji < z_coord[0]) continue;

                    // below it is assumed the wake function is stored using a
                    // quadratic grid
                    int iz = get_zindex_for_wake(
                        zji, wf.delta_z, wf.istart, wf.zstart);

                    // accumulate
                    sum_over_bunch(sum,
                                   iz,
                                   wf.size_wake,
                                   zji,
                                   bps.xmean(j_idx),
                                   bps.ymean(j_idx),
                                   bps.realnum(j_idx),
                                   wf.terms);
                }
            }

            wakes(z_grid * 0 + i) += sum[0];
            wakes(z_grid * 1 + i) += sum[1];
            wakes(z_grid * 2 + i) += sum[2];
            wakes(z_grid * 3 + i) += sum[3];
            wakes(z_grid * 4 + i) += sum[4];
        }
    };

    struct alg_turn_wake {
        typedef double value_type[];

        Bunch_params bp;
        Bunch_props bps;
        Wake_field wf;

        const double bunch_spacing;
        const double orbit_length;
        const bool full_machine;

        const int value_count = 5;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i, value_type sum) const
        {
            int turn = i + 1;
            int num_bunches = bps.num_bunches;
            double lsum[5] = {0, 0, 0, 0, 0};

            // current turn
            for (int j = 0; j < num_bunches; ++j) {
                int j_idx = bps.get_read_index(-turn, j);

                int j_bucket = bps.bucket_index[j_idx];
                if (turn == 1 && j_bucket >= bp.bucket) continue;

                double zji = bunch_spacing * (bp.bucket - j_bucket) +
                             orbit_length * turn +
                             (bps.zmean[j_idx] - bp.z_mean);

                // below it is assumed the wake function is stored using a
                // quadratic grid
                int iz =
                    get_zindex_for_wake(zji, wf.delta_z, wf.istart, wf.zstart);

                // accumulate
                sum_over_bunch(lsum,
                               iz,
                               wf.size_wake,
                               zji,
                               bps.xmean(j_idx),
                               bps.ymean(j_idx),
                               bps.realnum(j_idx),
                               wf.terms);
            }

            // full machine
            if (full_machine) {
                // TODO ...
            }

            sum[0] += lsum[0];
            sum[1] += lsum[1];
            sum[2] += lsum[2];
            sum[3] += lsum[3];
            sum[4] += lsum[4];
        }
    };

    struct alg_add_turn_wake {
        karray1d_dev wakes;
        karray1d_dev turn_wakes;
        const int z_grid;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            wakes(z_grid * 0 + i) += turn_wakes(0);
            wakes(z_grid * 1 + i) += turn_wakes(1);
            wakes(z_grid * 2 + i) += turn_wakes(2);
            wakes(z_grid * 3 + i) += turn_wakes(3);
            wakes(z_grid * 4 + i) += turn_wakes(4);
        }
    };

    struct alg_apply_kick {
        Particles p;
        ConstParticleMasks masks;

        karray1d_dev wf;
        const int z_grid;
        const double z_left;
        const double recip_h;
        const double wake_factor;

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            double* xw_lead = &wf(z_grid * 0);
            double* xw_trail = &wf(z_grid * 1);
            double* yw_lead = &wf(z_grid * 2);
            double* yw_trail = &wf(z_grid * 3);
            double* z_wake = &wf(z_grid * 4);

            if (masks(i)) {
                int bin = (p(i, 4) - z_left) * recip_h;
                if (bin < 0 || bin >= z_grid) return;

                double xkick = xw_lead[bin] + xw_trail[bin] * p(i, 0);
                double ykick = yw_lead[bin] + yw_trail[bin] * p(i, 2);
                double zkick = z_wake[bin];

                p(i, 1) += wake_factor * xkick;
                p(i, 3) += wake_factor * ykick;
                p(i, 5) += wake_factor * zkick;
            }
        }
    };

}

Impedance::Impedance(Impedance_options const& opts)
    : Collective_operator("impedance", 1.0)
    , opts(opts)
    , bunch_sim_id()
    , bps(1 /*num_bunches*/, opts.nstored_turns)
    , zbinning()
    , h_zbinning()
    , wakes()
    , h_wakes()
    , wake_field(opts.wake_file, opts.wake_type,
        opts.mwf_xlead, opts.mwf_xtrail, opts.mwf_ylead, opts.mwf_ytrail,
        opts.mwf_zwake)
{}

void
Impedance::apply_impl(Bunch_simulator& sim, double time_step, Logger& logger)
{
    if (sim[1].get_num_bunches())
        throw std::runtime_error(
            "Impedance cannot have bunches in secondary train");

    logger << "    Impedance\n";
    scoped_simple_timer timer("imp_total");

    // construct the workspace for a new bunch simulator
    if (bunch_sim_id != sim.id()) {
        construct_workspaces(sim);
        bunch_sim_id = sim.id();
    }

    // pre-work
    store_bunches_data(sim);

    // apply to bunches
    for (auto& train : sim.get_trains())
        for (auto& bunch : train.get_bunches())
            apply_bunch(bunch, time_step, logger);
}

void
Impedance::construct_workspaces(Bunch_simulator const& sim)
{
    zbinning = karray1d_dev("zbinning", opts.z_grid * 3);
    h_zbinning = Kokkos::create_mirror_view(zbinning);

    wakes = karray1d_dev("wakes", opts.z_grid * 5);
    h_wakes = Kokkos::create_mirror_view(wakes);

    int num_bunches = sim[0].get_num_bunches();
    if (num_bunches != bps.num_bunches)
        bps = Bunch_props(num_bunches, opts.nstored_turns);
}

void
Impedance::store_bunches_data(Bunch_simulator const& sim)
{
    scoped_simple_timer timer("imp_store_bunches_data");

    auto const& train = sim[0];
    auto num_bunches = train.get_num_bunches();
    auto num_local_bunches = train.get_num_local_bunches();

    // each bunch has 5 properties, x/y/z_mean, real_num, and bucket_idx
    karray1d_dev d_vbi_buf("vbi_buf", num_bunches * 5);
    karray1d_hst vbi_buf = Kokkos::create_mirror_view(d_vbi_buf);

    for (int i = 0; i < num_local_bunches; ++i) {
        auto const& bunch = train[i];
        auto means = Core_diagnostics::calculate_mean(bunch);

        int bucket_idx = bunch.get_bucket_index();
        int bunch_idx = bunch.get_bunch_index();

        if (opts.full_machine && (bucket_idx != bunch_idx))
            throw std::runtime_error(
                "for full_machine the buckets have to be occupied in order");
        if (bunch.get_comm().rank() == 0) {
            vbi_buf[bunch_idx * 5 + 0] = means[0];
            vbi_buf[bunch_idx * 5 + 1] = means[2];
            vbi_buf[bunch_idx * 5 + 2] = means[4];
            vbi_buf[bunch_idx * 5 + 3] = bunch.get_real_num();
            vbi_buf[bunch_idx * 5 + 4] = bucket_idx;
        }
    }

    int error = MPI_Allreduce(MPI_IN_PLACE,
                              (void*)vbi_buf.data(),
                              num_bunches * 5,
                              MPI_DOUBLE,
                              MPI_SUM,
                              train.get_comm());

    if (error != MPI_SUCCESS)
        throw std::runtime_error(
            "Impedance::store_bunches_data: MPI error in MPI_Allreduce");

    // copy the vbi_buf to device memory
    Kokkos::deep_copy(d_vbi_buf, vbi_buf);

    // copy the buffer to bps
    alg_write_bps write_bps{bps, d_vbi_buf};
    Kokkos::parallel_for(num_bunches, write_bps);

    // increment the registered turns in bps
    bps.increment_registered_turns();
}

void
Impedance::apply_bunch(Bunch& bunch, double time_step, Logger& logger)
{
    bunch.convert_to_fixed_t_lab();

    auto bp = calculate_moments_and_partitions(bunch);

    auto means = Core_diagnostics::calculate_mean(bunch);
    bp.z_mean = means[4];
    bp.N_factor = bunch.get_real_num() / bunch.get_total_num();
    bp.bucket = bunch.get_bucket_index();

    calculate_kicks(bunch, bp);

    /// N.B. the wakefiled file reads W/(Z_0*L), Z_0=1/(epsilon_0*c)
    double wake_factor = -4. * Kokkos::numbers::pi_v<double> * pconstants::rp;

    double gamma = bunch.get_reference_particle().get_gamma();
    double beta = bunch.get_reference_particle().get_beta();
    double w_f = wake_factor * time_step / (gamma * beta);

    apply_impedance_kick(bunch, bp, w_f);

    bunch.convert_to_fixed_z_lab();
}

Bunch_params
Impedance::calculate_moments_and_partitions(Bunch const& bunch)
{
    scoped_simple_timer timer("imp_moments_and_partitions");

    // output cell_size_z, xmom, ymom, zdensity
    Bunch_params bp;

    auto bunchmin = Core_diagnostics::calculate_min(bunch);
    bp.z_left = bunchmin[2];

    auto bunchmax = Core_diagnostics::calculate_max(bunch);
    double z_length = bunchmax[2] - bp.z_left;

    if (z_length <= 1.e-14) throw std::runtime_error("z_length too small ");

    // 1e-14 is to make sure the max-z particle falls in the last bin
    // bp.cell_size_z = z_length / double(opts.z_grid) + 1e-14;
    bp.cell_size_z = z_length / double(opts.z_grid);

    // double h = z_length/(opts.z_grid-1.0); // AM why have I done that???
    double h = bp.cell_size_z;

    // get binning results
    auto parts = bunch.get_local_particles();
    auto masks = bunch.get_local_particle_masks();

#if 0
    alg_z_binning alg(parts, masks, opts.z_grid, bp.z_left, h);
    Kokkos::parallel_reduce(bunch.size(), alg, h_zbinning.data());
    Kokkos::fence();
#endif

    // zero first
    zero_karray(zbinning);

    // z binning
    scatter_t scatter(zbinning);
    alg_z_binning_sv alg{
        parts, masks, scatter, opts.z_grid, bp.z_left, 1.0 / h};

    Kokkos::parallel_for(bunch.size(), alg);
    Kokkos::Experimental::contribute(zbinning, scatter);
    Kokkos::fence();

    // MPI reduction to get global z-binning results
    if (bunch.get_comm().size() > 1) {
        // copy to host
        Kokkos::deep_copy(h_zbinning, zbinning);

        // all reduce
        int error = MPI_Allreduce(MPI_IN_PLACE,
                                  h_zbinning.data(),
                                  opts.z_grid * 3,
                                  MPI_DOUBLE,
                                  MPI_SUM,
                                  bunch.get_comm());

        if (error != MPI_SUCCESS)
            throw std::runtime_error("MPI error in Impedance reduce z_binning");

        // copy back to device
        Kokkos::deep_copy(zbinning, h_zbinning);
    }

    // normalize
    alg_z_normalize alg2{zbinning, opts.z_grid};
    Kokkos::parallel_for(opts.z_grid, alg2);

#if 0
    Logger l;
    kt::print_arr_sum(l, zbinning, 0, opts.z_grid);
    kt::print_arr_sum(l, zbinning, opts.z_grid*1, opts.z_grid);
    kt::print_arr_sum(l, zbinning, opts.z_grid*2, opts.z_grid);
#endif

    return bp;
}

void
Impedance::calculate_kicks(Bunch const& bunch, Bunch_params const& bp)
{
    scoped_simple_timer timer("imp_calcualte_kicks");

    int num_trains = 0;
    int mean_bin = (int)((bp.z_mean - bp.z_left) / bp.cell_size_z);

    if (bps.registered_turns == 0)
        throw std::runtime_error(
            "registered_turns size cannot be zero, "
            "probably you propagate a bunch instead of a bunch_train");

    if (mean_bin < 0 || mean_bin >= opts.z_grid)
        throw std::runtime_error(
            "impedance: the index bin of beam min cannot be <0 or >z_grid, "
            "something is wrong");

    if (opts.full_machine) {
        /// num_trains is relevant only when the full machine option is
        /// considered a full machine consideres a num_train of bunches repeats
        /// with modulation wave wn[] all buckets are full, but only numbunches
        /// bunches properties are stored exemple: full_machine, all bunches
        /// identical, no wave across the machine: num_trains=num_buckets,
        /// wn=[0,0,0], it's a one bunch simulation example: full_machine, two
        /// bunch simulation, num_trains= num_buckets/2
        num_trains = int(opts.num_buckets / bps.num_bunches);

        if (std::abs(opts.num_buckets / float(bps.num_bunches) - num_trains) >
            1e-8)
            throw std::runtime_error(
                "full machine assumes repetitive numer of trains: "
                "num_buckets should be divisible to numbunches");

        if (opts.wn[0] < 0 || opts.wn[0] >= num_trains || opts.wn[1] < 0 ||
            opts.wn[1] >= num_trains || opts.wn[2] < 0 ||
            opts.wn[2] >= num_trains)
            throw std::runtime_error(
                "full machine wave number cannot be smaller than zero "
                "or larger than num_trains-1");
    }

    using Kokkos::TeamPolicy;
    using Kokkos::TeamThreadRange;

    // in-bunch z wake
    // zbinning: zdensity, xmom, ymom
    // wakes: xw_lead, xw_trail, yw_lead, yw_trail, zwake
    alg_z_wake ft_z_wake{bp, wake_field, zbinning, wakes};

    const int team_size_max =
        team_policy(opts.z_grid, 1)
            .team_size_max(ft_z_wake, Kokkos::ParallelForTag());

    Kokkos::parallel_for(TeamPolicy<>(opts.z_grid, team_size_max), ft_z_wake);

    // bunch-bunch wake
    // at the moment bucket 0 is in front of bucket 1,
    // which is in front of bucket 2, etc...
    alg_bunch_wake ft_bunch_wake(bp,
                                 wake_field,
                                 bps,
                                 wakes,
                                 opts.bunch_spacing,
                                 opts.orbit_length,
                                 opts.full_machine,
                                 opts.z_grid);
    Kokkos::parallel_for(opts.z_grid, ft_bunch_wake);

    // turn-turn wake
    if (bps.registered_turns > 1) {
        // calculate turn-turn wakes
        karray1d_dev turn_wakes("turn_wakes", 5);
        alg_turn_wake ft_turn_wake{bp,
                                   bps,
                                   wake_field,
                                   opts.bunch_spacing,
                                   opts.orbit_length,
                                   opts.full_machine};
        Kokkos::parallel_reduce(
            bps.registered_turns - 1, ft_turn_wake, turn_wakes);

        // add the turn-turn wakes to the final wakes
        alg_add_turn_wake ft_add_turn_wake{wakes, turn_wakes, opts.z_grid};
        Kokkos::parallel_for(opts.z_grid, ft_add_turn_wake);
    }

#if 0
    // prints
    Logger l;
    kt::print_arr_sum(l, wakes, 0, opts.z_grid);
    kt::print_arr_sum(l, wakes, opts.z_grid*1, opts.z_grid);
    kt::print_arr_sum(l, wakes, opts.z_grid*2, opts.z_grid);
    kt::print_arr_sum(l, wakes, opts.z_grid*3, opts.z_grid);
    kt::print_arr_sum(l, wakes, opts.z_grid*4, opts.z_grid);
#endif
}

void
Impedance::apply_impedance_kick(Bunch& bunch,
                                Bunch_params const& bp,
                                double wake_factor)
{
    scoped_simple_timer timer("imp_apply_kick");

    alg_apply_kick alg{bunch.get_local_particles(),
                       bunch.get_local_particle_masks(),
                       wakes,
                       opts.z_grid,
                       bp.z_left,
                       1.0 / bp.cell_size_z,
                       wake_factor};

    Kokkos::parallel_for(bunch.size(), alg);

#if 0
    Logger l(0, LoggerV::DEBUG);
    bunch.print_statistics(l);
#endif
}
