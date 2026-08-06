/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef HANDLE_WAKEFIELD_H
#define HANDLE_WAKEFIELD_H

#include "particles/ImpactXParticleContainer.H"
#include "elements/mixin/named.H"
#include "ChargeBinning.H"
#include "CSRBendElement.H"
#include "CSRKickModel.H"
#include "WakeConvolution.H"
#include "WakePush.H"

#include <AMReX_BLassert.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuParallelReduce.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <exception>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>


namespace impactx::particles::wakefields
{
    /** Function to handle CSR bending process including charge deposition,
     * wakefield generation (built-in analytic model or user-provided kick
     * model), and the momentum kick.
     *
     * Returns immediately if slice_ds is zero.
     *
     * @param[in] particle_container the particle species container
     * @param[in] element_variant variant type of the lattice element
     * @param[in] slice_ds slice spacing along s
     * @param[in] slice current slice index in the element (0-based, forward orientation)
     * @param[in] nslice number of slices in the element
     * @param[in] csr_kick_model optional user-provided CSR kick model;
     *            empty selects the built-in analytic model
     * @param[in] print_wakefield for debugging: print the wakefield to convolved_wakefield.txt
     */
    template <typename T_Element>
    void HandleWakefield (
        impactx::ImpactXParticleContainer& particle_container,
        T_Element const& element_variant,
        amrex::Real slice_ds,
        int slice,
        int nslice,
        CSRKickModel const & csr_kick_model = {},
        bool print_wakefield = false
    )
    {
        if (slice_ds == 0.0) { return; }

        BL_PROFILE("impactx::particles::wakefields::HandleWakefield")

        amrex::ParmParse pp_algo("algo");
        bool csr = false;
        pp_algo.queryAdd("csr", csr);

        // Call the CSR bend function
        auto const csr_bend = impactx::particles::wakefields::CSRBendElement(
            element_variant, particle_container.GetRefParticle());
        bool const element_has_csr = std::get<0>(csr_bend);
        amrex::Real const R = std::get<1>(csr_bend);
        amrex::Real const signed_R = std::get<2>(csr_bend);

        // Enter loop if lattice has bend element
        if (csr && element_has_csr)
        {
            bool const use_kick_model = static_cast<bool>(csr_kick_model);
#ifndef ImpactX_USE_FFT
            if (!use_kick_model) {
                throw std::runtime_error(
                    "algo.csr was requested but ImpactX was not compiled with FFT support. "
                    "Recompile with ImpactX_FFT=ON or provide a user-defined CSR kick model "
                    "(csr_kick_model)."
                );
            }
#endif

            int csr_bins = 150;
            pp_algo.queryAddWithParser("csr_bins", csr_bins);

            // Measure beam size, extract the min, max of particle positions
            [[maybe_unused]] auto const [x_min, y_min, t_min, x_max, y_max, t_max] =
                particle_container.MinAndMaxPositions();

            // Ensure we have a beam with a finite extent in t
            if (t_min == t_max)
            {
                auto const npart = particle_container.TotalNumberOfParticles();
                std::string error_msg = "Coherent Synchrotron Radiation (CSR) "
                                        "modeling was enabled, but ";
                if (npart == 1) {
                    error_msg.append("there is only one particle in the beam. ");
                } else {
                    error_msg.append("the beam is flat in t. ");
                }
                error_msg.append(
                    "CSR modeling is a collective effect and requires "
                    "a particle beam with a finite extent in t. "
                    "Please set the CSR option to false."
                );
                throw std::runtime_error(error_msg);
            }

            // Set parameters for charge deposition
            bool const is_unity_particle_weight = false; // Only true if w = 1
            bool const GetNumberDensity = true;

            int const num_bins = csr_bins;  // Set resolution
            amrex::Real const bin_min = t_min;
            amrex::Real const bin_max = t_max;
            amrex::Real const bin_size = (bin_max - bin_min) / (num_bins - 1);  // number of evaluation points

            // Allocate memory for the charge profile
            amrex::Gpu::DeviceVector<amrex::Real> charge_distribution(num_bins + 1, 0.0);

            // Call charge deposition function
            impactx::particles::wakefields::DepositCharge1D(particle_container, charge_distribution, bin_min, bin_size, is_unity_particle_weight);

            // Sum up all partial charge histograms to one MPI process to calculate the global wakefield.
            // Once calculated, we will distribute the kick per bin back to every MPI process.
            amrex::ParallelReduce::Sum(
                charge_distribution,
                amrex::ParallelDescriptor::IOProcessorNumber(),
                amrex::ParallelDescriptor::Communicator()
            );

            // Per-bin kick forces [N]. wake_pt is pre-sized on every rank to
            // the (global) wake length. On the IO processor it is overwritten
            // by the built-in convolution result or the kick model output of
            // the same length. The transverse components are optional and
            // only allocated if a user-provided kick model returns them.
            amrex::Gpu::DeviceVector<amrex::Real> wake_pt(num_bins, 0.0);
            amrex::Gpu::DeviceVector<amrex::Real> wake_px;
            amrex::Gpu::DeviceVector<amrex::Real> wake_py;

            if (use_kick_model)
            {
                // Per-bin weighted transverse sums, packed as
                // [sum_w | sum_w*x | sum_w*y] for a single MPI reduction.
                amrex::Gpu::DeviceVector<amrex::Real> packed_sums(3 * (num_bins + 1), 0.0);
                impactx::particles::wakefields::WeightedTransverseSums(
                    particle_container, packed_sums, bin_min, bin_size, is_unity_particle_weight);
                amrex::ParallelReduce::Sum(
                    packed_sums,
                    amrex::ParallelDescriptor::IOProcessorNumber(),
                    amrex::ParallelDescriptor::Communicator()
                );

                int const io_rank = amrex::ParallelDescriptor::IOProcessorNumber();
                int const nprocs = amrex::ParallelDescriptor::NProcs();
                int has_transverse = 0;  // bit 0: px is present, bit 1: py is present
                int model_succeeded = 1;
                std::string model_error;
                if (amrex::ParallelDescriptor::IOProcessor())
                {
                    auto invoke_model = [&] ()
                    {
                        // Assemble the binned beam profile
                        CSRProfile profile;
                        profile.charge = std::move(charge_distribution);
                        profile.mean_x.resize(num_bins + 1);
                        profile.mean_y.resize(num_bins + 1);
                        profile.bin_min = bin_min;
                        profile.bin_size = bin_size;
                        profile.num_bins = num_bins;

                        int const nb1 = num_bins + 1;
                        amrex::Real const * const dptr_sums = packed_sums.data();
                        amrex::Real * const dptr_mean_x = profile.mean_x.data();
                        amrex::Real * const dptr_mean_y = profile.mean_y.data();
                        amrex::ParallelFor(nb1, [=] AMREX_GPU_DEVICE (int i)
                        {
                            amrex::Real const sum_w = dptr_sums[i];
                            dptr_mean_x[i] = sum_w > 0
                                ? dptr_sums[nb1 + i] / sum_w
                                : amrex::Real(0.0);
                            dptr_mean_y[i] = sum_w > 0
                                ? dptr_sums[2 * nb1 + i] / sum_w
                                : amrex::Real(0.0);
                        });

                        // Assemble the element context
                        CSRElementContext context;
                        std::visit([&context] (auto const & element)
                        {
                            using T = std::decay_t<decltype(element)>;
                            context.element_type = T::type;
                            if constexpr (std::is_base_of_v<elements::mixin::Named, T>) {
                                if (element.has_name()) {
                                    context.element_name = element.name();
                                }
                            }
                            context.ds = element.ds();
                        }, element_variant);
                        context.rc = R;
                        context.signed_rc = signed_R;
                        context.nslice = nslice;
                        context.slice = slice;
                        context.s = static_cast<amrex::ParticleReal>(slice) * slice_ds;
                        context.slice_ds = slice_ds;
                        context.ref = particle_container.GetRefParticle();

                        // Ensure all deposition results are visible to the model
                        amrex::Gpu::streamSynchronize();

                        // Call the user-provided kick model
                        CSRKick kick = csr_kick_model(profile, context);

                        if (int(kick.pt.size()) != num_bins ||
                            !(kick.px.empty() || int(kick.px.size()) == num_bins) ||
                            !(kick.py.empty() || int(kick.py.size()) == num_bins))
                        {
                            throw std::runtime_error(
                                "HandleWakefield: the user-provided CSR kick model returned "
                                "arrays of wrong length for element type '" +
                                context.element_type + "' (name: '" + context.element_name +
                                "', slice " + std::to_string(slice) + "): expected " +
                                std::to_string(num_bins) + " values per component."
                            );
                        }

                        wake_pt = std::move(kick.pt);
                        wake_px = std::move(kick.px);
                        wake_py = std::move(kick.py);
                        has_transverse =
                            (wake_px.empty() ? 0 : 1) | (wake_py.empty() ? 0 : 2);
                    };

                    if (nprocs == 1)
                    {
                        // Preserve the original Python exception type and traceback in serial.
                        invoke_model();
                    }
                    else
                    {
                        try
                        {
                            invoke_model();
                        }
                        catch (std::exception const & error)
                        {
                            model_succeeded = 0;
                            model_error = error.what();
                        }
                        catch (...)
                        {
                            model_succeeded = 0;
                            model_error = "The CSR kick model raised an unknown exception.";
                        }
                    }
                }

                if (nprocs > 1)
                {
                    amrex::ParallelDescriptor::Bcast(&model_succeeded, 1, io_rank);
                    if (!model_succeeded)
                    {
                        int error_size = static_cast<int>(model_error.size());
                        amrex::ParallelDescriptor::Bcast(&error_size, 1, io_rank);
                        model_error.resize(error_size);
                        if (error_size > 0) {
                            amrex::ParallelDescriptor::Bcast(
                                model_error.data(), error_size, io_rank);
                        }
                        throw std::runtime_error(model_error);
                    }
                }

                // Broadcast which kick components are present and pre-size the
                // receiving ranks accordingly (the vectors are broadcast below).
                amrex::ParallelDescriptor::Bcast(&has_transverse, 1, io_rank);
                if (!amrex::ParallelDescriptor::IOProcessor())
                {
                    if (has_transverse & 1) { wake_px.resize(num_bins, 0.0); }
                    if (has_transverse & 2) { wake_py.resize(num_bins, 0.0); }
                }
            }
            else if (amrex::ParallelDescriptor::IOProcessor())
            {
                // Built-in analytic CSR model:
                // dN/ds convolved with the steady-state CSR wake function.

                // Call charge density derivative function
                amrex::Gpu::DeviceVector<amrex::Real> slopes(charge_distribution.size() - 1, 0.0);
                impactx::particles::wakefields::DerivativeCharge1D(charge_distribution, slopes, bin_size,
                                                                   GetNumberDensity); // Use number derivatives for convolution with CSR

                // Construct CSR wake function on 2N support
                amrex::Gpu::DeviceVector<amrex::Real> wake_function(num_bins*2, 0.0);
                amrex::Real *const dptr_wake_function = wake_function.data();
                auto const dR = R;  // for NVCC capture
                amrex::ParallelFor(num_bins*2, [=] AMREX_GPU_DEVICE(int i) {
                    if (i < num_bins) {
                        amrex::Real const s = static_cast<amrex::Real>(i) * bin_size;
                        dptr_wake_function[i] = impactx::particles::wakefields::w_l_csr(s, dR, bin_size);
                    }
                    else if (i > num_bins) {
                        amrex::Real const s = static_cast<amrex::Real>(i - 2*num_bins) * bin_size;
                        dptr_wake_function[i] = impactx::particles::wakefields::w_l_csr(s, dR, bin_size);
                    }
                });

                // Call convolution function
                wake_pt = impactx::particles::wakefields::convolve_fft(slopes, wake_function, bin_size);
            }

            // Broadcast the global kick per bin from the IO processor (which
            // alone computed it) to every MPI rank. Receivers are pre-sized above.
            amrex::ParallelDescriptor::Bcast(
                wake_pt,
                amrex::ParallelDescriptor::IOProcessorNumber(),
                amrex::ParallelDescriptor::Communicator()
            );
            if (!wake_px.empty())
            {
                amrex::ParallelDescriptor::Bcast(
                    wake_px,
                    amrex::ParallelDescriptor::IOProcessorNumber(),
                    amrex::ParallelDescriptor::Communicator()
                );
            }
            if (!wake_py.empty())
            {
                amrex::ParallelDescriptor::Bcast(
                    wake_py,
                    amrex::ParallelDescriptor::IOProcessorNumber(),
                    amrex::ParallelDescriptor::Communicator()
                );
            }

            // Check convolution output
            if (print_wakefield && amrex::ParallelDescriptor::IOProcessor())
            {
                // wake_pt is in (non-managed) device memory. Copy it
                // to the host before iterating it on the host for printing.
                amrex::Gpu::HostVector<amrex::Real> h_wakefield(wake_pt.size());
                amrex::Gpu::dtoh_memcpy(h_wakefield.data(), wake_pt.data(),
                                        wake_pt.size() * sizeof(amrex::Real));
                std::cout << "Convolved wakefield: ";
                std::ofstream outfile("convolved_wakefield.txt");
                for (amrex::Real const val : h_wakefield) {
                    std::cout << val << " ";
                    outfile << val << "\n";
                }
                std::cout << "\n";
                outfile.close();
            }

            // Call function to kick particles with wake
            impactx::particles::wakefields::WakePush(particle_container, wake_pt, wake_px, wake_py, slice_ds, bin_size, t_min);
        }
    }

} // namespace impactx::particles::wakefields

#endif // HANDLE_WAKEFIELD_H
