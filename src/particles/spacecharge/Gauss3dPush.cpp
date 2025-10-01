/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Ji Qiang
 * License: BSD-3-Clause-LBNL
 */
#include "Gauss3dPush.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>       // for Real
#include "diagnostics/ReducedBeamCharacteristics.H"


namespace impactx::particles::spacecharge
{
    /** Compute space-charge fields from a 3D Gaussian distribution.
     *
     * Compute integrals Eqs 45-47 used in the space-charge fields from a 3D Gaussian distribution.
     * "A two-and-a-half dimensional symplectic space-charge solver", LBNL Report Number: LBNL-2001674, 2025.
     * Input particle locations (x,y,z) and RMS sizes (sigx,sigy,sigz) and return the integrals for SC fields.
     *
     * @todo This function can be vectorized with AMREX_SIMD for better CPU performance.
     */
    AMREX_GPU_DEVICE
    void efldgauss (
        int nint,
        amrex::ParticleReal const xin,
        amrex::ParticleReal const yin,
        amrex::ParticleReal const zin,
        amrex::ParticleReal const sigx,
        amrex::ParticleReal const sigy,
        amrex::ParticleReal const sigz,
        amrex::ParticleReal const gamma,
        amrex::ParticleReal & pintex,
        amrex::ParticleReal & pintey,
        amrex::ParticleReal & pintez
    )
    {
        using namespace amrex::literals;

        // enforce an odd number of grid points
        if (nint % 2 == 0) nint += 1;

        amrex::ParticleReal const xmin = 0, xmax = 1;
        amrex::ParticleReal const h = (xmax - xmin) / (nint - 1);
        amrex::ParticleReal const xp = xin / sigx;
        amrex::ParticleReal const yp = yin / sigy;
        amrex::ParticleReal const zp = zin / sigz;
        amrex::ParticleReal const asp = sigx / sigy;
        amrex::ParticleReal const basp = sigx / (sigz * gamma);
        amrex::ParticleReal const xp2 = xp * xp;
        amrex::ParticleReal const yp2 = yp * yp;
        amrex::ParticleReal const zp2 = zp * zp;

        // integration results
        amrex::ParticleReal sum0ex = 0, sum0ey = 0, sum0ez = 0;

        amrex::ParticleReal fx = 0.0;
        amrex::ParticleReal fy = 0.0;
        amrex::ParticleReal fz = 0.0;

        for (int i = 0; i < nint; ++i)
        {
            amrex::ParticleReal const x = xmin + i * h;
            amrex::ParticleReal const t2 = asp * asp + (1_prt - asp * asp) * x * x;
            amrex::ParticleReal const t2sqrt = std::sqrt(t2);
            amrex::ParticleReal const bt2 = basp * basp + (1_prt - basp * basp) * x * x;
            amrex::ParticleReal const bt2sqrt = std::sqrt(bt2);
            amrex::ParticleReal const f2x = t2sqrt * bt2sqrt;
            amrex::ParticleReal const f2y = t2 * t2sqrt * bt2sqrt;
            amrex::ParticleReal const f2z = bt2 * t2sqrt * bt2sqrt;
            amrex::ParticleReal const f1 = x * x * std::exp(-x * x * (xp2 + yp2 / t2 + zp2 / bt2) / 2_prt);

            // Simpson's rule integration

            // Ex
            //amrex::ParticleReal const fx = f1 / f2x;
            fx = f1 / f2x;
            if (i%2 == 0)
                sum0ex += 2_prt * fx * h;
            else
                sum0ex += 4_prt * fx * h;

            // Ey
            //amrex::ParticleReal const fy = f1 / f2y;
            fy = f1 / f2y;
            if (i%2 == 0)
                sum0ey += 2_prt * fy * h;
            else
                sum0ey += 4_prt * fy * h;

            // Ez
            //amrex::ParticleReal const fz = f1 / f2z;
            fz = f1 / f2z;
            if (i%2 == 0)
                sum0ez += 2_prt * fz * h;
            else
                sum0ez += 4_prt * fz * h;
        }

        //end point correction
        sum0ex -= fx*h;
        sum0ey -= fy*h;
        sum0ez -= fz*h;

        pintex = sum0ex / 3_prt;
        pintey = sum0ey / 3_prt;
        pintez = sum0ez / 3_prt;
    }

    void Gauss3dPush (
        ImpactXParticleContainer & pc,
        amrex::ParticleReal const slice_ds
    )
    {
        BL_PROFILE("impactx::spacecharge::Gauss3dPush");

        using namespace amrex::literals;

        amrex::ParticleReal const charge = pc.GetRefParticle().charge;

        // get RMS sigma of beam
        // Standard deviation of the particle positions (speed of light times time delay for t) (unit: meter)
        std::unordered_map<std::string, amrex::ParticleReal> const rbc = diagnostics::reduced_beam_characteristics(pc);
        amrex::ParticleReal const sigx = rbc.at("sig_x");
        amrex::ParticleReal const sigy = rbc.at("sig_y");
        amrex::ParticleReal const sigz = rbc.at("sig_t");
        amrex::ParticleReal const bchchg = rbc.at("charge_C");

        // physical constants and reference quantities
        amrex::ParticleReal const c0_SI = 2.99792458e8;  // TODO move out
        amrex::ParticleReal const mc_SI = pc.GetRefParticle().mass * c0_SI;
        amrex::ParticleReal const pz_ref_SI = pc.GetRefParticle().beta_gamma() * mc_SI;
        amrex::ParticleReal const gamma = pc.GetRefParticle().gamma();
        amrex::ParticleReal const inv_gamma2 = 1.0_prt / (gamma * gamma);
        amrex::ParticleReal const rfpiepslon = c0_SI*c0_SI*1.0e-7;

        amrex::ParticleReal const dt = slice_ds / pc.GetRefParticle().beta() / c0_SI;

        // group together constants for the momentum
        using ablastr::constant::math::pi;
        amrex::ParticleReal const asp = sigx/sigy;
        amrex::ParticleReal const facx = 1.0/(sigx * sigx * sigz * std::sqrt(2_prt * pi));
        amrex::ParticleReal const facy = 1.0/(sigy * sigy * sigz * std::sqrt(2_prt * pi));
        amrex::ParticleReal const facz = 1.0/(sigz * sigz * sigz * std::sqrt(2_prt * pi));
        amrex::ParticleReal const push_consts = dt * charge * inv_gamma2 / pz_ref_SI
            * 2_prt * asp * bchchg * rfpiepslon;

        // loop over refinement levels
        int const nLevel = pc.finestLevel();
        for (int lev = 0; lev <= nLevel; ++lev)
        {
            // loop over all particle boxes
            using ParIt = ImpactXParticleContainer::iterator;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (ParIt pti(pc, lev); pti.isValid(); ++pti) {
                const int np = pti.numParticles();

                // preparing access to particle data: SoA of Reals
                auto& soa_real = pti.GetStructOfArrays().GetRealData();
                amrex::ParticleReal const * const AMREX_RESTRICT part_x = soa_real[RealSoA::x].dataPtr();
                amrex::ParticleReal const * const AMREX_RESTRICT part_y = soa_real[RealSoA::y].dataPtr();
                amrex::ParticleReal const * const AMREX_RESTRICT part_z = soa_real[RealSoA::z].dataPtr(); // note: currently for a fixed t
                amrex::ParticleReal* const AMREX_RESTRICT part_px = soa_real[RealSoA::px].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_py = soa_real[RealSoA::py].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT part_pz = soa_real[RealSoA::pz].dataPtr(); // note: currently for a fixed t

                // gather to each particle and push momentum
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) {
                    // access SoA Real data
                    amrex::ParticleReal const & AMREX_RESTRICT x = part_x[i];
                    amrex::ParticleReal const & AMREX_RESTRICT y = part_y[i];
                    amrex::ParticleReal const & AMREX_RESTRICT z = part_z[i];
                    amrex::ParticleReal & AMREX_RESTRICT px = part_px[i];
                    amrex::ParticleReal & AMREX_RESTRICT py = part_py[i];
                    amrex::ParticleReal & AMREX_RESTRICT pz = part_pz[i];

                    // field integrals from a 3D Gaussian bunch
                    int const nint = 101;  // TODO: should "nint" be user-configurable? Otherwise make it constexpr in efldgauss

                    amrex::ParticleReal eintx, einty, eintz;
                    efldgauss(nint,x,y,z,sigx,sigy,sigz,gamma,eintx,einty,eintz);

                    // push momentum
                    px += x * eintx * push_consts * facx;
                    py += y * einty * push_consts * facy;
                    pz += z * eintz * push_consts * facz;

                    // push position is done in the lattice elements
                });


            } // end loop over all particle boxes
        } // env mesh-refinement level loop
    }

}  // namespace impactx::particles::spacecharge
