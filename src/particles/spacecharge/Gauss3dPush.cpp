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
#include <AMReX_SPACE.H>      // for AMREX_D_DECL
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

        // TODO: GPU support, fuse loops (and avoid large, temporary arrays)
        std::vector<amrex::ParticleReal> xvec(nint), t2vec(nint), t2sqrt(nint), bt2vec(nint), bt2sqrt(nint);
        std::vector<amrex::ParticleReal> f2xvec(nint), f2yvec(nint), f2zvec(nint), f1vec(nint), fxvec(nint);

        for (int i = 0; i < nint; ++i)
            xvec[i] = xmin + i * h;

        for (int i = 0; i < nint; ++i)
            t2vec[i] = asp * asp + (1_prt - asp * asp) * xvec[i] * xvec[i];

        for (int i = 0; i < nint; ++i)
            t2sqrt[i] = std::sqrt(t2vec[i]);

        for (int i = 0; i < nint; ++i)
            bt2vec[i] = basp * basp + (1_prt - basp * basp) * xvec[i] * xvec[i];

        for (int i = 0; i < nint; ++i)
            bt2sqrt[i] = std::sqrt(bt2vec[i]);

        for (int i = 0; i < nint; ++i)
            f2xvec[i] = t2sqrt[i] * bt2sqrt[i];

        for (int i = 0; i < nint; ++i)
            f2yvec[i] = t2vec[i] * t2sqrt[i] * bt2sqrt[i];

        for (int i = 0; i < nint; ++i)
            f2zvec[i] = bt2vec[i] * t2sqrt[i] * bt2sqrt[i];

        amrex::ParticleReal xp2 = xp * xp, yp2 = yp * yp, zp2 = zp * zp;
        for (int i = 0; i < nint; ++i)
            f1vec[i] = xvec[i] * xvec[i] * std::exp(-xvec[i] * xvec[i] * (xp2 + yp2 / t2vec[i] + zp2 / bt2vec[i]) / 2_prt);

        // Simpson's rule integration
        amrex::ParticleReal sum0ex = 0, sum0ey = 0, sum0ez = 0;

        // Ex
        for (int i = 0; i < nint; ++i)
            fxvec[i] = f1vec[i] / f2xvec[i];
        for (int i = 0; i < nint; i += 2)
            sum0ex += 2_prt * fxvec[i] * h;
        for (int i = 1; i < nint; i += 2)
            sum0ex += 4_prt * fxvec[i] * h;
        pintex = sum0ex / 3_prt;

        // Ey
        for (int i = 0; i < nint; ++i)
            fxvec[i] = f1vec[i] / f2yvec[i];
        for (int i = 0; i < nint; i += 2)
            sum0ey += 2_prt * fxvec[i] * h;
        for (int i = 1; i < nint; i += 2)
            sum0ey += 4_prt * fxvec[i] * h;
        pintey = sum0ey / 3_prt;

        // Ez
        for (int i = 0; i < nint; ++i)
            fxvec[i] = f1vec[i] / f2zvec[i];
        for (int i = 0; i < nint; i += 2)
            sum0ez += 2_prt * fxvec[i] * h;
        for (int i = 1; i < nint; i += 2)
            sum0ez += 4_prt * fxvec[i] * h;
        pintez = sum0ez / 3_prt;
    }

    void Gauss3dPush (
        ImpactXParticleContainer & pc,
        amrex::ParticleReal const slice_ds
    )
    {
        BL_PROFILE("impactx::spacecharge::GatherAndPush");

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

        // group together constants for the momentum push
        amrex::ParticleReal const push_consts = dt * charge * inv_gamma2 / pz_ref_SI;

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
                    int const nint = 401;
                    amrex::ParticleReal eintx, einty, eintz;
                    efldgauss(nint,x,y,z,sigx,sigy,sigz,gamma,eintx,einty,eintz);

                    amrex::ParticleReal const asp = sigx/sigy;

                    // push momentum
                    using ablastr::constant::math::pi;
                    px += x*eintx*2_prt*asp/sigx/sigx/sigz/sqrt(2_prt*pi)*bchchg*rfpiepslon * push_consts;
                    py += y*einty*2_prt*asp/sigy/sigy/sigz/sqrt(2_prt*pi)*bchchg*rfpiepslon * push_consts;
                    pz += z*eintz*2_prt*asp/sigz/sigz/sigz/sqrt(2_prt*pi)*bchchg*rfpiepslon * push_consts;

                    // push position is done in the lattice elements
                });


            } // end loop over all particle boxes
        } // env mesh-refinement level loop
    }

}  // namespace impactx::particles::spacecharge
