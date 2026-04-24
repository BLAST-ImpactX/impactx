/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_MULTIPOLE_H
#define IMPACTX_MULTIPOLE_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/thin.H"
#include "mixin/lineartransport.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>
#include <stdexcept>


namespace impactx::elements
{
    struct Multipole
    : public mixin::Named,
      public mixin::BeamOptic<Multipole>,
      public mixin::LinearTransport<Multipole>,
      public mixin::Thin,
      public mixin::Alignment,
      public mixin::NoFinalize
#ifndef _WIN32
      // https://github.com/AMReX-Codes/amrex/issues/4609
      , public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
#endif
    {
        static constexpr auto type = "Multipole";
        using PType = ImpactXParticleContainer::ParticleType;

        /** compute factorial of multipole index
         *
         * Needs to be called whenever m_multipole is changed.
         */
        void
        compute_factorial ()
        {
            int const m = m_multipole - 1;
            m_mfactorial = 1;
            for( int n = 1; n < m + 1; n = n + 1 ) {
                m_mfactorial *= n;
            }
        }

        /** A general thin multipole element
         *
         * @param multipole index m (m=1 dipole, m=2 quadrupole, m=3 sextupole etc.)
         * @param K_normal Integrated normal multipole coefficient (1/meter^m)
         * @param K_skew Integrated skew multipole coefficient (1/meter^m)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param name a user defined and not necessarily unique name of the element
         */
        Multipole (
            int multipole,
            amrex::ParticleReal K_normal,
            amrex::ParticleReal K_skew,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            std::optional<std::string> name = std::nullopt
        )
        : Named(std::move(name)),
          Alignment(dx, dy, rotation_degree),
          m_multipole(multipole), m_Kn(K_normal), m_Ks(K_skew)
        {
            compute_factorial();
        }

        /** Reverse the element (negate strengths) */
        void reverse () { m_Kn = -m_Kn; m_Ks = -m_Ks; }

        /** Push all particles */
        using BeamOptic::operator();

        /** Compute and cache the constants for the push.
         *
         * In particular, used to pre-compute and cache variables that are
         * independent of the individually tracked particle.
         *
         * @param refpart reference particle
         */
        void compute_constants (RefPart const & refpart)
        {
            using namespace amrex::literals; // for _rt and _prt

            Alignment::compute_constants(refpart);

            // assign complex position and complex multipole strength
            m_alpha = {m_Kn, m_Ks};

            m_inv_factorial = 1.0_prt / static_cast<amrex::ParticleReal>(m_mfactorial);
        }

        /** This is a multipole functor, so that a variable of this type can be used like a
         *  multipole function.
         *
         * The @see compute_constants method must be called before pushing particles through this operator.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index (unused)
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            [[maybe_unused]] T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py,
            T_Real & AMREX_RESTRICT pt,
            [[maybe_unused]] T_IdCpu & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // a complex type with two T_Real
            using Complex = amrex::GpuComplex<T_Real>;

            // access reference particle values to find (beta*gamma)^2
            //amrex::ParticleReal const pt_ref = refpart.pt;
            //amrex::ParticleReal const betgam2 = powi<2>(pt_ref) - 1.0_prt;

            // initialize output values
            T_Real xout = x;
            T_Real yout = y;
            T_Real tout = t;
            T_Real pxout = px;
            T_Real pyout = py;
            T_Real ptout = pt;

            // assign complex position and complex multipole strength
            Complex const zeta(x, y);

            // compute complex momentum kick
            int const m = m_multipole - 1;
            Complex const kick = amrex::pow(zeta, m) * m_alpha;
            T_Real const dpx = -1.0_prt * kick.m_real * m_inv_factorial;
            T_Real const dpy = kick.m_imag * m_inv_factorial;

            // advance position and momentum
            // xout = x;
            pxout = px + dpx;

            // yout = y;
            pyout = py + dpy;

            // tout = t;
            ptout = pt;

            // assign updated values
            x = xout;
            y = yout;
            t = tout;
            px = pxout;
            py = pyout;
            pt = ptout;

            // undo shift due to alignment errors of the element
            shift_out(x, y, px, py);
        }

        /** This pushes the reference particle. */
        using Thin::operator();

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
         *
         * The thin multipole kick implemented by @see operator() is
         *   dpx - i*dpy = -(K_normal + i*K_skew) * (x + i*y)^(m-1) / (m-1)!
         * matching the MAD-X MULTIPOLE thin kick (MAD-X User Guide,
         * "Generalisation to normal and skew components"),
         *   dPx - i*dPy = -sum_n (KN_n + i*KS_n) * (x + i*y)^n / n!,
         * with the index mapping n = m_multipole - 1 and
         * K_normal <-> KNL[n], K_skew <-> KSL[n].
         *
         * The linear map is the Jacobian of this kick around the reference
         * trajectory, i.e. at (x, y) = (0, 0). Ideal behavior only: we
         * ignore alignment and rotation errors, and we ignore feed-down
         * from an off-axis reference orbit (the magnetic center and the
         * reference trajectory are assumed coincident).
         *
         * Three cases:
         *  - m_multipole == 1 (thin dipole kick): constant, yields the
         *    identity linear map.
         *  - m_multipole == 2 (thin quadrupole): the only nontrivial case;
         *    R(2,1) = -K_normal, R(2,3) = +K_skew,
         *    R(4,1) = +K_skew,   R(4,3) = +K_normal.
         *    (Wiedemann, Particle Accelerator Physics, 3rd ed., Sect. 6.3;
         *    Dragt, Lie Methods for Nonlinear Dynamics with Applications
         *    to Accelerator Physics, Sect. 7.3.)
         *  - m_multipole >= 3 (sextupole and higher): the kick vanishes at
         *    least quadratically at the origin, yielding the identity.
         *
         * Sign/units convention: ImpactX phase-space ordering is
         * (x, px, y, py, t, pt) with 1-based matrix indices, matching the
         * other linear elements in this directory (e.g. Drift, Quad).
         *
         * @param[in] refpart reference particle (unused)
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map ([[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt

            Map6x6 R = Map6x6::Identity();

            // Only m_multipole == 2 (thin quadrupole) has a nonzero Jacobian
            // at the reference trajectory; m_multipole == 1 is a constant
            // kick and m_multipole >= 3 vanishes at the origin.
            if (m_multipole == 2) {
                R(2,1) = -m_Kn;
                R(2,3) = +m_Ks;
                R(4,1) = +m_Ks;
                R(4,3) = +m_Kn;
            }

            return R;
        }

        int m_multipole; //! multipole index
        amrex::ParticleReal m_Kn; //! integrated normal multipole coefficient
        amrex::ParticleReal m_Ks; //! integrated skew multipole coefficient

      private:
        // constants that are independent of the individually tracked particle
        // see: compute_factorial() to refresh
        int m_mfactorial;  //! factorial of multipole index
        // see: compute_constants() to refresh
        amrex::GpuComplex<amrex::ParticleReal> m_alpha;
        amrex::ParticleReal m_inv_factorial;
    };

} // namespace impactx

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::Multipole)

#endif // IMPACTX_MULTIPOLE_H
