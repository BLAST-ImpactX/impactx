/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_CHRACC_H
#define IMPACTX_CHRACC_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/pipeaperture.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"
#include "mixin/lineartransport.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>
#include <stdexcept>


namespace impactx::elements
{
    struct ChrAcc
    : public mixin::Named,
      public mixin::BeamOptic<ChrAcc>,
      public mixin::LinearTransport<ChrAcc>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "ChrAcc";
        using PType = ImpactXParticleContainer::ParticleType;

        /** Acceleration in a uniform field Ez, with a uniform solenoidal field Bz.
         *
         * The Hamiltonian is expanded through second order in the
         * transverse variables (x,px,y,py), with the exact pt dependence retained.
         *
         * @param ds Segment length in m.
         * @param ez electric field strength in m^(-1)
         *           = (charge * electric field Ez in V/m) / (m*c^2)
         * @param bz magnetic field strength in m^(-1)
         *           = (charge * magnetic field Bz in T) / (m*c)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        ChrAcc (
            amrex::ParticleReal ds,
            amrex::ParticleReal ez,
            amrex::ParticleReal bz,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            amrex::ParticleReal aperture_x = 0,
            amrex::ParticleReal aperture_y = 0,
            int nslice = 1,
            std::optional<std::string> name = std::nullopt
        )
        : Named(std::move(name)),
          Thick(ds, nslice),
          Alignment(dx, dy, rotation_degree),
          PipeAperture(aperture_x, aperture_y),
          m_ez(ez), m_bz(bz)
        {
        }

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
            using amrex::Math::powi;

            Alignment::compute_constants(refpart);

            // length of the current slice
            m_slice_ds = m_ds / nslice();

            // access reference particle values (final, initial):
            m_ptf_ref = refpart.pt;
            m_pti_ref = m_ptf_ref + m_ez * m_slice_ds;
            m_bgf = std::sqrt(powi<2>(m_ptf_ref) - 1.0_prt);
            m_bgi = std::sqrt(powi<2>(m_pti_ref) - 1.0_prt);

            // compute focusing constant (1/m) and rotation angle (in rad)
            m_alpha = m_bz * 0.5_prt;
            m_alpha_iez = m_alpha / m_ez;
        }

        /** This is a chracc functor, so that a variable of this type can be used like a chracc function.
         *
         * The @see compute_constants method must be called before pushing particles through this operator.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py,
            T_Real & AMREX_RESTRICT pt,
            T_IdCpu & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;
            using namespace std; // for cmath(float)
            namespace stdx = amrex::simd::stdx;

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // initial conversion from static to dynamic units:
            px = px * m_bgi;
            py = py * m_bgi;
            pt = pt * m_bgi;

            // compute intermediate quantities related to acceleration
            T_Real const pti_tot = m_pti_ref + pt;
            T_Real const ptf_tot = m_ptf_ref + pt;
            T_Real const pti_tot2 = powi<2>(pti_tot);
            T_Real const ptf_tot2 = powi<2>(ptf_tot);

            // check whether particle lies within the domain of map definition
            auto const mask = pti_tot2 <= 1_prt || ptf_tot2 <= 1_prt;
            if (stdx::any_of(mask) == true) {
                amrex::ParticleIDWrapper<T_IdCpu>{idcpu}.make_invalid(mask);
            } else {
                T_Real const pzi_tot = sqrt(powi<2>(pti_tot) - 1_prt);
                T_Real const pzf_tot = sqrt(powi<2>(ptf_tot) - 1_prt);
                T_Real const pzi_ref = sqrt(powi<2>(m_pti_ref) - 1_prt);
                T_Real const pzf_ref = sqrt(powi<2>(m_ptf_ref) - 1_prt);

                T_Real const numer = -ptf_tot + pzf_tot;
                T_Real const denom = -pti_tot + pzi_tot;

                // compute focusing constant (1/m) and rotation angle (in rad)
                T_Real const theta = m_alpha_iez * log(numer / denom);
                auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);

                // initialize output values
                T_Real xout = x;
                T_Real yout = y;
                T_Real tout = t;
                T_Real pxout = px;
                T_Real pyout = py;
                T_Real ptout = pt;

                // advance positions and momenta using map for focusing
                xout = cos_theta * x + sin_theta / m_alpha * px;
                pxout = -m_alpha * sin_theta * x + cos_theta * px;

                yout = cos_theta * y + sin_theta / m_alpha * py;
                pyout = -m_alpha * sin_theta * y + cos_theta * py;

                // the correct symplectic update for t
                tout = t + (pzf_tot - pzf_ref - pzi_tot + pzi_ref) / m_ez;
                tout += (1_prt/pzi_tot - 1_prt/pzf_tot)
                      * (powi<2>(py - m_alpha * x) +
                         powi<2>(px + m_alpha * y))
                      / (2_prt * m_ez);
                ptout = pt;

                // assign intermediate momenta
                px = pxout;
                py = pyout;
                pt = ptout;

                // advance positions and momenta using map for rotation
                x = cos_theta * xout + sin_theta * yout;
                pxout = cos_theta * px + sin_theta * py;

                y = -sin_theta * xout + cos_theta * yout;
                pyout = -sin_theta * px + cos_theta * py;

                t = tout;
                ptout = pt;

                // assign updated momenta
                px = pxout;
                py = pyout;
                pt = ptout;

            }

            // final conversion from dynamic to static units:
            px = px / m_bgf;
            py = py / m_bgf;
            pt = pt / m_bgf;

            // apply transverse aperture
            apply_aperture(x, y, idcpu);

            // undo shift due to alignment errors of the element
            shift_out(x, y, px, py);
        }

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            // assign input reference particle values
            amrex::ParticleReal const x = refpart.x;
            amrex::ParticleReal const px = refpart.px;
            amrex::ParticleReal const y = refpart.y;
            amrex::ParticleReal const py = refpart.py;
            amrex::ParticleReal const z = refpart.z;
            amrex::ParticleReal const pz = refpart.pz;
            amrex::ParticleReal const t = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;
            amrex::ParticleReal const s = refpart.s;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // compute initial value of beta*gamma
            amrex::ParticleReal const bgi = std::sqrt(powi<2>(pt) - 1.0_prt);

            // advance pt (uniform acceleration)
            refpart.pt = pt - m_ez*slice_ds;

            // compute final value of beta*gamma
            amrex::ParticleReal const ptf = refpart.pt;
            amrex::ParticleReal const bgf = std::sqrt(powi<2>(ptf) - 1.0_prt);

            // update t
            refpart.t = t + (bgf - bgi)/m_ez;

            // advance position (x,y,z)
            refpart.x = x + slice_ds*px/bgi;
            refpart.y = y + slice_ds*py/bgi;
            refpart.z = z + slice_ds*pz/bgi;

            // advance momentum (px,py,pz)
            refpart.px = px*bgf/bgi;
            refpart.py = py*bgf/bgi;
            refpart.pz = pz*bgf/bgi;

            // advance integrated path length
            refpart.s = s + slice_ds;
        }

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
         *
         * @param[in] refpart reference particle
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map ([[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            throw std::runtime_error(std::string(type) + ": Envelope tracking is not yet implemented!");
            return Map6x6::Identity();
        }

        amrex::ParticleReal m_ez; //! electric field strength in 1/m
        amrex::ParticleReal m_bz; //! magnetic field strength in 1/m

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_slice_ds;  //! m_ds / nslice();
        amrex::ParticleReal m_ptf_ref, m_pti_ref, m_bgf, m_bgi;
        amrex::ParticleReal m_alpha, m_alpha_iez;
    };

} // namespace impactx

#endif // IMPACTX_CHRACC_H
