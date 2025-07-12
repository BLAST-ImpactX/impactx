/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_CHRPLASMALENS_H
#define IMPACTX_CHRPLASMALENS_H

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

#include <cmath>
#include <stdexcept>


namespace impactx::elements
{
    struct ChrPlasmaLens
    : public mixin::Named,
      public mixin::BeamOptic<ChrPlasmaLens>,
      public mixin::LinearTransport<ChrPlasmaLens>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::NoFinalize
    {
        static constexpr auto type = "ChrPlasmaLens";
        using PType = ImpactXParticleContainer::ParticleType;

        /** An active cylindrically symmetric plasma lens with chromatic focusing
         *
         * The Hamiltonian is expanded through second order in the
         * transverse variables (x,px,y,py), with the exact pt dependence retained.
         *
         * @param ds Segment length in m.
         * @param k  plasma lens focusing strength in m^(-2)
         *           = (azimuthal field gradient in T/m) / (rigidity in T-m)
         *        OR azimuthal magnetic field gradient in T/m (k > 0)
         * @param unit  Unit specification
         *           unit = 0 focusing strength in m^(-2)
         *           unit = 1 focusing strength in T/m
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        ChrPlasmaLens (
            amrex::ParticleReal ds,
            amrex::ParticleReal k,
            int unit,
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
          m_k(k), m_unit(unit)
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

            Alignment::compute_constants(refpart);

            // length of the current slice
            m_slice_ds = m_ds / nslice();

            // access reference particle values to find beta
            m_beta = refpart.beta();

            // normalize focusing strength units to MAD-X convention if needed
            m_g = m_unit == 1 ? m_k / refpart.rigidity_Tm() : m_k;
        }

        /** This is a chrplasmalens functor, so that a variable of this type can be used like a chrplasmalens function.
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
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            amrex::ParticleReal & AMREX_RESTRICT x,
            amrex::ParticleReal & AMREX_RESTRICT y,
            amrex::ParticleReal & AMREX_RESTRICT t,
            amrex::ParticleReal & AMREX_RESTRICT px,
            amrex::ParticleReal & AMREX_RESTRICT py,
            amrex::ParticleReal & AMREX_RESTRICT pt,
            uint64_t & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt
            using amrex::Math::powi;

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // compute particle momentum deviation delta + 1
            amrex::ParticleReal const delta1 = std::sqrt(1_prt - 2_prt*pt/m_beta + powi<2>(pt));
            amrex::ParticleReal const delta = delta1 - 1_prt;

            // compute phase advance per unit length in s (in rad/m)
            // chromatic dependence on delta is included
            amrex::ParticleReal const omega = std::sqrt(std::abs(m_g)/delta1);

            // initialize output values
            amrex::ParticleReal xout = x;
            amrex::ParticleReal yout = y;

            // initialize output values of momenta
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal const ptout = pt;

            // paceholder variables
            amrex::ParticleReal q1 = x;
            amrex::ParticleReal q2 = y;
            amrex::ParticleReal p1 = px;
            amrex::ParticleReal p2 = py;

            auto const [sin_ods, cos_ods] = amrex::Math::sincos(omega*m_slice_ds);

            // advance transverse position and momentum (focusing)
            xout = cos_ods * x + sin_ods / (omega * delta1) * px;
            pxout = -omega * delta1 * sin_ods * x + cos_ods * px;

            yout = cos_ods * y + sin_ods / (omega * delta1) * py;
            pyout = -omega * delta1 * sin_ods * y + cos_ods * py;

            // advance longitudinal position and momentum

            // the corresponding symplectic update to t
            amrex::ParticleReal const term = pt + delta/m_beta;
            amrex::ParticleReal const t0 = t - term*m_slice_ds/delta1;

            amrex::ParticleReal const w = omega*delta1;
            amrex::ParticleReal const term1 = -(powi<2>(p2)-powi<2>(q2) * powi<2>(w)) * std::sin(2_prt*m_slice_ds*omega);
            amrex::ParticleReal const term2 = -(powi<2>(p1)-powi<2>(q1) * powi<2>(w)) * std::sin(2_prt*m_slice_ds*omega);
            amrex::ParticleReal const term3 = -2_prt*q2*p2*w * std::cos(2_prt*m_slice_ds*omega);
            amrex::ParticleReal const term4 = -2_prt*q1*p1*w * std::cos(2_prt*m_slice_ds*omega);
            amrex::ParticleReal const term5 = 2_prt*omega*(q1*p1*delta1 + q2*p2*delta1
                                        -(powi<2>(p1) + powi<2>(p2))*m_slice_ds - (powi<2>(q1) + powi<2>(q2)) * powi<2>(w)*m_slice_ds);
            t = t0 + (-1_prt+m_beta*pt)/(8_prt*m_beta * powi<3>(delta1)*omega)
                     *(term1+term2+term3+term4+term5);

            // ptout = pt;

            // assign updated position & momenta
            x = xout;
            y = yout;
            px = pxout;
            py = pyout;
            pt = ptout;

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

            // assign intermediate parameter
            amrex::ParticleReal const step = slice_ds /std::sqrt(powi<2>(pt)-1.0_prt);

            // advance position and momentum (straight element)
            refpart.x = x + step*px;
            refpart.y = y + step*py;
            refpart.z = z + step*pz;
            refpart.t = t - step*pt;

            // advance integrated path length
            refpart.s = s + m_slice_ds;
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

        amrex::ParticleReal m_k; //! focusing strength in 1/m^2 (or T/m)
        int m_unit; //! unit specification for focusing strength

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_slice_ds;  //! m_ds / nslice();
        amrex::ParticleReal m_beta;
        amrex::ParticleReal m_g;
    };

} // namespace impactx

#endif // IMPACTX_CHRPLASMALENS_H
