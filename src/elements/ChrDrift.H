/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_CHRDRIFT_H
#define IMPACTX_CHRDRIFT_H

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
    struct ChrDrift
    : public mixin::Named,
      public mixin::BeamOptic<ChrDrift>,
      public mixin::LinearTransport<ChrDrift>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::NoFinalize
    {
        static constexpr auto type = "ChrDrift";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A drift with chromatic effects included
         *
         * The Hamiltonian is expanded through second order in the
         * transverse variables (x,px,y,py), with the exact pt dependence retained.
         *
         * @param ds Segment length in m
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        ChrDrift (
            amrex::ParticleReal ds,
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
          PipeAperture(aperture_x, aperture_y)
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

            // access reference particle values to find beta, gamma
            m_beta = refpart.beta();
            m_gamma = refpart.gamma();

            m_const1 = 1_prt / (2_prt * powi<3>(m_beta) * powi<2>(m_gamma));
        }

        /** This is a chrdrift functor, so that a variable of this type can be used like a chrdrift function.
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

            // access position data
            amrex::ParticleReal const xout = x;
            amrex::ParticleReal const yout = y;
            amrex::ParticleReal const tout = t;

            // initialize output values of momenta
            amrex::ParticleReal const pxout = px;
            amrex::ParticleReal const pyout = py;
            amrex::ParticleReal const ptout = pt;

            // compute particle momentum deviation delta + 1
            amrex::ParticleReal const idelta1 = 1_prt / std::sqrt(1_prt - 2_prt*pt/m_beta + powi<2>(pt));

            // advance transverse position and momentum (drift)
            x = xout + m_slice_ds * px * idelta1;
            // pxout = px;
            y = yout + m_slice_ds * py * idelta1;
            // pyout = py;

            // the corresponding symplectic update to t
            amrex::ParticleReal term = 2_prt * powi<2>(pt) + powi<2>(px) + powi<2>(py);
            term = 2_prt - 4_prt*m_beta*pt + powi<2>(m_beta)*term;
            term = -2_prt + powi<2>(m_gamma)*term;
            term = (-1_prt+m_beta*pt)*term;
            term = term * m_const1;
            t = tout - m_slice_ds * (1_prt / m_beta + term * powi<3>(idelta1));
            // ptout = pt;

            // assign updated momenta
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

            // advance position and momentum (drift)
            refpart.x = x + step*px;
            refpart.y = y + step*py;
            refpart.z = z + step*pz;
            refpart.t = t - step*pt;

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

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_slice_ds;  //! m_ds / nslice();
        amrex::ParticleReal m_beta;
        amrex::ParticleReal m_gamma;
        amrex::ParticleReal m_const1;
    };

} // namespace impactx

#endif // IMPACTX_CHRDRIFT_H
