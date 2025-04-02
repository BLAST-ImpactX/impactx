/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_EXACTQUAD_H
#define IMPACTX_EXACTQUAD_H

#include "particles/ImpactXParticleContainer.H"
#include "particles/integrators/Integrators.H"
#include "mixin/alignment.H"
#include "mixin/pipeaperture.H"
#include "mixin/beamoptic.H"
#include "mixin/thick.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"
#include "mixin/lineartransport.H"

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_Print.H>

#include <cmath>


namespace impactx::elements
{
    struct ExactQuad
    : public mixin::Named,
      public mixin::BeamOptic<ExactQuad>,
      public mixin::LinearTransport<ExactQuad>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::NoFinalize
    {
        static constexpr auto type = "ExactQuad";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A Quadrupole magnet treated using the exact Hamiltonian, which
         *  includes all nonlinear kinematic effects.  Like the element Quad,
         *  a hard-edge model of the magnetic field is used.  Symplectic
         *  integration is used based on a splitting of the Hamiltonian
         *  described in [to be cited].
         *
         * @param ds Segment length in m.
         * @param k  Quadrupole strength in m^(-2) (MADX convention)
         *           = (gradient in T/m) / (rigidity in T-m)
         *        OR Quadrupole strength in T/m (MaryLie convention)
         *           k > 0 horizontal focusing
         *           k < 0 horizontal defocusing
         * @param unit  Unit specification
         *           unit = 0 MADX convention (default)
         *           unit = 1 MaryLie convention
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param mapsteps number of integration steps per slice used for
         *        particle push in the applied fields
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        ExactQuad (
            amrex::ParticleReal ds,
            amrex::ParticleReal k,
            int unit,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            amrex::ParticleReal aperture_x = 0,
            amrex::ParticleReal aperture_y = 0,
            int mapsteps = 1,
            int nslice = 1,
            std::optional<std::string> name = std::nullopt
        )
        : Named(std::move(name)),
          Thick(ds, nslice),
          Alignment(dx, dy, rotation_degree),
          PipeAperture(aperture_x, aperture_y),
          m_k(k),m_unit(unit),m_mapsteps(mapsteps)
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

            // find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            m_betgam2 = std::pow(pt_ref, 2) - 1.0_prt;
            m_ibetgam2 = 1_prt / m_betgam2;
            m_beta = refpart.beta();
            m_ibeta = 1_prt / m_beta;

            // normalize quad units to MAD-X convention if needed
            m_g = m_unit == 1 ? m_k / refpart.rigidity_Tm() : m_k;

            // compute phase advance per unit length in s (in rad/m)
            m_omega = std::sqrt(std::abs(m_g));

        }

        /** This is a quad functor, so that a variable of this type can be used like a quad function.
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

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // coordinates within the quad (not needed)
            //amrex::ParticleReal const s = refpart.s;
            //amrex::ParticleReal const sedge = refpart.sedge;

            // numerical integration parameters
            amrex::ParticleReal const zin = 0_prt;
            amrex::ParticleReal const zout = m_slice_ds;
            int const nsteps = m_mapsteps;

            // initialize phase space 6-vector
            amrex::SmallVector<amrex::ParticleReal, 6, 1> particle{
                x, px, y, py, t, pt
            };

            // call integrator to advance through slice
            integrators::symp2_integrate_particle(particle,zin,zout,nsteps,*this);
//            integrators::symp4_integrate_particle(particle,zin,zout,nsteps,*this);

            // assign updated values
            x = particle(1);
            px = particle(2);
            y = particle(3);
            py = particle(4);
            t = particle(5);
            pt = particle(6);

            // apply transverse aperture
            apply_aperture(x, y, idcpu);

            // undo shift due to alignment errors of the element
            shift_out(x, y, px, py);
        }

        /** This pushes a particle in an ExactQuad element through
         *  the symplectic map associated with H_1 in the Hamiltonian
         *  splitting H = H_1 + H_2.  This is exactly the linear
         *  map of a quadrupole.
         *
         * @param tau Map step size in m
         * @param[in,out] particle particle phase space 6-vector
         * @param[in,out] zeval Longitudinal on-axis location in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void map1 (amrex::ParticleReal const tau,
                   amrex::SmallVector<amrex::ParticleReal, 6, 1> & particle,
                   amrex::ParticleReal & zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            amrex::ParticleReal const x = particle(1);
            amrex::ParticleReal const px = particle(2);
            amrex::ParticleReal const y = particle(3);
            amrex::ParticleReal const py = particle(4);
            amrex::ParticleReal const t = particle(5);
            amrex::ParticleReal const pt = particle(6);

            amrex::ParticleReal xout = x;
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal yout = y;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal tout = t;
            amrex::ParticleReal ptout = pt;

            amrex::ParticleReal const sin_omega_ds = std::sin(m_omega*tau);
            amrex::ParticleReal const cos_omega_ds = std::cos(m_omega*tau);
            amrex::ParticleReal const sinh_omega_ds = std::sinh(m_omega*tau);
            amrex::ParticleReal const cosh_omega_ds = std::cosh(m_omega*tau);
            amrex::ParticleReal const slice_bg = tau / m_betgam2;

            if (m_g > 0.0_prt)
            {
                // advance position and momentum (focusing quad)
                xout = cos_omega_ds*x + sin_omega_ds/m_omega*px;
                pxout = -m_omega*sin_omega_ds*x + cos_omega_ds*px;

                yout = cosh_omega_ds*y + sinh_omega_ds/m_omega*py;
                pyout = m_omega*sinh_omega_ds*y + cosh_omega_ds*py;

                tout = t + (slice_bg)*pt;
                // ptout = pt;
            } else if (m_g < 0.0_prt)
            {
                // advance position and momentum (defocusing quad)
                xout = cosh_omega_ds*x + sinh_omega_ds/m_omega*px;
                pxout = m_omega*sinh_omega_ds*x + cosh_omega_ds*px;

                yout = cos_omega_ds*y + sin_omega_ds/m_omega*py;
                pyout = -m_omega*sin_omega_ds*y + cos_omega_ds*py;

                tout = t + slice_bg*pt;
                // ptout = pt;
            } else {
                // advance position and momentum (zero strength = drift)
                xout = x + tau * px;
                // pxout = px;
                yout = y + tau * py;
                // pyout = py;
                tout = t + slice_bg * pt;
                // ptout = pt;
            }

            // push particle coordinates
            particle(1) = xout;
            particle(2) = pxout;
            particle(3) = yout;
            particle(4) = pyout;
            particle(5) = tout;
            particle(6) = ptout;

            zeval = zeval;
        }

        /** This pushes a particle in an ExactQuad element through
         *  the symplectic map associated with H_2 in the Hamiltonian
         *  splitting H = H_1 + H_2.
         *
         * @param tau Map step size in m
         * @param[in,out] particle particle phase space 6-vector
         * @param[in,out] zeval Longitudinal on-axis location in m
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void map2 (amrex::ParticleReal const tau,
                   amrex::SmallVector<amrex::ParticleReal, 6, 1> & particle,
                   amrex::ParticleReal & zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            amrex::ParticleReal const x = particle(1);
            amrex::ParticleReal const px = particle(2);
            amrex::ParticleReal const y = particle(3);
            amrex::ParticleReal const py = particle(4);
            amrex::ParticleReal const t = particle(5);
            amrex::ParticleReal const pt = particle(6);

            // compute the radical in the denominator (= pz):
            amrex::ParticleReal const inv_pzden = 1_prt / std::sqrt(
                std::pow(pt - m_ibeta, 2) -
                m_ibetgam2 -
                std::pow(px, 2) -
                std::pow(py, 2)
            );

            // The momenta are not affected.
            particle(1) = x + tau * px * (-1_prt + inv_pzden);
            particle(2) = px;
            particle(3) = y + tau * py * (-1_prt + inv_pzden);
            particle(4) = py;
            particle(5) = t + tau * (-m_ibeta - m_ibetgam2*pt + inv_pzden*(1_prt - pt*m_beta));
            particle(6) = pt;

            zeval += tau;
        }


        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const {  // TODO: update as well, but needs more careful placement of calc_constants

            using namespace amrex::literals; // for _rt and _prt

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
            amrex::ParticleReal const step = slice_ds / std::sqrt(std::pow(pt,2)-1.0_prt);

            // advance position and momentum (straight element)
            refpart.x = x + step*px;
            refpart.y = y + step*py;
            refpart.z = z + step*pz;
            refpart.t = t - step*pt;

            // advance integrated path length
            refpart.s = s + slice_ds;
        }

        /** This function returns the linear transport map.
         *
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map (RefPart const & AMREX_RESTRICT refpart) const  // TODO: update as well, but needs more careful placement of calc_constants
        {
            using namespace amrex::literals; // for _rt and _prt

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // access reference particle values to find beta*gamma^2
            amrex::ParticleReal const pt_ref = refpart.pt;
            amrex::ParticleReal const betgam2 = std::pow(pt_ref, 2) - 1.0_prt;

            // compute phase advance per unit length in s (in rad/m)
            amrex::ParticleReal const omega = std::sqrt(std::abs(m_k));

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();

            if (m_k > 0.0) {
                R(1,1) = std::cos(omega*slice_ds);
                R(1,2) = std::sin(omega*slice_ds)/omega;
                R(2,1) = -omega*std::sin(omega*slice_ds);
                R(2,2) = std::cos(omega*slice_ds);
                R(3,3) = std::cosh(omega*slice_ds);
                R(3,4) = std::sinh(omega*slice_ds)/omega;
                R(4,3) = omega*std::sinh(omega*slice_ds);
                R(4,4) = std::cosh(omega*slice_ds);
                R(5,6) = slice_ds/betgam2;
            } else if (m_k < 0.0) {
                R(1,1) = std::cosh(omega*slice_ds);
                R(1,2) = std::sinh(omega*slice_ds)/omega;
                R(2,1) = omega*std::sinh(omega*slice_ds);
                R(2,2) = std::cosh(omega*slice_ds);
                R(3,3) = std::cos(omega*slice_ds);
                R(3,4) = std::sin(omega*slice_ds)/omega;
                R(4,3) = -omega*std::sin(omega*slice_ds);
                R(4,4) = std::cos(omega*slice_ds);
                R(5,6) = slice_ds/betgam2;
            } else {
                R(1,2) = m_slice_ds;
                R(3,4) = m_slice_ds;
                R(5,6) = m_slice_ds / betgam2;
            }

            return R;
        }

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        amrex::ParticleReal m_k; //! quadrupole strength in 1/m^2 (or T/m)
        int m_unit; //! unit specification for quad strength

    private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_slice_ds;       //! m_ds / nslice();
        amrex::ParticleReal m_betgam2;        //! beta*gamma^2
        amrex::ParticleReal m_ibetgam2;       //! 1 / m_betgam2
        amrex::ParticleReal m_beta;           //! beta
        amrex::ParticleReal m_ibeta;          //! 1 / m_beta
        amrex::ParticleReal m_g;              //! quadrupole strength in 1/m^2
        amrex::ParticleReal m_omega;
        int m_mapsteps; //! number of integration steps per slice
    };

} // namespace impactx

#endif // IMPACTX_EXACTQUAD_H
