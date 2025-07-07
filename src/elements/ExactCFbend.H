/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_EXACTCFBEND_H
#define IMPACTX_EXACTCFBEND_H

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
#include <AMReX_GpuComplex.H>
#include <AMReX_Math.H>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace impactx::elements
{

    /** Dynamic data for the ExactCFbend elements
     *
     * Since we copy the element to the device, we cannot store this data on the element itself.
     * But we can store pointers to this data with the element and keep a lookup table here,
     * which we clean up in the end.
     */
    namespace ExactCFbendData
    {
        //! last used id for a created ExactCFbend
        inline int next_id = 0;

        //! host: normal multipole coefficients of the magnetic field
        inline std::map<int, std::vector<amrex::ParticleReal>> h_k_normal = {};
        //! host: skew multipole coefficients of the magnetic field
        inline std::map<int, std::vector<amrex::ParticleReal>> h_k_skew = {};

        //! device: normal multipole coefficients of the magnetic field
        inline std::map<int, amrex::Gpu::DeviceVector<amrex::ParticleReal>> d_k_normal = {};
        //! device: skew multipole coefficients of the magnetic field
        inline std::map<int, amrex::Gpu::DeviceVector<amrex::ParticleReal>> d_k_skew = {};

    } // namespace ExactCFbendData


    struct ExactCFbend
    : public mixin::Named,
      public mixin::BeamOptic<ExactCFbend>,
      public mixin::LinearTransport<ExactCFbend>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::NoFinalize
    {
        static constexpr auto type = "ExactCFbend";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A general (thick) combined-function sector dipole whose field consists
         *  of a superposition of ideal curvilinear hard-edge multipoles.
         *  The magnet is treated using the exact Hamiltonian, which
         *  includes all nonlinear kinematic effects.  This is a generalization
         *  of the element CFbend, which is a linear sector bend that
         *  supports quadrupole focusing.  Here, symplectic integration
         *  is used based on a splitting of the underlying Hamiltonian.
         *
         * @param ds Segment length in m.
         * @param k_normal Array of normal multipole coefficients (units below)
         * @param k_skew   Array of skew multipole coefficients (units below)
         * @param unit  Unit specification
         *           unit = 0 multipole strengths in units of 1/m^n (default)
         *           unit = 1 multipole strengths in units of T/m^(n-1)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param int_order order used for symplectic integration (2, 4, or 6)
         * @param mapsteps number of integration steps per slice used for
         *        particle push in the applied fields
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        ExactCFbend (
            amrex::ParticleReal ds,
            std::vector<amrex::ParticleReal> k_normal,
            std::vector<amrex::ParticleReal> k_skew,
            int unit,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            amrex::ParticleReal aperture_x = 0,
            amrex::ParticleReal aperture_y = 0,
            int int_order = 2,
            int mapsteps = 1,
            int nslice = 1,
            std::optional<std::string> name = std::nullopt
        )
        : Named(std::move(name)),
          Thick(ds, nslice),
          Alignment(dx, dy, rotation_degree),
          PipeAperture(aperture_x, aperture_y),
          m_unit(unit),m_int_order(int_order),m_mapsteps(mapsteps),m_id(ExactCFbendData::next_id)
        {
            // next created ExactCFbend has another id for its data
            ExactCFbendData::next_id++;

            // validate normal and skew coefficients are the same length
            m_ncoef = int(k_normal.size());
            if (m_ncoef != int(k_skew.size()))
                throw std::runtime_error("ExactMultipole: normal and skew coefficient arrays must have same length!");

            // host data
            ExactCFbendData::h_k_normal[m_id] = k_normal;
            ExactCFbendData::h_k_skew[m_id] = k_skew;
            m_k_normal_h_data = ExactCFbendData::h_k_normal[m_id].data();
            m_k_skew_h_data = ExactCFbendData::h_k_skew[m_id].data();

            // device data
            ExactCFbendData::d_k_normal.emplace(m_id, amrex::Gpu::DeviceVector<amrex::ParticleReal>(m_ncoef));
            ExactCFbendData::d_k_skew.emplace(m_id, amrex::Gpu::DeviceVector<amrex::ParticleReal>(m_ncoef));
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  k_normal.begin(), k_normal.end(),
                                  ExactCFbendData::d_k_normal[m_id].begin());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  k_skew.begin(), k_skew.end(),
                                  ExactCFbendData::d_k_skew[m_id].begin());
            amrex::Gpu::streamSynchronize();

            // low-level objects we can use on device
            m_k_normal_d_data = ExactCFbendData::d_k_normal[m_id].data();
            m_k_skew_d_data = ExactCFbendData::d_k_skew[m_id].data();

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

            // access reference particle values to find relativistic factors
            m_beta = refpart.beta();
            m_ibeta = 1_prt / m_beta;
            m_ibetgam2 = 1_prt / amrex::Math::powi<2>(refpart.beta_gamma());
            m_brho = refpart.rigidity_Tm();

#if AMREX_DEVICE_COMPILE
            amrex::ParticleReal* k_normal = m_k_normal_d_data;
#else
            amrex::ParticleReal* k_normal = m_k_normal_h_data;
#endif
            // find magnetic field and curvature; normalize units to MAD-X convention if needed
            amrex::ParticleReal curvature = m_unit == 1 ? k_normal[0] / m_brho : k_normal[0];
            m_rc = curvature == 0.0 ? 0.0_prt : 1.0_prt/curvature;

            // arc length and angle of the current slice
            m_slice_ds = m_ds / nslice();

        }

        /** This is an ExactCFbend functor, so that a variable of this type can be used like an ExactCFbend function.
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

            // coordinates within the multipole (not needed)
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

            // call integrator to advance through slice (int_order = 2 or 4, otherwise default 2)
            if (m_int_order == 2) {
                integrators::symp2_integrate_particle(particle,zin,zout,nsteps,*this);
            } else if (m_int_order == 4) {
                integrators::symp4_integrate_particle(particle,zin,zout,nsteps,*this);
            } else if (m_int_order == 6) {
                integrators::symp6_integrate_particle(particle,zin,zout,nsteps,*this);
            } else {
                integrators::symp2_integrate_particle(particle,zin,zout,nsteps,*this);
            }

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

        /** This pushes a particle in an ExactMultipole element through
         *  the symplectic map associated with H_1 in the Hamiltonian
         *  splitting H = H_1 + H_2.  This is the exact map for a sector bend.
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

            // treat the special case of zero field
            if (m_rc==0_prt) {
               // compute the radical in the denominator (= pz):
               amrex::ParticleReal const inv_pzden = 1_prt / std::sqrt(
                   amrex::Math::powi<2>(pt - m_ibeta) -
                   m_ibetgam2 -
                   amrex::Math::powi<2>(px) -
                   amrex::Math::powi<2>(py)
               );

               // advance position and momentum (exact drift)
               xout = x + tau * px * inv_pzden;
               yout = y + tau * py * inv_pzden;
               tout = t - tau * (m_ibeta + (pt - m_ibeta) * inv_pzden);

               // assign updated momenta
               pxout = px;
               pyout = py;
               ptout = pt;

            } else {

               // assign intermediate quantities
               amrex::ParticleReal const pperp2 = amrex::Math::powi<2>(pt)-2.0_prt * m_ibeta * pt - amrex::Math::powi<2>(py)+1.0_prt;
               amrex::ParticleReal const px2 = amrex::Math::powi<2>(px);

               // trigonometric evaluations
               amrex::ParticleReal const slice_phi = tau / m_rc;
               auto const [sin_phi, cos_phi] = amrex::Math::sincos(slice_phi);

               // determine if particle lies within the domain of map definition
               if (pperp2 > px2)
               {
                   amrex::ParticleReal const pperp = std::sqrt(pperp2);
                   amrex::ParticleReal const pzi = std::sqrt(amrex::Math::powi<2>(pperp) - amrex::Math::powi<2>(px));
                   amrex::ParticleReal const rho = m_rc + x;

                   // update momenta
                   pxout = px * cos_phi + (pzi - rho / m_rc) * sin_phi;
                   pyout = py;
                   ptout = pt;

                   // angle of momentum rotation
                   amrex::ParticleReal const px2f = amrex::Math::powi<2>(pxout);
                   // determine if particle lies within domain of map definition
                   if (pperp2 > px2f)
                   {
                       amrex::ParticleReal const pzf = std::sqrt(amrex::Math::powi<2>(pperp)-amrex::Math::powi<2>(pxout));
                       amrex::ParticleReal const theta = slice_phi + std::asin(px/pperp) - std::asin(pxout/pperp);

                       // update position coordinates
                       xout = -m_rc + rho*cos_phi + m_rc*(pzf + px*sin_phi - pzi*cos_phi);
                       yout = y + theta * m_rc * py;
                       tout = t - theta * m_rc * (pt - 1.0_prt * m_ibeta) - slice_phi * m_rc * m_ibeta;

                   } // else { amrex::Print() << "Warning: outside map domain." << "\n"; }
               } // else { amrex::Print() << "Warning: outside map domain." << "\n"; }
            }

            // push particle coordinates
            particle(1) = xout;
            particle(2) = pxout;
            particle(3) = yout;
            particle(4) = pyout;
            particle(5) = tout;
            particle(6) = ptout;

            zeval += 0_prt;
        }

        /** This pushes a particle in an ExactMultipole element through
         *  the symplectic map associated with H_2 in the Hamiltonian
         *  splitting H = H_1 + H_2.  This is the exact kick obtained
         *  in the ideal multipole field.
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

            amrex::ParticleReal xout = x;
            amrex::ParticleReal pxout = px;
            amrex::ParticleReal yout = y;
            amrex::ParticleReal pyout = py;
            amrex::ParticleReal tout = t;
            amrex::ParticleReal ptout = pt;

#if AMREX_DEVICE_COMPILE
            amrex::ParticleReal* k_normal = m_k_normal_d_data;
            amrex::ParticleReal* k_skew = m_k_skew_d_data;
#else
            amrex::ParticleReal* k_normal = m_k_normal_h_data;
            amrex::ParticleReal* k_skew = m_k_skew_h_data;
#endif

            // a complex type with two amrex::ParticleReal
            using Complex = amrex::GpuComplex<amrex::ParticleReal>;

            // assign complex position and complex multipole strength
            Complex const zeta(x, y);

            // initialize the momentum kick
            amrex::ParticleReal dpx = 0_prt;
            amrex::ParticleReal dpy = 0_prt;

            int factorial = 1;

            // loop over array of multipole coefficients
            for (int m = 1; m < m_ncoef; m++) {

               // normalize units to MAD-X convention if needed
               amrex::ParticleReal knormal_value = m_unit == 1 ? k_normal[m] / m_brho : k_normal[m];
               amrex::ParticleReal kskew_value = m_unit == 1 ? k_skew[m] / m_brho : k_skew[m];

               // multipole contribution to complex momentum kick
               Complex alpha(knormal_value, kskew_value);
               Complex kick = amrex::pow(zeta, m) * alpha;
               dpx += -1.0_prt * kick.m_real / factorial;
               dpy += kick.m_imag / factorial;

               factorial *= m+1;

            }

            // advance momentum by a step of length tau
            pxout = pxout + tau * dpx;
            pyout = pyout + tau * dpy;

            //update the particle coordinates
            particle(1) = xout;
            particle(2) = pxout;
            particle(3) = yout;
            particle(4) = pyout;
            particle(5) = tout;
            particle(6) = ptout;

            zeval += tau;
        }


        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const
        {
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

            // access reference particle values to find relativistic factors
            //m_beta = refpart.beta();
            //m_ibeta = 1_prt / m_beta;
            //m_ibetgam2 = 1_prt / amrex::Math::powi<2>(refpart.beta_gamma());
            amrex::ParticleReal const brho = refpart.rigidity_Tm();   
                
#if AMREX_DEVICE_COMPILE
            amrex::ParticleReal* k_normal = m_k_normal_d_data;
#else
            amrex::ParticleReal* k_normal = m_k_normal_h_data;
#endif
            // find magnetic field and curvature; normalize units to MAD-X convention if needed
            amrex::ParticleReal const curvature = m_unit == 1 ? k_normal[0] / brho : k_normal[0];
            amrex::ParticleReal const rc = curvature == 0.0 ? 0.0_prt : 1.0_prt/curvature;
            
            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // special case of zero field = an exact drift
            if (rc==0_prt) {
               // advance position and momentum (drift)
               amrex::ParticleReal const step = slice_ds /std::sqrt(amrex::Math::powi<2>(pt)-1.0_prt);
               refpart.x = x + step*px;
               refpart.y = y + step*py;
               refpart.z = z + step*pz;
               refpart.t = t - step*pt;

            } else {

               // assign intermediate parameters
               amrex::ParticleReal const B = refpart.beta_gamma() /rc;
               amrex::ParticleReal const theta = slice_ds / rc;

               // calculate expensive terms once
               auto const [sin_theta, cos_theta] = amrex::Math::sincos(theta);

               // advance position and momentum (bend)
               refpart.px = px*cos_theta - pz*sin_theta;
               refpart.py = py;
               refpart.pz = pz*cos_theta + px*sin_theta;
               refpart.pt = pt;

               refpart.x = x + (refpart.pz - pz)/B;
               refpart.y = y + (theta/B)*py;
               refpart.z = z - (refpart.px - px)/B;
               refpart.t = t - (theta/B)*pt;

            }

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
            amrex::ParticleReal const betgam2 = amrex::Math::powi<2>(pt_ref) - 1.0_prt;

            // pick the right data depending if we are on the host side
            // (reference particle push) or device side (particles):
#if AMREX_DEVICE_COMPILE
            amrex::ParticleReal* k_normal = m_k_normal_d_data;
            amrex::ParticleReal* k_skew = m_k_skew_d_data;
#else
            amrex::ParticleReal* k_normal = m_k_normal_h_data;
            amrex::ParticleReal* k_skew = m_k_skew_h_data;
#endif

            // normalize quad units to MAD-X convention if needed
            amrex::ParticleReal const kn = m_unit == 1 ? k_normal[1] / refpart.rigidity_Tm() : k_normal[1];
            amrex::ParticleReal const ks = m_unit == 1 ? k_skew[1] / refpart.rigidity_Tm() : k_skew[1];

            amrex::ParticleReal const kabs = std::sqrt(amrex::Math::powi<2>(kn) + amrex::Math::powi<2>(ks));
            amrex::ParticleReal const kplus = kn + kabs;
            amrex::ParticleReal const kminus = kn - kabs;

            // compute phase advance per unit length in s (in rad/m)
            amrex::ParticleReal const omega = std::sqrt(kabs);
            amrex::ParticleReal const sin_omega_ds = std::sin(omega*slice_ds);
            amrex::ParticleReal const cos_omega_ds = std::cos(omega*slice_ds);
            amrex::ParticleReal const sinh_omega_ds = std::sinh(omega*slice_ds);
            amrex::ParticleReal const cosh_omega_ds = std::cosh(omega*slice_ds);

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();

            // linear matrix for a quadrupole with both normal and skew coefficients
            R(1,1) = (kplus * cos_omega_ds - kminus * cosh_omega_ds) / (2_prt * kabs);
            R(1,2) = (kplus * sin_omega_ds - kminus * sinh_omega_ds) / (2_prt * kabs * omega);
            R(1,3) = ks * (-cos_omega_ds + cosh_omega_ds) / (2_prt * kabs);
            R(1,4) = ks * (-sin_omega_ds + sinh_omega_ds) / (2_prt * kabs * omega);
            R(2,1) = (-kplus * sin_omega_ds - kminus * sinh_omega_ds) / (2_prt * omega);
            R(2,2) = (kplus * cos_omega_ds - kminus * cosh_omega_ds) / (2_prt * kabs);
            R(2,3) = ks * (sin_omega_ds + sinh_omega_ds) / (2_prt * omega);
            R(2,4) = ks * (-cos_omega_ds + cosh_omega_ds) / (2_prt * kabs);
            R(3,1) = ks * (-cos_omega_ds + cosh_omega_ds) / (2_prt * kabs);
            R(3,2) = ks * (-sin_omega_ds + sinh_omega_ds) / (2_prt * kabs * omega);
            R(3,3) = (-kminus * cos_omega_ds + kplus * cosh_omega_ds) / (2_prt * kabs);
            R(3,4) = (-kminus * sin_omega_ds + kplus * sinh_omega_ds) / (2_prt * kabs * omega);
            R(4,1) = ks * (sin_omega_ds + sinh_omega_ds) / (2_prt * omega);
            R(4,2) = ks * (-cos_omega_ds + cosh_omega_ds) / (2_prt * kabs);
            R(4,3) = (kminus * sin_omega_ds + kplus * sinh_omega_ds) / (2_prt * omega);
            R(4,4) = (-kminus * cos_omega_ds + kplus * cosh_omega_ds) / (2_prt * kabs);
            R(5,6) = slice_ds/betgam2;

            for (int j=1; j<=6; j++) {
               for (int k=1; k<=6; k++) {
                  amrex::Print() << "R(" << j << "," << k << ") = " << R(j,k) << "\n";
               }
            }

            return R;
        }

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();


        /** Close and deallocate all data and handles.
         */
        void
        finalize ()
        {
            // remove from unique data map
            if (ExactCFbendData::h_k_normal.count(m_id) != 0u)
                ExactCFbendData::h_k_normal.erase(m_id);
            if (ExactCFbendData::h_k_skew.count(m_id) != 0u)
                ExactCFbendData::h_k_skew.erase(m_id);

            if (ExactCFbendData::d_k_normal.count(m_id) != 0u)
                ExactCFbendData::d_k_normal.erase(m_id);
            if (ExactCFbendData::d_k_skew.count(m_id) != 0u)
                ExactCFbendData::d_k_skew.erase(m_id);
        }

        int m_unit; //! unit specification for Multipole strength
        int m_int_order; //! order used for the symplectic integration (2 or 4)
        int m_mapsteps; //! number of integration steps per slice
        int m_id; //! unique ExactMultipole id used for data lookup map

        int m_ncoef = 0; //! number of Fourier coefficients
        amrex::ParticleReal* m_k_normal_h_data = nullptr; //! non-owning pointer to host cosine coefficients
        amrex::ParticleReal* m_k_skew_h_data = nullptr; //! non-owning pointer to host sine coefficients
        amrex::ParticleReal* m_k_normal_d_data = nullptr; //! non-owning pointer to device cosine coefficients
        amrex::ParticleReal* m_k_skew_d_data = nullptr; //! non-owning pointer to device sine coefficients

    private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_slice_ds;       //! m_ds / nslice();
        amrex::ParticleReal m_ibetgam2;       //! 1 / m_betgam2
        amrex::ParticleReal m_beta;           //! beta
        amrex::ParticleReal m_ibeta;          //! 1 / m_beta
        amrex::ParticleReal m_brho;           //! magnetic ridigity in T-m
        amrex::ParticleReal m_rc;             //! magnet radius of curvature
    };

} // namespace impactx

#endif // IMPACTX_EXACTCFBEND_H
