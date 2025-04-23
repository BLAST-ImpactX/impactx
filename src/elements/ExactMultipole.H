/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_EXACTMULTIPOLE_H
#define IMPACTX_EXACTMULTIPOLE_H

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
#include <stdexcept>
#include <vector>

namespace impactx::elements
{

    /** Dynamic data for the ExactMultipole elements
     *
     * Since we copy the element to the device, we cannot store this data on the element itself.
     * But we can store pointers to this data with the element and keep a lookup table here,
     * which we clean up in the end.
     */
    namespace MultipoleData
    {           
        //! last used id for a created ExactMultipole
        inline int next_id = 0;
         
        //! host: normal multipole coefficients of the magnetic field
        inline std::map<int, std::vector<amrex::ParticleReal>> h_k_normal = {};
        //! host: skew multipole coefficients of the magnetic field
        inline std::map<int, std::vector<amrex::ParticleReal>> h_k_skew = {};
         
        //! device: normal multipole coefficients of the magnetic field
        inline std::map<int, amrex::Gpu::DeviceVector<amrex::ParticleReal>> d_k_normal = {};
        //! device: skew multipole coefficients of the magnetic field
        inline std::map<int, amrex::Gpu::DeviceVector<amrex::ParticleReal>> d_k_skew = {};
         
    } // namespace MultipoleData


    struct ExactMultipole
    : public mixin::Named,
      public mixin::BeamOptic<ExactMultipole>,
      public mixin::LinearTransport<ExactMultipole>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture,
      public mixin::NoFinalize
    {
        static constexpr auto type = "ExactMultipole";
        using PType = ImpactXParticleContainer::ParticleType;

        /** A general (thick) straight-axis magnetic element whose field consists
         *  of a superposition of ideal hard-edge multipoles.  Like the element
         *  ExactQuad, the magnet is treated using the exact Hamiltonian, which 
         *  includes all nonlinear kinematic effects.  This is a generalization
         *  of the (thin) element Multipole, which applies the kick associated
         *  with a single ideal multipole.  Here, symplectic integration
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
         * @param int_order order used for symplectic integration (2 or 4)
         * @param mapsteps number of integration steps per slice used for
         *        particle push in the applied fields
         * @param nslice number of slices used for the application of space charge
         * @param name a user defined and not necessarily unique name of the element
         */
        ExactMultipole (
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
          m_unit(unit),m_int_order(int_order),m_mapsteps(mapsteps),m_id(MultipoleData::next_id)
        {
            // next created ExactMultipole has another id for its data  
            MultipoleData::next_id++;
         
            // validate normal and skew coefficients are the same length
            m_ncoef = int(k_normal.size());
            if (m_ncoef != int(k_skew.size()))
                throw std::runtime_error("ExactMultipole: normal and skew coefficient arrays must have same length!");
            
            // host data
            MultipoleData::h_k_normal[m_id] = k_normal;
            MultipoleData::h_k_skew[m_id] = k_skew;
            m_k_normal_h_data = MultipoleData::h_k_normal[m_id].data();
            m_k_skew_h_data = MultipoleData::h_k_skew[m_id].data();
            
            // device data
            MultipoleData::d_k_normal.emplace(m_id, amrex::Gpu::DeviceVector<amrex::ParticleReal>(m_ncoef));
            MultipoleData::d_k_skew.emplace(m_id, amrex::Gpu::DeviceVector<amrex::ParticleReal>(m_ncoef));
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  k_normal.begin(), k_normal.end(),
                                  MultipoleData::d_k_normal[m_id].begin());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  k_skew.begin(), k_skew.end(),
                                  MultipoleData::d_k_skew[m_id].begin());
            amrex::Gpu::streamSynchronize();
            
            // low-level objects we can use on device
            m_k_normal_d_data = MultipoleData::d_k_normal[m_id].data();
            m_k_skew_d_data = MultipoleData::d_k_skew[m_id].data();

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

            // print coefficients
            // pick the right data depending if we are on the host side
            // (reference particle push) or device side (particles):
#if AMREX_DEVICE_COMPILE
            amrex::ParticleReal* k_normal = m_k_normal_d_data;
//            amrex::ParticleReal* k_skew = m_k_skew_d_data;
#else
            amrex::ParticleReal* k_normal = m_k_normal_h_data;
//            amrex::ParticleReal* k_skew = m_k_skew_h_data;
#endif

            // normalize quad units to MAD-X convention if needed
            m_g = m_unit == 1 ? k_normal[2] / refpart.rigidity_Tm() : k_normal[2];
            
            // compute phase advance per unit length in s (in rad/m)
            m_omega = std::sqrt(std::abs(m_g));
        }

        /** This is a Multipole functor, so that a variable of this type can be used like a Multipole function.
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
                // advance position and momentum (focusing Multipole)
                xout = cos_omega_ds*x + sin_omega_ds/m_omega*px;
                pxout = -m_omega*sin_omega_ds*x + cos_omega_ds*px;

                yout = cosh_omega_ds*y + sinh_omega_ds/m_omega*py;
                pyout = m_omega*sinh_omega_ds*y + cosh_omega_ds*py;

                tout = t + (slice_bg)*pt;
                // ptout = pt;
            } else if (m_g < 0.0_prt)
            {
                // advance position and momentum (defocusing Multipole)
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

            zeval += 0_prt;
        }

        /** This pushes a particle in an ExactMultipole element through
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
            amrex::ParticleReal const omega = std::sqrt(std::abs(m_g));

            // initialize linear map matrix elements
            Map6x6 R = Map6x6::Identity();

            if (m_g > 0.0) {
                R(1,1) = std::cos(omega*slice_ds);
                R(1,2) = std::sin(omega*slice_ds)/omega;
                R(2,1) = -omega*std::sin(omega*slice_ds);
                R(2,2) = std::cos(omega*slice_ds);
                R(3,3) = std::cosh(omega*slice_ds);
                R(3,4) = std::sinh(omega*slice_ds)/omega;
                R(4,3) = omega*std::sinh(omega*slice_ds);
                R(4,4) = std::cosh(omega*slice_ds);
                R(5,6) = slice_ds/betgam2;
            } else if (m_g < 0.0) {
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


        /** Close and deallocate all data and handles.
         */ 
        void
        finalize ()
        {
            // remove from unique data map
            if (MultipoleData::h_k_normal.count(m_id) != 0u)
                MultipoleData::h_k_normal.erase(m_id);
            if (MultipoleData::h_k_skew.count(m_id) != 0u)
                MultipoleData::h_k_skew.erase(m_id);
         
            if (MultipoleData::d_k_normal.count(m_id) != 0u)
                MultipoleData::d_k_normal.erase(m_id);
            if (MultipoleData::d_k_skew.count(m_id) != 0u)
                MultipoleData::d_k_skew.erase(m_id); 
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
        amrex::ParticleReal m_betgam2;        //! beta*gamma^2
        amrex::ParticleReal m_ibetgam2;       //! 1 / m_betgam2
        amrex::ParticleReal m_beta;           //! beta
        amrex::ParticleReal m_ibeta;          //! 1 / m_beta
        amrex::ParticleReal m_g;              //! quadrupole strength in 1/m^2
        std::vector<amrex::ParticleReal> g_normal;
        std::vector<amrex::ParticleReal> g_skew;
        amrex::ParticleReal m_omega;
    };

} // namespace impactx

#endif // IMPACTX_EXACTMULTIPOLE_H
