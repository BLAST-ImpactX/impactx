/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_MAGNETOSTATIC_VECTOR_POTENTIAL_H
#define IMPACTX_MAGNETOSTATIC_VECTOR_POTENTIAL_H

#include "particles/ImpactXParticleContainer.H"
#include "mixin/gpuelement.H"
#include "particles/integrators/Integrators.H"
#include "helper/TrackedParser.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/lineartransport.H"
#include "mixin/pipeaperture.H"
#include "mixin/thick.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>
#include <AMReX_SmallMatrix.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>


namespace impactx::elements
{
    /** The prepared, resource-free pusher for a user-defined magnetostatic vector potential.
     *
     * This is what a particle-push kernel is handed. It owns nothing: the user-defined
     * formulas reach it as compiled, non-owning executors, prepared by
     * @see MagnetostaticVectorPotentialGPU::make_particle_pusher for the lifetime of one
     * element slice.
     */
    struct MagnetostaticVectorPotentialPusher
    : public mixin::ParticlePush<MagnetostaticVectorPotentialPusher>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture
    {
        static constexpr auto type = "MagnetostaticVectorPotential";
        using PType = ImpactXParticleContainer::ParticleType;

        static constexpr int DEFAULT_unit = 0;
        static constexpr int DEFAULT_int_order = 2;
        static constexpr int DEFAULT_mapsteps = 10;

        /** The particle-facing configuration of a user-defined magnetostatic vector potential
         *
         * @param ds Segment length in m.
         * @param unit  Unit specification for the vector potential
         *           unit = 0 dimensionless (default)
         *           unit = 1 SI units of T-m
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param int_order order used for symplectic integration (2, 4 or 6)
         * @param mapsteps number of integration steps per slice used for
         *        particle push in the applied fields
         * @param nslice number of slices used for the application of space charge
         */
        MagnetostaticVectorPotentialPusher (
            amrex::ParticleReal ds,
            int unit,
            amrex::ParticleReal dx,
            amrex::ParticleReal dy,
            amrex::ParticleReal rotation_degree,
            amrex::ParticleReal aperture_x,
            amrex::ParticleReal aperture_y,
            int int_order,
            int mapsteps,
            int nslice
        )
          : Thick(ds, nslice),
            Alignment(dx, dy, rotation_degree),
            PipeAperture(aperture_x, aperture_y),
            m_unit(unit), m_int_order(int_order), m_mapsteps(mapsteps)
        {
        }

        /** Reverse the element (negate ds)
         *
         * The user-supplied potential A(x,y,z) is defined in the element's
         * local frame. Reversing only negates the path-length integration
         * direction; the user is responsible for the field's z-symmetry.
         */
        void reverse () { Thick::reverse(); }

        /** Push all particles on a tile */
        using mixin::ParticlePush<MagnetostaticVectorPotentialPusher>::operator();

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

            m_ibeta = 1_prt / refpart.beta();

            // normalize the vector potential if SI units are specified
            m_scale = m_unit == 1 ? 1_prt / refpart.rigidity_Tm() : 1_prt;
        }

        /** Push a single particle through this element.
         *
         * The pusher must be prepared through ``make_particle_pusher`` on the owning element
         * before particles are pushed through this operator.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param idcpu particle global index
         * @param refpart reference particle
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
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // longitudinal coordinate z within the element, measured from its center:
            // the reference particle was already advanced through this slice
            amrex::ParticleReal const s = refpart.s - m_slice_ds;
            amrex::ParticleReal const zin = s - refpart.sedge - m_ds/2_prt;
            amrex::ParticleReal const zout = zin + m_slice_ds;
            int const nsteps = m_mapsteps;

            // initialize phase space 6-vector
            amrex::SmallVector<amrex::ParticleReal, 6, 1> particle{
                x, px, y, py, t, pt
            };

            // call integrator to advance through slice; a particle it cannot advance is lost
            bool const converged = integrators::symp_integrate_semiexplicit(
                particle, zin, zout, m_int_order, nsteps, *this);
            amrex::ParticleIDWrapper<uint64_t>{idcpu}.make_invalid(!converged);

            // assign updated values
            x  = particle(1);
            px = particle(2);
            y  = particle(3);
            py = particle(4);
            t  = particle(5);
            pt = particle(6);
        }

        /** This pushes two copies of the particle coordinates in the
         *  extended (doubled) phase space using the derivatives of
         *  the specified Hamiltonian.
         *
         * @param tau Map step size in m
         * @param[in,out] particle1 particle phase space 6-vector, copy 1
         * @param[in,out] particle2 particle phase space 6-vector, copy 2
         * @param[in,out] zeval1 lattice position in m for particle1
         * @param[in,out] zeval2 lattice position in m for particle2
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void map (amrex::ParticleReal const tau,
                  amrex::SmallVector<amrex::ParticleReal, 6, 1> & particle1,
                  amrex::SmallVector<amrex::ParticleReal, 6, 1> & particle2,
                  amrex::ParticleReal & zeval1,
                  amrex::ParticleReal & zeval2) const
        {
            amrex::ParticleReal const x1  = particle1(1);
            amrex::ParticleReal const px1 = particle1(2);
            amrex::ParticleReal const y1  = particle1(3);
            amrex::ParticleReal const py1 = particle1(4);
            amrex::ParticleReal const t1  = particle1(5);
            amrex::ParticleReal const pt1 = particle1(6);

            amrex::ParticleReal const x2  = particle2(1);
            amrex::ParticleReal const px2 = particle2(2);
            amrex::ParticleReal const y2  = particle2(3);
            amrex::ParticleReal const py2 = particle2(4);
            amrex::ParticleReal const t2  = particle2(5);
            amrex::ParticleReal const pt2 = particle2(6);

            // Evaluate the derivatives of the Hamiltonian.
            auto [dHdx, dHdpx, dHdy, dHdpy, dHdt, dHdpt] =
                Hamiltonian_Gradient(x1, px2, y1, py2, t1, pt2, zeval1);

            // Update particle coordinates and momenta.
            particle1(1) = x1;
            particle1(2) = px1 - tau * dHdx;
            particle1(3) = y1;
            particle1(4) = py1 - tau * dHdy;
            particle1(5) = t1;
            particle1(6) = pt1 - tau * dHdt;

            particle2(1) = x2 + tau * dHdpx;
            particle2(2) = px2;
            particle2(3) = y2 + tau * dHdpy;
            particle2(4) = py2;
            particle2(5) = t2 + tau * dHdpt;
            particle2(6) = pt2;

            zeval2 += tau;
        }

        /** Evaluate the gradient of the Hamiltonian.
         *
         * @param x horizontal coordinate in m
         * @param px horizontal momentum
         * @param y vertical coordinate in m
         * @param py vertical momentum
         * @param t temporal coordinate in m (unused)
         * @param pt temporal momentum (energy deviation)
         * @param zeval longitudinal position in m, measured from the element center
         * @returns tuple containing derivatives of the Hamiltonian
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        std::tuple<amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
                   amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal>
        Hamiltonian_Gradient (amrex::ParticleReal x, amrex::ParticleReal px,
                              amrex::ParticleReal y, amrex::ParticleReal py,
                              [[maybe_unused]] amrex::ParticleReal t, amrex::ParticleReal pt,
                              amrex::ParticleReal zeval) const
        {
            using namespace amrex::literals; // for _rt and _prt

            amrex::ParticleReal const z = zeval;

            // Evaluate the vector potential and its derivatives in the local frame.
            // Not inlined: one shared parser instance instead of eight at every stage
            // of the integrator, which keeps the GPU compile time and stack small.
            auto const ax    = amrex::callNoinline(m_ax, x, y, z)    * m_scale;
            auto const ay    = amrex::callNoinline(m_ay, x, y, z)    * m_scale;
            auto const daxdx = amrex::callNoinline(m_daxdx, x, y, z) * m_scale;
            auto const daxdy = amrex::callNoinline(m_daxdy, x, y, z) * m_scale;
            auto const daydx = amrex::callNoinline(m_daydx, x, y, z) * m_scale;
            auto const daydy = amrex::callNoinline(m_daydy, x, y, z) * m_scale;
            auto const dazdx = amrex::callNoinline(m_dazdx, x, y, z) * m_scale;
            auto const dazdy = amrex::callNoinline(m_dazdy, x, y, z) * m_scale;

            // Rectilinear magnetostatic element: pz is the radical appearing
            // in the Hamiltonian.
            amrex::ParticleReal const pxmech = px - ax;
            amrex::ParticleReal const pymech = py - ay;
            amrex::ParticleReal const radical = std::sqrt(
                1_prt + pt*pt - pxmech*pxmech - pymech*pymech - 2_prt*pt*m_ibeta
            );

            amrex::ParticleReal const dHdx = -(pxmech * daxdx + pymech * daydx) / radical - dazdx;
            amrex::ParticleReal const dHdy = -(pxmech * daxdy + pymech * daydy) / radical - dazdy;
            amrex::ParticleReal const dHdt = 0_prt;

            amrex::ParticleReal const dHdpx = pxmech / radical;
            amrex::ParticleReal const dHdpy = pymech / radical;
            amrex::ParticleReal const dHdpt = -m_ibeta + (m_ibeta - pt) / radical;

            return std::make_tuple(dHdx, dHdpx, dHdy, dHdpy, dHdt, dHdpt);
        }

        int m_unit; //! unit specification for the vector potential
        int m_int_order; //! order used for the symplectic integration (2, 4 or 6)
        int m_mapsteps; //! number of integration steps per slice

        /** Formulas compiled for this pusher's execution space; non-owning
         *
         * Valid for the slice they were prepared for.
         */
        amrex::ParserExecutor<3> m_ax;     //! A_x(x,y,z)
        amrex::ParserExecutor<3> m_ay;     //! A_y(x,y,z)
        amrex::ParserExecutor<3> m_daxdx;  //! dA_x/dx(x,y,z)
        amrex::ParserExecutor<3> m_daxdy;  //! dA_x/dy(x,y,z)
        amrex::ParserExecutor<3> m_daydx;  //! dA_y/dx(x,y,z)
        amrex::ParserExecutor<3> m_daydy;  //! dA_y/dy(x,y,z)
        amrex::ParserExecutor<3> m_dazdx;  //! dA_z/dx(x,y,z)
        amrex::ParserExecutor<3> m_dazdy;  //! dA_z/dy(x,y,z)

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_slice_ds; //! m_ds / nslice()
        amrex::ParticleReal m_ibeta;    //! 1 / beta
        amrex::ParticleReal m_scale;    //! scale factor defining vector-potential units
    };

    /** User-defined formulas of the vector potential and its transverse derivatives
     *
     * Functions of (x, y, z) in the element-local frame, where z is measured from
     * the element center.
     */
    struct VectorPotentialFormulas
    {
        TrackedParser<3> ax;     //!< A_x(x,y,z)
        TrackedParser<3> ay;     //!< A_y(x,y,z)
        TrackedParser<3> daxdx;  //!< dA_x/dx(x,y,z)
        TrackedParser<3> daxdy;  //!< dA_x/dy(x,y,z)
        TrackedParser<3> daydx;  //!< dA_y/dx(x,y,z)
        TrackedParser<3> daydy;  //!< dA_y/dy(x,y,z)
        TrackedParser<3> dazdx;  //!< dA_z/dx(x,y,z)
        TrackedParser<3> dazdy;  //!< dA_z/dy(x,y,z)
    };

    /** A user-defined magnetostatic vector potential: the canonical element that owns its formulas.
     *
     * The reference-particle push, serialization and finalization run here on the host.
     * Particle pushes run against @see MagnetostaticVectorPotentialPusher.
     */
    struct MagnetostaticVectorPotentialGPU
    : public mixin::BeamOptic<MagnetostaticVectorPotentialGPU>,
      public mixin::LinearTransport<MagnetostaticVectorPotentialGPU, false>,
      public MagnetostaticVectorPotentialPusher
    {
        /** Symplectic integration in a user-defined magnetostatic vector
         *  potential, using the exact Hamiltonian, which includes all
         *  nonlinear kinematic effects.  Integration is performed with
         *  respect to a Cartesian coordinate system local to the body of
         *  the element.  A symmetric, semi-explicit symplectic integration
         *  scheme is used, based on:
         *
         *  B. Jayawardana and T. Ohsawa, ``Semiexplicit symplectic
         *  integrators for non-separable Hamiltonian systems,"
         *  Math. Comput. 92, pp. 251-281 (2022),
         *  https://doi.org/10.1090/mcom/3778
         *
         * The formulas are functions of (x, y, z) in the element-local frame,
         * where z is measured from the element center.
         *
         * @param ds Segment length in m.
         * @param unit  Unit specification for the vector potential
         *           unit = 0 dimensionless (default)
         *           unit = 1 SI units of T-m
         * @param str_ax function string for A_x(x,y,z)
         * @param str_ay function string for A_y(x,y,z)
         * @param str_daxdx function string for dA_x/dx(x,y,z)
         * @param str_daxdy function string for dA_x/dy(x,y,z)
         * @param str_daydx function string for dA_y/dx(x,y,z)
         * @param str_daydy function string for dA_y/dy(x,y,z)
         * @param str_dazdx function string for dA_z/dx(x,y,z)
         * @param str_dazdy function string for dA_z/dy(x,y,z)
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param aperture_x horizontal half-aperture in m
         * @param aperture_y vertical half-aperture in m
         * @param int_order order used for symplectic integration (2, 4 or 6)
         * @param mapsteps number of integration steps per slice used for
         *        particle push in the applied fields
         * @param nslice number of slices used for the application of space charge
         */
        MagnetostaticVectorPotentialGPU (
            amrex::ParticleReal ds,
            int unit = DEFAULT_unit,
            std::string str_ax = "0",
            std::string str_ay = "0",
            std::string str_daxdx = "0",
            std::string str_daxdy = "0",
            std::string str_daydx = "0",
            std::string str_daydy = "0",
            std::string str_dazdx = "0",
            std::string str_dazdy = "0",
            amrex::ParticleReal dx = DEFAULT_dx,
            amrex::ParticleReal dy = DEFAULT_dy,
            amrex::ParticleReal rotation_degree = DEFAULT_rotation_degree,
            amrex::ParticleReal aperture_x = DEFAULT_aperture_x,
            amrex::ParticleReal aperture_y = DEFAULT_aperture_y,
            int int_order = DEFAULT_int_order,
            int mapsteps = DEFAULT_mapsteps,
            int nslice = DEFAULT_nslice
        )
          : MagnetostaticVectorPotentialPusher(ds, unit, dx, dy, rotation_degree,
                                               aperture_x, aperture_y, int_order, mapsteps, nslice),
            m_potential{
                {std::move(str_ax),    variables()},
                {std::move(str_ay),    variables()},
                {std::move(str_daxdx), variables()},
                {std::move(str_daxdy), variables()},
                {std::move(str_daydx), variables()},
                {std::move(str_daydy), variables()},
                {std::move(str_dazdx), variables()},
                {std::move(str_dazdy), variables()}
            }
        {
            check(m_unit, m_int_order, m_mapsteps);
        }

        /** Reject an integration setup the element cannot run with
         *
         * @param unit unit specification for the vector potential (0 or 1)
         * @param int_order order used for symplectic integration (2, 4 or 6)
         * @param mapsteps number of integration steps per slice (> 0)
         */
        static void check (int unit, int int_order, int mapsteps)
        {
            if (unit != 0 && unit != 1)
                throw std::invalid_argument("MagnetostaticVectorPotential: unit must be 0 or 1.");
            if (int_order != 2 && int_order != 4 && int_order != 6)
                throw std::invalid_argument("MagnetostaticVectorPotential: The order used for symplectic integration must be 2, 4 or 6.");
            if (mapsteps <= 0)
                throw std::invalid_argument("MagnetostaticVectorPotential: mapsteps must be > 0!");
        }

        /** The variables of the user-defined formulas */
        static amrex::Vector<std::string> variables () { return {"x", "y", "z"}; }

        /** Push all particles */
        using mixin::BeamOptic<MagnetostaticVectorPotentialGPU>::operator();

        /** This pushes the covariance matrix. */
        using mixin::LinearTransport<MagnetostaticVectorPotentialGPU, false>::operator();

        /** This pushes the reference particle.
         *
         * @param[in,out] refpart reference particle
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        void operator() (RefPart & AMREX_RESTRICT refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt

            // assign input reference particle values
            amrex::ParticleReal const x  = refpart.x;
            amrex::ParticleReal const px = refpart.px;
            amrex::ParticleReal const y  = refpart.y;
            amrex::ParticleReal const py = refpart.py;
            amrex::ParticleReal const z  = refpart.z;
            amrex::ParticleReal const pz = refpart.pz;
            amrex::ParticleReal const t  = refpart.t;
            amrex::ParticleReal const pt = refpart.pt;
            amrex::ParticleReal const s  = refpart.s;

            // length of the current slice
            amrex::ParticleReal const slice_ds = m_ds / nslice();

            // assign intermediate parameter
            amrex::ParticleReal const step =
                slice_ds / std::sqrt(amrex::Math::powi<2>(pt) - 1.0_prt);

            // advance position and momentum (straight element)
            refpart.x = x + step*px;
            refpart.y = y + step*py;
            refpart.z = z + step*pz;
            refpart.t = t - step*pt;

            // advance integrated path length
            refpart.s = s + slice_ds;
        }

        /** Prepare the pusher for one element slice
         *
         * Compiles the formulas on first use in an AMReX session.
         *
         * @param refpart reference particle, already advanced for this slice
         * @return a resource-free pusher holding compiled views of the formulas
         */
        [[nodiscard]] MagnetostaticVectorPotentialPusher make_particle_pusher (RefPart const & refpart) const
        {
            MagnetostaticVectorPotentialPusher pusher = *this;
            pusher.compute_constants(refpart);
            pusher.m_ax    = m_potential.ax.executor();
            pusher.m_ay    = m_potential.ay.executor();
            pusher.m_daxdx = m_potential.daxdx.executor();
            pusher.m_daxdy = m_potential.daxdy.executor();
            pusher.m_daydx = m_potential.daydx.executor();
            pusher.m_daydy = m_potential.daydy.executor();
            pusher.m_dazdx = m_potential.dazdx.executor();
            pusher.m_dazdy = m_potential.dazdy.executor();
            return pusher;
        }

        /** Release the compiled formulas; the formulas themselves survive for reuse */
        void finalize ()
        {
            m_potential.ax.release();
            m_potential.ay.release();
            m_potential.daxdx.release();
            m_potential.daxdy.release();
            m_potential.daydx.release();
            m_potential.daydy.release();
            m_potential.dazdx.release();
            m_potential.dazdy.release();
        }

        //! the formulas this element owns, compiled on demand
        VectorPotentialFormulas m_potential;
    };

    /** MagnetostaticVectorPotential: GPU state plus common host-only metadata */
    using MagnetostaticVectorPotential = Element<MagnetostaticVectorPotentialGPU>;

} // namespace impactx::elements

IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::MagnetostaticVectorPotentialGPU)

#endif // IMPACTX_MAGNETOSTATIC_VECTOR_POTENTIAL_H
