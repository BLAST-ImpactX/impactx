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
#include "particles/integrators/Integrators.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/dynamicdata.H"
#include "mixin/lineartransport.H"
#include "mixin/named.H"
#include "mixin/pipeaperture.H"
#include "mixin/thick.H"

#include <AMReX_Extension.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>
#include <AMReX_SmallMatrix.H>
#include <AMReX_Vector.H>

#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>


namespace impactx::elements
{
    /** Host-side parsers and device-callable executors for the user-defined
     *  vector-potential components and their transverse derivatives.
     *
     * The amrex::Parser is host-only; the corresponding amrex::ParserExecutor
     * is device-callable. We keep both: host-side parsers live in this struct
     * and the executors are obtained on demand in compute_constants().
     *
     * Host data persists across simulations; the executors are released when
     * AMReX shuts down.
     */
    struct MagnetostaticVectorPotentialParsers
    {
        amrex::Parser ax;
        amrex::Parser ay;
        amrex::Parser daxdx;
        amrex::Parser daxdy;
        amrex::Parser daydx;
        amrex::Parser daydy;
        amrex::Parser dazdx;
        amrex::Parser dazdy;
    };

    struct MagnetostaticVectorPotential
    : public mixin::Named,
      public mixin::BeamOptic<MagnetostaticVectorPotential>,
      public mixin::LinearTransport<MagnetostaticVectorPotential>,
      public mixin::Thick,
      public mixin::Alignment,
      public mixin::PipeAperture
    {
        static constexpr auto type = "MagnetostaticVectorPotential";
        using PType = ImpactXParticleContainer::ParticleType;

        using DynamicData = mixin::GPUDataRegistry<MagnetostaticVectorPotentialParsers>;

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
         * @param name a user defined and not necessarily unique name of the element
         */
        MagnetostaticVectorPotential (
            amrex::ParticleReal ds,
            int unit,
            std::string const & str_ax = "0",
            std::string const & str_ay = "0",
            std::string const & str_daxdx = "0",
            std::string const & str_daxdy = "0",
            std::string const & str_daydx = "0",
            std::string const & str_daydy = "0",
            std::string const & str_dazdx = "0",
            std::string const & str_dazdy = "0",
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            amrex::ParticleReal aperture_x = 0,
            amrex::ParticleReal aperture_y = 0,
            int int_order = 2,
            int mapsteps = 10,
            int nslice = 1,
            std::optional<std::string> name = std::nullopt
        )
        : Named(std::move(name)),
          Thick(ds, nslice),
          Alignment(dx, dy, rotation_degree),
          PipeAperture(aperture_x, aperture_y),
          m_unit(unit), m_int_order(int_order), m_mapsteps(mapsteps),
          m_id(DynamicData::allocate_id())
        {
            // variables of the user-defined functions
            amrex::Vector<std::string> const varnames = {"x", "y", "z"};

            // build host parsers
            amrex::ParmParse const pp;
            auto & parsers = DynamicData::emplace(m_id);
            parsers.ax    = pp.makeParser(str_ax,    varnames);
            parsers.ay    = pp.makeParser(str_ay,    varnames);
            parsers.daxdx = pp.makeParser(str_daxdx, varnames);
            parsers.daxdy = pp.makeParser(str_daxdy, varnames);
            parsers.daydx = pp.makeParser(str_daydx, varnames);
            parsers.daydy = pp.makeParser(str_daydy, varnames);
            parsers.dazdx = pp.makeParser(str_dazdx, varnames);
            parsers.dazdy = pp.makeParser(str_dazdy, varnames);
        }

        /** Reverse the element (negate ds)
         *
         * The user-supplied potential A(x,y,z) is defined in the element's
         * local frame. Reversing only negates the path-length integration
         * direction; the user is responsible for the field's z-symmetry.
         */
        void reverse () { Thick::reverse(); }

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

            m_ibeta = 1_prt / refpart.beta();

            // normalize the vector potential if SI units are specified
            m_scale = m_unit == 1 ? 1_prt / refpart.rigidity_Tm() : 1_prt;

            // compile device-callable executors (cached across calls)
            constexpr int num_vars = 3;
            auto const & parsers = *DynamicData::get(m_id);
            m_dfunc_ax    = parsers.ax.compile<num_vars>();
            m_dfunc_ay    = parsers.ay.compile<num_vars>();
            m_dfunc_daxdx = parsers.daxdx.compile<num_vars>();
            m_dfunc_daxdy = parsers.daxdy.compile<num_vars>();
            m_dfunc_daydx = parsers.daydx.compile<num_vars>();
            m_dfunc_daydy = parsers.daydy.compile<num_vars>();
            m_dfunc_dazdx = parsers.dazdx.compile<num_vars>();
            m_dfunc_dazdy = parsers.dazdy.compile<num_vars>();
        }

        /** Push a single particle through this element.
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

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

            // coordinates within the element
            amrex::ParticleReal const s = refpart.s - m_slice_ds;
            amrex::ParticleReal const sedge = refpart.sedge;

            // numerical integration parameters
            amrex::ParticleReal const zin = s - sedge - m_ds/2_prt;
            amrex::ParticleReal const zout = zin + m_slice_ds;
            int const nsteps = m_mapsteps;

            // initialize phase space 6-vector
            amrex::SmallVector<amrex::ParticleReal, 6, 1> particle{
                x, px, y, py, t, pt
            };

            // call integrator to advance through slice (int_order = 2, 4 or 6, otherwise default 2)
            integrators::symp_integrate_semiexplicit(particle, zin, zout, m_int_order, nsteps, *this);

            // assign updated values
            x  = particle(1);
            px = particle(2);
            y  = particle(3);
            py = particle(4);
            t  = particle(5);
            pt = particle(6);

            // apply transverse aperture
            apply_aperture(x, y, idcpu);

            // undo shift due to alignment errors of the element
            shift_out(x, y, px, py);
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
         * @param t temporal coordinate in m
         * @param pt temporal momentum (energy deviation)
         * @param zeval lattice position in m
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

            // Evaluate the vector potential and its derivatives in the local frame
            auto const ax    = m_dfunc_ax(x, y, z)    * m_scale;
            auto const ay    = m_dfunc_ay(x, y, z)    * m_scale;
            auto const daxdx = m_dfunc_daxdx(x, y, z) * m_scale;
            auto const daxdy = m_dfunc_daxdy(x, y, z) * m_scale;
            auto const daydx = m_dfunc_daydx(x, y, z) * m_scale;
            auto const daydy = m_dfunc_daydy(x, y, z) * m_scale;
            auto const dazdx = m_dfunc_dazdx(x, y, z) * m_scale;
            auto const dazdy = m_dfunc_dazdy(x, y, z) * m_scale;

            // Rectilinear magnetostatic element: pz is the radical appearing
            // in the Hamiltonian.
            amrex::ParticleReal const pxmech = px - ax;
            amrex::ParticleReal const pymech = py - ay;
            amrex::ParticleReal const radical = std::sqrt(
                1_prt + pt*pt - px*px - py*py - 2_prt*pt*m_ibeta
            );

            amrex::ParticleReal const dHdx = -(pxmech * daxdx + pymech * daydx) / radical - dazdx;
            amrex::ParticleReal const dHdy = -(pxmech * daxdy + pymech * daydy) / radical - dazdy;
            amrex::ParticleReal const dHdt = 0_prt;

            amrex::ParticleReal const dHdpx = pxmech / radical;
            amrex::ParticleReal const dHdpy = pymech / radical;
            amrex::ParticleReal const dHdpt = -m_ibeta + (m_ibeta - pt) / radical;

            return std::make_tuple(dHdx, dHdpx, dHdy, dHdpy, dHdt, dHdpt);
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
            amrex::ParticleReal const step = slice_ds / std::sqrt(std::pow(pt, 2) - 1.0_prt);

            // advance position and momentum (straight element)
            refpart.x = x + step*px;
            refpart.y = y + step*py;
            refpart.z = z + step*pz;
            refpart.t = t - step*pt;

            // advance integrated path length
            refpart.s = s + slice_ds;
        }

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** Close and deallocate all parser handles.
         *
         * The amrex::Parser instances we own depend on AMReX runtime state, so
         * we must drop them before AMReX::Finalize() tears that state down.
         */
        void
        finalize ()
        {
            DynamicData::registry().erase(m_id);
        }

        /** This function returns the linear transport map.
         *
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map ([[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            throw std::runtime_error(std::string(type) + ": Envelope tracking is not yet implemented!");
            return Map6x6::Identity();
        }

        int m_unit; //! unit specification for the vector potential
        int m_int_order; //! order used for the symplectic integration (2, 4 or 6)
        int m_mapsteps; //! number of integration steps per slice
        int m_id; //! unique element id used for the parser data lookup

        // device-callable executors (cached after first compile_constants() call)
        amrex::ParserExecutor<3> m_dfunc_ax;     //! A_x(x,y,z)
        amrex::ParserExecutor<3> m_dfunc_ay;     //! A_y(x,y,z)
        amrex::ParserExecutor<3> m_dfunc_daxdx;  //! dA_x/dx(x,y,z)
        amrex::ParserExecutor<3> m_dfunc_daxdy;  //! dA_x/dy(x,y,z)
        amrex::ParserExecutor<3> m_dfunc_daydx;  //! dA_y/dx(x,y,z)
        amrex::ParserExecutor<3> m_dfunc_daydy;  //! dA_y/dy(x,y,z)
        amrex::ParserExecutor<3> m_dfunc_dazdx;  //! dA_z/dx(x,y,z)
        amrex::ParserExecutor<3> m_dfunc_dazdy;  //! dA_z/dy(x,y,z)

      private:
        // constants that are independent of the individually tracked particle,
        // see: compute_constants() to refresh
        amrex::ParticleReal m_slice_ds; //! m_ds / nslice()
        amrex::ParticleReal m_ibeta;    //! 1 / beta
        amrex::ParticleReal m_scale;    //! scale factor defining vector-potential units
    };

} // namespace impactx::elements

IMPACTX_GPUDATA_EXTERN(impactx::elements::MagnetostaticVectorPotential)
IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::MagnetostaticVectorPotential)

#endif // IMPACTX_MAGNETOSTATIC_VECTOR_POTENTIAL_H
