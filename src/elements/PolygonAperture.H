/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Eric G. Stern, Chad Mitchell, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_POLYGONAPERTURE_H
#define IMPACTX_POLYGONAPERTURE_H


#include "particles/ImpactXParticleContainer.H"
#include "mixin/alignment.H"
#include "mixin/beamoptic.H"
#include "mixin/dynamicdata.H"
#include "mixin/lineartransport.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"
#include "mixin/spintransport.H"
#include "mixin/thin.H"
#include "mixin/TrackedVector.H"

#include <ablastr/constant.H>

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>

#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>


namespace impactx::elements
{

    /** Polygon vertices for aperture boundary definition.
     *
     * Host data persists across simulations; device data is uploaded
     * on-demand "lazily" on first use and released when AMReX shuts down.
     */
    struct PolygonVertices
    {
        mixin::TrackedVector<amrex::ParticleReal> x;
        mixin::TrackedVector<amrex::ParticleReal> y;
    };

    struct PolygonAperture
    : public mixin::Named,
      public mixin::BeamOptic<PolygonAperture>,
      public mixin::LinearTransport<PolygonAperture>,
      public mixin::Thin,
      public mixin::Alignment,
      public mixin::SpinTransport,
      public mixin::NoFinalize
      // public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "PolygonAperture";
        using PType = ImpactXParticleContainer::ParticleType;

        using DynamicData = mixin::GPUDataRegistry<PolygonVertices>;

        enum Action
        {
            transmit,
            absorb
        };

        static std::string
        action_name (Action const & action);

        /** A thin collimator element that applies a transverse aperture boundary defined
         *  by points of a polygon vertices.
         *  Particles outside the boundary are considered lost.
         *
         * @param action specify action of domain (transmit/absorb)
         * @param vertices_x array of coordinates of the horizontal vertex positions (m)
         * @param vertices_y array of coordinates of the vertical vertex positions (m)
         * @param min_radius2 a minimum radius**2 below which all particles are transmitted
         * @param repeat_x horizontal period for repeated masking, optional (m)
         * @param repeat_y vertical period for repeated masking, optional (m)
         * @param shift_odd_x for hexagonal/triangular mask patterns: horizontal shift of every 2nd (odd) vertical period by repeat_x / 2. Use alignment offsets dx,dy to move whole mask as needed.
         * @param dx horizontal translation error in m
         * @param dy vertical translation error in m
         * @param rotation_degree rotation error in the transverse plane [degrees]
         * @param name a user defined and not necessarily unique name of the element
         */
        PolygonAperture (
            std::vector< amrex::ParticleReal> vertices_x,
            std::vector< amrex::ParticleReal> vertices_y,
            amrex::ParticleReal min_radius2,
            amrex::ParticleReal repeat_x,
            amrex::ParticleReal repeat_y,
            bool shift_odd_x,
            Action action,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            std::optional<std::string> name = std::nullopt
        )
        : Named(std::move(name)),
          Alignment(dx, dy, rotation_degree),
          m_action(action), m_repeat_x(repeat_x), m_repeat_y(repeat_y), m_shift_odd_x(shift_odd_x),
          m_id(DynamicData::allocate_id()),
          m_min_radius2(min_radius2)
        {
            using namespace amrex::literals; // for _rt and _prt

            if (m_repeat_y == 0_prt && m_shift_odd_x) {
                throw std::runtime_error("PolygonAperture: shift_odd_x can only be used if repeat_y is set, too!");
            }

            m_nvert = vertices_x.size();
            if (m_nvert != vertices_y.size()) {
                throw std::runtime_error("PolygonAperture: horizontal and vertical vertices arrays must have same length!");
            }

            if ((vertices_x[0] != vertices_x[m_nvert-1]) ||
                (vertices_y[0] != vertices_y[m_nvert-1])) {
                throw std::runtime_error("PolygonAperture: first and last vertex must be exactly equal!");
            }

            auto& vert = DynamicData::emplace(
                m_id,
                std::move(vertices_x),
                std::move(vertices_y)
            );
            m_vertices_x_h_data = vert.x.host_const().data();
            m_vertices_y_h_data = vert.y.host_const().data();
        }

        /** Reverse the element (no-op: aperture is not invertible for absorbed particles) */
        void reverse () { /* no-op */ }

        /** Push all particles */
        using BeamOptic::operator();

        /** Compute and cache the constants for the push.
         *
         * In particular, used to pre-compute and cache variables that are
         * independent of the individually tracked particle.
         *
         * @param refpart reference particle (unused)
         */
        void compute_constants ([[maybe_unused]] RefPart const & refpart)
        {
            Alignment::compute_constants(refpart);

            // Ensure dynamic vertex data is on GPU and we have fresh pointers to it
            auto const & vert = *DynamicData::get(m_id);
            m_vertices_x_d_data = vert.x.device_const().data();
            m_vertices_y_d_data = vert.y.device_const().data();
        }

        /** This is an aperture functor, so that a variable of this type can be used like an
         *  aperture function.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t (unused)
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t (unused)
         * @param idcpu particle global index
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator() (
            T_Real const & AMREX_RESTRICT x,
            T_Real const & AMREX_RESTRICT y,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT t,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT px,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT py,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT pt,
            T_IdCpu & AMREX_RESTRICT idcpu,
            [[maybe_unused]] RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            // a complex type with two T_Real
            using Complex = amrex::GpuComplex<T_Real>;

            using namespace amrex::literals; // for _rt and _prt
            namespace stdx = amrex::simd::stdx;
            using amrex::Math::powi;
            using std::fmod;
            using std::abs;
            using amrex::arg;
            using VBool = decltype(x>0_prt);

            // define intermediate variables
            amrex::ParticleReal const dx = m_repeat_x * 0.5_prt;
            amrex::ParticleReal const dy = m_repeat_y * 0.5_prt;

            // shift every 2nd row by half of m_repeat_x
            T_Real sx = x;
            T_Real const sy = y;
            VBool const shift_x = m_shift_odd_x == false ? VBool{false} : fmod(floor((y+dy)/m_repeat_y), 2) != T_Real{0};
#ifdef AMREX_USE_SIMD
            amrex::simd::stdx::where(shift_x, sx) += dx;
#else
            sx += shift_x ? dx : 0_prt;
#endif

            // if the aperture is periodic, shift sx,sy coordinates to the fundamental domain
            T_Real const u = (m_repeat_x == 0_prt) ? sx : (fmod(abs(sx)+dx, m_repeat_x)-dx);
            T_Real const v = (m_repeat_y == 0_prt) ? sy : (fmod(abs(sy)+dy, m_repeat_y)-dy);

            // compare against the aperture boundary
            VBool inside_aperture;

            // calculate whether point is inside polygon or not
            T_Real const r2 = powi<2>(u) + powi<2>(v);
            // are we inside the minimum radius?
            if (r2 < m_min_radius2) {
                inside_aperture = VBool(true); // yes!
            } else{

#if AMREX_DEVICE_COMPILE
                amrex::ParticleReal const * vertices_x = m_vertices_x_d_data;
                amrex::ParticleReal const * vertices_y = m_vertices_y_d_data;
#else
                amrex::ParticleReal const * vertices_x = m_vertices_x_h_data;
                amrex::ParticleReal const * vertices_y = m_vertices_y_h_data;
#endif

                // no, we have to do the polygon calculation

                /*
                 * The algorithms works as follows:
                 *
                 * You have a polygon (in 2D) described by points V_1, ..., V_n arranged
                 * counter-clockwise. You have a particle with coordinates w = (x,y).
                 * Get the difference vector from the particle to each of
                 * the vertices in turn expressed as a complex number z_k, k=1,n. Now the angle
                 * between the particle and vertex i and vertex j is
                 *
                 * theta_{ij} = arg(\bar{z_i} z_j)
                 * This is easily demonstrated by expressing the complex numbers in exponential
                 * polar notation and also the angle is positive in the counter-clockwise direction.
                 *
                 * Take the sum of all the angles. The angle sum for particle inside the polygon will
                 * be 2 pi. The sum of angles of particle outside of the polygon will be 0.
                */
                T_Real thetasum = 0.0;
                Complex v1conj, v2;

                for (size_t i = 0; i < m_nvert-1; ++i) {
                    v1conj = Complex(u - vertices_x[i], -(v-vertices_y[i]));
                    v2 = Complex(u - vertices_x[i+1], v-vertices_y[i+1]);
                    thetasum += arg(v2 * v1conj);
                }
                inside_aperture = abs(thetasum/(2*ablastr::constant::math::pi)) > 1.0e-12_prt;
            }

            // For transmit (default): invalidate if outside aperture
            // For absorb: invalidate if inside aperture
            auto const invalid_mask = m_action == Action::transmit ? !inside_aperture : inside_aperture;
            amrex::ParticleIDWrapper<T_IdCpu>{idcpu}.make_invalid(invalid_mask);
        }

        /** This pushes the reference particle. */
        using Thin::operator();

        /** This operator pushes the particle spin.
         *
         * An aperture is field-free, so the spin is unaffected.
         *
         * @param x particle position in x
         * @param y particle position in y
         * @param t particle position in t
         * @param px particle momentum in x
         * @param py particle momentum in y
         * @param pt particle momentum in t
         * @param sx particle spin x-component (unused)
         * @param sy particle spin y-component (unused)
         * @param sz particle spin z-component (unused)
         * @param idcpu particle global index
         * @param refpart reference particle (unused)
         */
        template<typename T_Real=amrex::ParticleReal, typename T_IdCpu=uint64_t>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void spin_and_phasespace_push (
            T_Real const & AMREX_RESTRICT x,
            T_Real const & AMREX_RESTRICT y,
            T_Real const & AMREX_RESTRICT t,
            T_Real const & AMREX_RESTRICT px,
            T_Real const & AMREX_RESTRICT py,
            T_Real const & AMREX_RESTRICT pt,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT sx,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT sy,
            [[maybe_unused]] T_Real const & AMREX_RESTRICT sz,
            T_IdCpu & AMREX_RESTRICT idcpu,
            RefPart const & AMREX_RESTRICT refpart
        ) const
        {
            // spin is unaffected

            // phase space push
            (*this)(x, y, t, px, py, pt, idcpu, refpart);
        }

        /** This pushes the covariance matrix. */
        using LinearTransport::operator();

        /** This function returns the linear transport map.
         *
         * A polygon aperture (potentially) removes particles; it does not
         * deflect or drift the beam, so the linear transport map is the
         * identity.
         *
         * @param[in] refpart reference particle (unused)
         * @returns 6x6 transport matrix
         */
        AMREX_GPU_HOST AMREX_FORCE_INLINE
        Map6x6
        transport_map ([[maybe_unused]] RefPart const & AMREX_RESTRICT refpart) const
        {
            return Map6x6::Identity();
        }

        Action m_action; //! action type (transmit, absorb)
        amrex::ParticleReal m_repeat_x; //! horizontal period for repeated masking
        amrex::ParticleReal m_repeat_y; //! vertical period for repeated masking
        bool m_shift_odd_x; //! for hexagonal/triangular mask patterns: horizontal shift of every 2nd (odd) vertical period by m_repeat_x / 2. Use Alignment offsets to move whole mask as needed.

        int m_id; //! unique PolygonAperture id used for data lookup map

        amrex::ParticleReal m_min_radius2; // minimum radius in which all pass squared

        std::vector<double>::size_type m_nvert = 0; //! number of vertices
        amrex::ParticleReal const * m_vertices_x_h_data = nullptr; //! non-owning pointer to host x vertices
        amrex::ParticleReal const * m_vertices_y_h_data = nullptr; //! non-owning pointer to host vertices
        amrex::ParticleReal const * m_vertices_x_d_data = nullptr; //! non-owning pointer to device x vertices
        amrex::ParticleReal const * m_vertices_y_d_data = nullptr; //! non-owning pointer to device y vertices

    };

} // namespace impactx

IMPACTX_GPUDATA_EXTERN(impactx::elements::PolygonAperture)
IMPACTX_PUSH_EXTERN_TEMPLATE(impactx::elements::PolygonAperture)

#endif // IMPACTX_POLYGONAPERTURE_H
