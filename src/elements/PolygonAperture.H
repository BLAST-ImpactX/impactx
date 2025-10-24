/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
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
#include "mixin/thin.H"
#include "mixin/lineartransport.H"
#include "mixin/named.H"
#include "mixin/nofinalize.H"

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_SIMD.H>
#include <ablastr/constant.H>
#include <cmath>
#include <stdexcept>


namespace impactx::elements
{

    /** Dynamic data for the PolygonAperture elements
     *
     * Since we copy the element to the device, we cannot store this data on the element itself.
     * But we can store pointers to this data with the element and keep a lookup table here,
     * which we clean up in the end.
     */
    namespace PolygonApertureData
    {
        //! last used id for a created ExactMultipole
        inline int next_id = 0;

        //! host: normal multipole coefficients of the magnetic field
        inline std::map<int, std::vector<amrex::ParticleReal>> h_vertices_x = {};
        //! host: skew multipole coefficients of the magnetic field
        inline std::map<int, std::vector<amrex::ParticleReal>> h_vertices_y = {};

        //! device: normal multipole coefficients of the magnetic field
        inline std::map<int, amrex::Gpu::DeviceVector<amrex::ParticleReal>> d_vertices_x = {};
        //! device: skew multipole coefficients of the magnetic field
        inline std::map<int, amrex::Gpu::DeviceVector<amrex::ParticleReal>> d_vertices_y = {};

    } // namespace PolygonApertureData

    struct PolygonAperture
    : public mixin::Named,
      public mixin::BeamOptic<PolygonAperture>,
      public mixin::LinearTransport<PolygonAperture>,
      public mixin::Thin,
      public mixin::Alignment,
      public mixin::NoFinalize,
      public amrex::simd::Vectorized<amrex::simd::native_simd_size_particlereal>
    {
        static constexpr auto type = "PolygonAperture";
        using PType = ImpactXParticleContainer::ParticleType;

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
            std::vector< amrex::ParticleReal> vertices_xx,
            std::vector< amrex::ParticleReal> vertices_yy,
            amrex::ParticleReal min_radius2x,
            amrex::ParticleReal repeat_x,
            amrex::ParticleReal repeat_y,
            bool shift_odd_x,
            Action action,
            amrex::ParticleReal dx = 0,
            amrex::ParticleReal dy = 0,
            amrex::ParticleReal rotation_degree = 0,
            std::optional<std::string> name = std::nullopt
        )
    //
        : Named(std::move(name)),
          Alignment(dx, dy, rotation_degree),
          m_action(action), m_repeat_x(repeat_x), m_repeat_y(repeat_y), m_shift_odd_x(shift_odd_x),
      vertices_x(vertices_xx), vertices_y(vertices_yy), min_radius2(min_radius2x)
        {
            using namespace amrex::literals; // for _rt and _prt

            if (m_repeat_y == 0_prt && m_shift_odd_x) {
                throw std::runtime_error("Aperture: shift_odd_x can only be used if repeat_y is set, too!");
            }

            // next created MulPolygonAperture has another id for its data
            PolygonApertureData::next_id++;

            m_nvert = vertices_x.size();
            // there must be the same number of horizontal and vertical vertices
            if (m_nvert != vertices_y.size()) {
                throw std::runtime_error("PolygonAperture: horizontal and vertical vertices arrays must have same length!");
            }

            // host data
            PolygonApertureData::h_vertices_x[m_id] = vertices_x;
            PolygonApertureData::h_vertices_y[m_id] = vertices_y;
            m_vertices_x_h_data = PolygonApertureData::h_vertices_x[m_id].data();
            m_vertices_y_h_data = PolygonApertureData::h_vertices_y[m_id].data();

            // device data
            PolygonApertureData::d_vertices_x.emplace(m_id, amrex::Gpu::DeviceVector<amrex::ParticleReal>(m_nvert));
            PolygonApertureData::d_vertices_y.emplace(m_id, amrex::Gpu::DeviceVector<amrex::ParticleReal>(m_nvert));
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  vertices_x.begin(), vertices_x.end(),
                                  PolygonApertureData::d_vertices_x[m_id].begin());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  vertices_y.begin(), vertices_y.end(),
                                  PolygonApertureData::d_vertices_y[m_id].begin());
            amrex::Gpu::streamSynchronize();

            // low-level objects we can use on device
            m_vertices_x_d_data = PolygonApertureData::d_vertices_x[m_id].data();
            m_vertices_y_d_data = PolygonApertureData::d_vertices_y[m_id].data();

        }

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
            T_Real & AMREX_RESTRICT x,
            T_Real & AMREX_RESTRICT y,
            [[maybe_unused]] T_Real & AMREX_RESTRICT t,
            T_Real & AMREX_RESTRICT px,
            T_Real & AMREX_RESTRICT py,
            [[maybe_unused]] T_Real & AMREX_RESTRICT pt,
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

            // shift due to alignment errors of the element
            shift_in(x, y, px, py);

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
            T_Real u = (m_repeat_x == 0_prt) ? sx : (fmod(abs(sx)+dx, m_repeat_x)-dx);
            T_Real v = (m_repeat_y == 0_prt) ? sy : (fmod(abs(sy)+dy, m_repeat_y)-dy);

            // compare against the aperture boundary
            VBool inside_aperture;

            // calculate whether point is inside polygon or not
            T_Real r2 = powi<2>(u) + powi<2>(v);
            // are we inside the minimum radius?
            if (r2 < min_radius2) {
                inside_aperture = VBool(true); // yes!
            } else{
                // no, we have to do the polygon calculation
                T_Real thetasum = 0.0;
                Complex v1conj, v2;
                // vertices_x and vertices_y must be the same size
                for (size_t i = 0; i <vertices_x.size()-1; ++i) {
                    v1conj = Complex(u - vertices_x[i], -(v-vertices_y[i]));
                    v2 = Complex(u - vertices_x[i+1], v-vertices_y[i+1]);
                    thetasum += arg(v2 * v1conj);
                }
                inside_aperture = abs(thetasum/(2*ablastr::constant::math::pi)) < 1.0e-12_prt;
            }

            // For transmit (default): invalidate if outside aperture
            // For absorb: invalidate if inside aperture
            auto const invalid_mask = m_action == Action::transmit ? !inside_aperture : inside_aperture;
            amrex::ParticleIDWrapper<T_IdCpu>{idcpu}.make_invalid(invalid_mask);

            // undo shift due to alignment errors of the element
            shift_out(x, y, px, py);
        }

        /** This pushes the reference particle. */
        using Thin::operator();

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

        ///////////////////////////////////////////////////////////////
        // class data

        Action m_action; //! action type (transmit, absorb)
        amrex::ParticleReal m_repeat_x; //! horizontal period for repeated masking
        amrex::ParticleReal m_repeat_y; //! vertical period for repeated masking
        bool m_shift_odd_x; //! for hexagonal/triangular mask patterns: horizontal shift of every 2nd (odd) vertical period by m_repeat_x / 2. Use Alignment offsets to move whole mask as needed.

        int m_id; //! unique PolygonAperture id used for data lookup map

        std::vector< amrex::ParticleReal> vertices_x;
        std::vector< amrex::ParticleReal> vertices_y;
        amrex::ParticleReal min_radius2; // minimum radius in which all pass squared

        std::vector<double>::size_type m_nvert = 0; //! number of vertices
        amrex::ParticleReal* m_vertices_x_h_data = nullptr; //! non-owning pointer to host cosine coefficients
        amrex::ParticleReal* m_vertices_y_h_data = nullptr; //! non-owning pointer to host sine coefficients
        amrex::ParticleReal* m_vertices_x_d_data = nullptr; //! non-owning pointer to device cosine coefficients
        amrex::ParticleReal* m_vertices_y_d_data = nullptr; //! non-owning pointer to device sine coefficients

    };

} // namespace impactx

#endif // IMPACTX_POLYGONAPERTURE_H
