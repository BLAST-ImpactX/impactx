/* Copyright 2022-2026 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_ELEMENTS_MIXIN_BEND_GEOMETRY_H
#define IMPACTX_ELEMENTS_MIXIN_BEND_GEOMETRY_H

#include "particles/ImpactXParticleContainer.H"
#include "thick.H"

#include <ablastr/constant.H>

#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>


namespace impactx::elements::mixin
{
    /** This is a helper class for elements that bend the reference orbit: the parameters
     *  that specify the bend
     *
     * A thick bend (@see ThickBend) is specified by its segment length ``ds`` plus exactly one of
     *
     * - ``rc``: the radius of curvature in m,
     * - ``phi``: the bend angle in degrees, so that ``rc = ds / phi``,
     * - ``B``: the magnetic field in T, so that ``rc = rigidity / B``,
     *
     * or by ``phi`` together with ``B``: a bend of fixed angle whose radius follows from the
     * field. Its arc length ``rc * phi`` then depends on the beam energy, and ``ds`` only
     * advances the integrated path length ``s``.
     *
     * A combined-function bend can also leave all of them unset; its dipole field then comes
     * from its multipole coefficients. A thin bend is specified by its bend angle, called
     * ``theta``, together with ``rc`` or ``B``. The edge of a bend is specified by ``rc`` or ``B``.
     *
     * A radius of ``rc = 0`` or a field of ``B = 0`` is an infinite radius, and a bend with an
     * infinite radius or a bend angle of 0 is straight: its reference orbit is not bent and
     * ``signed_rc`` is infinite. Because the rigidity carries the sign of the charge, a bend
     * specified by ``B`` bends particles of opposite charge the opposite way.
     *
     * The center of curvature lies at ``x = -signed_rc``.
     */
    struct BendGeometry
    {
        /** Which kind of bend element this is, which sets how it can be specified */
        enum class Kind : int
        {
            thick,    //!< ds plus one of rc, phi or B, or phi together with B
            combined, //!< as thick, or none of them (the dipole field is a multipole coefficient)
            thin,     //!< the bend angle theta together with rc or B
            edge      //!< one of rc or B
        };

        /** How the geometry of the bend is specified */
        enum class Spec : int
        {
            rc,     //!< radius of curvature
            phi,    //!< bend angle
            B,      //!< magnetic field
            phi_B,  //!< bend angle and magnetic field
            phi_rc, //!< bend angle and radius of curvature (thin bends)
            none    //!< not specified here (combined-function bends)
        };

        /** The name of the bend angle parameter: ``theta`` for thin bends, ``phi`` otherwise
         *
         * @param kind the kind of bend
         */
        static constexpr char const *
        angle_name (Kind kind)
        {
            return kind == Kind::thin ? "theta" : "phi";
        }

        /** Is this the legacy way of leaving ``B`` unset: a thick bend with an angle and ``B = 0``?
         *
         * Callers accept it as a bend angle alone and emit a deprecation warning in their
         * own channel (Python warning, input-file warning).
         *
         * @param kind the kind of bend
         * @param rc radius of curvature in m, if given
         * @param angle_degree bend angle in degrees, if given
         * @param B magnetic field in T, if given
         */
        static bool
        is_legacy_unset_B (
            Kind kind,
            std::optional<amrex::ParticleReal> rc,
            std::optional<amrex::ParticleReal> angle_degree,
            std::optional<amrex::ParticleReal> B
        )
        {
            return kind == Kind::thick &&
                   !rc.has_value() && angle_degree.has_value() && B.has_value() && *B == 0;
        }

        /** The parameters that specify a bend
         *
         * @param kind the kind of bend
         * @param element_type the element type, used in error messages
         * @param rc radius of curvature in m, if given
         * @param angle_degree bend angle in degrees, if given
         * @param B magnetic field in T, if given
         * @throws std::invalid_argument if the combination does not specify one geometry
         */
        BendGeometry (
            Kind kind,
            char const * element_type,
            std::optional<amrex::ParticleReal> rc,
            std::optional<amrex::ParticleReal> angle_degree,
            std::optional<amrex::ParticleReal> B
        )
        : m_kind(kind)
        {
            set_geometry(element_type, rc, angle_degree, B);
        }

        /** Change how the geometry of the bend is specified
         *
         * @param element_type the element type, used in error messages
         * @param rc radius of curvature in m, if given
         * @param angle_degree bend angle in degrees, if given
         * @param B magnetic field in T, if given
         * @throws std::invalid_argument if the combination does not specify one geometry;
         *         the element is left unchanged then
         */
        void set_geometry (
            char const * element_type,
            std::optional<amrex::ParticleReal> rc,
            std::optional<amrex::ParticleReal> angle_degree,
            std::optional<amrex::ParticleReal> B
        )
        {
            using namespace amrex::literals; // for _rt and _prt

            std::string const type(element_type);
            std::string const angle(angle_name(m_kind));
            bool const has_rc = rc.has_value();
            bool const has_angle = angle_degree.has_value();
            // the legacy B = 0 with an angle leaves B unset
            bool const has_B = B.has_value() && !is_legacy_unset_B(m_kind, rc, angle_degree, B);

            if (has_rc && has_B) {
                throw std::invalid_argument(type + ": rc cannot be combined with B. " + allowed());
            }

            switch (m_kind) {
                case Kind::thick:
                case Kind::combined:
                    if (has_rc && has_angle) {
                        throw std::invalid_argument(type + ": rc cannot be combined with phi. " + allowed());
                    }
                    if (!has_rc && !has_angle && !has_B) {
                        if (m_kind == Kind::combined) { m_spec = Spec::none; break; }
                        throw std::invalid_argument(type + ": the bend needs " + allowed(false));
                    }
                    m_spec = has_rc ? Spec::rc
                           : has_angle && has_B ? Spec::phi_B
                           : has_angle ? Spec::phi
                           : Spec::B;
                    break;
                case Kind::thin:
                    if (!has_angle || (!has_rc && !has_B)) {
                        throw std::invalid_argument(type + ": the bend needs " + allowed(false));
                    }
                    m_spec = has_rc ? Spec::phi_rc : Spec::phi_B;
                    break;
                case Kind::edge:
                    if (has_angle) {
                        throw std::invalid_argument(type + ": the edge has no " + angle + ". " + allowed());
                    }
                    if (!has_rc && !has_B) {
                        throw std::invalid_argument(type + ": the edge needs " + allowed(false));
                    }
                    m_spec = has_rc ? Spec::rc : Spec::B;
                    break;
            }

            m_rc = has_rc ? *rc : 0_prt;
            m_phi = (has_angle ? *angle_degree : 0_prt) * bend_degree2rad;
            m_B = has_B ? *B : 0_prt;
        }

        /** The radius of curvature in m, if the bend is specified by it */
        [[nodiscard]] std::optional<amrex::ParticleReal>
        rc_parameter () const
        {
            if (m_spec == Spec::rc || m_spec == Spec::phi_rc) { return m_rc; }
            return std::nullopt;
        }

        /** The bend angle in degrees, if the bend is specified by it */
        [[nodiscard]] std::optional<amrex::ParticleReal>
        angle_parameter () const
        {
            if (m_spec == Spec::phi || m_spec == Spec::phi_B || m_spec == Spec::phi_rc) {
                return m_phi / bend_degree2rad;
            }
            return std::nullopt;
        }

        /** The magnetic field in T, if the bend is specified by it */
        [[nodiscard]] std::optional<amrex::ParticleReal>
        B_parameter () const
        {
            if (m_spec == Spec::B || m_spec == Spec::phi_B) { return m_B; }
            return std::nullopt;
        }

        /** Change one parameter of the bend, keeping how its geometry is specified
         *
         * @param element_type the element type, used in error messages
         * @param name the parameter: "rc", the angle (@see angle_name, in degrees) or "B"
         * @param value the new value
         * @throws std::invalid_argument if the bend is not specified by this parameter;
         *         @see set_geometry changes how it is specified
         */
        void set_parameter (
            char const * element_type,
            std::string const & name,
            amrex::ParticleReal value
        )
        {
            bool const is_rc = name == "rc" && rc_parameter().has_value();
            bool const is_angle = name == angle_name(m_kind) && angle_parameter().has_value();
            bool const is_B = name == "B" && B_parameter().has_value();
            if (!is_rc && !is_angle && !is_B) {
                throw std::invalid_argument(
                    std::string(element_type) + ": this bend is not specified by " + name +
                    ". Use set_geometry() to change how the bend is specified.");
            }

            if (is_rc) { m_rc = value; }
            if (is_angle) { m_phi = value * bend_degree2rad; }
            if (is_B) { m_B = value; }
        }

        /** Is the geometry specified by rc, an angle or B? False only for a combined-function
         *  bend whose dipole field is a multipole coefficient.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool
        has_geometry () const
        {
            return m_spec != Spec::none;
        }

        /** Is the reference orbit a straight line? (false if not specified here) */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool
        is_straight () const
        {
            using namespace amrex::literals; // for _rt and _prt

            switch (m_spec) {
                case Spec::rc: return m_rc == 0_prt;
                case Spec::phi: return m_phi == 0_prt;
                case Spec::B: return m_B == 0_prt;
                case Spec::phi_rc: return m_phi == 0_prt || m_rc == 0_prt;
                case Spec::none: return false;
                default: return m_phi == 0_prt || m_B == 0_prt;  // Spec::phi_B
            }
        }

        /** Signed radius of curvature of the reference orbit in m, from rc or B
         *
         * A straight bend returns +infinity, as does a bend not specified here. The radius of a
         * thick bend specified by its angle alone depends on its length, @see ThickBend::signed_rc
         *
         * @param refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal
        signed_rc (RefPart const & refpart) const
        {
            if (is_straight() || m_spec == Spec::phi || m_spec == Spec::none) {
                return std::numeric_limits<amrex::ParticleReal>::infinity();
            }
            if (m_spec == Spec::rc || m_spec == Spec::phi_rc) { return m_rc; }
            return refpart.rigidity_Tm() / m_B;  // Spec::B, Spec::phi_B
        }

        /** Reverse the bend angle
         *
         * Together with the reversal of the element's length, this reverses the bend.
         */
        void reverse ()
        {
            m_phi = -m_phi;
        }

      protected:
        //! not named degree2rad, which the bends inherit from Alignment
        static constexpr amrex::ParticleReal bend_degree2rad =
            ablastr::constant::math::pi_v<amrex::ParticleReal> / amrex::ParticleReal{180};

        Kind m_kind;               //! the kind of bend, which sets how it can be specified
        Spec m_spec;               //! how the geometry is specified
        amrex::ParticleReal m_rc;  //! radius of curvature in m (Spec::rc, Spec::phi_rc)
        amrex::ParticleReal m_phi; //! bend angle in rad (Spec::phi, Spec::phi_B, Spec::phi_rc)
        amrex::ParticleReal m_B;   //! magnetic field in T (Spec::B, Spec::phi_B)

      private:
        /** The allowed combinations of parameters, for error messages
         *
         * @param prefix start with "Specify the bend by"
         */
        [[nodiscard]] std::string
        allowed (bool prefix = true) const
        {
            std::string const by = prefix ? "Specify the bend by " : "";
            switch (m_kind) {
                case Kind::thin: return by + "theta together with one of rc or B.";
                case Kind::edge: return by + "one of rc or B.";
                case Kind::combined:
                    return by + "exactly one of rc, phi or B, or phi together with B, "
                           "or by none of them.";
                default: return by + "exactly one of rc, phi or B, or phi together with B.";
            }
        }
    };

    /** This is a helper class for thick elements that bend the reference orbit
     *
     * It derives the geometry of the reference orbit from the parameters that specify the
     * bend, @see BendGeometry, and from the element's length and number of slices.
     *
     * The center of curvature lies at ``x = -signed_rc``.
     */
    struct ThickBend
    : public Thick,
      public BendGeometry
    {
        static constexpr Kind bend_kind = Kind::thick;

        /** A thick element that bends the reference orbit
         *
         * @param element_type the element type, used in error messages
         * @param ds Segment length in m.
         * @param nslice number of slices used for the application of space charge
         * @param rc radius of curvature in m, if given
         * @param phi_degree bend angle in degrees, if given
         * @param B magnetic field in T, if given
         * @param kind thick, or combined for a bend whose dipole field may be a multipole
         * @throws std::invalid_argument if the combination does not specify one geometry
         */
        ThickBend (
            char const * element_type,
            amrex::ParticleReal ds,
            int nslice,
            std::optional<amrex::ParticleReal> rc,
            std::optional<amrex::ParticleReal> phi_degree,
            std::optional<amrex::ParticleReal> B,
            Kind kind = Kind::thick
        )
        : Thick(ds, nslice),
          BendGeometry(kind, element_type, rc, phi_degree, B)
        {
        }

        /** Signed radius of curvature of the reference orbit in m
         *
         * The center of curvature lies at ``x = -signed_rc``; a straight bend returns +infinity.
         *
         * @param refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal
        signed_rc (RefPart const & refpart) const
        {
            if (m_spec == Spec::phi && !is_straight()) { return m_ds / m_phi; }
            return BendGeometry::signed_rc(refpart);
        }

        /** Signed bend angle of the reference orbit in one slice, in rad
         *
         * @param refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal
        slice_angle (RefPart const & refpart) const
        {
            using namespace amrex::literals; // for _rt and _prt

            if (is_straight()) { return 0_prt; }

            switch (m_spec) {
                case Spec::phi: return m_phi / nslice();
                case Spec::phi_B: {
                    // the angle is fixed; its direction follows the field (and charge) and ds
                    amrex::ParticleReal const ds = m_ds;
                    amrex::ParticleReal const rc = signed_rc(refpart);
                    amrex::ParticleReal const bend_sign = ds * rc / std::abs(ds * rc);
                    return bend_sign * std::abs(m_phi) / nslice();
                }
                default: return (m_ds / nslice()) / signed_rc(refpart);  // Spec::rc, Spec::B
            }
        }

        /** Arc length of the reference orbit in one slice, in m
         *
         * This equals ``ds / nslice``, except for a bend specified by ``phi`` together
         * with ``B``, whose arc length follows from its angle and radius.
         *
         * @param refpart reference particle
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal
        slice_arc_length (RefPart const & refpart) const
        {
            if (m_spec == Spec::phi_B && !is_straight()) {
                return signed_rc(refpart) * slice_angle(refpart);
            }
            return m_ds / nslice();
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

            if (is_straight()) {
               // advance position and momentum (drift)
               amrex::ParticleReal const step = slice_ds / std::sqrt(powi<2>(pt) - 1.0_prt);
               refpart.x = x + step*px;
               refpart.y = y + step*py;
               refpart.z = z + step*pz;
               refpart.t = t - step*pt;

            } else {

               // assign intermediate parameters
               amrex::ParticleReal const rc = signed_rc(refpart);
               amrex::ParticleReal const theta = slice_angle(refpart);
               amrex::ParticleReal const B = refpart.beta_gamma() / rc;

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

        /** Reverse the element (negate ds and the bend angle) */
        void reverse ()
        {
            Thick::reverse();
            BendGeometry::reverse();
        }
    };

} // namespace impactx::elements::mixin

#endif // IMPACTX_ELEMENTS_MIXIN_BEND_GEOMETRY_H
