#ifndef FF_DIPEDGE_H
#define FF_DIPEDGE_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/utils/simple_timer.h"


// Calculations of dipedge taken from MAD-X SUBROUTINE tmfrng from
// file twiss.f90.

// Comments from that file:
/*!----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for fringe field of a dipole.                      *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     h         (double)  curvature of magnet body.                    *
  !     sk1       (double)  quadrupole strength in magnet body.          *
  !     edge      (double)  edge focussing angle.                        *
  !     he        (double)  curvature of pole face.                      *
  !     sig       (double)  sign: +1 for entry, -1 for exit.             *
  !     corr      (double)  correction factor according to slac 75.      *
  !     Output:                                                          *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second order terms.                          *
*/

/* MAD-X Documentation on dipedge:

label:
DIPEDGE, H=real, E1=real, FINT=real,
HGAP=real, TILT=real;

A DIPEDGE  has zero length and five attributes.
H  Is angle/length or 1/rho (default: 0 m−1 - for the default the dipedge element
has no effect). (must be equal to that of the associated SBEND)

E1 The rotation angle for the pole face. The sign convention is as for a SBEND
bending magnet. Note that it is different for an entrance and an exit. (default:
0 rad).

FINT  field integral as for SBEND sector bend. Note that each DIPEDGE has its own
FINT, so that specifying FINTX is no longer necessary.

HGAP  half gap height of the associated SBEND bending magnet.

TILT  The roll angle about the longitudinal axis (default: 0 rad, i.e. a horizontal
bend). A positive angle represents a clockwise rotation.
 */

/*
 * Synergia extension for dipedge
  *
  *  If second_order attribute is nonzero, then the dipedge propagation
  *  will take second order terms into account. Normally, this will default to 0
  *  because these second order terms are not symplectic, but they can be
  *  enabled for the calculation of chromaticity which is affected by second terms.
*/

// p [Gev/c] = -- * B*rho [ Tesla meters ]
#define PH_CNV_brho_to_p   (1.0e-9 * pconstants::c)

namespace dipedge_impl
{
    struct DipedgeParams
    {

        double e1;
        double h;
        double fint;
        double hgap;
        double tilt;

        double scale;

        double re_2_1;
        double re_4_3;
        double te[10];

        double pref_b;
        double m_b;
        double ref_cdt;

        double ce1, se1;
        double tanedge;
        double secedge;
        double sn;
        double psip;

    };

    template<class BP> // BP is bunch particles?
    struct PropDipedge
    {
        typename BP::parts_t p;
        typename BP::const_masks_t masks;
        const DipedgeParams dipedge_params;

        PropDipedge(BP & bp, DipedgeParams const& dipedge_params)
            : p(bp.parts), masks(bp.masks), dipedge_params(dipedge_params)
            { }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            if (masks(i))
            {

                // bend
                FF_algorithm::dipedge_unit(
                        p(i,0), p(i,1), p(i,2), p(i,3),
                        dipedge_params.re_2_1, dipedge_params.re_4_3,
                        dipedge_params.te );

           }
        }
    };

    template<class BP>
    struct PropDipedgeSimd
    {
        using gsv_t = typename BP::gsv_t;
        using parts_t = typename BP::parts_t;

        typename BP::parts_t p;
        typename BP::const_masks_t masks;
        const DipedgeParams dipedge_params;

        PropDipedgeSimd(BP & bp, DipedgeParams const& dipedge_params)
            : p(bp.parts), masks(bp.masks), dipedge_params(dipedge_params)
        { }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int idx) const
        {
            int i = idx * gsv_t::size();

            int m = 0;
            for(int x=i; x<i+gsv_t::size(); ++x) m |= masks(x);

            if (m)
            {
                gsv_t p0(&p(i, 0));
                gsv_t p1(&p(i, 1));
                gsv_t p2(&p(i, 2));
                gsv_t p3(&p(i, 3));

 
                // bend
                FF_algorithm::dipedge_unit<gsv_t>(
                        p0, p1, p2, p3,
                        dipedge_params.re_2_1,
                        dipedge_params.re_4_3,
                        dipedge_params.te );


                p0.store(&p(i, 0));
                p1.store(&p(i, 1));
                p2.store(&p(i, 2));
                p3.store(&p(i, 3));
            }
        }
    };

    inline void prop_reference(
            Reference_particle & ref_l, 
            DipedgeParams & dipedge_params )
    {
        double pref_l = ref_l.get_momentum();
        double    m_l = ref_l.get_mass();

        // propagate the reference particle, and set the edge kick strength 
        // from the reference particle
        double    x_l = ref_l.get_state()[Bunch::x];
        double   xp_l = ref_l.get_state()[Bunch::xp];
        double    y_l = ref_l.get_state()[Bunch::y];
        double   yp_l = ref_l.get_state()[Bunch::yp];
        double  cdt_l = 0.0;
        double dpop_l = ref_l.get_state()[Bunch::dpop];


        FF_algorithm::dipedge_unit(
                x_l, xp_l, y_l, yp_l,
                dipedge_params.re_2_1,
                dipedge_params.re_4_3,
                dipedge_params.te);
        
        ref_l.set_state(x_l, xp_l, y_l, yp_l,
            cdt_l, dpop_l);

    }
}

namespace FF_dipedge {

    template <class BunchT>
    inline void
    apply(Lattice_element_slice const& slice, BunchT& bunch)
    {
        using namespace dipedge_impl;

        scoped_simple_timer timer("libFF_dipedge");

        auto const& ele = slice.get_lattice_element();

        DipedgeParams dipedge_params;

        /*------------------------
           double tanedg =   tan(edge);
           double secedg =   1.0 / cos(edge);
           double sn     =   sin(edge);
           double psip   =   edge - corr * secedg * (1.0 + sn * sn);
           double re_2_1 =   h * tanedg;
           double re_4_3 = - h * tan(psip);
       //--------------------------*/

        const bool fsec = false;
        const double sig = 1.0;
        const double he = 0.0;
        const double sk1 = 0.0;

        const double edge = ele.get_double_attribute("e1", 0.0);
        const double h = ele.get_double_attribute("h", 0.0);
        const double fint = ele.get_double_attribute("fint", 0.0);
        const double hgap = ele.get_double_attribute("hgap", 0.0);
        const double tilt = ele.get_double_attribute("tilt", 0.0);
        const double corr = (h + h) * hgap * fint;

        double tanedg = tan(edge);
        double secedg = 1.0 / cos(edge);
        double sn = sin(edge);
        double psip = edge - corr * secedg * (1.0 + sn * sn);

        dipedge_params.re_2_1 = h * tanedg;
        dipedge_params.re_4_3 = -h * tan(psip);

        const bool second_order =
            (ele.get_double_attribute("second_order", 0.0) != 0.0);
        if (!second_order) {
            dipedge_params.te[0] = dipedge_params.te[1] = dipedge_params.te[2] =
                dipedge_params.te[3] = dipedge_params.te[4] = dipedge_params.te[5] =
                dipedge_params.te[6] = dipedge_params.te[7] = dipedge_params.te[8] =
                dipedge_params.te[9] = 0.0;
        } else {
            double hh = sig * (h / 2.0);

            dipedge_params.te[0] = -hh * tanedg * tanedg;
            dipedge_params.te[1] = hh * secedg * secedg;

            dipedge_params.te[2] =
                (h / 2.0) * he * (secedg * secedg * secedg) + sk1 * tanedg;
            dipedge_params.te[3] = -dipedge_params.te[0];
            dipedge_params.te[4] =
                hh * h * (tanedg * tanedg * tanedg) - dipedge_params.te[2];
            dipedge_params.te[5] = dipedge_params.te[0];

            dipedge_params.te[6] = -dipedge_params.te[0];

            dipedge_params.te[7] = -dipedge_params.te[2];
            dipedge_params.te[8] = dipedge_params.te[0];
            dipedge_params.te[9] = -dipedge_params.te[1];
            if (sig > 0) {
                dipedge_params.te[4] +=
                    (h * secedg) * (h * secedg) * (tanedg / 2.0);
            } else {
                dipedge_params.te[2] -=
                    (h * secedg) * (h * secedg) * (tanedg / 2.0);
                dipedge_params.te[7] +=
                    (h * secedg) * (h * secedg) * (tanedg / 2.0);
            }
        }

        Reference_particle& ref_l = bunch.get_design_reference_particle();
        Reference_particle const& ref_b = bunch.get_reference_particle();

        dipedge_params.scale =
            ref_l.get_momentum() /
            (ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]));

        // lattice reference
        double pref_l = ref_l.get_momentum();
        double brho_l = pref_l / PH_CNV_brho_to_p;
        double m_l = ref_l.get_mass();
        int charge_l = ref_l.get_charge();

        // the particle dpop is with respect to this momentum which goes
        // with the bunch
        double pref_b = ref_b.get_momentum();
        double brho_b = pref_b / PH_CNV_brho_to_p;
        double m_b = bunch.get_mass();

        // common
        dipedge_params.pref_b = pref_b;
        dipedge_params.m_b = m_b;

        dipedge_params.ce1 = cos(-edge);
        dipedge_params.se1 = sin(-edge);

        using namespace Kokkos;
        using exec = typename BunchT::exec_space;

        auto apply = [&](ParticleGroup pg) {
            auto bp = bunch.get_bunch_particles(pg);
            if (!bp.num_valid()) return;

            // dipedge is a 0 length element with no multipole moments
            // (currently)
#if LIBFF_USE_GSV
            prop_reference(ref_l, dipedge_params);

            auto range = RangePolicy<exec>(0, bp.size_in_gsv());
            PropDipedgeSimd<typename BunchT::bp_t> dipedge(bp, dipedge_params);
            Kokkos::parallel_for(range, dipedge);

#else  // LIBFF_USE_GSV
            auto range = RangePolicy<exec>(0, bp.size());
            PropDipedge<typename BunchT::bp_t> dipedge(bp, dipedge_params);
            Kokkos::parallel_for(range, dipedge);
#endif // LIBFF_USE_GSV
        };

        apply(ParticleGroup::regular);
        apply(ParticleGroup::spectator);
        // dipedge has 0 length so there is no need to increment trajectory here
        // bunch.get_reference_particle().increment_trajectory(0);
        // Also no need to increment bunch abs_time
        Kokkos::fence();
    }
}

#endif // FF_DIPEDGE_H

