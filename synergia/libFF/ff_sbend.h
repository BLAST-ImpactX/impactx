#ifndef FF_SBEND_H
#define FF_SBEND_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/foundation/physical_constants.h"


// p [Gev/c] = -- * B*rho [ Tesla meters ]
#define PH_CNV_brho_to_p   (1.0e-9 * pconstants::c)


namespace sbend_impl
{
    struct SbendParams
    {
        bool ledge;
        bool redge;
        bool ledgekick;
        bool redgekick;

        double length;
        double angle;
        double r0;

        int    steps;

        double k_l[10];  // up to decapole
        double scale;

        double strength;
        double pref_b;   // bunch reference momentum
        double m_b;
        double ref_cdt;

        double e1;  // usFaceAngle,  usAngle
        double e2;  // dsFaceAngle, -dsAngle

        double ce1;
        double se1;
        double ce2;
        double se2;

        double us_edge_k;
        double ds_edge_k;

        double us_edge_k_p;
        double ds_edge_k_p;

        double us_edge_k_x;
        double us_edge_k_y;

        double ds_edge_k_x;
        double ds_edge_k_y;
    };

    template<class BP>
    struct PropSbend
    {
        typename BP::parts_t p;
        typename BP::const_masks_t masks;
        const SbendParams sp;

        const double dphi;
        const Kokkos::complex<double> phase;
        const Kokkos::complex<double> term;

        PropSbend(BP & bp, SbendParams const& sp)
            : p(bp.parts), masks(bp.masks), sp(sp)
            , dphi( -(sp.angle - (sp.e1 + sp.e2)) )   // -psi
            , phase(std::exp(std::complex<double>(0.0, -dphi)))
            , term( std::complex<double>(0.0, sp.length / sp.angle) *
                    std::complex<double>(1.0 - cos(sp.angle), - sin(sp.angle)) *
                    std::complex<double>(cos(sp.e2), -sin(sp.e2)) )
        { }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            if (masks(i))
            {
                if (sp.ledge)
                {
                    // slot
                    FF_algorithm::slot_unit( 
                            p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                            sp.ce1, sp.se1, sp.pref_b, sp.m_b );

                    // edge:
                    // 1. chef fixed angle (only kicks yp)
                    // FF_algorithm::edge_unit(
                    //        p(i,2), p(i,3), us_edge_k );
                    //
                    // 2. chef per-particle angle
                    // FF_algorithm::edge_unit(
                    //        p(i,2), p(i,1), p(i,3), dpop, us_edge_k_p);
                    //
                    // 3. ref particle angle (kicks both xp and yp)
                    if (sp.ledgekick)
                    {
                        FF_algorithm::edge_unit(
                            p(i,2), p(i,1), p(i,3), 
                            sp.us_edge_k_x, sp.us_edge_k_y, 0); 
                    } 
                }

                // bend
                FF_algorithm::bend_unit(
                        p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                        dphi, sp.strength, sp.pref_b, sp.m_b, sp.ref_cdt, 
                        phase, term );

                if (sp.redge)
                {
                    // edge
                    // FF_algorithm::edge_unit(y, yp, ds_edge_k);
                    // FF_algorithm::edge_unit(y, xp, yp, dpop, ds_edge_k_p);
                    if (sp.redgekick)
                    {
                        FF_algorithm::edge_unit(
                            p(i,2), p(i,1), p(i,3), 
                            sp.ds_edge_k_x, sp.ds_edge_k_y, 0);
                    }

                    // slot
                    FF_algorithm::slot_unit(
                            p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                            sp.ce2, sp.se2, sp.pref_b, sp.m_b);
                }
            }
        }
    };

    template<class BP>
    struct PropSbendSimd
    {
        using gsv_t = typename BP::gsv_t;
        using parts_t = typename BP::parts_t;

        typename BP::parts_t p;
        typename BP::const_masks_t masks;
        const SbendParams sp;

        const double dphi;
        const Kokkos::complex<double> phase;
        const Kokkos::complex<double> term;

        PropSbendSimd(BP & bp, SbendParams const& sp)
            : p(bp.parts), masks(bp.masks), sp(sp)
            , dphi( -(sp.angle - (sp.e1 + sp.e2)) )   // -psi
            , phase(std::exp(std::complex<double>(0.0, -dphi)))
            , term( std::complex<double>(0.0, sp.length / sp.angle) *
                    std::complex<double>(1.0 - cos(sp.angle), - sin(sp.angle)) *
                    std::complex<double>(cos(sp.e2), -sin(sp.e2)) )
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
                gsv_t p4(&p(i, 4));
                gsv_t p5(&p(i, 5));

                if (sp.ledge)
                {
                    // slot
                    FF_algorithm::slot_unit<gsv_t>( 
                            p0, p1, p2, p3, p4, p5,
                            sp.ce1, sp.se1, sp.pref_b, sp.m_b );

                    // edge:
                    // 1. chef fixed angle (only kicks yp)
                    // FF_algorithm::edge_unit<gsv_t>(
                    //        p2, p3, sp.us_edge_k );
                    //
                    // 2. chef per-particle angle
                    // FF_algorithm::edge_unit<gsv_t>(
                    //        p2, p1, p3, p5, sp.us_edge_k_p);
                    //
                    // 3. ref particle angle (kicks both xp and yp)
                    if (sp.ledgekick)
                    {
                        FF_algorithm::edge_unit<gsv_t>(p2, p1, p3, sp.us_edge_k_x, sp.us_edge_k_y, 0);
                    }
                }

                // bend
                FF_algorithm::bend_unit<gsv_t>(
                        p0, p1, p2, p3, p4, p5,
                        dphi, sp.strength, sp.pref_b, sp.m_b, sp.ref_cdt, 
                        phase, term );

                if (sp.redge)
                {
                    // edge
                    // FF_algorithm::edge_unit<gsv_t>(p2, p3, sp.ds_edge_k);
                    // FF_algorithm::edge_unit<gsv_t>(p2, p1, p3, p5, sp.ds_edge_k_p);
                    if (sp.redgekick)
                    {
                        FF_algorithm::edge_unit<gsv_t>( 
                            p2, p1, p3, sp.ds_edge_k_x, sp.ds_edge_k_y, 0); 
		            }

                    // slot
                    FF_algorithm::slot_unit<gsv_t>(
                            p0, p1, p2, p3, p4, p5,
                            sp.ce2, sp.se2, sp.pref_b, sp.m_b);
                }

                p0.store(&p(i, 0));
                p1.store(&p(i, 1));
                p2.store(&p(i, 2));
                p3.store(&p(i, 3));
                p4.store(&p(i, 4));
                p5.store(&p(i, 5));
            }
        }
    };

    namespace fa  = ::FF_algorithm;
    namespace by6 = ::FF_algorithm::yoshida_constants::bend_yoshida6;

    template<class BP>
    struct PropSbendCF
    {
        using gsv_t = typename BP::gsv_t;
        using parts_t = typename BP::parts_t;
        using const_masks_t = typename BP::const_masks_t;

        parts_t p;
        const_masks_t masks;
        const SbendParams sp;

        const double step_angle;
        const double step_ref_cdt;

        const double step_kl[10];

        const Kokkos::complex<double> phase_e1;
        const Kokkos::complex<double> phase_e2;

        const Kokkos::complex<double> step_phase[4];
        const Kokkos::complex<double> step_term[4];
        const double step_dphi[4];

        PropSbendCF( BP & bp, SbendParams const& sp )
            : p(bp.parts), masks(bp.masks), sp(sp)

            , step_angle(sp.angle/sp.steps)
            , step_ref_cdt(sp.ref_cdt/sp.steps)

            , step_kl{  // k_b[i] = k_l[i] * scale
                sp.k_l[0] * sp.scale * sp.length / sp.steps,
                sp.k_l[1] * sp.scale * sp.length / sp.steps,
                sp.k_l[2] * sp.scale * sp.length / sp.steps,
                sp.k_l[3] * sp.scale * sp.length / sp.steps,
                sp.k_l[4] * sp.scale * sp.length / sp.steps,
                sp.k_l[5] * sp.scale * sp.length / sp.steps,
                sp.k_l[6] * sp.scale * sp.length / sp.steps,
                sp.k_l[7] * sp.scale * sp.length / sp.steps,
                sp.k_l[8] * sp.scale * sp.length / sp.steps,
                sp.k_l[9] * sp.scale * sp.length / sp.steps }

            , phase_e1(fa::bend_edge_phase(sp.e1))
            , phase_e2(fa::bend_edge_phase(sp.e2))

            , step_phase{
                fa::sbend_unit_phase(by6::c1, step_angle),
                fa::sbend_unit_phase(by6::c2, step_angle),
                fa::sbend_unit_phase(by6::c3, step_angle),
                fa::sbend_unit_phase(by6::c4, step_angle) }

            , step_term{
                fa::sbend_unit_term(by6::c1, step_angle, sp.r0),
                fa::sbend_unit_term(by6::c2, step_angle, sp.r0),
                fa::sbend_unit_term(by6::c3, step_angle, sp.r0),
                fa::sbend_unit_term(by6::c4, step_angle, sp.r0) }

            , step_dphi{
                fa::sbend_dphi(by6::c1, step_angle),
                fa::sbend_dphi(by6::c2, step_angle),
                fa::sbend_dphi(by6::c3, step_angle),
                fa::sbend_dphi(by6::c4, step_angle) }
        { }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            if (masks(i))
            {
                if (sp.ledge)
                {
                    // slot
                    FF_algorithm::slot_unit(
                            p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                            sp.ce1, sp.se1, sp.pref_b, sp.m_b );

                    // edge
                    //FF_algorithm::edge_unit(y, yp, us_edge_k);
                    //FF_algorithm::edge_unit(y, xp, yp, dpop, us_edge_k_p);
                    if (sp.ledgekick)
                    {
                        FF_algorithm::edge_unit(
                            p(i,2), p(i,1), p(i,3), 
                            sp.us_edge_k_x, sp.us_edge_k_y, 0); 
                    }

                    // bend edge (thin, but with face angle)
                    FF_algorithm::bend_edge(
                            p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                            sp.e1, phase_e1, sp.strength, sp.pref_b, sp.m_b);
                }

                // bend body
                FF_algorithm::bend_yoshida6<typename BP::part_t,
                    fa::thin_cf_kick_5<typename BP::part_t>, 5> ( 
                            p(i,0), p(i,1), p(i,2), 
                            p(i,3), p(i,4), p(i,5),
                            sp.pref_b, sp.m_b, step_ref_cdt,
                            step_kl, step_dphi, step_phase, step_term,
                            sp.r0, sp.strength, sp.steps );

                if (sp.redge)
                {
                    // bend edge (thin, but with face angle)
                    FF_algorithm::bend_edge(
                            p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                            sp.e2, phase_e2, sp.strength, sp.pref_b, sp.m_b);

                    // edge
                    //FF_algorithm::edge_unit(y, yp, ds_edge_k);
                    //FF_algorithm::edge_unit(y, xp, yp, dpop, ds_edge_k_p);
                    if (sp.redgekick)
                    {
                        FF_algorithm::edge_unit(
                            p(i,2), p(i,1), p(i,3), 
                            sp.ds_edge_k_x, sp.ds_edge_k_y, 0); 
                    }

                    // slot
                    FF_algorithm::slot_unit(
                            p(i,0), p(i,1), p(i,2), p(i,3), p(i,4), p(i,5),
                            sp.ce2, sp.se2, sp.pref_b, sp.m_b);
                }
            }
        }
    };

    template<class BP>
    struct PropSbendCFSimd
    {
        using gsv_t = typename BP::gsv_t;
        using parts_t = typename BP::parts_t;
        using const_masks_t = typename BP::const_masks_t;

        parts_t p;
        const_masks_t masks;
        const SbendParams sp;

        const double step_angle;
        const double step_ref_cdt;

        const double step_kl[10];

        const Kokkos::complex<double> phase_e1;
        const Kokkos::complex<double> phase_e2;

        const Kokkos::complex<double> step_phase[4];
        const Kokkos::complex<double> step_term[4];
        const double step_dphi[4];

        PropSbendCFSimd( BP & bp, SbendParams const& sp )
            : p(bp.parts), masks(bp.masks), sp(sp)

            , step_angle(sp.angle/sp.steps)
            , step_ref_cdt(sp.ref_cdt/sp.steps)

            , step_kl{  // k_b[i] = k_l[i] * scale
                sp.k_l[0] * sp.scale * sp.length / sp.steps,
                sp.k_l[1] * sp.scale * sp.length / sp.steps,
                sp.k_l[2] * sp.scale * sp.length / sp.steps,
                sp.k_l[3] * sp.scale * sp.length / sp.steps,
                sp.k_l[4] * sp.scale * sp.length / sp.steps,
                sp.k_l[5] * sp.scale * sp.length / sp.steps,
                sp.k_l[6] * sp.scale * sp.length / sp.steps,
                sp.k_l[7] * sp.scale * sp.length / sp.steps,
                sp.k_l[8] * sp.scale * sp.length / sp.steps,
                sp.k_l[9] * sp.scale * sp.length / sp.steps }

            , phase_e1(fa::bend_edge_phase(sp.e1))
            , phase_e2(fa::bend_edge_phase(sp.e2))

            , step_phase{
                fa::sbend_unit_phase(by6::c1, step_angle),
                fa::sbend_unit_phase(by6::c2, step_angle),
                fa::sbend_unit_phase(by6::c3, step_angle),
                fa::sbend_unit_phase(by6::c4, step_angle) }

            , step_term{
                fa::sbend_unit_term(by6::c1, step_angle, sp.r0),
                fa::sbend_unit_term(by6::c2, step_angle, sp.r0),
                fa::sbend_unit_term(by6::c3, step_angle, sp.r0),
                fa::sbend_unit_term(by6::c4, step_angle, sp.r0) }

            , step_dphi{
                fa::sbend_dphi(by6::c1, step_angle),
                fa::sbend_dphi(by6::c2, step_angle),
                fa::sbend_dphi(by6::c3, step_angle),
                fa::sbend_dphi(by6::c4, step_angle) }
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
                gsv_t p4(&p(i, 4));
                gsv_t p5(&p(i, 5));

                if (sp.ledge)
                {
                    // slot
                    FF_algorithm::slot_unit<gsv_t>(
                            p0, p1, p2, p3, p4, p5,
                            sp.ce1, sp.se1, sp.pref_b, sp.m_b );

                    // edge
                    //FF_algorithm::edge_unit(y, yp, us_edge_k);
                    //FF_algorithm::edge_unit(y, xp, yp, dpop, us_edge_k_p);
                    if (sp.ledgekick)
                    {
                        FF_algorithm::edge_unit<gsv_t>(
                                p2, p1, p3, 
                                sp.us_edge_k_x, sp.us_edge_k_y, 0); 
                    }

                    // bend edge (thin, but with face angle)
                    FF_algorithm::bend_edge<gsv_t>(
                            p0, p1, p2, p3, p4, p5,
                            sp.e1, phase_e1, sp.strength, sp.pref_b, sp.m_b);
                }

                // bend body
                FF_algorithm::bend_yoshida6<gsv_t, 
                    fa::thin_cf_kick_5<gsv_t>, 5> ( 
                            p0, p1, p2, p3, p4, p5,
                            sp.pref_b, sp.m_b, step_ref_cdt,
                            step_kl, step_dphi, step_phase, step_term,
                            sp.r0, sp.strength, sp.steps );

                if (sp.redge)
                {
                    // bend edge (thin, but with face angle)
                    FF_algorithm::bend_edge<gsv_t>(
                            p0, p1, p2, p3, p4, p5,
                            sp.e2, phase_e2, sp.strength, sp.pref_b, sp.m_b);

                    // edge
                    //FF_algorithm::edge_unit(y, yp, ds_edge_k);
                    //FF_algorithm::edge_unit(y, xp, yp, dpop, ds_edge_k_p);
                    if (sp.redgekick)
                    {
                        FF_algorithm::edge_unit<gsv_t>(
                            p2, p1, p3, 
                            sp.ds_edge_k_x, sp.ds_edge_k_y, 0); 
                    }

                    // slot
                    FF_algorithm::slot_unit<gsv_t>(
                            p0, p1, p2, p3, p4, p5,
                            sp.ce2, sp.se2, sp.pref_b, sp.m_b);
                }

                p0.store(&p(i, 0));
                p1.store(&p(i, 1));
                p2.store(&p(i, 2));
                p3.store(&p(i, 3));
                p4.store(&p(i, 4));
                p5.store(&p(i, 5));
            }
        }
    };

    inline void prop_reference(
            Reference_particle & ref_l, 
            SbendParams & sp )
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
        double dpop_l = 0.0;

        //for reference particle start with dp/p==0
        //double dpop_l = ref_l.get_state()[Bunch::dpop];

        double dphi =  -(sp.angle - (sp.e1 + sp.e2));

        Kokkos::complex<double> phase = 
            std::exp(std::complex<double>(0.0, -dphi));

        Kokkos::complex<double> term  = 
            std::complex<double>(0.0, sp.length / sp.angle) *
            std::complex<double>(1.0 - cos(sp.angle), - sin(sp.angle)) *
            std::complex<double>(cos(sp.e2), -sin(sp.e2));

        if (sp.ledge)
        {
            // slot
            FF_algorithm::slot_unit(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.ce1, sp.se1, pref_l, m_l);

            double p_l = 1.0 + dpop_l;
            double zp_l = sqrt(p_l*p_l - xp_l*xp_l - yp_l*yp_l);

            sp.us_edge_k_x = sp.us_edge_k_p * (xp_l/zp_l);
            sp.us_edge_k_y = sp.us_edge_k_p * (yp_l/zp_l);

            // edge kick that is in accordance with 
            // method 1 (chef fixed angle)
            // double brho_l = pref_l / PH_CNV_brho_to_p;
            // sp.us_edge_k = sp.strength * tan(atan2(xp_l, zp_l)) / brho_l;
            // FF_algorithm::edge_unit(y_l, yp_l, sp.us_edge_k );

            // edge kick strenth are scaled to bunch. so need to div by "scale" to scale
            // it to the lattice reference
            // in accordance with method 3
            if (sp.ledgekick)
            {
                FF_algorithm::edge_unit(y_l, xp_l, yp_l, 
                    sp.us_edge_k_x/sp.scale, sp.us_edge_k_y/sp.scale, 0);
            }

        }

        FF_algorithm::bend_complete(
                x_l, xp_l, y_l, yp_l, cdt_l, dpop_l,
                dphi, sp.strength, pref_l, m_l, 0.0/*ref cdt*/, 
                phase, term);

        if (sp.redge)
        {
            double p_l = 1.0 + dpop_l;
            double zp_l = sqrt(p_l*p_l - xp_l*xp_l - yp_l*yp_l);

            sp.ds_edge_k_x = sp.ds_edge_k_p * (xp_l/zp_l);
            sp.ds_edge_k_y = sp.ds_edge_k_p * (yp_l/zp_l);

            // edge kick that is in accordance with 
            // method 1 (chef fixed angle)
            // double brho_l = pref_l / PH_CNV_brho_to_p;
            // sp.ds_edge_k = - sp.strength * tan(atan2(xp_l, zp_l)) / brho_l;
            // FF_algorithm::edge_unit(y_l, yp_l, sp.ds_edge_k );

            // edge kick strenth are scaled to bunch. so need to div by "scale" to scale
            // it to the lattice reference
            if (sp.redgekick)
            {
                FF_algorithm::edge_unit(y_l, xp_l, yp_l, 
                    sp.ds_edge_k_x/sp.scale, sp.ds_edge_k_y/sp.scale, 0);
            }

            // slot
            FF_algorithm::slot_unit(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.ce2, sp.se2, pref_l, m_l);
        }

        ref_l.set_state(x_l, xp_l, y_l, yp_l, cdt_l, dpop_l);
        sp.ref_cdt = cdt_l;
    }

    inline void prop_reference_cf(
            Reference_particle & ref_l, 
            SbendParams & sp )
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
        // for reference cdt, start with dp/p==0
        // double dpop_l = ref_l.get_state()[Bunch::dpop];
        double dpop_l = 0.0;

        double step_length = sp.length / sp.steps;
        double step_angle = sp.angle/sp.steps;

        // step strength for reference
        double step_kl[10];
        for (int i=0; i<10; ++i) step_kl[i] = sp.k_l[i] * step_length;

        Kokkos::complex<double> phase_e1 = FF_algorithm::bend_edge_phase(sp.e1);
        Kokkos::complex<double> phase_e2 = FF_algorithm::bend_edge_phase(sp.e2);

        Kokkos::complex<double> step_phase[4] = { 
            fa::sbend_unit_phase(by6::c1, step_angle),
            fa::sbend_unit_phase(by6::c2, step_angle),
            fa::sbend_unit_phase(by6::c3, step_angle),
            fa::sbend_unit_phase(by6::c4, step_angle) 
        };

        Kokkos::complex<double> step_term[4] = {
            fa::sbend_unit_term(by6::c1, step_angle, sp.r0),
            fa::sbend_unit_term(by6::c2, step_angle, sp.r0),
            fa::sbend_unit_term(by6::c3, step_angle, sp.r0),
            fa::sbend_unit_term(by6::c4, step_angle, sp.r0)
        };

        double step_dphi[4] = {
            fa::sbend_dphi(by6::c1, step_angle),
            fa::sbend_dphi(by6::c2, step_angle),
            fa::sbend_dphi(by6::c3, step_angle),
            fa::sbend_dphi(by6::c4, step_angle),
        };

        if (sp.ledge)
        {
            // slot
            FF_algorithm::slot_unit(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.ce1, sp.se1, pref_l, m_l);

            double p_l = 1.0 + dpop_l;
            double zp_l = sqrt(p_l*p_l - xp_l*xp_l - yp_l*yp_l);

            sp.us_edge_k_x = sp.us_edge_k_p * (xp_l/zp_l);
            sp.us_edge_k_y = sp.us_edge_k_p * (yp_l/zp_l);

            // edge
            //FF_algorithm::edge_unit(y_l, yp_l, us_edge_k/scale);
            //FF_algorithm::edge_unit(y_l, xp_l, yp_l, dpop_l, us_edge_k_p/scale);
            if (sp.ledgekick)
            {
                FF_algorithm::edge_unit(y_l, xp_l, yp_l, 
                    sp.us_edge_k_x/sp.scale, sp.us_edge_k_y/sp.scale, 0);
            }

            // bend edge (thin)
            FF_algorithm::bend_edge(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.e1, phase_e1, sp.strength, pref_l, m_l);
        }

        FF_algorithm::bend_yoshida6<double, 
            FF_algorithm::thin_cf_kick_5<double>, 5> ( 
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l,
                    pref_l, m_l, 0.0 /* step ref_cdt */,
                    step_kl, step_dphi, step_phase, step_term,
                    sp.r0, sp.strength, sp.steps);

        if (sp.redge)
        {
            // bend edge (thin)
            if (sp.redgekick)
            {
                FF_algorithm::bend_edge(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.e2, phase_e2, sp.strength, pref_l, m_l);
            }

            double p_l = 1.0 + dpop_l;
            double zp_l = sqrt(p_l*p_l - xp_l*xp_l - yp_l*yp_l);

            sp.ds_edge_k_x = sp.ds_edge_k_p * (xp_l/zp_l);
            sp.ds_edge_k_y = sp.ds_edge_k_p * (yp_l/zp_l);

            // edge
            //FF_algorithm::edge_unit(y_l, yp_l, ds_edge_k);
            //FF_algorithm::edge_unit(y_l, xp_l, yp_l, dpop_l, ds_edge_k_p/scale);
            FF_algorithm::edge_unit(y_l, xp_l, yp_l, 
                    sp.ds_edge_k_x/sp.scale, sp.ds_edge_k_y/sp.scale, 0);

            // slot
            FF_algorithm::slot_unit(
                    x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, 
                    sp.ce2, sp.se2, pref_l, m_l);
        }

        ref_l.set_state(x_l, xp_l, y_l, yp_l, cdt_l, dpop_l);
        sp.ref_cdt = cdt_l;
    }

    // max multipole moments
    constexpr const int max_mp_order = 5;

    template<class BP>
    struct PropThinPoleSimd
    {
        using gsv_t = typename BP::gsv_t;

        typename BP::parts_t p;
        typename BP::const_masks_t masks;
        kt::arr_t<double, 2*max_mp_order> k;

        KOKKOS_INLINE_FUNCTION
        void operator()(const int idx) const
        {
            int i = idx * gsv_t::size();

            int m = 0;
            for(int x=i; x<i+gsv_t::size(); ++x) m |= masks(x);

            if (m)
            {
                gsv_t  x(&p(i, 0));
                gsv_t xp(&p(i, 1));
                gsv_t  y(&p(i, 2));
                gsv_t yp(&p(i, 3));

                if (k[2] || k[3])   
                    FF_algorithm::thin_quadrupole_unit(x, xp, y, yp, k.data+2);

                if (k[4] || k[5])   
                    FF_algorithm::thin_sextupole_unit(x, xp, y, yp, k.data+4);

                if (k[6] || k[7])   
                    FF_algorithm::thin_octupole_unit(x, xp, y, yp, k.data+6);

                if (k[8] || k[9])   
                    FF_algorithm::thin_magnet_unit(x, xp, y, yp, k.data+8, 5);

#if 0
                if (k[10] || k[11]) 
                    FF_algorithm::thin_magnet_unit(x, xp, y, yp, k.data+10, 6);

                if (k[12] || k[13]) 
                    FF_algorithm::thin_magnet_unit(x, xp, y, yp, k.data+12, 7);

                if (k[14] || k[15]) 
                    FF_algorithm::thin_magnet_unit(x, xp, y, yp, k.data+14, 8);
#endif

                xp.store(&p(i, 1));
                yp.store(&p(i, 3));
            }
        }
    };


}

namespace FF_sbend
{

template<class BunchT>
inline void apply(Lattice_element_slice const& slice, BunchT & bunch)
{
    using namespace sbend_impl;

    scoped_simple_timer timer("libFF_sbend");

    auto const& ele = slice.get_lattice_element();

    SbendParams sp;

    double a = ele.get_double_attribute("angle");
    double l = ele.get_double_attribute("l");

    sp.length = slice.get_right() - slice.get_left();
    sp.angle  = ( sp.length / l ) * a;
    sp.r0     = l / a;
    sp.steps  = FF_algorithm::default_steps;

    sp.ledge = slice.has_left_edge();
    sp.redge = slice.has_right_edge();
    sp.ledgekick = ele.get_double_attribute("kill_entry_kick", 0.0) == 0.0;
    sp.redgekick = ele.get_double_attribute("kill_exit_kick", 0.0 ) == 0.0;

    sp.e1 = ele.get_double_attribute("e1", 0.0);
    sp.e2 = ele.get_double_attribute("e2", 0.0);

    Reference_particle       & ref_l = bunch.get_design_reference_particle();
    Reference_particle const & ref_b = bunch.get_reference_particle();

    // quad, sextupole, and octupole components
    sp.k_l[0] = ele.get_double_attribute("k1", 0.0);
    sp.k_l[2] = ele.get_double_attribute("k2", 0.0);
    sp.k_l[4] = ele.get_double_attribute("k3", 0.0);
    sp.k_l[6] = ele.get_double_attribute("k4", 0.0);
    sp.k_l[8] = ele.get_double_attribute("k5", 0.0);

    sp.k_l[1] = ele.get_double_attribute("k1s", 0.0);
    sp.k_l[3] = ele.get_double_attribute("k2s", 0.0);
    sp.k_l[5] = ele.get_double_attribute("k3s", 0.0);
    sp.k_l[7] = ele.get_double_attribute("k4s", 0.0);
    sp.k_l[9] = ele.get_double_attribute("k5s", 0.0);

    int cf = 0;  // combined function

    if ((sp.k_l[0] != 0.0) || (sp.k_l[1] != 0.0)) cf = 1;
    if ((sp.k_l[2] != 0.0) || (sp.k_l[3] != 0.0)) cf = 2;
    if ((sp.k_l[4] != 0.0) || (sp.k_l[5] != 0.0)) cf = 3;
    if ((sp.k_l[6] != 0.0) || (sp.k_l[7] != 0.0)) cf = 4;
    if ((sp.k_l[8] != 0.0) || (sp.k_l[9] != 0.0)) cf = 5;

    // error multipole terms
    kt::arr_t<double, 2*max_mp_order> kn;

    // multipole moments 
    // kn[2] ... kn[2*max_mp_order-1]
    bool has_mp = false;

    std::string a_attr = "a0";
    std::string b_attr = "b0";

    for(int i=1; i<max_mp_order; ++i)
    {
        a_attr[1] = '0' + i;
        b_attr[1] = '0' + i;

        kn[i*2+0] = ele.get_double_attribute(b_attr, 0.0);
        kn[i*2+1] = ele.get_double_attribute(a_attr, 0.0);

        if (kn[i*2+0] || kn[i*2+1]) 
        {
            has_mp = true;

            double coeff = FF_algorithm::factorial(i);

            kn[i*2+0] *= a * coeff;
            kn[i*2+1] *= a * coeff;

            if (sp.k_l[(i-1)*2+0] || sp.k_l[(i-1)*2+1])
            {
                throw std::runtime_error(
                        "having both kn/ks and a/b coefficents are forbidden in sbends");
            }

            sp.k_l[(i-1)*2+0] = kn[i*2+0];
            sp.k_l[(i-1)*2+1] = kn[i*2+1];
        }
    }


    sp.scale = ref_l.get_momentum() / 
        (ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]));

    if (!sp.redge) sp.e2 = 0;
    if (!sp.ledge) sp.e1 = 0;

    // lattice reference
    double pref_l = ref_l.get_momentum();
    double brho_l = pref_l / PH_CNV_brho_to_p;
    double    m_l = ref_l.get_mass();
    int  charge_l = ref_l.get_charge();

    // the particle dpop is with respect to this momentum which goes with the bunch
    double pref_b = ref_b.get_momentum();
    double brho_b = pref_b / PH_CNV_brho_to_p;
    double    m_b = bunch.get_mass();

    // common
    sp.strength = brho_l * a / l;
    sp.pref_b   = pref_b;
    sp.m_b      = m_b;

    sp.ce1 = cos(-sp.e1);
    sp.se1 = sin(-sp.e1);
    sp.ce2 = cos(-sp.e2);
    sp.se2 = sin(-sp.e2);

    // edge kick strength (scaled to particle)
    sp.us_edge_k =   ((charge_l > 0) ? 1.0 : -1.0) * sp.strength * tan( sp.e1) / brho_b;
    sp.ds_edge_k = - ((charge_l > 0) ? 1.0 : -1.0) * sp.strength * tan(-sp.e2) / brho_b;

    // edge kick (per) particle (full strength, angles are calculated at each particle)
    sp.us_edge_k_p =   ((charge_l > 0) ? 1.0 : -1.0) * sp.strength / brho_b;
    sp.ds_edge_k_p = - ((charge_l > 0) ? 1.0 : -1.0) * sp.strength / brho_b;

    // edge kick x/y
    sp.us_edge_k_x = 0.0;
    sp.us_edge_k_y = 0.0;

    sp.ds_edge_k_x = 0.0;
    sp.ds_edge_k_y = 0.0;

    using namespace Kokkos;
    using exec = typename BunchT::exec_space;

    if (cf == 0)
    {
        auto apply = [&](ParticleGroup pg) {
            auto bp = bunch.get_bunch_particles(pg);
            if (!bp.num_valid()) return;

#if LIBFF_USE_GSV
            if (has_mp)
            {
#if 0
                // bend - multipole kick - bend
                auto range = RangePolicy<exec>(0, bp.size_in_gsv());

                auto sp1 = sp;
                auto sp2 = sp;

                sp1.redge = false;
                sp1.length *= 0.5;
                sp1.angle *= 0.5;
                sp1.e2 = 0.0;

                sp2.ledge = false;
                sp2.length *= 0.5;
                sp2.angle *= 0.5;
                sp2.e1 = 0.0;

                prop_reference(ref_l, sp1);
                prop_reference(ref_l, sp2);

                PropSbendSimd<typename BunchT::bp_t> sbend1(bp, sp1);
                Kokkos::parallel_for(range, sbend1);

                PropThinPoleSimd<typename BunchT::bp_t> ptp{
                    bp.parts, bp.masks, kn};
                Kokkos::parallel_for(range, ptp);

                PropSbendSimd<typename BunchT::bp_t> sbend2(bp, sp2);
                Kokkos::parallel_for(range, sbend2);
#else

                // Yoshida (same as CF sbends)
                prop_reference_cf(ref_l, sp);  // has a_n/b_n multipoles

                auto range = RangePolicy<exec>(0, bp.size_in_gsv());
                PropSbendCFSimd<typename BunchT::bp_t> sbend(bp, sp);
                Kokkos::parallel_for(range, sbend);
#endif
            }
            else
            {
                prop_reference(ref_l, sp);  // no multipoles (not combined-function)

                auto range = RangePolicy<exec>(0, bp.size_in_gsv());
                PropSbendSimd<typename BunchT::bp_t> sbend(bp, sp);
                Kokkos::parallel_for(range, sbend);
            }
#else
            auto range = RangePolicy<exec>(0, bp.size());
            PropSbend<typename BunchT::bp_t> sbend(bp, sp);
            Kokkos::parallel_for(range, sbend);
#endif
        };

        apply(ParticleGroup::regular);
        apply(ParticleGroup::spectator);
    }
    else
    {
        // propagate reference
        prop_reference_cf(ref_l, sp); // has k1, k2 attributes

        auto apply = [&](ParticleGroup pg) {
            auto bp = bunch.get_bunch_particles(pg);
            if (!bp.num_valid()) return;

#if LIBFF_USE_GSV
            auto range = RangePolicy<exec>(0, bp.size_in_gsv());
            PropSbendCFSimd<typename BunchT::bp_t> sbend(bp, sp);
            Kokkos::parallel_for(range, sbend);
#else
            auto range = RangePolicy<exec>(0, bp.size());
            PropSbendCF<typename BunchT::bp_t> sbend(bp, sp);
            Kokkos::parallel_for(range, sbend);
#endif
        };

        apply(ParticleGroup::regular);
        apply(ParticleGroup::spectator);
    }

    bunch.get_reference_particle().increment_trajectory(sp.length);
    // absolute time
    double const velocity = bunch.get_reference_particle().get_beta()*pconstants::c;
    bunch.get_reference_particle().increment_bunch_abs_time(sp.length/velocity);

    Kokkos::fence();
}

}

#endif // FF_SBEND_H
