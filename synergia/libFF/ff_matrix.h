#ifndef FF_MATRIX_H
#define FF_MATRIX_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/utils/simple_timer.h"

// Implementation of MAD-X arbitrary matrix element

/* MAD-X Documentation on matrix:

Arbitrary Matrix Element

label: MATRIX, TYPE=string, L=real,
    KICK1=real,..., KICK6=real,
    RM11=real, ..., RM66=real,
   TM111=real, ..., TM666=real;

The MATRIX element allows the definition of an arbitrary transfer matrix. It has
four real array attributes: LLength of the element, which may be zero.

KICKi  Defines the kick of the element acting on the six phase space
coordinates.

RMik   Defines the linear transfer matrix (6 × 6 terms) of the element.

TMikl  Defines the second-order terms (6 × 6 × 6 terms) of the element.

Data values that are not explicitly entered are taken from the identity
transformation for the RMik matrix elements, and taken as zero for the KICKi
kick factors and the TMikl second order terms. In the thin-lens tracking module
a non-zero length for an arbitrary matrix is accepted,98 CHAPTER 11. ELEMENT
TYPES however no non-zero second order terms are allowed to avoid non symplectic
tracking runs. In the latter case the tracking run is aborted.label:

 */

// p [Gev/c] = -- * B*rho [ Tesla meters ]
#define PH_CNV_brho_to_p (1.0e-9 * pconstants::c)

namespace matrix_impl {
    struct MatrixParams {
        double l;           // length
        double kick[6];     // KICKn parameters
        double rm[6][6];    // RM linear matrix
        double tm[6][6][6]; // TM 2nd order tensor
        // these are not currently used but may be in the future if
        // we decide that scale matrix elements by momentum makes sense (iffy)
        double scale, m_b;
        double pref_b;
    };

    template <class BP> // BP is bunch particles?
    struct PropMatrix {
        typename BP::parts_t p;
        typename BP::const_masks_t masks;
        const MatrixParams mp;

        PropMatrix(BP& bp, MatrixParams const& mp)
            : p(bp.parts), masks(bp.masks), mp(mp)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int i) const
        {
            if (masks(i)) {

                // bend
                FF_algorithm::matrix_unit(p(i, 0),
                                          p(i, 1),
                                          p(i, 2),
                                          p(i, 3),
                                          p(i, 4),
                                          p(i, 5),
                                          mp.kick,
                                          mp.rm,
                                          mp.tm);
            }
        }
    };

    template <class BP>
    struct PropMatrixSimd {
        using gsv_t = typename BP::gsv_t;
        using parts_t = typename BP::parts_t;

        typename BP::parts_t p;
        typename BP::const_masks_t masks;
        const MatrixParams mp;

        PropMatrixSimd(BP& bp, MatrixParams const& mp)
            : p(bp.parts), masks(bp.masks), mp(mp)
        {}

        KOKKOS_INLINE_FUNCTION
        void
        operator()(const int idx) const
        {
            int i = idx * gsv_t::size();

            int m = 0;
            for (int x = i; x < i + gsv_t::size(); ++x)
                m |= masks(x);

            if (m) {
                gsv_t p0(&p(i, 0));
                gsv_t p1(&p(i, 1));
                gsv_t p2(&p(i, 2));
                gsv_t p3(&p(i, 3));
                gsv_t p4(&p(i, 4));
                gsv_t p5(&p(i, 5));

                // matrix
                FF_algorithm::matrix_unit<gsv_t>(
                    p0, p1, p2, p3, p4, p5, mp.kick, mp.rm, mp.tm);

                p0.store(&p(i, 0));
                p1.store(&p(i, 1));
                p2.store(&p(i, 2));
                p3.store(&p(i, 3));
                p4.store(&p(i, 4));
                p5.store(&p(i, 5));
            }
        }
    };

    // mstrix must have 0 length
    inline void
    prop_reference(Reference_particle& ref_l, MatrixParams& mp)
    {
        if (mp.l != 0) {
            throw std::runtime_error(
                "error, matrix element has non-zero length");
        }

        double pref_l = ref_l.get_momentum();
        double m_l = ref_l.get_mass();

        // propagate the reference particle, and set the edge kick strength
        // from the reference particle
        double x_l = ref_l.get_state()[Bunch::x];
        double xp_l = ref_l.get_state()[Bunch::xp];
        double y_l = ref_l.get_state()[Bunch::y];
        double yp_l = ref_l.get_state()[Bunch::yp];
        double cdt_l = 0.0;
        double dpop_l = ref_l.get_state()[Bunch::dpop];

        FF_algorithm::matrix_unit(
            x_l, xp_l, y_l, yp_l, cdt_l, dpop_l, mp.kick, mp.rm, mp.tm);
    }
}

namespace FF_matrix {

    template <class BunchT>
    inline void
    apply(Lattice_element_slice const& slice, BunchT& bunch)
    {
        using namespace matrix_impl;

        scoped_simple_timer timer("libFF_matrix");

        auto const& ele = slice.get_lattice_element();

        MatrixParams mp;

        mp.l = ele.get_double_attribute("l", 0.0);
        // offset pointer so I can use 1-based indexing to match MAD-X
        // convention
        mp.kick[0] = ele.get_double_attribute("kick1", 0.0);
        mp.kick[1] = ele.get_double_attribute("kick2", 0.0);
        mp.kick[2] = ele.get_double_attribute("kick3", 0.0);
        mp.kick[3] = ele.get_double_attribute("kick4", 0.0);
        mp.kick[4] = ele.get_double_attribute("kick5", 0.0);
        mp.kick[5] = ele.get_double_attribute("kick6", 0.0);

        // diagonal elements default to 1, all others default to 0
        mp.rm[0][0] = ele.get_double_attribute("rm11", 1.0); // diagonal
        mp.rm[0][1] = ele.get_double_attribute("rm12", 0.0);
        mp.rm[0][2] = ele.get_double_attribute("rm13", 0.0);
        mp.rm[0][3] = ele.get_double_attribute("rm14", 0.0);
        mp.rm[0][4] = ele.get_double_attribute("rm15", 0.0);
        mp.rm[0][5] = ele.get_double_attribute("rm16", 0.0);

        mp.rm[1][0] = ele.get_double_attribute("rm21", 0.0);
        mp.rm[1][1] = ele.get_double_attribute("rm22", 1.0); // diagonal
        mp.rm[1][2] = ele.get_double_attribute("rm23", 0.0);
        mp.rm[1][3] = ele.get_double_attribute("rm24", 0.0);
        mp.rm[1][4] = ele.get_double_attribute("rm25", 0.0);
        mp.rm[1][5] = ele.get_double_attribute("rm26", 0.0);

        mp.rm[2][0] = ele.get_double_attribute("rm31", 0.0);
        mp.rm[2][1] = ele.get_double_attribute("rm32", 0.0);
        mp.rm[2][2] = ele.get_double_attribute("rm33", 1.0); // diagonal
        mp.rm[2][3] = ele.get_double_attribute("rm34", 0.0);
        mp.rm[2][4] = ele.get_double_attribute("rm35", 0.0);
        mp.rm[2][5] = ele.get_double_attribute("rm36", 0.0);

        mp.rm[3][0] = ele.get_double_attribute("rm41", 0.0);
        mp.rm[3][1] = ele.get_double_attribute("rm42", 0.0);
        mp.rm[3][2] = ele.get_double_attribute("rm43", 0.0);
        mp.rm[3][3] = ele.get_double_attribute("rm44", 1.0); // diagonal
        mp.rm[3][4] = ele.get_double_attribute("rm45", 0.0);
        mp.rm[3][5] = ele.get_double_attribute("rm46", 0.0);

        mp.rm[4][0] = ele.get_double_attribute("rm51", 0.0);
        mp.rm[4][1] = ele.get_double_attribute("rm52", 0.0);
        mp.rm[4][2] = ele.get_double_attribute("rm53", 0.0);
        mp.rm[4][3] = ele.get_double_attribute("rm54", 0.0);
        mp.rm[4][4] = ele.get_double_attribute("rm55", 1.0); // diagonal
        mp.rm[4][5] = ele.get_double_attribute("rm56", 0.0);

        mp.rm[5][0] = ele.get_double_attribute("rm61", 0.0);
        mp.rm[5][1] = ele.get_double_attribute("rm62", 0.0);
        mp.rm[5][2] = ele.get_double_attribute("rm63", 0.0);
        mp.rm[5][3] = ele.get_double_attribute("rm64", 0.0);
        mp.rm[5][4] = ele.get_double_attribute("rm65", 0.0);
        mp.rm[5][5] = ele.get_double_attribute("rm66", 1.0); // diagonal

        mp.tm[0][0][0] = ele.get_double_attribute("tm111", 0.0);
        mp.tm[0][0][1] = ele.get_double_attribute("tm112", 0.0);
        mp.tm[0][0][2] = ele.get_double_attribute("tm113", 0.0);
        mp.tm[0][0][3] = ele.get_double_attribute("tm114", 0.0);
        mp.tm[0][0][4] = ele.get_double_attribute("tm115", 0.0);
        mp.tm[0][0][5] = ele.get_double_attribute("tm116", 0.0);

        mp.tm[0][1][0] = ele.get_double_attribute("tm121", 0.0);
        mp.tm[0][1][1] = ele.get_double_attribute("tm122", 0.0);
        mp.tm[0][1][2] = ele.get_double_attribute("tm123", 0.0);
        mp.tm[0][1][3] = ele.get_double_attribute("tm124", 0.0);
        mp.tm[0][1][4] = ele.get_double_attribute("tm125", 0.0);
        mp.tm[0][1][5] = ele.get_double_attribute("tm126", 0.0);

        mp.tm[0][2][0] = ele.get_double_attribute("tm131", 0.0);
        mp.tm[0][2][1] = ele.get_double_attribute("tm132", 0.0);
        mp.tm[0][2][2] = ele.get_double_attribute("tm133", 0.0);
        mp.tm[0][2][3] = ele.get_double_attribute("tm134", 0.0);
        mp.tm[0][2][4] = ele.get_double_attribute("tm135", 0.0);
        mp.tm[0][2][5] = ele.get_double_attribute("tm136", 0.0);

        mp.tm[0][3][0] = ele.get_double_attribute("tm141", 0.0);
        mp.tm[0][3][1] = ele.get_double_attribute("tm142", 0.0);
        mp.tm[0][3][2] = ele.get_double_attribute("tm143", 0.0);
        mp.tm[0][3][3] = ele.get_double_attribute("tm144", 0.0);
        mp.tm[0][3][4] = ele.get_double_attribute("tm145", 0.0);
        mp.tm[0][3][5] = ele.get_double_attribute("tm146", 0.0);

        mp.tm[0][4][0] = ele.get_double_attribute("tm151", 0.0);
        mp.tm[0][4][1] = ele.get_double_attribute("tm152", 0.0);
        mp.tm[0][4][2] = ele.get_double_attribute("tm153", 0.0);
        mp.tm[0][4][3] = ele.get_double_attribute("tm154", 0.0);
        mp.tm[0][4][4] = ele.get_double_attribute("tm155", 0.0);
        mp.tm[0][4][5] = ele.get_double_attribute("tm156", 0.0);

        mp.tm[0][5][0] = ele.get_double_attribute("tm161", 0.0);
        mp.tm[0][5][1] = ele.get_double_attribute("tm162", 0.0);
        mp.tm[0][5][2] = ele.get_double_attribute("tm163", 0.0);
        mp.tm[0][5][3] = ele.get_double_attribute("tm164", 0.0);
        mp.tm[0][5][4] = ele.get_double_attribute("tm165", 0.0);
        mp.tm[0][5][5] = ele.get_double_attribute("tm166", 0.0);

        mp.tm[1][0][0] = ele.get_double_attribute("tm211", 0.0);
        mp.tm[1][0][1] = ele.get_double_attribute("tm212", 0.0);
        mp.tm[1][0][2] = ele.get_double_attribute("tm213", 0.0);
        mp.tm[1][0][3] = ele.get_double_attribute("tm214", 0.0);
        mp.tm[1][0][4] = ele.get_double_attribute("tm215", 0.0);
        mp.tm[1][0][5] = ele.get_double_attribute("tm216", 0.0);

        mp.tm[1][1][0] = ele.get_double_attribute("tm221", 0.0);
        mp.tm[1][1][1] = ele.get_double_attribute("tm222", 0.0);
        mp.tm[1][1][2] = ele.get_double_attribute("tm223", 0.0);
        mp.tm[1][1][3] = ele.get_double_attribute("tm224", 0.0);
        mp.tm[1][1][4] = ele.get_double_attribute("tm225", 0.0);
        mp.tm[1][1][5] = ele.get_double_attribute("tm226", 0.0);

        mp.tm[1][2][0] = ele.get_double_attribute("tm231", 0.0);
        mp.tm[1][2][1] = ele.get_double_attribute("tm232", 0.0);
        mp.tm[1][2][2] = ele.get_double_attribute("tm233", 0.0);
        mp.tm[1][2][3] = ele.get_double_attribute("tm234", 0.0);
        mp.tm[1][2][4] = ele.get_double_attribute("tm235", 0.0);
        mp.tm[1][2][5] = ele.get_double_attribute("tm236", 0.0);

        mp.tm[1][3][0] = ele.get_double_attribute("tm241", 0.0);
        mp.tm[1][3][1] = ele.get_double_attribute("tm242", 0.0);
        mp.tm[1][3][2] = ele.get_double_attribute("tm243", 0.0);
        mp.tm[1][3][3] = ele.get_double_attribute("tm244", 0.0);
        mp.tm[1][3][4] = ele.get_double_attribute("tm245", 0.0);
        mp.tm[1][3][5] = ele.get_double_attribute("tm246", 0.0);

        mp.tm[1][4][0] = ele.get_double_attribute("tm251", 0.0);
        mp.tm[1][4][1] = ele.get_double_attribute("tm252", 0.0);
        mp.tm[1][4][2] = ele.get_double_attribute("tm253", 0.0);
        mp.tm[1][4][3] = ele.get_double_attribute("tm254", 0.0);
        mp.tm[1][4][4] = ele.get_double_attribute("tm255", 0.0);
        mp.tm[1][4][5] = ele.get_double_attribute("tm256", 0.0);

        mp.tm[1][5][0] = ele.get_double_attribute("tm261", 0.0);
        mp.tm[1][5][1] = ele.get_double_attribute("tm262", 0.0);
        mp.tm[1][5][2] = ele.get_double_attribute("tm263", 0.0);
        mp.tm[1][5][3] = ele.get_double_attribute("tm264", 0.0);
        mp.tm[1][5][4] = ele.get_double_attribute("tm265", 0.0);
        mp.tm[1][5][5] = ele.get_double_attribute("tm266", 0.0);

        mp.tm[2][0][0] = ele.get_double_attribute("tm311", 0.0);
        mp.tm[2][0][1] = ele.get_double_attribute("tm312", 0.0);
        mp.tm[2][0][2] = ele.get_double_attribute("tm313", 0.0);
        mp.tm[2][0][3] = ele.get_double_attribute("tm314", 0.0);
        mp.tm[2][0][4] = ele.get_double_attribute("tm315", 0.0);
        mp.tm[2][0][5] = ele.get_double_attribute("tm316", 0.0);

        mp.tm[2][1][0] = ele.get_double_attribute("tm321", 0.0);
        mp.tm[2][1][1] = ele.get_double_attribute("tm322", 0.0);
        mp.tm[2][1][2] = ele.get_double_attribute("tm323", 0.0);
        mp.tm[2][1][3] = ele.get_double_attribute("tm324", 0.0);
        mp.tm[2][1][4] = ele.get_double_attribute("tm325", 0.0);
        mp.tm[2][1][5] = ele.get_double_attribute("tm326", 0.0);

        mp.tm[2][2][0] = ele.get_double_attribute("tm331", 0.0);
        mp.tm[2][2][1] = ele.get_double_attribute("tm332", 0.0);
        mp.tm[2][2][2] = ele.get_double_attribute("tm333", 0.0);
        mp.tm[2][2][3] = ele.get_double_attribute("tm334", 0.0);
        mp.tm[2][2][4] = ele.get_double_attribute("tm335", 0.0);
        mp.tm[2][2][5] = ele.get_double_attribute("tm336", 0.0);

        mp.tm[2][3][0] = ele.get_double_attribute("tm341", 0.0);
        mp.tm[2][3][1] = ele.get_double_attribute("tm342", 0.0);
        mp.tm[2][3][2] = ele.get_double_attribute("tm343", 0.0);
        mp.tm[2][3][3] = ele.get_double_attribute("tm344", 0.0);
        mp.tm[2][3][4] = ele.get_double_attribute("tm345", 0.0);
        mp.tm[2][3][5] = ele.get_double_attribute("tm346", 0.0);

        mp.tm[2][4][0] = ele.get_double_attribute("tm351", 0.0);
        mp.tm[2][4][1] = ele.get_double_attribute("tm352", 0.0);
        mp.tm[2][4][2] = ele.get_double_attribute("tm353", 0.0);
        mp.tm[2][4][3] = ele.get_double_attribute("tm354", 0.0);
        mp.tm[2][4][4] = ele.get_double_attribute("tm355", 0.0);
        mp.tm[2][4][5] = ele.get_double_attribute("tm356", 0.0);

        mp.tm[2][5][0] = ele.get_double_attribute("tm361", 0.0);
        mp.tm[2][5][1] = ele.get_double_attribute("tm362", 0.0);
        mp.tm[2][5][2] = ele.get_double_attribute("tm363", 0.0);
        mp.tm[2][5][3] = ele.get_double_attribute("tm364", 0.0);
        mp.tm[2][5][4] = ele.get_double_attribute("tm365", 0.0);
        mp.tm[2][5][5] = ele.get_double_attribute("tm366", 0.0);

        mp.tm[3][0][0] = ele.get_double_attribute("tm411", 0.0);
        mp.tm[3][0][1] = ele.get_double_attribute("tm412", 0.0);
        mp.tm[3][0][2] = ele.get_double_attribute("tm413", 0.0);
        mp.tm[3][0][3] = ele.get_double_attribute("tm414", 0.0);
        mp.tm[3][0][4] = ele.get_double_attribute("tm415", 0.0);
        mp.tm[3][0][5] = ele.get_double_attribute("tm416", 0.0);

        mp.tm[3][1][0] = ele.get_double_attribute("tm421", 0.0);
        mp.tm[3][1][1] = ele.get_double_attribute("tm422", 0.0);
        mp.tm[3][1][2] = ele.get_double_attribute("tm423", 0.0);
        mp.tm[3][1][3] = ele.get_double_attribute("tm424", 0.0);
        mp.tm[3][1][4] = ele.get_double_attribute("tm425", 0.0);
        mp.tm[3][1][5] = ele.get_double_attribute("tm426", 0.0);

        mp.tm[3][2][0] = ele.get_double_attribute("tm431", 0.0);
        mp.tm[3][2][1] = ele.get_double_attribute("tm432", 0.0);
        mp.tm[3][2][2] = ele.get_double_attribute("tm433", 0.0);
        mp.tm[3][2][3] = ele.get_double_attribute("tm434", 0.0);
        mp.tm[3][2][4] = ele.get_double_attribute("tm435", 0.0);
        mp.tm[3][2][5] = ele.get_double_attribute("tm436", 0.0);

        mp.tm[3][3][0] = ele.get_double_attribute("tm441", 0.0);
        mp.tm[3][3][1] = ele.get_double_attribute("tm442", 0.0);
        mp.tm[3][3][2] = ele.get_double_attribute("tm443", 0.0);
        mp.tm[3][3][3] = ele.get_double_attribute("tm444", 0.0);
        mp.tm[3][3][4] = ele.get_double_attribute("tm445", 0.0);
        mp.tm[3][3][5] = ele.get_double_attribute("tm446", 0.0);

        mp.tm[3][4][0] = ele.get_double_attribute("tm451", 0.0);
        mp.tm[3][4][1] = ele.get_double_attribute("tm452", 0.0);
        mp.tm[3][4][2] = ele.get_double_attribute("tm453", 0.0);
        mp.tm[3][4][3] = ele.get_double_attribute("tm454", 0.0);
        mp.tm[3][4][4] = ele.get_double_attribute("tm455", 0.0);
        mp.tm[3][4][5] = ele.get_double_attribute("tm456", 0.0);

        mp.tm[3][5][0] = ele.get_double_attribute("tm461", 0.0);
        mp.tm[3][5][1] = ele.get_double_attribute("tm462", 0.0);
        mp.tm[3][5][2] = ele.get_double_attribute("tm463", 0.0);
        mp.tm[3][5][3] = ele.get_double_attribute("tm464", 0.0);
        mp.tm[3][5][4] = ele.get_double_attribute("tm465", 0.0);
        mp.tm[3][5][5] = ele.get_double_attribute("tm466", 0.0);

        mp.tm[4][0][0] = ele.get_double_attribute("tm511", 0.0);
        mp.tm[4][0][1] = ele.get_double_attribute("tm512", 0.0);
        mp.tm[4][0][2] = ele.get_double_attribute("tm513", 0.0);
        mp.tm[4][0][3] = ele.get_double_attribute("tm514", 0.0);
        mp.tm[4][0][4] = ele.get_double_attribute("tm515", 0.0);
        mp.tm[4][0][5] = ele.get_double_attribute("tm516", 0.0);

        mp.tm[4][1][0] = ele.get_double_attribute("tm521", 0.0);
        mp.tm[4][1][1] = ele.get_double_attribute("tm522", 0.0);
        mp.tm[4][1][2] = ele.get_double_attribute("tm523", 0.0);
        mp.tm[4][1][3] = ele.get_double_attribute("tm524", 0.0);
        mp.tm[4][1][4] = ele.get_double_attribute("tm525", 0.0);
        mp.tm[4][1][5] = ele.get_double_attribute("tm526", 0.0);

        mp.tm[4][2][0] = ele.get_double_attribute("tm531", 0.0);
        mp.tm[4][2][1] = ele.get_double_attribute("tm532", 0.0);
        mp.tm[4][2][2] = ele.get_double_attribute("tm533", 0.0);
        mp.tm[4][2][3] = ele.get_double_attribute("tm534", 0.0);
        mp.tm[4][2][4] = ele.get_double_attribute("tm535", 0.0);
        mp.tm[4][2][5] = ele.get_double_attribute("tm536", 0.0);

        mp.tm[4][3][0] = ele.get_double_attribute("tm541", 0.0);
        mp.tm[4][3][1] = ele.get_double_attribute("tm542", 0.0);
        mp.tm[4][3][2] = ele.get_double_attribute("tm543", 0.0);
        mp.tm[4][3][3] = ele.get_double_attribute("tm544", 0.0);
        mp.tm[4][3][4] = ele.get_double_attribute("tm545", 0.0);
        mp.tm[4][3][5] = ele.get_double_attribute("tm546", 0.0);

        mp.tm[4][4][0] = ele.get_double_attribute("tm551", 0.0);
        mp.tm[4][4][1] = ele.get_double_attribute("tm552", 0.0);
        mp.tm[4][4][2] = ele.get_double_attribute("tm553", 0.0);
        mp.tm[4][4][3] = ele.get_double_attribute("tm554", 0.0);
        mp.tm[4][4][4] = ele.get_double_attribute("tm555", 0.0);
        mp.tm[4][4][5] = ele.get_double_attribute("tm556", 0.0);

        mp.tm[4][5][0] = ele.get_double_attribute("tm561", 0.0);
        mp.tm[4][5][1] = ele.get_double_attribute("tm562", 0.0);
        mp.tm[4][5][2] = ele.get_double_attribute("tm563", 0.0);
        mp.tm[4][5][3] = ele.get_double_attribute("tm564", 0.0);
        mp.tm[4][5][4] = ele.get_double_attribute("tm565", 0.0);
        mp.tm[4][5][5] = ele.get_double_attribute("tm566", 0.0);

        mp.tm[5][0][0] = ele.get_double_attribute("tm611", 0.0);
        mp.tm[5][0][1] = ele.get_double_attribute("tm612", 0.0);
        mp.tm[5][0][2] = ele.get_double_attribute("tm613", 0.0);
        mp.tm[5][0][3] = ele.get_double_attribute("tm614", 0.0);
        mp.tm[5][0][4] = ele.get_double_attribute("tm615", 0.0);
        mp.tm[5][0][5] = ele.get_double_attribute("tm616", 0.0);

        mp.tm[5][1][0] = ele.get_double_attribute("tm621", 0.0);
        mp.tm[5][1][1] = ele.get_double_attribute("tm622", 0.0);
        mp.tm[5][1][2] = ele.get_double_attribute("tm623", 0.0);
        mp.tm[5][1][3] = ele.get_double_attribute("tm624", 0.0);
        mp.tm[5][1][4] = ele.get_double_attribute("tm625", 0.0);
        mp.tm[5][1][5] = ele.get_double_attribute("tm626", 0.0);

        mp.tm[5][2][0] = ele.get_double_attribute("tm631", 0.0);
        mp.tm[5][2][1] = ele.get_double_attribute("tm632", 0.0);
        mp.tm[5][2][2] = ele.get_double_attribute("tm633", 0.0);
        mp.tm[5][2][3] = ele.get_double_attribute("tm634", 0.0);
        mp.tm[5][2][4] = ele.get_double_attribute("tm635", 0.0);
        mp.tm[5][2][5] = ele.get_double_attribute("tm636", 0.0);

        mp.tm[5][3][0] = ele.get_double_attribute("tm641", 0.0);
        mp.tm[5][3][1] = ele.get_double_attribute("tm642", 0.0);
        mp.tm[5][3][2] = ele.get_double_attribute("tm643", 0.0);
        mp.tm[5][3][3] = ele.get_double_attribute("tm644", 0.0);
        mp.tm[5][3][4] = ele.get_double_attribute("tm645", 0.0);
        mp.tm[5][3][5] = ele.get_double_attribute("tm646", 0.0);

        mp.tm[5][4][0] = ele.get_double_attribute("tm651", 0.0);
        mp.tm[5][4][1] = ele.get_double_attribute("tm652", 0.0);
        mp.tm[5][4][2] = ele.get_double_attribute("tm653", 0.0);
        mp.tm[5][4][3] = ele.get_double_attribute("tm654", 0.0);
        mp.tm[5][4][4] = ele.get_double_attribute("tm655", 0.0);
        mp.tm[5][4][5] = ele.get_double_attribute("tm656", 0.0);

        mp.tm[5][5][0] = ele.get_double_attribute("tm661", 0.0);
        mp.tm[5][5][1] = ele.get_double_attribute("tm662", 0.0);
        mp.tm[5][5][2] = ele.get_double_attribute("tm663", 0.0);
        mp.tm[5][5][3] = ele.get_double_attribute("tm664", 0.0);
        mp.tm[5][5][4] = ele.get_double_attribute("tm665", 0.0);
        mp.tm[5][5][5] = ele.get_double_attribute("tm666", 0.0);

        // jury is still out whether to scale matrix by momentum
        Reference_particle& ref_l = bunch.get_design_reference_particle();
        Reference_particle const& ref_b = bunch.get_reference_particle();

        mp.scale =
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
        mp.pref_b = pref_b;
        mp.m_b = m_b;

        using namespace Kokkos;
        using exec = typename BunchT::exec_space;

        auto apply = [&](ParticleGroup pg) {
            auto bp = bunch.get_bunch_particles(pg);
            if (!bp.num_valid()) return;

                // matrix is a 0 length element with no multipole moments

#if LIBFF_USE_GSV
            prop_reference(ref_l, mp);

            auto range = RangePolicy<exec>(0, bp.size_in_gsv());
            PropMatrixSimd<typename BunchT::bp_t> matrix(bp, mp);
            Kokkos::parallel_for(range, matrix);

#else  // LIBFF_USE_GSV
            auto range = RangePolicy<exec>(0, bp.size());
            PropMatrix<typename BunchT::bp_t> matrix(bp, mp);
            Kokkos::parallel_for(range, matrix);
#endif // LIBFF_USE_GSV
        };

        apply(ParticleGroup::regular);
        apply(ParticleGroup::spectator);
        // matrix has 0 length so there is no need to increment trajectory here
        // bunch.get_reference_particle().increment_trajectory(0);
        Kokkos::fence();
    }
}

#endif // FF_MATRIX_H
