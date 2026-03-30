#ifndef FF_ELEMENT_H
#define FF_ELEMENT_H

#include "synergia/lattice/lattice_element_slice.h"
#include "synergia/lattice/lattice_element.h"

#include "synergia/libFF/ff_drift.h"
#include "synergia/libFF/ff_sbend.h"
#include "synergia/libFF/ff_quadrupole.h"
#include "synergia/libFF/ff_multipole.h"
#include "synergia/libFF/ff_sextupole.h"
#include "synergia/libFF/ff_octupole.h"
#include "synergia/libFF/ff_kicker.h"
#include "synergia/libFF/ff_solenoid.h"
#include "synergia/libFF/ff_rfcavity.h"
#include "synergia/libFF/ff_elens.h"
#include "synergia/libFF/ff_nllens.h"
#include "synergia/libFF/ff_foil.h"
#include "synergia/libFF/ff_dipedge.h"
#include "synergia/libFF/ff_matrix.h"

namespace FF_element
{
    template<class BUNCH>
    void apply(Lattice_element_slice const& slice, BUNCH & b)
    {
        auto const& elm = slice.get_lattice_element();
        auto t = elm.get_type();

        switch(t)
        {
        case element_type::drift:      FF_drift::apply(slice, b); break; 
        case element_type::sbend:      FF_sbend::apply(slice, b); break;
        case element_type::quadrupole: FF_quadrupole::apply(slice, b); break;

        case element_type::multipole:  FF_multipole::apply(slice, b); break;
        case element_type::sextupole:  FF_sextupole::apply(slice, b); break;
        case element_type::octupole:   FF_octupole::apply(slice, b); break;

        case element_type::hkicker:    FF_kicker::apply(slice, b); break;
        case element_type::vkicker:    FF_kicker::apply(slice, b); break;
        case element_type::kicker:     FF_kicker::apply(slice, b); break;

        case element_type::solenoid:   FF_solenoid::apply(slice, b); break;
        case element_type::rfcavity:   FF_rfcavity::apply(slice, b); break;
        case element_type::elens:      FF_elens::apply(slice, b); break;
        case element_type::nllens:     FF_nllens::apply(slice, b); break;

        case element_type::monitor:    FF_drift::apply(slice, b); break;
        case element_type::hmonitor:   FF_drift::apply(slice, b); break;
        case element_type::vmonitor:   FF_drift::apply(slice, b); break;
        case element_type::marker:     FF_drift::apply(slice, b); break;
        case element_type::instrument: FF_drift::apply(slice, b); break;
        case element_type::rcollimator:FF_drift::apply(slice, b); break;

        case element_type::foil:       FF_foil::apply(slice, b); break;
        case element_type::dipedge:    FF_dipedge::apply(slice, b); break;
        case element_type::matrix:     FF_matrix::apply(slice, b); break;

        default: 
            throw std::runtime_error(
                    "FF_element::apply() unknown element type = " + 
                        elm.get_type_name() +
                    ", element name = " + 
                        elm.get_name() );
        }
    }

    template<class BUNCH>
    void apply(Lattice_element const& element, BUNCH & b)
    { apply(Lattice_element_slice(element), b); }
};



#if 0
class JetParticle;
class FF_element
{
public:
    FF_element()
    : steps(default_steps), order(default_order) { };

    virtual ~FF_element() = default;

#if 0
    template<size_t I>
    virtual void apply(Lattice_element_slice const& slice, 
            Trigon_particle<I> & trigon)
    { }
#endif

    virtual void apply(Lattice_element_slice const& slice, 
            Bunch & bunch)
    { }

#if 0
    Reference_particle &
        get_ref_particle_from_slice(Lattice_element_slice & slice) const
    { return slice.get_lattice_element().get_lattice().get_reference_particle(); }
#endif

    Reference_particle const &
        get_ref_particle_from_slice(Lattice_element_slice const & slice) const
    { return slice.get_lattice_element().get_lattice().get_reference_particle(); }

    void set_yoshida_steps(int s)
    { steps = s; }

    void set_yoshida_order(int o)
    { order = o; }

    int get_yoshida_steps() const
    { return steps; }

    int get_yoshida_order() const
    { return order; }


    // static members
    static void set_default_yoshida_steps(int s)
    { default_steps = s; }

    static void set_default_yoshida_order(int o)
    { default_order = o; }

    static int get_default_yoshdia_steps()
    { return default_steps; }

    static int get_default_yoshida_order()
    { return default_order; }

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);

protected:

    static int default_steps;
    static int default_order;

    int steps;
    int order;
};
#endif

#endif // FF_ELEMENT_H
