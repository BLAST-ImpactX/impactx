
#ifndef LATTICE_ELEMENT_PROCESSOR_H
#define LATTICE_ELEMENT_PROCESSOR_H

#include "synergia/lattice/lattice_element.h"


// Pre-process the lattice element so it has necessary
// attributes in order to be used in the progation.
//
// this processing needs to add all the default attributes
// to the element as previously constructed in the element
// adaptors, as well as to derive the internal attributes
// such as the "arclength" of a rbend element.
//
// the derivation of external attributes (e.g., rf frequency)
// is postponed since the lattice might not be completed
// at the moment of calling the process() method.

class Lattice_element_processor
{
public:

    static Lattice_element
        process(Lattice_element const & element);

private:

    static void drift(Lattice_element & element);
    static void sbend(Lattice_element & element);
    static void quadrupole(Lattice_element & element);
    static void multipole(Lattice_element & element);
    static void rfcavity(Lattice_element & element);

    static void hkicker(Lattice_element & element);
    static void vkicker(Lattice_element & element);
    static void kicker(Lattice_element & element);
};

#endif
