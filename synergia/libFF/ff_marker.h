#ifndef FF_MARKER_H
#define FF_MARKER_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/libFF/ff_element.h"

class FF_marker : public FF_element
{
private:
    double get_reference_cdt(double length, Reference_particle & reference_particle);
public:
    FF_marker();
    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);
    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_marker();
};

typedef boost::shared_ptr<FF_marker > FF_marker_sptr;

#endif // FF_MARKER_H
