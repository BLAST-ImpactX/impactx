#ifndef FF_CONSTFOC_H
#define FF_CONSTFOC_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/libFF/ff_element.h"

class FF_constfoc : public FF_element
{
private:
    double get_reference_cdt(
            double length, double csl, double snl, double BL, double iBL, 
            Reference_particle & reference_particle);
public:
    FF_constfoc();
    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);
    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_constfoc();
};

typedef boost::shared_ptr<FF_constfoc > FF_constfoc_sptr;

#endif // FF_CONSTFOC_H
