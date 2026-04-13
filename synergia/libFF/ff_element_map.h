#ifndef FF_ELEMENT_MAP_H
#define FF_ELEMENT_MAP_H

#include <map>
#include <string>

#include "synergia/libFF/ff_element.h"

class FF_element_map
{
private:
    std::map<std::string, FF_element_sptr > element_map;

public:
    FF_element_map();
    bool has_element_type(std::string const& type) const;
    void set_element_type(std::string const& type, FF_element_sptr element_sptr);
    FF_element_sptr get_element_type(std::string const& type) const;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    ~FF_element_map();
};

void construct_big_giant_global_ff_element_map();
extern FF_element_map the_big_giant_global_ff_element_map;
#endif // FF_ELEMENT_MAP_H
