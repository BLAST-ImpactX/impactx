#!/usr/bin/env python
# -*- coding: utf-8 -*-

import synergia

no_simplify_attribute = 'no_simplify'
extractor_type_attribute = 'extractor_type'

def no_simplify(element):
    retval = False
    if element.has_string_attribute(no_simplify_attribute):
        if not element.get_string_attribute(no_simplify_attribute).lower == 'false':
            retval = True
    if element.has_string_attribute(extractor_type_attribute):
        retval = True
    return retval

def _copy_attributes_ancestors(element, new_element):
    new_element.copy_attributes_from(element)
    for ancestor in element.get_ancestors():
        new_element.add_ancestor(ancestor)

def eliminate_markers(lattice):
    new_lattice = synergia.lattice.Lattice(lattice.get_name())
    new_lattice.set_reference_particle(lattice.get_reference_particle())
    for element in lattice.get_elements():
            if (element.get_type_name() == 'marker') and \
                (not no_simplify(element)):
                pass
            else:
                new_lattice.append(element)
    return new_lattice

def convert_monitors(lattice):
    new_lattice = synergia.lattice.Lattice(lattice.get_name())
    new_lattice.set_reference_particle(lattice.get_reference_particle())
    for element in lattice.get_elements():
            if ((element.get_type_name() == 'monitor') or \
                (element.get_type_name() == 'hmonitor') or \
                (element.get_type_name() == 'vmonitor')) and \
                (not no_simplify(element)):
                new_element = synergia.lattice.Lattice_element('drift', element.get_name())
                _copy_attributes_ancestors(element, new_element)
                new_lattice.append(new_element)
            else:
                new_lattice.append(element)
    return new_lattice

def _is_magnet(type):
    suffix = 'pole'
    return type.find(suffix) == (len(type) - len(suffix))

def convert_magnets(lattice):
    new_lattice = synergia.lattice.Lattice(lattice.get_name())
    new_lattice.set_reference_particle(lattice.get_reference_particle())
    for element in lattice.get_elements():
            if _is_magnet(element.get_type_name()) and \
                (not no_simplify(element)):
                has_nonzero_strength = False
                for attr in element.get_double_attribute_names():
                    if (len(attr) == 1) and (attr == 'k'):
                        if element.get_double_attribute(attr) != 0.0:
                            has_nonzero_strength = True
                    else:
                        if (attr[0] == 'k') and attr[1].isdigit():
                            if element.get_double_attribute(attr) != 0.0:
                                has_nonzero_strength = True
                if has_nonzero_strength:
                    new_lattice.append(element)
                else:                    
                    new_element = synergia.lattice.Lattice_element('drift', element.get_name())
                    _copy_attributes_ancestors(element, new_element)
                    new_lattice.append(new_element)
            else:
                new_lattice.append(element)
    return new_lattice

def convert_kickers(lattice):
    new_lattice = synergia.lattice.Lattice(lattice.get_name())
    new_lattice.set_reference_particle(lattice.get_reference_particle())
    for element in lattice.get_elements():
            if ((element.get_type_name() == 'hicker') or \
                (element.get_type_name() == 'vkicker')) and \
                (not no_simplify(element)):
                if element.get_double_attribute('kick') == 0.0:
                    new_element = synergia.lattice.Lattice_element('drift', element.get_name())
                    _copy_attributes_ancestors(element, new_element)
                    new_lattice.append(new_element)
                else:
                    new_lattice.append(element)
            elif (element.get_type_name() == 'kicker') and \
                (not no_simplify(element)):
                if (element.get_double_attribute('hkick') == 0.0) and \
                    (element.get_double_attribute('vkick') == 0.0):
                    new_element = synergia.lattice.Lattice_element('drift', element.get_name())
                    _copy_attributes_ancestors(element, new_element)
                    new_lattice.append(new_element)                    
                else:
                    new_lattice.append(element)
            else:
                new_lattice.append(element)
    return new_lattice

def combine_drifts(lattice):
    new_lattice = synergia.lattice.Lattice(lattice.get_name())
    new_lattice.set_reference_particle(lattice.get_reference_particle())
    last_element = None
    for element in lattice.get_elements():
        if last_element:
            last_l = last_element.get_length()
            ancestors = last_element.get_ancestors()
            if (element.get_type_name() == 'drift') and \
                (not no_simplify(element)):
                concat_name = last_element.get_name()
                concat_name += '+'
                concat_name += element.get_name()
                this_l = element.get_length()
                last_element = synergia.lattice.Lattice_element('drift', concat_name)
                last_element.set_double_attribute('l', last_l + this_l)
                for ancestor in ancestors:
                    last_element.add_ancestor(ancestor)
            else:
                new_lattice.append(last_element)
                new_lattice.append(element)
                last_element = None
        else:
            if element.get_type_name() == 'drift':
                last_element = element
            else:
                new_lattice.append(element)
                last_element = None
    if last_element:
        new_lattice.append(last_element)
    return new_lattice

def simplify_all(lattice):
    return \
        combine_drifts(\
            convert_kickers(\
                convert_magnets(\
                    convert_monitors(\
                        eliminate_markers(lattice)))))

