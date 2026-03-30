#!/usr/bin/env python
import pytest
import synergia
from synergia.lattice import Lattice, Lattice_element
from synergia.foundation import Reference_particle

def test_set_get_remove_double_attribute():
    le = Lattice_element("drift", "d")
    le.set_double_attribute('foo', 2.0)
    assert le.has_double_attribute('foo')
    assert le.get_double_attribute('foo') == 2.0
    le.remove_double_attribute('foo')
    assert not le.has_double_attribute('foo')
    pass

def test_set_get_remove_string_attribute():
    le = Lattice_element("drift", "d")
    le.set_string_attribute('foo', 'gerbil')
    assert le.has_string_attribute('foo')
    assert le.get_string_attribute('foo') == 'gerbil'
    le.remove_string_attribute('foo')
    assert not le.has_string_attribute('foo')
    pass

def test_set_get_remove_vector_attribute():
    le = Lattice_element("drift", "d")
    le.set_vector_attribute('foov', [1.0, 2.0, 3.0])
    assert le.has_vector_attribute('foov')
    vattr = le.get_vector_attribute('foov')
    assert len(vattr) == 3
    assert vattr[0] == 1.0
    assert vattr[1] == 2.0
    assert vattr[2] == 3.0
    le.remove_vector_attribute('foov')
    assert not le.has_vector_attribute('foov')
    pass
