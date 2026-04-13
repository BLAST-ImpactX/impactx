#!/usr/bin/env python
import synergia
import pytest

lattice_empty_matrix = """
m: matrix, l=0;
s: sequence, l=0;
m, at=0.0;
endsequence;
"""


lattice_matrix = """
m: matrix, l=0, rm11=0.5, rm12=-.25, rm64=3.53125, tm253=-42.21875;
s: sequence, l=0;
m, at=0.0;
endsequence;
"""


def test_parse_empty_matrix():
    reader = synergia.lattice.MadX_reader()
    # can we parse the empty matrix
    reader.parse(lattice_empty_matrix)
    # can we get the lattice
    lattice = reader.get_lattice("s")

    # Is the element present in the lattice?
    found_name = False
    found_type_name = False
    found_type = False

    # print(lattice)

    for elem in lattice.get_elements():
        if elem.get_name() == "m":
            found_name = True
            print(elem)

            # check the element is what we think it is
            tn = "matrix"
            if elem.get_type_name() == tn:
                found_type_name = True

            if elem.get_type() == synergia.lattice.element_type.matrix:
                found_type = True

    assert found_name
    assert found_type_name
    assert found_type


def test_parse_arb_matrix():
    reader = synergia.lattice.MadX_reader()
    # can we parse the empty matrix
    reader.parse(lattice_matrix)
    # can we get the lattice
    lattice = reader.get_lattice("s")
    print(lattice)

    # Find the element in the lattice?

    for elem in lattice.get_elements():
        if elem.get_name() == "m":
            # check the attributes
            assert elem.get_double_attribute("rm11") == 0.5
            assert elem.get_double_attribute("rm12") == -0.25
            assert elem.get_double_attribute("rm64") == 3.53125
            assert elem.get_double_attribute("tm253") == -42.21875


if __name__ == "__main__":
    # test_parse_empty_matrix()
    test_parse_arb_matrix()
