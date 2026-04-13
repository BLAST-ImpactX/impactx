#!/usr/bin/env python

from synergia.foundation import Reference_particle, Four_momentum
import pytest
import numpy

charge = 1
mass = 2.2
total_energy = 3.0
total_energy2 = 5.1
four_momentum = Four_momentum(mass)
four_momentum.set_total_energy(total_energy)
state = numpy.array([1.1, 2.2, 3.3, 4.4, 5.5, 6.6], numpy.double)
step_length = 1.234
steps = 17
turns = 7


def test_construct():
    r = Reference_particle(charge, mass, total_energy)


def test_construct2():
    r = Reference_particle(charge, four_momentum)


def test_construct3():
    r = Reference_particle(charge, four_momentum, state)


def test_get_charge():
    r = Reference_particle(charge, four_momentum)
    assert r.get_charge() == charge


def test_set_get_four_momentum():
    r = Reference_particle(charge, four_momentum)
    # get_four_momentum() returns a const& in C++, and a copy of object
    # in Python. Therefore, set_total_energy() would not change the value
    # in the original reference particle
    r.get_four_momentum().set_total_energy(total_energy2)
    assert r.get_four_momentum().get_total_energy() == pytest.approx(total_energy)


def test_set_get_four_momentum2():
    r = Reference_particle(charge, four_momentum)
    f2 = Four_momentum(mass)
    f2.set_total_energy(total_energy2)
    assert r.get_four_momentum().get_total_energy() == pytest.approx(total_energy)
    r.set_four_momentum(f2)
    assert r.get_four_momentum().get_total_energy() == pytest.approx(total_energy2)


def test_set_get_state():
    r = Reference_particle(charge, four_momentum)
    r.set_state(state)
    new_state = r.get_state()
    for i in range(0, 6):
        assert new_state[i] == pytest.approx(state[i])


def test_set_get_total_energy():
    r = Reference_particle(charge, four_momentum)
    r.set_total_energy(total_energy2)
    assert r.get_total_energy() == pytest.approx(total_energy2)


def test_get_beta():
    r = Reference_particle(charge, four_momentum)
    assert r.get_beta() == pytest.approx(four_momentum.get_beta())


def test_get_gamma():
    r = Reference_particle(charge, four_momentum)
    assert r.get_gamma() == pytest.approx(four_momentum.get_gamma())


def test_get_momentum():
    r = Reference_particle(charge, four_momentum)
    assert r.get_momentum() == pytest.approx(four_momentum.get_momentum())


def test_set_get_trajectory():
    r = Reference_particle(charge, four_momentum)
    partial_s = 3 * step_length
    r.set_trajectory(turns, steps * step_length, partial_s)
    assert r.get_repetition() == turns
    assert r.get_repetition_length() == pytest.approx(steps * step_length)
    assert r.get_s_n() == pytest.approx(partial_s)
    assert r.get_s() == pytest.approx(turns * steps * step_length + partial_s)


def test_set_get_increment_abs_time():
    r = Reference_particle(charge, four_momentum)
    abstime = 1.0 / 1024.0  # exactly representable
    incrtime = 1 / 32768.0  # exactly representable
    r.set_bunch_abs_time(abstime)
    assert r.get_bunch_abs_time() == abstime
    r.increment_bunch_abs_time(incrtime)
    assert r.get_bunch_abs_time() == (abstime + incrtime)


def test_set_get_abs_offset():
    r = Reference_particle(charge, four_momentum)
    absoffset = 1 / 4096.0  # exactly representable
    r.set_bunch_abs_offset(absoffset)
    assert r.get_bunch_abs_offset() == absoffset
