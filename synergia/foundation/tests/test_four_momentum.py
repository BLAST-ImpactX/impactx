#!/usr/bin/env python

from synergia.foundation import Four_momentum
import pytest

mass = 3.0
total_energy = 27.4
value = 3.1415
beta_value = 0.678


def test_construct():
    f = Four_momentum(mass)


def test_construct2():
    f = Four_momentum(mass, total_energy)


def test_get_mass():
    f = Four_momentum(mass)
    assert f.get_mass() == pytest.approx(mass)


def test_set_get_total_energy():
    f = Four_momentum(mass)
    f.set_total_energy(value)
    assert f.get_total_energy() == pytest.approx(value)


def test_set_get_kinetic_energy():
    f = Four_momentum(mass)
    f.set_kinetic_energy(value)
    assert f.get_kinetic_energy() == pytest.approx(value)


def test_set_get_total_energy():
    f = Four_momentum(mass)
    f.set_total_energy(value)
    assert f.get_total_energy() == pytest.approx(value)


def test_set_get_beta():
    f = Four_momentum(mass)
    f.set_beta(beta_value)
    assert f.get_beta() == pytest.approx(beta_value)


def test_equal():
    tolerance = 1.0e-12
    f = Four_momentum(mass)
    f.set_beta(beta_value)
    f_same = Four_momentum(mass)
    f_same.set_beta(beta_value)
    f_different = Four_momentum(mass)
    f_different.set_beta(beta_value * 1.1)
    assert f.equal(f_same, tolerance)
    assert not f.equal(f_different, tolerance)
