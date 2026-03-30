#!/usr/bin/env python
import sys
import os
import numpy as np
import synergia
import pytest


def test_map_to_twiss():
    # list of alpha, beta, qs to test
    twiss_list = [
        (0.1, 22.0, 0.1257),
        (-0.1, 24.0, 0.87323),
        (0.04, 933.0, 0.005),
        (0.04017325428125798, 974.397979726345, 0.023655638652905025),
        (-0.57, 974.397979726345, 0.023655638652905025),
    ]

    for t in twiss_list:
        a = t[0]
        b = t[1]
        mu = 2.0 * np.pi * t[2]
        map = np.array(
            [
                [np.cos(mu) + a * np.sin(mu), b * np.sin(mu)],
                [-(1 + a**2) * np.sin(mu) / b, np.cos(mu) - a * np.sin(mu)],
            ]
        )

        # print(f'alpha: {a}, beta: {b}, Qs: {t[2]}')
        # print('map: ', map)
        # check determinant
        det = map[0, 0] * map[1, 1] - map[0, 1] * map[1, 0]
        assert det == pytest.approx(1.0, abs=1e-3)

        atest, btest, qtest = synergia.simulation.Lattice_simulator.map_to_twiss(map)
        # print(f'determined values from map: alpha: {atest}, beta: {btest}, Qs: {qtest}')
        assert a == pytest.approx(atest)
        assert b == pytest.approx(btest)
        assert qtest == pytest.approx(t[2])


def test_booster_map():
    fullmap_list = []
    map_names = []

    map_names.append("booster-00000")
    fullmap_list.append(
        np.array(
            [
                [
                    -2.35859774e-01,
                    -3.26696182e01,
                    0.00000000e00,
                    0.00000000e00,
                    2.45055034e-03,
                    3.71103891e00,
                ],
                [
                    2.87168479e-02,
                    -2.61017827e-01,
                    0.00000000e00,
                    0.00000000e00,
                    1.40666320e-06,
                    -9.05527355e-02,
                ],
                [
                    0.00000000e00,
                    0.00000000e00,
                    3.62462342e-01,
                    -4.88366765e00,
                    0.00000000e00,
                    0.00000000e00,
                ],
                [
                    0.00000000e00,
                    0.00000000e00,
                    1.77021951e-01,
                    3.73786762e-01,
                    0.00000000e00,
                    0.00000000e00,
                ],
                [
                    -1.09309492e-01,
                    4.61051782e00,
                    0.00000000e00,
                    0.00000000e00,
                    9.74429598e-01,
                    -1.35798796e02,
                ],
                [
                    -3.96372273e-06,
                    3.05883524e-03,
                    0.00000000e00,
                    0.00000000e00,
                    7.98388448e-04,
                    9.14701148e-01,
                ],
            ]
        )
    )

    map_names.append("booster-00001")
    fullmap_list.append(
        np.array(
            [
                [
                    -2.37127010e-01,
                    -3.26549541e01,
                    0.00000000e00,
                    0.00000000e00,
                    2.43070643e-03,
                    3.71815633e00,
                ],
                [
                    2.87077309e-02,
                    -2.62185732e-01,
                    0.00000000e00,
                    0.00000000e00,
                    1.38661631e-06,
                    -9.05361380e-02,
                ],
                [
                    0.00000000e00,
                    0.00000000e00,
                    3.61251662e-01,
                    -4.88662790e00,
                    0.00000000e00,
                    0.00000000e00,
                ],
                [
                    0.00000000e00,
                    0.00000000e00,
                    1.77087341e-01,
                    3.72390550e-01,
                    0.00000000e00,
                    0.00000000e00,
                ],
                [
                    -1.09318281e-01,
                    4.61679399e00,
                    0.00000000e00,
                    0.00000000e00,
                    9.74654822e-01,
                    -1.35807288e02,
                ],
                [
                    -3.92571435e-06,
                    3.03291015e-03,
                    0.00000000e00,
                    0.00000000e00,
                    7.91548892e-04,
                    9.15289396e-01,
                ],
            ]
        )
    )

    map_names.append("booster-01000")
    fullmap_list.append(
        np.array(
            [
                [
                    -2.36525139e-01,
                    -3.26619064e01,
                    0.00000000e00,
                    0.00000000e00,
                    2.44514656e-03,
                    3.71423486e00,
                ],
                [
                    2.87120666e-02,
                    -2.61630993e-01,
                    0.00000000e00,
                    0.00000000e00,
                    1.39802933e-06,
                    -9.05449942e-02,
                ],
                [
                    0.00000000e00,
                    0.00000000e00,
                    3.61826369e-01,
                    -4.88522393e00,
                    0.00000000e00,
                    0.00000000e00,
                ],
                [
                    0.00000000e00,
                    0.00000000e00,
                    1.77056355e-01,
                    3.73053252e-01,
                    0.00000000e00,
                    0.00000000e00,
                ],
                [
                    -1.09313464e-01,
                    4.61363722e00,
                    0.00000000e00,
                    0.00000000e00,
                    9.74494860e-01,
                    -1.35798153e02,
                ],
                [
                    -3.95190205e-06,
                    3.05146867e-03,
                    0.00000000e00,
                    0.00000000e00,
                    7.96430995e-04,
                    9.14835680e-01,
                ],
            ]
        )
    )

    # print(f'len(fullmap_list): {len(fullmap_list)}')
    # print(type(fullmap_list))

    for fullmap, nm in zip(fullmap_list, map_names):
        print("testing map ", nm)
        fullmap[5, 4] = -fullmap[5, 4]
        fullmap[4, 5] = -fullmap[4, 5]
        # print(f'alpha: {a}, beta: {b}, Qs: {t[2]}')
        # print('map: ', map)
        # check determinant
        det = fullmap[4, 4] * fullmap[4, 4] - fullmap[4, 5] * fullmap[5, 4]
        # print('det: ', det)
        assert det == pytest.approx(1.0, abs=6e-2)

        atest, btest, qtest = synergia.simulation.Lattice_simulator.map_to_twiss(
            fullmap[4:6, 4:6]
        )
        # print(f'determined values from map: alpha: {atest}, beta: {btest}, Qs: {qtest}')
        assert atest == pytest.approx(0.0909, 5.0e-3)
        assert btest == pytest.approx(413.0, 1.0e-2)
        qtest == pytest.approx(0.053, 1.0e-3)


def test_map_to_twiss_slice():
    betaz = 974.0
    qs = 0.0244
    mu = 2.0 * np.pi * qs
    map1 = np.array(
        [[np.cos(mu), betaz * np.sin(mu)], [-np.sin(mu) / betaz, np.cos(mu)]]
    )

    a, b, c = synergia.simulation.Lattice_simulator.map_to_twiss(map1)
    assert a == pytest.approx(0.0)
    assert b == pytest.approx(betaz)
    assert c == pytest.approx(qs)

    map2 = np.zeros((6, 6))
    for i in range(6):
        for j in range(6):
            map2[i, j] = 10 * (i + 1) + j + 1

    map2a = np.array(map2)
    map2a[0:2, 0:2] = map1[:, :]
    print(np.array2string(map2a, separator=","))
    a, b, c = synergia.simulation.Lattice_simulator.map_to_twiss(map2a[0:2, 0:2])
    assert a == pytest.approx(0.0)
    assert b == pytest.approx(betaz)
    assert c == pytest.approx(qs)


if __name__ == "__main__":
    test_map_to_twiss_slice()
    test_map_to_twiss()
    test_booster_map()
