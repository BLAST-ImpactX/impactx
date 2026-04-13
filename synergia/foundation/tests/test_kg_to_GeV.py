#!/usr/bin/env python

import synergia.foundation
import pytest


def test_kg_to_Gev():
    kg_to_GeV_derived = (
        1.0e-9 * synergia.foundation.pconstants.c**2 / synergia.foundation.pconstants.e
    )
    print("kg_to_GeV_derived: ", kg_to_GeV_derived)
    print("pconstants.kg_to_GeV: ", synergia.foundation.pconstants.kg_to_GeV)
    assert kg_to_GeV_derived == pytest.approx(
        synergia.foundation.pconstants.kg_to_GeV, rel=1.0e-15
    )


if __name__ == "__main__":
    test_kg_to_Gev()
