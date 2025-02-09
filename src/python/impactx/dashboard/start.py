"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from . import setup_server
from .Input.defaults import DashboardDefaults

server, state, ctrl = setup_server()


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def initialize_states():
    """
    Initializes all dashboard state values upon call.

    The issue as of now is it initialize all at once instead of by section.
    """
    for name, value in DashboardDefaults.DEFAULT_VALUES.items():
        setattr(state, name, value)


def main():
    """
    Launches Trame application server
    """
    initialize_states()
    server.start()
    return 0
