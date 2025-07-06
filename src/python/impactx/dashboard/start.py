"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from . import server, state
from .app import application
from .Input.defaults import DashboardDefaults
from .Toolbar.sim_history.ui import load_my_js

# -----------------------------------------------------------------------------
# Core setup logic
# -----------------------------------------------------------------------------


def initialize_states():
    """
    Initializes all states with default values.
    """
    for name, value in DashboardDefaults.DEFAULT_VALUES.items():
        setattr(state, name, value)


def setup_dashboard():
    initialize_states()
    load_my_js(server)
    return application()


# -----------------------------------------------------------------------------
# Application classes
# -----------------------------------------------------------------------------


class DashboardApp:
    """
    Full ImpactX Dashboard app.
    """

    def start(self):
        setup_dashboard()
        server.start()
        return 0


class JupyterApp:
    """
    Jupyter-compatible version of the dashboard.
    """

    def __init__(self):
        self.ui = setup_dashboard()
