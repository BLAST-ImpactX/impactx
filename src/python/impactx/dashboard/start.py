"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

import asyncio

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

        # Ensure an event loop exists for Python 3.10+ (required on macOS)
        # This prevents RuntimeError when server.start() tries to get_event_loop()
        # In Python 3.10+, get_event_loop() raises RuntimeError if no loop exists
        # We create and set one if needed
        try:
            asyncio.get_event_loop()
        except RuntimeError:
            # No event loop exists in this thread, create and set one
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)

        server.start()
        return 0


class JupyterApp:
    """
    Jupyter-compatible version of the dashboard.
    """

    def __init__(self):
        self.ui = setup_dashboard()
