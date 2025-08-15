from trame.widgets import html
from trame.widgets import vuetify3 as vuetify

# isort: off

from .server import setup_server

# Create single server instance for the entire dashboard
server, state, ctrl = setup_server()
from .Toolbar.general import GeneralToolbar

from .Analyze.ui import AnalyzeSimulation
from .Input.csr import CSRConfiguration
from .Input.isr import ISRConfiguration
from .Input.distribution import DistributionConfiguration
from .Input.simulation_parameters import SimulationParameters
from .Input.lattice import LatticeConfiguration
from .Input.components.navigation import NavigationComponents
from .Input.space_charge import SpaceChargeConfiguration

from .start import JupyterApp
# isort: on


__all__ = [
    "html",
    "JupyterApp",
    "setup_server",
    "server",
    "state",
    "ctrl",
    "html",
    "vuetify",
    "AnalyzeSimulation",
    "NavigationComponents",
    "CSRConfiguration",
    "ISRConfiguration",
    "DistributionConfiguration",
    "SimulationParameters",
    "LatticeConfiguration",
    "SpaceChargeConfiguration",
    "GeneralToolbar",
]
