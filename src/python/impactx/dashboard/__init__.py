from trame.widgets import html
from trame.widgets import vuetify3 as vuetify

# isort: off

from .server import setup_server

# Create single server instance for the entire dashboard
server, state, ctrl = setup_server()
from .Toolbar.general import GeneralToolbar

from .Analyze.ui import AnalyzeSimulation
from .Input.csr.ui import csrConfiguration
from .Input.isr.ui import isrConfiguration
from .Input.distribution.ui import DistributionParameters
from .Input.inputParameters.inputMain import InputParameters
from .Input.lattice.ui import LatticeConfiguration
from .Input.components.navigation import NavigationComponents
from .Input.space_charge_configuration.spaceChargeMain import SpaceChargeConfiguration

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
    "csrConfiguration",
    "isrConfiguration",
    "DistributionParameters",
    "InputParameters",
    "LatticeConfiguration",
    "SpaceChargeConfiguration",
    "GeneralToolbar",
]
