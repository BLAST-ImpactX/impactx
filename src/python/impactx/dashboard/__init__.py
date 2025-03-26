from trame.widgets import html
from trame.widgets import vuetify3 as vuetify

# isort: off

from .trame_setup import setup_server
from .Toolbar.controls import GeneralToolbar

from .Analyze.plotsMain import AnalyzeSimulation
from .Input.csrConfiguration.csrMain import csrConfiguration
from .Input.distributionParameters.distributionMain import DistributionParameters
from .Input.inputParameters.inputMain import InputParameters
from .Input.latticeConfiguration.latticeMain import LatticeConfiguration
from .Input.components.navigation import NavigationComponents
from .Input.space_charge_configuration.spaceChargeMain import SpaceChargeConfiguration

from .jupyterApplication import JupyterMainApplication as JupyterApp
# isort: on


__all__ = [
    "html",
    "JupyterApp",
    "setup_server",
    "html",
    "vuetify",
    "AnalyzeSimulation",
    "NavigationComponents",
    "csrConfiguration",
    "DistributionParameters",
    "InputParameters",
    "LatticeConfiguration",
    "SpaceChargeConfiguration",
    "GeneralToolbar",
]
