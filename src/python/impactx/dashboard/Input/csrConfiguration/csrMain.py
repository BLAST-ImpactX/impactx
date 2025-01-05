from trame.widgets import vuetify

from ...Input.trameFunctions import TrameFunctions
from ...trame_setup import setup_server
from ..generalFunctions import generalFunctions

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Default State Variables
# -----------------------------------------------------------------------------

state.csr_bins = generalFunctions.get_default("csr_bins", "default_values")
state.csr_bins_error_message = ""

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------


@state.change("csr_bins")
def on_csr_bins_change(csr_bins, **kwargs):
    error_message = generalFunctions.validate_against(csr_bins, "int", ["positive"])
    state.csr_bins_error_message = error_message
    generalFunctions.update_simulation_validation_status()


# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------


class csrConfiguration:
    @staticmethod
    def card():
        """
        Creates UI content for CSR.
        """

        with vuetify.VCard(v_show="csr", style="width: 170px;"):
            with vuetify.VCardTitle("CSR"):
                vuetify.VSpacer()
                vuetify.VIcon(
                    "mdi-information",
                    classes="ml-2",
                    click=lambda: generalFunctions.documentation("CSR"),
                    style="color: #00313C;",
                )
            vuetify.VDivider()
            with vuetify.VCardText():
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(classes="py-0"):
                        TrameFunctions.select(
                            label="Particle Shape",
                            v_model_name="particle_shape",
                        )
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(classes="py-0"):
                        TrameFunctions.text_field(
                            label="CSR Bins", v_model_name="csr_bins"
                        )
