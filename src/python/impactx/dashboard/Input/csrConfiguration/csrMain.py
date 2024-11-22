from trame.widgets import html, vuetify

from ...trame_setup import setup_server
from ..generalFunctions import generalFunctions

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Default State Variables
# -----------------------------------------------------------------------------

state.csr_bins = 150
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
                with vuetify.VTooltip(bottom=True, nudge_top="10"):
                    with vuetify.Template(v_slot_activator="{ on, attrs }"):
                        vuetify.VSelect(
                            label="Particle Shape",
                            v_model=("particle_shape",),
                            items=([1, 2, 3],),
                            dense=True,
                            v_on="on",
                            v_bind="attrs",
                        )
                    html.Span("{{ parameter_tooltips.particle_shape }}")
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(classes="py-0"):
                        with vuetify.VTooltip(bottom=True, nudge_top="10"):
                            with vuetify.Template(v_slot_activator="{ on, attrs }"):
                                vuetify.VTextField(
                                    label="CSR Bins",
                                    v_model=("csr_bins",),
                                    error_messages=("csr_bins_error_message",),
                                    type="number",
                                    dense=True,
                                    v_on="on",
                                    v_bind="attrs",
                                )
                            html.Span("{{ parameter_tooltips.csr_bins }}")
