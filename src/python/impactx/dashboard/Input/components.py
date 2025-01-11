from . import generalFunctions, vuetify


class CardComponents:
    """
    Class contains staticmethods to build
    card components using Vuetify's VCard.
    """

    @staticmethod
    def input_header(section_name: str, additional_components=None) -> None:
        """
        Creates a standardized header look for inputs.

        :param section_name: The name for the input section.
        """

        documentation_name = section_name.lower().replace(" ", "_")
        with vuetify.VCardTitle(section_name):
            vuetify.VSpacer()
            if additional_components:
                additional_components()
            CardComponents.refresh_icon(documentation_name)
            CardComponents.documentation_icon(documentation_name)
        vuetify.VDivider()

    @staticmethod
    def documentation_icon(section_name: str) -> vuetify.VIcon:
        """
        Takes user to input section's documentation.

        :param section_name: The name for the input section.
        """

        return vuetify.VIcon(
            "mdi-information",
            style="color: #00313C;",
            click=lambda: generalFunctions.documentation(section_name),
        )

    @staticmethod
    def refresh_icon(section_name: str) -> vuetify.VIcon:
        """
        Resets input values to default.

        :param section_name: The name for the input section.
        """

        return vuetify.VIcon(
            "mdi-refresh",
            style="color: #00313C;",
            click=lambda: generalFunctions.reset_inputs(section_name),
        )
