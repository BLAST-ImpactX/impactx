from trame.widgets import html

from ... import setup_server, vuetify
from .. import (
    CardComponents,
    InputComponents,
    NavigationComponents,
    generalFunctions,
)

server, state, ctrl = setup_server()


init_value = ""
state.variables = [
    {"name": init_value, "value": init_value, "error_message": init_value}
]
state.is_only_variable = len(state.variables) == 1


class LatticeVariableHandler:
    """
    Stores all functionality for dashboard variable referencing.
    """

    @state.change("variables")
    def on_variables_list_change(variables, **kwargs):
        for lattice_element in state.lattice_params_bound_or_pending_variable.values():
            ctrl.updateLatticeElementParameters(
                lattice_element["index"],
                lattice_element["parameter_name"],
                lattice_element["ui_value"],
                lattice_element["parameter_type"],
            )

    # -----------------------------------------------------------------------------
    # Controllers
    # -----------------------------------------------------------------------------

    @ctrl.add("add_variable")
    def on_add_change() -> None:
        """
        Adds a new variable to the dashboard's variable
        with empty values and updates UI.
        Stored in a state which contains a list with dictionaries.
        """

        new_variable = {key: "" for key in state.variables[0]}
        state.variables.append(new_variable)
        state.dirty("variables")
        LatticeVariableHandler.update_delete_availability()

    @ctrl.add("delete_variable")
    def on_delete_change(index: int) -> None:
        """
        Deleted the variable defined by the user
        provided the index

        :param index: The index of the variable
        """

        state.variables.pop(index)
        state.dirty("variables")
        LatticeVariableHandler.update_delete_availability()

    @ctrl.add("update_variable")
    def on_variable_change(key_name: str, index: int, event) -> None:
        """
        Called when a variable name or value changes.
        Validates the value and updates it's stored value.

        :param key_name: The name of the variable.
        :param index: The index of the variable.
        :param event: Either the variable's new name or value.
        """

        if key_name == "name":
            LatticeVariableHandler.validate_variable_name(event, index)
            state.variables[index]["name"] = event
        else:
            state.variables[index][key_name] = event
        state.dirty("variables")

    @ctrl.add("reset_variables")
    def on_reset_variables() -> None:
        """
        Resets the dashboard's variables to default.
        """

        state.variables = [
            {"name": init_value, "value": init_value, "error_message": init_value}
        ]
        state.dirty("variables")
        LatticeVariableHandler.update_delete_availability()

    # -----------------------------------------------------------------------------
    # Methods
    # -----------------------------------------------------------------------------

    @staticmethod
    def update_delete_availability() -> None:
        """
        Updates the state flag that controls whether the delete variable
        functionality should be disabled. The delete functionality is disabled
        when there is only one variable in the list.
        """

        state.is_only_variable = True if len(state.variables) == 1 else False
        state.dirty("is_only_variable")

    @staticmethod
    def get_duplicate_indexes(new_name: str, current_index: int) -> list:
        """
        Returns the indexes of duplicate variable names.

        :param new_name: The name of the variable.
        :current_index: The index of the variable.
        """

        duplicates = [
            index
            for index, var in enumerate(state.variables)
            if var["name"] == new_name and index != current_index
        ]

        if duplicates:
            duplicates.append(current_index)
        return duplicates

    @staticmethod
    def element_potentially_using_element(ui_value, _dummy_index=None):
        return bool(ui_value) and ui_value[0].isalpha()

    @staticmethod
    def validate_variable_name(new_name, index) -> None:
        """
        Validates the variable name and outputs an error message if any.

        :param new_name: The name of the variable.
        :index: The index of the variable.
        """

        def set_var_error_message(message):
            state.variables[index]["error_message"] = message
            state.dirty("variables")

        alpha = new_name and new_name[0].isalpha()
        duplicate_indexes = LatticeVariableHandler.get_duplicate_indexes(
            new_name, index
        )
        send_error = set_var_error_message("error")

        if not alpha:
            send_error
            generalFunctions.update_simulation_validation_status()  # need to optimize function later
            state.dirty("variables")
            return True
        elif duplicate_indexes:
            for index in duplicate_indexes:
                send_error
            generalFunctions.update_simulation_validation_status()  # need to optimize function later
            state.dirty("variables")
            return True
        else:
            set_var_error_message("")
            generalFunctions.update_simulation_validation_status()  # need to optimize function later

    @staticmethod
    def determine_if_variable(var_name):
        """
        Determines if var_name is already a variable in the current
        list of variables.

        :param: var_name: The name of the variable.
        :return: A bool and [if found] the index of the variable.
        """

        found_index = next(
            (i for i, var in enumerate(state.variables) if var["name"] == var_name),
            None,
        )
        return (found_index is not None, found_index)

    # -----------------------------------------------------------------------------
    # UI
    # -----------------------------------------------------------------------------

    @staticmethod
    def dialog_settings():
        dialog_name = "lattice_configuration_dialog_tab_settings"

        with NavigationComponents.create_dialog_tabs(
            dialog_name, 2, ["Variables", "Defaults"]
        ):
            with vuetify.VTabsWindow(v_model=(dialog_name, 0)):
                with vuetify.VTabsWindowItem():
                    with vuetify.VCardText():
                        with vuetify.VContainer(fluid=True):
                            with vuetify.VRow(
                                v_for="(variable, index) in variables",
                                classes="align-center justify-center py-0",
                            ):
                                with vuetify.VCol(cols=5, classes="pr-0"):
                                    vuetify.VTextField(
                                        placeholder="Name",
                                        v_model=("variable.name",),
                                        variant="outlined",
                                        density="compact",
                                        background_color="grey lighten-4",
                                        update_modelValue=(
                                            ctrl.update_variable,
                                            "['name', index, $event]",
                                        ),
                                        error_messages=("variable.error_message", []),
                                        hide_details=True,
                                        clearable=True,
                                    )
                                with vuetify.VCol(cols=1, classes="px-0 text-center"):
                                    html.Span("=", classes="mx-0")
                                with vuetify.VCol(cols=4, classes="pl-0"):
                                    vuetify.VTextField(
                                        placeholder="Value",
                                        v_model=("variable.value",),
                                        variant="outlined",
                                        density="compact",
                                        type="number",
                                        background_color="grey lighten-4",
                                        update_modelValue=(
                                            ctrl.update_variable,
                                            "['value', index, $event]",
                                        ),
                                        hide_details=True,
                                        clearable=True,
                                    )
                                with vuetify.VCol(cols=2, classes="d-flex"):
                                    with html.Div(classes="mr-2"):
                                        CardComponents.card_button(
                                            "mdi-plus",
                                            color="primary",
                                            description="Add Variable",
                                            click=ctrl.add_variable,
                                            v_show="index === variables.length - 1",
                                            density="default",
                                            size="x-small",
                                            variant="elevated",
                                        )
                                    with html.Div():
                                        CardComponents.card_button(
                                            "mdi-delete",
                                            color="secondary",
                                            description="Delete Variable",
                                            click=(ctrl.delete_variable, "[index]"),
                                            disabled=("is_only_variable",),
                                            density="default",
                                            size="x-small",
                                            variant="elevated",
                                        )
                            with vuetify.VRow(classes="mt-2"):
                                with vuetify.VCol():
                                    vuetify.VBtn(
                                        "Reset Variables",
                                        color="accent",
                                        click=ctrl.reset_variables,
                                        block=True,
                                    )
                with vuetify.VTabsWindowItem():
                    with vuetify.VCardText():
                        with vuetify.VRow():
                            with vuetify.VCol(cols=3):
                                InputComponents.text_field(
                                    label="nslice",
                                    v_model_name="nslice",
                                    change=(
                                        ctrl.nsliceDefaultChange,
                                        "['nslice', $event]",
                                    ),
                                )
