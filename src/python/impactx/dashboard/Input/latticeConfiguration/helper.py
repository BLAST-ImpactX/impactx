from trame.widgets import html


from .. import generalFunctions, InputComponents, NavigationComponents, setup_server, vuetify

server, state, ctrl = setup_server()


init_value = ""
state.variables = [{"name": init_value, "value": init_value, "error_message": init_value}]
state.is_only_variable = len(state.variables) == 1
state.toolbar_settings = False


@ctrl.add("add_variable")
def on_add_change():
    new_variable = { key: "" for key in state.variables[0] }
    state.variables.append(new_variable)
    state.dirty("variables")
    LatticeVariableHandler.update_delete_availability()

@ctrl.add("delete_variable")
def on_delete_change(index) -> None:
    """
    Deleted the variable defined by the user
    provided the index

    :param index: The index of the variable
    """
    print(f"index is {index}")
    state.variables.pop(index)
    state.dirty("variables")
    LatticeVariableHandler.update_delete_availability()
    print(f"Deleted variable at index {index}. Updated list: {state.variables}")

@ctrl.add("update_variable")
def on_variable_change(key_name: str, index: int, event):
    if key_name == "name":
        if LatticeVariableHandler.validate_variable_name(event, index):
            return
    state.variables[index][key_name] = event
    print(state.variables)

@ctrl.add("reset_variables")
def on_reset_variables():
    """
    Resets the variables list.
    """
    state.variables = [{"name": init_value, "value": init_value, "error_message": init_value}]
    state.dirty("variables")
    LatticeVariableHandler.update_delete_availability()

class LatticeVariableHandler:
    """
    Stores all functionality for dashboard variable referencing.
    """

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
    def get_duplicate_indexes(new_name, current_index):
        duplicates = [
            index
            for index, var in enumerate(state.variables)
            if var["name"] == new_name and index != current_index
        ]

        if duplicates:
            duplicates.append(current_index)
        return duplicates


    @staticmethod
    def validate_variable_name(new_name, index):

        def set_var_error_message(message):
            state.variables[index]["error_message"] = message
            state.dirty("variables")

        alpha = new_name and new_name[0].isalpha()
        duplicate_indexes = LatticeVariableHandler.get_duplicate_indexes(new_name, index)
        send_error = set_var_error_message("error")

        if not alpha: 
            send_error
            return True
        elif duplicate_indexes:
            for index in duplicate_indexes:
                send_error
            return True
        else:
           set_var_error_message("")

    @staticmethod
    def determine_if_variable(var_name):
        found_index = next((i for i, var in enumerate(state.variables) if var["name"] == var_name), None)
        return (found_index is not None, found_index)

    @staticmethod
    def variable_btn(**kwargs) -> vuetify.VBtn:
        """
        Creates a templated button.
        """
        return vuetify.VBtn(
            icon=True,
            small=True,
            elevation=2,
            **kwargs,
        )

    @staticmethod
    def variable_btn_icon(mdi_name: str) -> vuetify.VIcon:
        """
        Creates a templated icon for the button.
        """
        return vuetify.VIcon(mdi_name,small=True)

    @staticmethod
    def dialog_settings():
        dialog_name = "lattice_configuration_dialog"

        NavigationComponents.create_dialog_tabs(dialog_name, 2, ["Variables", "Defaults"])

        with vuetify.VTabsItems(v_model=(dialog_name, 0)):
            with vuetify.VTabItem():
                with vuetify.VCardText():
                    with vuetify.VContainer(fluid=True):
                        with vuetify.VRow(
                            v_for="(variable, index) in variables",
                            classes="align-center justify-center py-0",
                        ):
                            with vuetify.VCol(cols=5, classes="pr-0"):
                                vuetify.VTextField(
                                    placeholder="Name",
                                    outlined=True,
                                    dense=True,
                                    background_color="grey lighten-4",
                                    input=(ctrl.update_variable, "['name', index, $event]"),
                                    error_messages=("variable.error_message", []),
                                    hide_details=True,
                                )
                            with vuetify.VCol(cols=1, classes="px-0 text-center"):
                                html.Span("=", classes="mx-0")
                            with vuetify.VCol(cols=4, classes="pl-0"):
                                vuetify.VTextField(
                                    placeholder="Value",
                                    outlined=True,
                                    dense=True,
                                    type="number",
                                    background_color="grey lighten-4",
                                    change=(ctrl.update_variable, "['value', index, $event]"),
                                    hide_details=True,
                                )
                            with vuetify.VCol(cols=2, classes="d-flex"):
                                with html.Div(classes="mr-2"):
                                    with LatticeVariableHandler.variable_btn(
                                        color="primary",
                                        click=ctrl.add_variable,
                                        v_show="index === variables.length - 1"
                                    ):
                                        LatticeVariableHandler.variable_btn_icon("mdi-plus")
                                with html.Div():
                                    with LatticeVariableHandler.variable_btn(
                                        color="secondary",
                                        click=(ctrl.delete_variable, "[index]"),
                                        disabled=("is_only_variable",)
                                    ):
                                        LatticeVariableHandler.variable_btn_icon("mdi-delete")
                        with vuetify.VRow(classes="mt-2"):
                            with vuetify.VCol():
                                vuetify.VBtn(
                                    "Reset Variables",
                                    color="accent",
                                    click=ctrl.reset_variables,
                                    block=True,
                                )
            with vuetify.VTabItem():
                with vuetify.VCardText():
                    with vuetify.VRow():
                        with vuetify.VCol(cols=3):
                            InputComponents.text_field(
                                label="nslice",
                                v_model_name="nslice",
                                change=(ctrl.nsliceDefaultChange, "['nslice', $event]"),
                            )