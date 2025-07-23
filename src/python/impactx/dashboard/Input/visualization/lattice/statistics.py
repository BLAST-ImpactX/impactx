"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .... import html, state, vuetify
from ...utils import GeneralFunctions

state.total_elements = 0
state.total_length = 0
state.max_length = 0
state.min_length = 0
state.avg_length = 0
state.total_steps = 0
state.element_counts = {}
state.length_stats_content = ""
state.lattice_is_empty = len(state.selected_lattice_list) == 0


class LatticeVisualizerStatisticUtils:
    @staticmethod
    def _extract_parameter_values(parameter_name: str, value_type=float):
        """
        Helper function to extract parameter values from the lattice list.

        :param parameter_name: Name of the parameter to extract (case-insensitive)
        :param value_type: Type to convert values to (float, int, etc.)
        :return: List of extracted values
        """
        values = []

        for element in state.selected_lattice_list:
            for param in element.get("parameters", []):
                if param.get("parameter_name", "").lower() == parameter_name.lower():
                    try:
                        values.append(value_type(param.get("sim_input", 0)))
                    except (ValueError, TypeError):
                        pass

        return values

    @staticmethod
    def update_length_statistics() -> None:
        """
        Computes and return the total, min, max, and average length of the
        lattice configuration. Sums all elements' 'ds' (length) parameters.
        """
        lengths = LatticeVisualizerStatisticUtils._extract_parameter_values("ds", float)
        if lengths:
            state.total_length = f"{round(sum(lengths), 2)}m"
            state.min_length = f"{round(min(lengths), 3)}m"
            state.max_length = f"{round(max(lengths), 3)}m"
            state.avg_length = f"{round(sum(lengths) / len(lengths), 3)}m"
            state.length_stats_content = [
                f"Longest: {state.max_length}",
                f"Shortest: {state.min_length}",
                f"Average: {state.avg_length}",
            ]
        else:
            state.total_length = None
            state.min_length = None
            state.max_length = None
            state.avg_length = None
            state.length_stats_content = []

    @staticmethod
    def update_element_counts() -> dict[str, int]:
        """
        Computes the element counts in the lattice list.

        :return: Dictionary of element names and their counts, sorted by count descending.
        """
        counts = {}
        for element in state.selected_lattice_list:
            key = str(element["name"]).lower()
            # can't do += 1 because key is not already initialized
            counts[key] = counts.get(key, 0) + 1

        state.lattice_is_empty = len(counts) == 0
        # sort from desc. so we see top elements left to right
        sorted_counts = dict(
            sorted(counts.items(), key=lambda item: item[1], reverse=True)
        )

        state.element_counts = sorted_counts
        return sorted_counts

    @staticmethod
    def update_total_steps() -> int:
        """
        Computes the total number of steps by summing 'nslice'
        across all lattice elements.

        :return: Total number of slices.
        """
        steps = LatticeVisualizerStatisticUtils._extract_parameter_values("nslice", int)
        return sum(steps)


class LatticeVisualizerStatisticComponents:
    @staticmethod
    def _stat(title: str) -> None:
        """
        Displays a statistic block for the statistics section
        in the lattice visualizer.

        :param title: The statistic name
        """
        title_state_name = GeneralFunctions.normalize_for_v_model(title)
        is_stat_length = "length" in title.lower()

        vuetify.VCardSubtitle(title, classes="pb-0 mb-0")

        with vuetify.VCardTitle(
            f"{{{{ {title_state_name} || '-' }}}}",
            classes="d-flex align-center justify-center my-0 py-0",
        ):
            if is_stat_length:
                LatticeVisualizerStatisticComponents._additional_length_stats()

    @staticmethod
    def _additional_length_stats():
        with vuetify.VTooltip(
            location="bottom",
        ):
            with vuetify.Template(v_slot_activator="{ props }"):
                vuetify.VIcon(
                    "mdi-information",
                    size="x-small",
                    v_bind="props",
                    disabled=("lattice_is_empty",),
                    classes="ml-2",
                )
            with vuetify.Template(v_for="line in length_stats_content"):
                html.Div("{{ line }}")

    @staticmethod
    def statistics():
        with vuetify.VCardText():
            # row 1: numerical stats
            with vuetify.VRow(classes="text-center"):
                with vuetify.VCol():
                    LatticeVisualizerStatisticComponents._stat("Total Elements")
                with vuetify.VCol():
                    LatticeVisualizerStatisticComponents._stat("Total Length")
                with vuetify.VCol():
                    LatticeVisualizerStatisticComponents._stat("Total Steps")
                with vuetify.VCol():
                    LatticeVisualizerStatisticComponents._stat("Periods")

            # row 2: element breakdown
            with vuetify.VRow(classes="pt-0 mt-0"):
                with vuetify.VCol(cols=12):
                    vuetify.VCardSubtitle("Element Breakdown")
                    with vuetify.Template(v_if="lattice_is_empty"):
                        vuetify.VCardTitle("Lattice list is empty.")
                    with vuetify.Template(v_else=True):
                        with vuetify.VChipGroup():
                            with vuetify.Template(
                                v_for="(count, name) in element_counts", key="name"
                            ):
                                vuetify.VChip(
                                    "{{ name.charAt(0).toUpperCase() + name.slice(1) }}: {{ count }}",
                                    style="font-size: 0.75rem;",
                                )
