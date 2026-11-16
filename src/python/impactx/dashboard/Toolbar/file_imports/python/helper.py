"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import ast
import functools
import re
import tokenize
from typing import Any

from impactx import RefPart
from impactx import distribution as impactx_distribution

from ....Input.defaults import DashboardDefaults
from ....Input.defaults_helper import InputDefaultsHelper
from .lattice_helper import (
    DEFINITIONS,
    TRACKING_MODES,
    _definition_time_nodes,
    _is_block,
    _tracking_functions,
    first_tracking_call,
    module_bindings,
    parse_script,
    top_level_statements,
    tracking_method,
    walk_outside_definitions,
    walk_with_definitions,
)

# calls that take the reference particle of the beam, with its keyword
REFERENCE_USES = {
    "init_envelope": "ref",
    "set_ref_particle": "refpart",
    "track_reference": None,
}

# reference particle setters, applied in the order of the script like in ImpactX
REFERENCE_SETTERS = ("set_species", "set_mass_MeV", "set_charge_qe")

# the setters of RefPart, which return the reference particle for chaining
REFPART_SETTERS = {name for name in dir(RefPart) if name.startswith("set_")}

# the reference particle of the beam, as found by _reference_origin, and of the
# envelope, which sim.init_envelope() sets
BEAM_REFERENCE = "beam"
ENVELOPE_REFERENCE = "envelope"

# expressions that evaluate some of their parts never, or more than once
COMPREHENSIONS = (ast.DictComp, ast.GeneratorExp, ast.ListComp, ast.SetComp)


class DashboardParserHelper:
    """
    Helper functions for building dashboard parser.
    """

    @staticmethod
    def import_file_content(file: str, state) -> dict:
        """
        Retrieves and prints the content of the uploaded simulation file.

        :param content: The content of the ImpactX simulation file.
        """
        if file:
            content = file["content"].decode("utf-8-sig")  # also with a BOM
            return content
        else:
            state.file_content = ""
            return ""

    @staticmethod
    def parse_single_inputs(content: str) -> dict:
        """
        Parses individual simulation parameters from the simulation file content.

        :param content: The content of the ImpactX simulation file.
        """
        reference_dictionary = DashboardDefaults.DEFAULT_VALUES.copy()
        tree = parse_script(content)
        unresolved_intensities = {}

        parsing_patterns = [
            r"\b{}\s*=\s*([^#\n]+)",  # (param = value)
            r"set_{}\(([^)]+)\)",  # (set_param(value))
        ]

        found = set()
        for parameter_name in reference_dictionary.keys():
            if parameter_name.endswith("_list"):
                continue
            if parameter_name in ("mass_MeV", "charge_qe"):
                continue  # see parse_reference_particle

            for pattern in parsing_patterns:
                pattern_match = re.search(pattern.format(parameter_name), content)
                if pattern_match:
                    text = pattern_match.group(1).strip()
                    line = content.count("\n", 0, pattern_match.start()) + 1
                    try:
                        value = DashboardParserHelper.literal_value(
                            tree, ast.parse(text, mode="eval").body, line
                        )
                    except (ValueError, SyntaxError) as error:
                        if parameter_name not in ("bunch_charge_C", "beam_current_A"):
                            raise
                        # only the intensity of the selected model must resolve
                        unresolved_intensities[parameter_name] = error
                        break

                    if parameter_name == "space_charge" and isinstance(value, bool):
                        # ImpactX takes True as the deprecated name of "3D"
                        value = "3D" if value else "false"

                    reference_dictionary[parameter_name] = value
                    found.add(parameter_name)
                    break

        reference_dictionary["tracking_mode"] = DashboardParserHelper._tracking_mode(
            tree, reference_dictionary["tracking_mode"]
        )
        # 2D space charge takes the beam current in A as the intensity
        intensity = (
            "beam_current_A"
            if reference_dictionary["space_charge"] == "2D"
            else "bunch_charge_C"
        )
        tracks_beam = reference_dictionary["tracking_mode"] != "Reference Tracking"
        if tracks_beam and intensity in unresolved_intensities:
            raise unresolved_intensities[intensity]
        assigned = found | set(unresolved_intensities)
        if (
            tracks_beam
            and intensity == "beam_current_A"
            and "beam_current_A" not in assigned
            and "bunch_charge_C" in assigned
        ):
            raise ValueError(
                "Cannot import the beam intensity: with 2D space charge, ImpactX "
                "takes a beam current, which the dashboard imports from "
                "beam_current_A, but the script sets only bunch_charge_C."
            )

        # the input field shows the kinetic energy that the simulation uses
        if "kin_energy_MeV" in found:
            reference_dictionary["kin_energy_on_ui"] = reference_dictionary[
                "kin_energy_MeV"
            ]

        reference_dictionary.update(
            DashboardParserHelper.parse_reference_particle(tree)
        )

        return reference_dictionary

    @staticmethod
    def _tracking_mode(tree: ast.Module, default: str) -> str:
        """
        Returns the dashboard tracking mode of the first tracking call.

        :param tree: The parsed ImpactX simulation file.
        :param default: The tracking mode without a tracking call.
        """
        method = tracking_method(tree)
        if method is None:
            return default
        return TRACKING_MODES[method]

    @staticmethod
    def parse_reference_particle(tree: ast.Module) -> dict:
        """
        Parses the mass and charge of the reference particle that the first
        tracking call uses.

        The script is followed like ImpactX up to that call. set_species() sets
        the mass and charge of a reference particle, set_mass_MeV() and
        set_charge_qe() one each, and a later call overrides an earlier one.
        The particles are the one of the beam, e.g., sim.beam.ref or
        pc.ref_particle(), and new ones, e.g., RefPart() or ref.copy().
        pc.set_ref_particle(p) copies p to the one of the beam, and
        sim.init_envelope(p, ...) to the one of the envelope.
        sim.track_reference(p) tracks p, sim.track_envelope() the one of the
        envelope, and sim.track_particles() the one of the beam.

        Code in function definitions runs only when they are called (@see
        check_script_functions). Calls on the tracked particle that may run
        never or more than once, e.g., in a loop, a condition or a conditional
        expression, raise a ValueError, as do calls on a particle the dashboard
        cannot tell. In a part of a script, a name it does not assign, e.g.,
        ref, is the one of the beam.

        :param tree: The parsed ImpactX simulation file.
        :return: Dictionary with the parsed 'mass_MeV' and 'charge_qe'.
        """
        tracking_call = first_tracking_call(tree, outside_definitions=True)
        context = {
            "parents": {
                child: node
                for node in ast.walk(tree)
                for child in ast.iter_child_nodes(node)
            },
            "bindings": module_bindings(tree),
            "top_level": {id(node) for node in top_level_statements(tree.body)},
            # names of the class, e.g., from impactx import RefPart as Particle
            "refpart_names": {"RefPart"}
            | {
                alias.asname
                for node in ast.walk(tree)
                if isinstance(node, ast.ImportFrom)
                for alias in node.names
                if alias.name == "RefPart" and alias.asname
            },
        }
        parents = context["parents"]
        reference_methods = set(REFERENCE_SETTERS) | set(REFERENCE_USES)

        def before_tracking(node: ast.AST) -> bool:
            # the arguments of the tracking call run before it
            return tracking_call is None or (node.end_lineno, node.end_col_offset) <= (
                tracking_call.end_lineno,
                tracking_call.end_col_offset,
            )

        def particle(node: ast.expr, at: ast.AST):
            """The key of the reference particle that an expression is."""
            origin = DashboardParserHelper._reference_origin(node, at, context)
            if origin is None:
                raise ValueError(
                    f"Cannot import the reference particle: the dashboard cannot "
                    f"tell which one '{ast.unparse(node)}' in line {node.lineno} "
                    f"is."
                )
            return origin if isinstance(origin, str) else id(origin)

        def conditional(node: ast.AST, statement: ast.stmt) -> bool:
            """Checks if a node may run never, or more than once."""
            if _is_block(statement):
                return True
            child = node
            while child is not statement:
                parent = parents[child]
                if (
                    (isinstance(parent, ast.IfExp) and child is not parent.test)
                    or (
                        isinstance(parent, ast.BoolOp) and child is not parent.values[0]
                    )
                    or (
                        isinstance(parent, COMPREHENSIONS)
                        and child is not parent.generators[0]
                    )
                    or (
                        isinstance(parent, ast.comprehension)
                        and not (
                            child is parent.iter
                            and parents[parent].generators[0] is parent
                        )
                    )
                    or isinstance(parent, ast.Assert)
                ):
                    return True
                child = parent
            return False

        def changed_particles(call: ast.Call) -> list:
            """
            The particle of ref.reset() or ref.load_file(...): the receiver, or
            the first argument of RefPart.reset(ref), also as ref=ref.
            """
            return (
                [call.func.value]
                + call.args[:1]
                + [item.value for item in call.keywords if item.arg == "ref"]
            )

        def particle_origin(node: ast.expr, at: ast.AST) -> bool:
            """Checks for a reference particle, also of a name the script binds."""
            if isinstance(node, ast.Name) and node.id not in context["bindings"]:
                return False  # e.g., a built-in function
            origin = DashboardParserHelper._reference_origin(node, at, context)
            return origin is not None

        statements = []
        for statement in top_level_statements(tree.body):
            if not before_tracking(statement):
                break
            if not isinstance(statement, DEFINITIONS):
                statements.append(statement)
                continue
            # parts of a definition that run where it is, e.g., a class body
            if isinstance(statement, ast.ClassDef):
                parts = [walk_with_definitions(statement, class_bodies=True)]
            else:
                parts = map(walk_with_definitions, _definition_time_nodes(statement))
            for part in parts:
                for node in part:
                    if (
                        isinstance(node, ast.Call)
                        and isinstance(node.func, ast.Attribute)
                        and node.func.attr in reference_methods
                    ):
                        raise ValueError(
                            f"Cannot import '{ast.unparse(node)}' in line "
                            f"{node.lineno}: it runs where '{statement.name}' is "
                            f"defined."
                        )

        # the calls on reference particles, in the order they run, also where a
        # class or function inside a block is defined
        events = []
        for statement in statements:
            for node in walk_with_definitions(statement, class_bodies=True):
                if not hasattr(node, "end_lineno") or not before_tracking(node):
                    continue
                # a method called later, or by its name, cannot be followed
                if (
                    isinstance(node, ast.Attribute)
                    and node.attr in reference_methods
                    and not (
                        isinstance(parents.get(node), ast.Call)
                        and parents[node].func is node
                    )
                ) or (
                    isinstance(node, ast.Constant) and node.value in reference_methods
                ):
                    raise ValueError(
                        f"Cannot import the reference particle: "
                        f"'{ast.unparse(node)}' in line {node.lineno} is not "
                        f"called directly."
                    )
                # e.g., ref.mass = 1.67e-27, in SI units, which the import does not follow
                if (
                    isinstance(node, ast.Attribute)
                    and node.attr in ("mass", "charge")
                    and not isinstance(node.ctx, ast.Load)
                    and DashboardParserHelper._reference_origin(
                        node.value, node, context
                    )
                    is not None
                ):
                    raise ValueError(
                        f"Cannot import the reference particle: '{ast.unparse(node)}' "
                        f"in line {node.lineno} sets an attribute; use set_species(), "
                        f"set_mass_MeV() or set_charge_qe()."
                    )
                if not (
                    isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                ):
                    continue
                method = node.func.attr
                if (
                    method in REFERENCE_SETTERS
                    or method in ("init_envelope", "set_ref_particle")
                    # e.g., ref.reset() or ref.load_file("beam.madx"), also as
                    # RefPart.load_file(ref, "beam.madx")
                    or (
                        method in ("reset", "load_file")
                        and any(
                            particle_origin(item, node)
                            for item in changed_particles(node)
                        )
                    )
                    or (
                        method == "copy"
                        and not node.args
                        and not node.keywords
                        and DashboardParserHelper._reference_origin(
                            node.func.value, node, context
                        )
                        is not None
                    )
                ):
                    events.append((node, statement))
        events.sort(key=lambda event: (event[0].end_lineno, event[0].end_col_offset))

        # the setter calls of each particle, or why it cannot be followed
        states = {}

        def state(key) -> dict:
            return states.setdefault(key, {"calls": [], "unknown": None})

        for node, statement in events:
            method = node.func.attr
            reason = None
            if conditional(node, statement):
                reason = (
                    f"'{ast.unparse(node)}' in line {node.lineno} is inside a "
                    f"loop, condition or other block, or a conditional expression"
                )
            if method in ("reset", "load_file"):
                for item in changed_particles(node):
                    if particle_origin(item, node):
                        target = state(particle(item, node))
                        target["unknown"] = target["unknown"] or (
                            f"'{ast.unparse(node)}' in line {node.lineno} changes "
                            f"it in a way the dashboard does not follow"
                        )
                continue
            if method in REFERENCE_SETTERS:
                target = state(particle(node.func.value, node))
                if reason:
                    target["unknown"] = target["unknown"] or reason
                else:
                    target["calls"].append(node)
                continue

            if method == "copy":
                source, destination = particle(node.func.value, node), id(node)
            else:
                keyword = REFERENCE_USES[method]
                arguments = node.args[:1] + [
                    item.value for item in node.keywords if item.arg == keyword
                ]
                if not arguments:
                    raise ValueError(f"Cannot import '{ast.unparse(node)}'.")
                source = particle(arguments[0], node)
                destination = (
                    ENVELOPE_REFERENCE if method == "init_envelope" else BEAM_REFERENCE
                )
            copied = state(source)
            if reason:
                state(destination)["unknown"] = reason
            else:
                states[destination] = {
                    "calls": list(copied["calls"]),
                    "unknown": copied["unknown"],
                }

        method = tracking_method(tree)
        tracked = BEAM_REFERENCE
        if method == "track_envelope" and ENVELOPE_REFERENCE in states:
            tracked = ENVELOPE_REFERENCE
        elif method == "track_reference" and tracking_call is not None:
            argument, at = DashboardParserHelper._tracked_reference(tree, tracking_call)
            tracked = particle(argument, at)

        tracked_state = state(tracked)
        if tracked_state["unknown"]:
            raise ValueError(
                f"Cannot import the reference particle: {tracked_state['unknown']}."
            )

        values = {}
        for call in tracked_state["calls"]:
            if len(call.args) != 1 or call.keywords:
                raise ValueError(f"Cannot import '{ast.unparse(call)}'.")
            value = DashboardParserHelper.literal_value(tree, call.args[0])

            if call.func.attr == "set_species":
                try:
                    reference = RefPart().set_species(value)
                except (RuntimeError, TypeError, ValueError) as error:
                    raise ValueError(
                        f"Cannot import '{ast.unparse(call)}' in line "
                        f"{call.lineno}: {error}"
                    ) from None
                values["mass_MeV"] = reference.mass_MeV
                values["charge_qe"] = reference.charge_qe
            else:
                values[call.func.attr.removeprefix("set_")] = value

        # the dashboard input for the charge is an integer
        charge_qe = values.get("charge_qe")
        if isinstance(charge_qe, float) and charge_qe.is_integer():
            values["charge_qe"] = int(charge_qe)

        return values

    @staticmethod
    def _tracked_reference(tree: ast.Module, tracking_call: ast.Call):
        """
        Returns the reference particle that the first tracking call, a
        track_reference(), tracks: the expression, and the node where it is
        evaluated. In a function of the script that starts with the tracking
        call, e.g., def run(p): sim.track_reference(p), a parameter is the
        argument of the call of the function, and other names are read when
        it is called.
        """
        where = f"in line {tracking_call.lineno}"
        if not isinstance(tracking_call.func, ast.Name):
            if len(tracking_call.args) != 1 or tracking_call.keywords:
                raise ValueError(
                    f"Cannot import '{ast.unparse(tracking_call)}' {where}."
                )
            return tracking_call.args[0], tracking_call

        name = tracking_call.func.id
        definition = next(
            node
            for node in top_level_statements(tree.body)
            if isinstance(node, ast.FunctionDef) and node.name == name
        )
        inner = _tracking_functions(tree)[name]
        if len(inner.args) != 1 or inner.keywords:
            raise ValueError(
                f"Cannot import '{ast.unparse(inner)}' in line {inner.lineno}."
            )
        argument = inner.args[0]
        cannot_follow = ValueError(
            f"Cannot import the reference particle: the dashboard cannot tell "
            f"which one '{ast.unparse(argument)}' in line {argument.lineno} is, "
            f"when '{name}' is called {where}."
        )

        arguments = definition.args
        positional = arguments.posonlyargs + arguments.args
        parameters = [item.arg for item in positional + arguments.kwonlyargs]
        if isinstance(argument, ast.Name) and argument.id in parameters:
            if any(isinstance(item, ast.Starred) for item in tracking_call.args) or any(
                item.arg is None for item in tracking_call.keywords
            ):
                raise cannot_follow  # e.g., run(*particles)
            index = parameters.index(argument.id)
            if index < len(positional) and index < len(tracking_call.args):
                return tracking_call.args[index], tracking_call
            for item in tracking_call.keywords:
                if item.arg == argument.id:
                    return item.value, tracking_call
            # the default value, evaluated where the function is defined
            if index < len(positional):
                offset = index - (len(positional) - len(arguments.defaults))
                if offset >= 0:
                    return arguments.defaults[offset], definition
            else:
                default = arguments.kw_defaults[index - len(positional)]
                if default is not None:
                    return default, definition
            raise cannot_follow

        # names of the function, which are no names of the module
        local_names = set(parameters) | {
            node.id
            for node in ast.walk(definition)
            if isinstance(node, ast.Name) and not isinstance(node.ctx, ast.Load)
        }
        for item in (arguments.vararg, arguments.kwarg):
            if item is not None:
                local_names.add(item.arg)
        if any(
            isinstance(node, ast.Name) and node.id in local_names
            for node in ast.walk(argument)
        ) or any(isinstance(node, ast.Call) for node in ast.walk(argument)):
            raise cannot_follow
        return argument, tracking_call

    @staticmethod
    def _reference_origin(node: ast.expr, at: ast.AST, context: dict, depth: int = 0):
        """
        Returns which reference particle an expression is, at the position of
        a node: BEAM_REFERENCE for the one of the beam, e.g., sim.beam.ref or
        pc.ref_particle(), the call that creates a new one, e.g., RefPart() or
        ref.copy(), or None if the dashboard cannot tell. A name counts if all
        its assignments are top-level statements, e.g., ref = sim.beam.ref.

        :param context: The parents of the nodes, the bindings of the
            module-level names, the ids of the top-level statements, and the
            names of RefPart in the script.
        """
        # the setters of RefPart return it, e.g., RefPart().set_species("proton")
        while (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr in REFPART_SETTERS
        ):
            node = node.func.value

        if (isinstance(node, ast.Attribute) and node.attr == "ref") or (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "ref_particle"
            and not node.args
            and not node.keywords
        ):
            return BEAM_REFERENCE
        if isinstance(node, ast.Call) and not node.args and not node.keywords:
            function = node.func
            if (
                isinstance(function, ast.Name)
                and function.id in context["refpart_names"]
            ) or (isinstance(function, ast.Attribute) and function.attr == "RefPart"):
                return node  # a new reference particle
            if isinstance(function, ast.Attribute) and function.attr == "copy":
                copied = DashboardParserHelper._reference_origin(
                    function.value, at, context, depth + 1
                )
                return node if copied is not None else None  # a new one
        if isinstance(node, ast.Name) and node.id not in context["bindings"]:
            # a part of a script, e.g., ref.set_species("proton") without ref = ...
            return BEAM_REFERENCE
        if isinstance(node, ast.Name) and depth < 100:
            assignment = DashboardParserHelper._assignment_at(node.id, at, context)
            if assignment is not None:
                return DashboardParserHelper._reference_origin(
                    assignment.value, assignment, context, depth + 1
                )
        return None

    @staticmethod
    def _assignment_at(name: str, at: ast.AST, context: dict):
        """
        Returns the assignment of a name in effect where a node starts, if all
        assignments of the name are top-level statements, else None.
        """
        assignments = []
        for node in context["bindings"].get(name, []):
            parent = context["parents"].get(node)
            simple = (
                isinstance(parent, ast.Assign)
                and any(target is node for target in parent.targets)
            ) or (
                isinstance(parent, ast.AnnAssign)
                and parent.target is node
                and parent.value is not None
            )
            if not simple or id(parent) not in context["top_level"]:
                return None
            assignments.append(parent)
        before = [
            assignment
            for assignment in assignments
            if (assignment.end_lineno, assignment.end_col_offset)
            <= (at.lineno, at.col_offset)
        ]
        if not before:
            return None
        return max(before, key=lambda node: (node.end_lineno, node.end_col_offset))

    @staticmethod
    def _simple_assignments(tree: ast.Module, name: str):
        """
        Returns the top-level assignments of a name, e.g., Ib = 0.5, if these
        are all its bindings in the script, else None: e.g., Ib *= 2, a tuple
        assignment, an assignment in a block, or a global in a function can
        change it in ways the dashboard does not follow.
        """
        simple = [
            statement
            for statement in top_level_statements(tree.body)
            if isinstance(statement, ast.Assign)
            and len(statement.targets) == 1
            and isinstance(statement.targets[0], ast.Name)
            and statement.targets[0].id == name
        ]
        simple_targets = {id(statement.targets[0]) for statement in simple}

        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)):
                if any(
                    isinstance(child, (ast.Global, ast.Nonlocal))
                    and name in child.names
                    for child in ast.walk(node)
                ):
                    return None
        for node in walk_outside_definitions(tree):
            if (
                isinstance(node, ast.Name)
                and node.id == name
                and not isinstance(node.ctx, ast.Load)
                and id(node) not in simple_targets
            ):
                return None
            if isinstance(node, ast.alias) and (node.asname or node.name) == name:
                return None
            if getattr(node, "name", None) == name and isinstance(
                node, (ast.ExceptHandler, ast.MatchAs, ast.MatchStar)
            ):
                return None
        return simple

    @staticmethod
    def literal_value(tree: ast.Module, node: ast.expr, line: int = None) -> Any:
        """
        Evaluates a literal, e.g., "proton" or -1, or a variable that is
        assigned a literal at the top level of the script before its use.

        :param tree: The parsed ImpactX simulation file.
        :param node: The expression to evaluate.
        :param line: The line of the expression, if node is not from tree.
        :return: The value of the expression.
        """
        try:
            return ast.literal_eval(node)
        except ValueError:
            pass

        line = line or node.lineno
        inside_definition = any(
            isinstance(
                definition, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
            )
            and definition.body[0].lineno <= line <= definition.end_lineno
            for definition in ast.walk(tree)
        )
        if isinstance(node, ast.Name) and not inside_definition:
            # a name in a function or class may be one of its own, e.g., a parameter
            assignments = DashboardParserHelper._simple_assignments(tree, node.id)
            earlier = [a for a in assignments or [] if a.lineno < line]
            if earlier:
                return DashboardParserHelper.literal_value(tree, earlier[-1].value)

        raise ValueError(
            f"Cannot import '{ast.unparse(node)}': the dashboard only evaluates "
            f"literal values and variables that are assigned one."
        )

    @staticmethod
    def parse_list_inputs(content: str) -> dict:
        """
        Parses list-based simulation parameters from the simulation file content.

        Some of the dashboard inputs are lists, such as n_cell = [16, 16, 16].
        We have to look inside of the brackets to retrieve the individual values.

        :param content: The content of the ImpactX simulation file.
        """
        dictionary = {}
        list_inputs = [
            "n_cell",
            "prob_relative",
            "blocking_factor_x",
            "blocking_factor_y",
            "blocking_factor_z",
        ]

        # The below regex will capture everything inside of the brackets
        # it assumes that the list will be placed in a single line (not multiple).
        # Example: n_cell = [16, 16, 16] -->  [16, 16, 16]
        list_parsing = r"{}\s*=\s*(\[[^\]]*\])"

        for input_name in list_inputs:
            match = re.search(list_parsing.format(input_name), content)
            if match:
                values = ast.literal_eval(match.group(1).strip())

                if input_name == "n_cell":
                    for i, dim in enumerate(["x", "y", "z"]):
                        dictionary[f"n_cell_{dim}"] = values[i]

                elif input_name == "prob_relative":
                    dictionary["prob_relative"] = values

                elif input_name.startswith("blocking_factor_"):
                    # Convert [16] -> 16
                    dictionary[input_name] = values[0]

        return dictionary

    @staticmethod
    def parse_distribution(content: str) -> dict:
        """
        Parses distribution section from the simulation file content.

        :param content: The content of the ImpactX simulation file.
        """

        dictionary = {"distribution": {"name": "", "type": "", "parameters": {}}}

        # e.g., a parenthesis in a comment does not end the arguments
        content = DashboardParserHelper._without_comments(content)
        distribution_name = re.search(r"distribution\.(\w+)\(", content)
        distribution_type_twiss = re.search(r"twiss\((.*?)\)", content, re.DOTALL)
        distribution_type_quadratic = re.search(
            r"distribution\.\w+\((.*?)\)", content, re.DOTALL
        )
        parameters = {}

        def extract_parameters(distribution_type):
            # the value runs up to a comma, parenthesis or line end,
            # EX: alpha_x=-0.5, emitt_x=1.0e-6
            parameter_pairs = re.findall(
                r"(\w+)\s*=\s*([^,\)\n]+)", distribution_type.group(1)
            )
            return {name: value.strip() for name, value in parameter_pairs}

        if distribution_name:
            dictionary["distribution"]["name"] = distribution_name.group(1)

        if distribution_type_twiss:
            dictionary["distribution"]["type"] = "Twiss"
            parameters = extract_parameters(distribution_type_twiss)
        elif distribution_type_quadratic:
            dictionary["distribution"]["type"] = "Quadratic"
            parameters = extract_parameters(distribution_type_quadratic)

        if distribution_name:
            positional = DashboardParserHelper._positional_arguments(
                content, distribution_name
            )
            for name in positional:
                if name in parameters:
                    raise ValueError(
                        f"Cannot import distribution.{distribution_name.group(1)}"
                        f"(...): multiple values for '{name}'."
                    )
            parameters = {**positional, **parameters}

        dictionary["distribution"]["parameters"] = parameters

        return dictionary

    @staticmethod
    def _without_comments(content: str) -> str:
        """
        Returns the content with spaces in place of its comments, so that all
        other code keeps its position.

        :param content: The content of the ImpactX simulation file.
        """
        # the line breaks of Python: a form feed does not end a line
        lines = re.findall(r"[^\r\n]*(?:\r\n|\r|\n)|[^\r\n]+$", content)
        tokens = tokenize.generate_tokens(functools.partial(next, iter(lines), ""))
        for token in tokens:
            if token.type == tokenize.COMMENT:
                (row, start), (_, end) = token.start, token.end
                line = lines[row - 1]
                lines[row - 1] = line[:start] + " " * (end - start) + line[end:]
        return "".join(lines)

    @staticmethod
    def _positional_arguments(content: str, match: re.Match) -> dict:
        """
        Parses the positional arguments of the distribution call at a regex
        match of ``distribution.<Name>(``, e.g.,
        ``distribution.SpinvMF(0.4, 0.9, 0.1)``, named after its signature.

        :param content: The content of the ImpactX simulation file.
        :param match: The match of the distribution call.
        :return: Dictionary of the parameter names and the source code of their values.
        """
        name = match.group(1)
        line = content.count("\n", 0, match.start()) + 1
        line_start = content.rfind("\n", 0, match.start()) + 1
        # end of "distribution", also for impactx.distribution.<Name>(
        end_column = len(content[line_start : match.start()].encode("utf-8")) + len(
            "distribution"
        )

        calls = [
            node
            for node in ast.walk(ast.parse(content))
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == name
            and (node.func.value.end_lineno, node.func.value.end_col_offset)
            == (line, end_column)
        ]
        if not calls or not calls[0].args:
            return {}  # e.g., a call in a comment
        call = calls[0]

        distribution_class = getattr(impactx_distribution, name, None)
        positional = (
            InputDefaultsHelper.positional_parameters(distribution_class)
            if distribution_class is not None
            else []
        )
        if len(call.args) > len(positional) or any(
            isinstance(argument, ast.Starred) for argument in call.args
        ):
            raise ValueError(
                f"Cannot import distribution.{name}(...): it takes at most "
                f"{len(positional)} positional arguments."
            )
        return {
            parameter_name: ast.get_source_segment(content, argument)
            or ast.unparse(argument)
            for parameter_name, argument in zip(positional, call.args)
        }

    @staticmethod
    def parse_variables(content: str, used_vars: set) -> dict:
        r"""
        Parses variable definitions from the simulation file content.

        The function looks for lines that define variables in similar format of:
            variable_name = value

        variables_regex breakdown:
            ^\s* - Remove leading whitespace before variable name
            (\w+) => Retrieve the variable name
            \s*=\s*  => Allow the line to have spaces around '=' (ie. 'a= 10' or 'a =10')
            ([^#\n]+) => Captures everything after '=' up to a comment '#' or end of line

        EX when calling this function:
            content = '''
                a = 10
                b = 2.5  # some comment
                c = "some string"
            '''
            used_vars = {"a", "b"}
            result = parse_variables(content, used_vars)
            # result: { "a": 10, "b": 2.5}

        :param content: The content of the ImpactX simulation file.
        :param used_vars: Set of variable names to filter by.
        :return: Dictionary of parsed variables.
        """

        variables = {}
        variables_regex = r"^\s*(\w+)\s*=\s*([^#\n]+)"
        variables_found = re.findall(variables_regex, content, re.MULTILINE)

        for var_name, var_value in variables_found:
            if var_name not in used_vars:
                continue
            try:
                value = ast.literal_eval(var_value.strip())
                # only allowing numbers at the moment
                if isinstance(value, (int, float)):
                    variables[var_name] = value
            except Exception:
                continue
        return variables
