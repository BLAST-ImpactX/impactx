import ast
import inspect
import re
from typing import Any, Dict, List, Optional

from impactx import elements as impactx_elements

from ....Input.defaults_helper import InputDefaultsHelper

# methods of sim.lattice: followed by the dashboard, without effect on the
# lattice, and everything else, which the dashboard cannot follow
LATTICE_SUPPORTED_METHODS = {"append", "clear", "extend", "pop_back"}
LATTICE_READ_METHODS = {
    "count",
    "count_by_kind",
    "generation",
    "get_kinds",
    "has_kind",
    "index",
    "is_empty",
    "isclose",
    "map_trace",
    "plot_survey",
    "select",
    "size",
    "to_dicts",
    "to_py",
    "transfer_map",
}

# methods without effect on a selection, sim.lattice.select(...), or an element
SELECTION_READ_METHODS = {"count_by_kind", "get_kinds", "has_kind", "isclose", "select"}
ELEMENT_READ_METHODS = {"copy", "finalize", "has_name", "isclose", "push", "to_dict"}
ELEMENT_READ_METHODS |= {"rc", "transfer_map"}  # e.g., Sbend.rc(ref)
ELEMENT_READ_METHODS |= {  # e.g., quad.__repr__() or quad.__eq__
    "__dir__",
    "__eq__",
    "__format__",
    "__ge__",
    "__getstate__",
    "__gt__",
    "__hash__",
    "__le__",
    "__lt__",
    "__ne__",
    "__reduce__",
    "__reduce_ex__",
    "__repr__",
    "__sizeof__",
    "__str__",
}

# list methods without effect on the list
LIST_READ_METHODS = {"copy", "count", "index"}

# built-in functions that neither change their arguments nor return a container of them
SAFE_FUNCTIONS = {
    "abs",
    "all",
    "any",
    "bool",
    "float",
    "format",
    "id",
    "int",
    "isinstance",
    "len",
    "print",
    "range",
    "repr",
    "round",
    "str",
    "type",
}

# built-in functions that run code from text or reach any variable by its name
DYNAMIC_FUNCTIONS = {"__import__", "eval", "exec", "globals", "locals", "vars"}

# attributes that reach any variable by its name, e.g., sys.modules
DYNAMIC_ATTRIBUTES = {"__dict__", "__globals__", "f_globals", "f_locals", "modules"}

# functions that import a module, e.g., the script itself by __name__
IMPORT_FUNCTIONS = {"__import__", "import_module", "reload", "run_module", "run_path"}

# decorators that do not call the function where it is defined
DEFINING_DECORATORS = {
    "abstractmethod",
    "cached_property",
    "classmethod",
    "deleter",
    "getter",
    "property",
    "setter",
    "staticmethod",
}

# built-in functions whose results run no code of the script in operators
INERT_FUNCTIONS = {
    "all",
    "any",
    "bool",
    "float",
    "format",
    "id",
    "int",
    "isinstance",
    "len",
    "print",
    "range",
    "repr",
    "str",
}

# values of expressions whose operators run no code of the script, e.g., cell + [q]
KNOWN_OPERANDS = (bool, bytes, complex, float, int, list, str, tuple, type(None))

# upper limits for the elements of a repetition, e.g., cell * n, which the
# dashboard shows one by one, and of a lattice, which is kept in memory
MAX_REPEATED_ELEMENTS = 1000
MAX_LATTICE_ELEMENTS = 100_000

# simulation methods that track the beam, and the dashboard tracking mode of each
TRACKING_MODES = {
    "track_particles": "Particle Tracking",
    "track_envelope": "Envelope Tracking",
    "track_reference": "Reference Tracking",
    "evolve": "Particle Tracking",  # deprecated
}

DEFINITIONS = (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Lambda)


def _is_block(statement: ast.stmt) -> bool:
    """
    Checks for compound statements, e.g., loops, conditions, function definitions.
    """
    return hasattr(statement, "body") or hasattr(statement, "cases")


def _is_main_guard(statement: ast.stmt) -> bool:
    """
    Checks for ``if __name__ == "__main__":``.
    """
    if not isinstance(statement, ast.If):
        return False
    test = statement.test
    return (
        isinstance(test, ast.Compare)
        and isinstance(test.left, ast.Name)
        and test.left.id == "__name__"
        and len(test.ops) == 1
        and isinstance(test.ops[0], ast.Eq)
        and isinstance(test.comparators[0], ast.Constant)
        and test.comparators[0].value == "__main__"
    )


def top_level_statements(body: List[ast.stmt]):
    """
    Yields the top-level statements of a script, including those in the body
    of ``if __name__ == "__main__":``.
    """
    for statement in body:
        if _is_main_guard(statement):
            yield from top_level_statements(statement.body)
        else:
            yield statement


def walk_outside_definitions(node: ast.AST):
    """
    Like ast.walk, but does not enter function and class definitions, whose
    code does not run where it is defined.
    """
    yield node
    for child in ast.iter_child_nodes(node):
        if not isinstance(child, DEFINITIONS):
            yield from walk_outside_definitions(child)


def _is_tracking_call(node: ast.AST) -> bool:
    """
    Checks for a call that tracks the beam, e.g., sim.track_particles().
    """
    return (
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr in TRACKING_MODES
    )


def _position(node: ast.AST):
    return (node.lineno, node.col_offset)


def parse_script(content: str) -> ast.Module:
    """
    Parses an ImpactX simulation file.

    A script that runs in a function, e.g., def main(): ..., called by its
    last statement, main(), also inside if __name__ == "__main__":, runs the
    body of the function there, once. The body is moved there, to the top
    level, if the function has no parameters, decorators, return, yield or
    global statements, and is defined right before the call and used
    nowhere else. A first tracking call in a loop, e.g., for run in
    range(3):, is followed in the first iteration (@see
    _inline_first_iteration).

    :param content: The content of the ImpactX simulation file.
    :return: The parsed file.
    """
    tree = ast.parse(content)
    _inline_main_function(tree)
    _inline_first_iteration(tree)
    return tree


def _inline_main_function(tree: ast.Module) -> ast.Module:
    """
    Moves the body of a main function to its call, @see parse_script.
    """
    body = tree.body
    if len(body) < 2:
        return tree
    container, call = body, body[-1]
    if _is_main_guard(call) and not call.orelse and len(call.body) == 1:
        container, call = call.body, call.body[0]
    definition = body[-2]
    if not (
        isinstance(call, ast.Expr)
        and isinstance(call.value, ast.Call)
        and isinstance(call.value.func, ast.Name)
        and not call.value.args
        and not call.value.keywords
        and isinstance(definition, ast.FunctionDef)
        and definition.name == call.value.func.id
        and not definition.decorator_list
        and definition.returns is None
    ):
        return tree
    arguments = definition.args
    if (
        arguments.posonlyargs
        or arguments.args
        or arguments.kwonlyargs
        or arguments.vararg
        or arguments.kwarg
    ):
        return tree
    # the body runs once, from its start to its end
    if any(
        isinstance(
            node, (ast.Return, ast.Yield, ast.YieldFrom, ast.Global, ast.Nonlocal)
        )
        for node in walk_outside_definitions(definition)
    ):
        return tree
    name = definition.name
    if module_bindings(tree).get(name) != [definition] or [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Name) and node.id == name
    ] != [call.value.func]:
        return tree
    body.remove(definition)
    container[-1:] = definition.body
    return tree


def _first_item(iterable: ast.expr, tree: ast.Module) -> Optional[ast.expr]:
    """
    Returns the first item of a range of ints, e.g., range(3) or range(runs)
    with runs = 3, or of a list or tuple literal, as an expression, or None.
    A name in the range must be assigned an int once, at the top level, before
    the loop.
    """
    if isinstance(iterable, (ast.List, ast.Tuple)):
        # all items are evaluated before the loop: none may run code
        if iterable.elts and not any(
            isinstance(node, (ast.Call, ast.Lambda, ast.NamedExpr, ast.Starred))
            for node in ast.walk(iterable)
        ):
            return iterable.elts[0]
        return None
    if not (
        isinstance(iterable, ast.Call)
        and isinstance(iterable.func, ast.Name)
        and iterable.func.id == "range"
        and not iterable.keywords
        and 1 <= len(iterable.args) <= 3
    ):
        return None
    bindings = module_bindings(tree)
    top_level = {id(statement) for statement in top_level_statements(tree.body)}

    def int_value(argument: ast.expr) -> Any:
        if isinstance(argument, ast.Name):
            nodes = bindings.get(argument.id, [])
            assignment = next(
                (
                    statement
                    for statement in top_level_statements(tree.body)
                    if isinstance(statement, ast.Assign)
                    and len(statement.targets) == 1
                    and nodes == [statement.targets[0]]
                ),
                None,
            )
            if (
                assignment is None
                or id(assignment) not in top_level
                or _position(assignment) >= _position(iterable)
            ):
                raise ValueError(f"'{argument.id}' is not an int literal")
            argument = assignment.value
        return ast.literal_eval(argument)

    try:
        bounds = [int_value(argument) for argument in iterable.args]
        values = range(*bounds)
    except (ValueError, TypeError, SyntaxError):
        return None
    if not values:
        return None
    return ast.copy_location(ast.Constant(value=values[0]), iterable)


def _inline_first_iteration(tree: ast.Module) -> None:
    """
    A first tracking call that is a statement of a loop over a range or a list
    literal, e.g., for run in range(3): ...; sim.track_particles(), runs in the
    first iteration: the loop is replaced by an assignment of the first item
    and its body, to the top level, if nothing before the call in the body
    can end the iteration, e.g., break or continue.
    """
    tracking_call = first_tracking_call(tree, outside_definitions=True)
    if tracking_call is None:
        return
    containers = [tree.body] + [
        statement.body for statement in tree.body if _is_main_guard(statement)
    ]
    for container in containers:
        for index, loop in enumerate(container):
            if not isinstance(loop, ast.For):
                continue
            position = next(
                (
                    i
                    for i, statement in enumerate(loop.body)
                    if isinstance(statement, (ast.Expr, ast.Assign))
                    and statement.value is tracking_call
                ),
                None,
            )
            if position is None or not isinstance(loop.target, ast.Name):
                continue
            if any(
                isinstance(node, (ast.Break, ast.Continue))
                for statement in loop.body[:position]
                for node in walk_outside_definitions(statement)
            ):
                return
            first = _first_item(loop.iter, tree)
            if first is None:
                return
            target = ast.copy_location(
                ast.Name(id=loop.target.id, ctx=ast.Store()), loop.target
            )
            assignment = ast.Assign(targets=[target], value=first)
            assignment.lineno, assignment.col_offset = target.lineno, target.col_offset
            assignment.end_lineno = loop.iter.end_lineno
            assignment.end_col_offset = loop.iter.end_col_offset
            container[index : index + 1] = [assignment] + loop.body
            return


def _definition_time_nodes(definition: ast.AST) -> List[ast.AST]:
    """
    Returns the parts of a function, class or lambda definition that run where
    it is defined: decorators, default values, annotations and base classes.
    """
    nodes = list(getattr(definition, "decorator_list", []))
    nodes += getattr(definition, "bases", []) + getattr(definition, "keywords", [])
    arguments = getattr(definition, "args", None)
    if isinstance(arguments, ast.arguments):
        nodes += arguments.defaults + [d for d in arguments.kw_defaults if d]
        for argument in (
            arguments.posonlyargs
            + arguments.args
            + arguments.kwonlyargs
            + [arguments.vararg, arguments.kwarg]
        ):
            if argument is not None and argument.annotation is not None:
                nodes.append(argument.annotation)
    if getattr(definition, "returns", None) is not None:
        nodes.append(definition.returns)
    return nodes


def walk_with_definitions(node: ast.AST, class_bodies: bool = False):
    """
    Like walk_outside_definitions, but also yields the function, class and
    lambda definitions, and the parts of them that run where they are defined,
    e.g., decorators and default values, without their bodies. With
    class_bodies, also the bodies of classes, which run there, too.
    """
    yield node
    yield from _walk_children(ast.iter_child_nodes(node), class_bodies)


def _walk_children(children, class_bodies: bool):
    """
    Walks nodes for walk_with_definitions.
    """
    for child in children:
        if isinstance(child, DEFINITIONS):
            yield child
            for part in _definition_time_nodes(child):
                yield from walk_with_definitions(part, class_bodies)
            if class_bodies and isinstance(child, ast.ClassDef):
                yield from _walk_children(child.body, class_bodies)
        else:
            yield from walk_with_definitions(child, class_bodies)


def _free_names(definition: ast.AST) -> set:
    """
    Returns the names a function, class or lambda reads from the module, or
    may read: names it loads but does not bind, and names it declares global.
    Its definition-time parts, e.g., default values, run outside of it.
    """
    loads, locals_, declared = set(), set(), set()
    nested = []
    definition_time = {id(node) for node in _definition_time_nodes(definition)}

    def visit(node: ast.AST) -> None:
        for child in ast.iter_child_nodes(node):
            if id(child) in definition_time:
                continue
            if isinstance(child, DEFINITIONS):
                nested.append(child)
                if hasattr(child, "name"):
                    locals_.add(child.name)
                for part in _definition_time_nodes(child):
                    visit_part(part)
                continue
            if isinstance(child, ast.comprehension):
                # its target is local to the comprehension: count only loads
                visit_part(child.iter)
                for condition in child.ifs:
                    visit_part(condition)
                continue
            if isinstance(child, ast.Name):
                (loads if isinstance(child.ctx, ast.Load) else locals_).add(child.id)
            elif isinstance(child, ast.arg):
                locals_.add(child.arg)
            elif isinstance(child, ast.alias):
                locals_.add(child.asname or child.name.split(".")[0])
            elif isinstance(child, (ast.ExceptHandler, ast.MatchAs, ast.MatchStar)):
                if child.name:
                    locals_.add(child.name)
            elif isinstance(child, ast.MatchMapping) and child.rest:
                locals_.add(child.rest)
            elif isinstance(child, (ast.Global, ast.Nonlocal)):
                declared.update(child.names)
            visit(child)

    def visit_part(node: ast.AST) -> None:
        if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load):
            loads.add(node.id)
        visit(node)

    visit(definition)
    if isinstance(definition, ast.ClassDef):
        # names of the class body are not visible in its methods
        free = loads | set().union(*(map(_free_names, nested)))
    else:
        own = locals_ - declared
        free = (loads - own) | (set().union(*map(_free_names, nested)) - own)
    return free | declared


def module_bindings(tree: ast.Module) -> Dict[str, List[ast.AST]]:
    """
    Returns the nodes that bind each module-level name: assignments, imports,
    definitions and global declarations.
    """
    bindings: Dict[str, List[ast.AST]] = {}
    for node in walk_with_definitions(tree):
        name = None
        if isinstance(node, ast.Name) and not isinstance(node.ctx, ast.Load):
            name = node.id
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            name = node.name
        elif isinstance(node, (ast.ExceptHandler, ast.MatchAs, ast.MatchStar)):
            name = node.name
        elif isinstance(node, ast.MatchMapping):
            name = node.rest  # e.g., case {**rest}
        if name:
            bindings.setdefault(name, []).append(node)
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            for alias in node.names:
                name = alias.asname or alias.name.split(".")[0]
                bindings.setdefault(name, []).append(alias)
    for node in ast.walk(tree):
        if isinstance(node, ast.Global):
            for name in node.names:
                bindings.setdefault(name, []).append(node)
    return bindings


def _tracking_functions(tree: ast.Module) -> Dict[str, ast.Call]:
    """
    Returns the top-level functions of the script that track the beam, with the
    first tracking call in each. A name that is also bound otherwise is left out.
    """
    bindings = module_bindings(tree)
    functions = {}
    for node in top_level_statements(tree.body):
        if isinstance(node, ast.FunctionDef):  # calling an async one does not run it
            calls = [child for child in ast.walk(node) if _is_tracking_call(child)]
            if calls and bindings.get(node.name) == [node]:
                functions[node.name] = min(calls, key=_position)
    return functions


# attributes in the code of a function that make it relevant for the import
RELEVANT_ATTRIBUTES = {
    "beam",
    "init_envelope",
    "lattice",
    "ref",
    "ref_particle",
    "set_charge_qe",
    "set_kin_energy_MeV",
    "set_mass_MeV",
    "set_ref_particle",
    "set_species",
} | set(TRACKING_MODES)


def check_script_functions(tree: ast.Module, source) -> None:
    """
    Raises if the script uses a function, class or lambda of its own before the
    first tracking call, and its code uses module-level variables, the lattice,
    the reference particle or tracking: the dashboard does not follow where
    and when such code runs. A direct call of a top-level function that is the
    first tracking call is allowed. A hook, e.g., sim.hook["before_period"] = f,
    raises: its code runs while the beam is tracked and may change the lattice.

    :param tree: The parsed ImpactX simulation file.
    :param source: Function that returns the source code of a node.
    """
    tracking_call = first_tracking_call(tree, outside_definitions=True)
    end = _position(tracking_call) if tracking_call is not None else None
    bindings = module_bindings(tree)
    variables = {
        name
        for name, nodes in bindings.items()
        if any(
            not isinstance(
                node, (ast.alias, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
            )
            for node in nodes
        )
    }

    definitions = [node for node in ast.walk(tree) if isinstance(node, DEFINITIONS)]
    free_names = {id(definition): _free_names(definition) for definition in definitions}

    def dynamic_access(node: ast.AST) -> bool:
        """
        Checks for code that runs text or reaches any variable by its name,
        e.g., exec(...), globals() or sys.modules.
        """
        return (
            (
                isinstance(node, ast.Name)
                and node.id in DYNAMIC_FUNCTIONS | {"__builtins__"}
                and node.id not in bindings
            )
            or (isinstance(node, ast.Attribute) and node.attr in DYNAMIC_ATTRIBUTES)
            or (
                isinstance(node, ast.Attribute)
                and node.attr in DYNAMIC_FUNCTIONS
                and isinstance(node.value, ast.Name)
                and node.value.id in ("builtins", "__builtins__")
            )
            # the module of the script itself, e.g., from __main__ import cell, or
            # of the built-in functions, e.g., import builtins
            or (isinstance(node, ast.alias) and node.name in ("__main__", "builtins"))
            or (
                isinstance(node, ast.ImportFrom)
                and node.module in ("__main__", "builtins")
            )
            or (
                isinstance(node, ast.Call)
                and (
                    getattr(node.func, "id", None) in IMPORT_FUNCTIONS
                    or getattr(node.func, "attr", None) in IMPORT_FUNCTIONS
                )
                and any(
                    (
                        isinstance(argument, ast.Constant)
                        and argument.value == "__main__"
                    )
                    or (isinstance(argument, ast.Name) and argument.id == "__name__")
                    for argument in node.args
                )
            )
        )

    relevant = set()
    for definition in definitions:
        if any(
            (
                isinstance(node, ast.Attribute)
                and (node.attr in RELEVANT_ATTRIBUTES or "lattice" in node.attr)
            )
            or dynamic_access(node)
            for node in ast.walk(definition)
        ) or (free_names[id(definition)] & variables):
            relevant.add(definition)
    names = {getattr(d, "name", None) for d in relevant} - {None}
    changed = True
    while changed:
        # a definition that uses a relevant function, or contains a relevant
        # definition, e.g., a method, is relevant, too
        changed = False
        for definition in definitions:
            if definition not in relevant and (
                free_names[id(definition)] & names
                or any(
                    node in relevant
                    for node in ast.walk(definition)
                    if node is not definition
                )
            ):
                relevant.add(definition)
                names.add(getattr(definition, "name", None))
                changed = True

    outside_definitions = {id(node) for node in walk_outside_definitions(tree)}

    def defines_only(definition: ast.AST) -> bool:
        """Checks for decorators that do not call it, e.g., @staticmethod."""
        return all(
            getattr(decorator, "id", None) in DEFINING_DECORATORS
            or getattr(decorator, "attr", None) in DEFINING_DECORATORS
            for decorator in getattr(definition, "decorator_list", [])
        )

    # the receiver and arguments of the tracking call run before it
    inside_tracking = set()
    if tracking_call is not None:
        inside_tracking = {id(node) for node in ast.walk(tracking_call)} - {
            id(tracking_call),
            id(tracking_call.func),
        }

    def before_end(node: ast.AST) -> bool:
        return hasattr(node, "lineno") and (
            end is None or _position(node) < end or id(node) in inside_tracking
        )

    # e.g., sim.hook["before_period"] = f, whose code runs while the beam is tracked
    for node in walk_outside_definitions(tree):
        if isinstance(node, ast.Attribute) and node.attr == "hook" and before_end(node):
            raise ValueError(
                f"Cannot import the lattice: '{source(node)}' in line "
                f"{node.lineno} uses a hook, whose code runs while the beam is "
                f"tracked, which the dashboard does not import."
            )

    for node in walk_with_definitions(tree, class_bodies=True):
        if not before_end(node):
            continue
        used = None
        if node in relevant and (
            isinstance(node, ast.Lambda) or not defines_only(node)
        ):
            used = node  # a lambda, or a decorated definition, may be stored
        elif _is_tracking_call(node) and id(node) not in outside_definitions:
            raise ValueError(
                f"Cannot import the lattice: '{source(node)}' in line "
                f"{node.lineno} runs where a class or function is defined, "
                f"which the dashboard does not follow."
            )
        elif dynamic_access(node):
            raise ValueError(
                f"Cannot import the lattice: '{source(node)}' in line "
                f"{node.lineno} may run code or change variables the dashboard "
                f"does not follow."
            )
        elif (
            isinstance(node, ast.Name)
            and isinstance(node.ctx, ast.Load)
            and node.id in names
            and not (tracking_call is not None and node is tracking_call.func)
        ):
            used = node
        if used is not None:
            raise ValueError(
                f"Cannot import the lattice: '{source(used)}' in line "
                f"{used.lineno} uses a function of the script that uses the "
                f"simulation or its variables, which the dashboard does not follow."
            )

    # a function of the script that tracks must do so first: the dashboard does
    # not follow code that runs in it before
    if tracking_call is not None and isinstance(tracking_call.func, ast.Name):
        definition = next(
            node
            for node in top_level_statements(tree.body)
            if isinstance(node, ast.FunctionDef) and node.name == tracking_call.func.id
        )
        body = definition.body
        if (
            isinstance(body[0], ast.Expr)
            and isinstance(body[0].value, ast.Constant)
            and isinstance(body[0].value.value, str)
        ):
            body = body[1:]  # the docstring
        first = body[0] if body else None
        if not (
            isinstance(first, (ast.Expr, ast.Assign))
            and _is_tracking_call(first.value)
            # its receiver and arguments run no code, e.g., sim.track_particles()
            and not any(
                isinstance(
                    node,
                    (
                        ast.Call,
                        ast.Lambda,
                        ast.NamedExpr,
                        ast.Starred,
                        ast.DictComp,
                        ast.GeneratorExp,
                        ast.ListComp,
                        ast.SetComp,
                    ),
                )
                for node in ast.walk(first.value)
                if node is not first.value
            )
        ):
            raise ValueError(
                f"Cannot import the lattice: the function '{definition.name}' "
                f"called in line {tracking_call.lineno} runs other code before "
                f"it tracks the beam, which the dashboard does not follow."
            )

    # the first simulation must start with a statement of its own
    if tracking_call is not None:
        for statement in top_level_statements(tree.body):
            if any(node is tracking_call for node in ast.walk(statement)):
                plain = (
                    isinstance(statement, (ast.Expr, ast.Assign))
                    and statement.value is tracking_call
                )
                if not plain:
                    raise ValueError(
                        f"Cannot import the lattice: the tracking call in line "
                        f"{tracking_call.lineno} must be a statement of its own, "
                        f"outside of loops, conditions and other blocks."
                    )
                break


def first_tracking_call(
    tree: ast.Module, outside_definitions: bool = False
) -> Optional[ast.Call]:
    """
    Returns the call that starts the first simulation, in the order of the
    source code: a tracking call, e.g., sim.track_particles(), or a call of a
    function of the script that tracks. The dashboard imports this simulation.

    Calls outside of function and class definitions take precedence, because
    they run where they are. With outside_definitions, only those count.
    """
    functions = _tracking_functions(tree)

    def starts_tracking(node: ast.AST) -> bool:
        return _is_tracking_call(node) or (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id in functions
        )

    walks = [walk_outside_definitions(tree)]
    if not outside_definitions:
        walks.append(ast.walk(tree))

    for nodes in walks:
        calls = [node for node in nodes if starts_tracking(node)]
        if calls:
            return min(calls, key=_position)
    return None


def tracking_method(tree: ast.Module) -> Optional[str]:
    """
    Returns the method of the first tracking call, e.g., "track_particles".
    """
    call = first_tracking_call(tree)
    if call is not None and not _is_tracking_call(call):
        call = _tracking_functions(tree)[call.func.id]
    return None if call is None else call.func.attr


class _Element:
    """
    An element constructor call, e.g., elements.Drift(ds=1.0).
    """

    def __init__(self, call: ast.Call):
        self.call = call
        self.unresolved: Optional[str] = None  # reason, if it may have changed


class _Unresolved:
    """
    A value that the dashboard cannot evaluate, with the reason.
    """

    def __init__(self, reason: str):
        self.reason = reason


class _Inert(_Unresolved):
    """
    A value whose operators run no code of the script, e.g., a str constant
    or the result of len().
    """


class _Lattice:
    """
    The value of sim.lattice.
    """


class _Selection:
    """
    The value of sim.lattice.select(...).
    """


class _Simulation:
    """
    The value of ImpactX(), the owner of sim.lattice.
    """


class _SimulationClass:
    """
    The class ImpactX under another name, e.g., from impactx import ImpactX as IX.
    """


class _Function:
    """
    A function, lambda or class of the script: calling it runs its code.
    """

    def __init__(self, node: ast.AST):
        self.node = node


class _Instance:
    """
    An instance of a class of the script: its methods run the code of the class.
    """

    def __init__(self, function: _Function):
        self.function = function


LATTICE = _Lattice()


class DashboardLatticeParser:
    """
    Helper class to parse the lattice configuration from Impactx .py-compatible simulation files.

    The parser runs the top-level statements of the script in order, up to the
    first tracking call, e.g., sim.track_particles(), and follows like Python:
        - element constructor calls, e.g., elements.Drift(ds=1.0) or
          elements.BeamMonitor("monitor", backend="h5")
        - lists and tuples of elements, also nested, with * unpacking, + and
          repetition, e.g., cell * 3, indexing and slicing with integers
        - list methods, e.g., cell.append(), cell.reverse(), and assignments to
          items and slices
        - sim.lattice.append(), .extend(), .clear() and .pop_back()

    The body of ``if __name__ == "__main__":`` counts as top level.
    A list or element that a loop, condition, function or other block, or a
    function call, may change, cannot be evaluated afterwards. A lattice that
    uses anything the dashboard cannot evaluate, or other changes of
    sim.lattice, raise a ValueError instead of importing a wrong lattice. So
    do changes of element attributes after construction, e.g., quad.k = 2.0,
    for an element of the lattice.
    """

    def __init__(self, content: str):
        self._content = content
        self._lines: Optional[List[str]] = None  # the lines of content, with their ends
        self._lattice: List[_Element] = []
        self._lattice_unresolved: Optional[str] = None  # reason, if it may have changed
        self._variables: Dict[str, Any] = {}
        # lists and tuples that may have changed, by id, with the reason
        self._unresolved_sequences: Dict[int, tuple] = {}
        self._running_functions: set = set()
        self._simulation: Optional[_Simulation] = None  # the one of the lattice
        self._script_names: set = set()  # the module-level names the script binds
        self._module_names: set = set()  # names bound only by imports

    def parse(self) -> Dict[str, Any]:
        """
        Extracts the lattice configuration from provided ImpactX simulation file.

        Example return:
        {
            "lattice_elements": [
                {"name": "Drift", "parameters": {"ds": "0.5"}},
                {"name": "Quad", "parameters": {"k": "quad_strength", "ds": "0.3"}}
            ]
        }

        :return: Parsed dictionary containing 'lattice_elements'.
        """
        lattice_elements = [
            self.parse_element(element) for element in self.collect_lattice_elements()
        ]
        return {"lattice_elements": lattice_elements}

    def parse_element(self, element: _Element) -> Dict[str, Any]:
        """
        Parses the name and the arguments of an element constructor call.
        Positional arguments get their names from the signature of the element.

        EX:
            elements.Drift(name="d1", ds=1.0, nslice=ns)

            Results in:
            {"name": "Drift", "parameters": {"name": "d1", "ds": "1.0", "nslice": "ns"}}

        :param element: The element.
        :return: A dictionary with the element name and the source code of each argument.
        """
        call = element.call
        name = call.func.attr
        parameters = {}

        # e.g., elements.BeamMonitor("monitor", backend="h5")
        element_class = getattr(impactx_elements, name, None)
        if call.args and element_class is not None:
            positional = InputDefaultsHelper.positional_parameters(element_class)
            if len(call.args) > len(positional) or any(
                isinstance(argument, ast.Starred) for argument in call.args
            ):
                raise ValueError(
                    f"Cannot import '{self._source(call)}': {name} takes at most "
                    f"{len(positional)} positional arguments."
                )
            for parameter_name, argument in zip(positional, call.args):
                parameters[parameter_name] = self._argument_text(argument)

        for keyword in call.keywords:
            if keyword.arg is None:
                raise ValueError(
                    f"Cannot import '{self._source(call)}': unpacking of "
                    f"arguments with ** is not supported."
                )
            if keyword.arg in parameters:
                raise ValueError(
                    f"Cannot import '{self._source(call)}': multiple values for "
                    f"'{keyword.arg}'."
                )
            parameters[keyword.arg] = self._argument_text(keyword.value)

        return {"name": name, "parameters": parameters}

    def collect_lattice_elements(self) -> List[_Element]:
        """
        Runs the top-level statements of the script that build the lattice.

        EX:
            drift1 = elements.Drift(ds=1.0)
            cell = [drift1, elements.Quad(ds=0.5, k=1.0)]
            sim.lattice.extend(cell * 2)

            Results in the elements of:
            [
                elements.Drift(ds=1.0),
                elements.Quad(ds=0.5, k=1.0),
                elements.Drift(ds=1.0),
                elements.Quad(ds=0.5, k=1.0),
            ]

        :return: Each element in sim.lattice, in order.
        """
        tree = parse_script(self._content)
        tracking_call = first_tracking_call(tree, outside_definitions=True)
        check_script_functions(tree, self._source)

        # a tracking method called later under another name, e.g.,
        # track = sim.track_particles, could end the lattice anywhere
        called = {
            id(node.func) for node in ast.walk(tree) if isinstance(node, ast.Call)
        }
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Attribute)
                and node.attr in TRACKING_MODES
                and id(node) not in called
            ):
                raise ValueError(
                    f"Cannot import the lattice: '{self._source(node)}' in line "
                    f"{node.lineno} is not called directly."
                )

        self._lattice = []
        self._lattice_unresolved = None
        self._variables = {}
        self._unresolved_sequences = {}
        self._running_functions = set()
        self._simulation = None
        bindings = module_bindings(tree)
        self._script_names = set(bindings)
        self._module_names = {
            name
            for name, nodes in bindings.items()
            if all(isinstance(node, ast.alias) for node in nodes)
        }

        for statement in top_level_statements(tree.body):
            if (
                isinstance(statement, (ast.FunctionDef, ast.AsyncFunctionDef))
                and not statement.decorator_list
            ):
                self._define_function(statement)
            elif _is_block(statement):
                self._run_block(statement)
            if any(node is tracking_call for node in ast.walk(statement)):
                # its receiver and arguments run before it, later changes of
                # the lattice belong to another simulation
                if isinstance(tracking_call.func, ast.Attribute):
                    self._evaluate(tracking_call.func.value)
                for argument in tracking_call.args:
                    self._evaluate_argument(argument)
                for keyword in tracking_call.keywords:
                    self._evaluate(keyword.value)
                break
            if not _is_block(statement):
                self._run(statement)

        if self._lattice_unresolved:
            raise ValueError(f"Cannot import the lattice: {self._lattice_unresolved}")
        if not self._lattice:
            # e.g., a module of functions that build the lattice, none called
            for definition in ast.walk(tree):
                if not isinstance(definition, DEFINITIONS):
                    continue
                for node in ast.walk(definition):
                    if (
                        isinstance(node, ast.Call)
                        and isinstance(node.func, ast.Attribute)
                        and isinstance(node.func.value, ast.Attribute)
                        and node.func.value.attr == "lattice"
                        and node.func.attr not in LATTICE_READ_METHODS
                    ):
                        raise ValueError(
                            f"Cannot import the lattice: it is built only in "
                            f"'{self._source(node)}' in line {node.lineno}, "
                            f"inside a function or class that does not run "
                            f"before the first tracking call."
                        )
        for element in self._lattice:
            if element.unresolved:
                raise ValueError(f"Cannot import the lattice: {element.unresolved}")

        return list(self._lattice)

    def extract_lattice_inputs(self, parsed_lattice: Dict[str, Any]) -> List[str]:
        """
        Extracts all parameter values from parsed lattice data.
        """
        used_variables = set()

        for element in parsed_lattice.get("lattice_elements", []):
            for value in element["parameters"].values():
                used_variables.add(str(value).strip())

        return list(used_variables)

    # -----------------------------------------------------------------------------
    # Statements
    # -----------------------------------------------------------------------------

    def _run(self, statement: ast.stmt) -> None:
        """
        Runs one top-level statement that is not a block.
        """
        if isinstance(statement, ast.Expr):
            self._evaluate(statement.value)
        elif isinstance(statement, ast.Assign):
            value = self._evaluate(statement.value)
            for target in statement.targets:
                self._assign(target, value)
        elif isinstance(statement, ast.AnnAssign):
            if statement.value is not None:
                value = self._evaluate(statement.value)
                self._assign(statement.target, value)
        elif isinstance(statement, ast.AugAssign):
            self._augmented_assign(statement)
        elif isinstance(statement, ast.Delete):
            for target in statement.targets:
                self._delete(target)
        elif isinstance(statement, ast.Assert):
            self._evaluate(statement.test)  # the message runs only on failure
        elif isinstance(statement, (ast.Import, ast.ImportFrom)):
            for alias in statement.names:
                name = alias.asname or alias.name.split(".")[0]
                if alias.name == "ImpactX" and alias.asname:
                    self._variables[name] = _SimulationClass()  # e.g., as IX
                else:
                    self._variables[name] = _Unresolved(f"'{name}' is imported")
        elif not isinstance(statement, (ast.Global, ast.Nonlocal, ast.Pass)):
            self._run_block(statement)

    def _run_block(self, statement: ast.stmt) -> None:
        """
        Follows a block, e.g., a loop or a function definition, or a statement
        the dashboard does not run: lists and elements it uses may change, so
        they cannot be evaluated afterwards. Changes of sim.lattice raise.
        """
        where = f"in line {statement.lineno}"
        self._raise_if_changes_lattice(statement, where)

        definition = isinstance(statement, DEFINITIONS)
        parameters = set()
        if definition and not isinstance(statement, ast.ClassDef):
            arguments = statement.args
            for argument in (
                arguments.posonlyargs
                + arguments.args
                + arguments.kwonlyargs
                + [arguments.vararg, arguments.kwarg]
            ):
                if argument is not None:
                    parameters.add(argument.arg)

        reason = f"a list or element it uses may be changed by the code {where}"
        changed = f"may be changed by the code {where}"
        for node in self._walk_without_lambdas(statement):
            self._mark_used(node, reason, parameters)
            if isinstance(node, ast.Name) and node.id not in parameters:
                if not isinstance(node.ctx, ast.Load) and not definition:
                    self._variables[node.id] = _Unresolved(f"'{node.id}' {changed}")
            elif isinstance(node, (ast.Global, ast.Nonlocal)):
                for name in node.names:
                    self._mark_unresolved(self._variables.get(name), reason)
                    self._variables[name] = _Unresolved(f"'{name}' {changed}")
            elif isinstance(node, (ast.Import, ast.ImportFrom)) and not definition:
                for alias in node.names:
                    name = alias.asname or alias.name.split(".")[0]
                    self._variables[name] = _Unresolved(f"'{name}' {changed}")
            elif isinstance(node, ast.ExceptHandler) and node.name:
                self._variables[node.name] = _Unresolved(f"'{node.name}' {changed}")
            elif isinstance(node, (ast.MatchAs, ast.MatchStar, ast.MatchMapping)):
                # e.g., case [*cell]: or case {**cell}:
                name = node.rest if isinstance(node, ast.MatchMapping) else node.name
                if name and name not in parameters and not definition:
                    self._variables[name] = _Unresolved(f"'{name}' {changed}")
            elif self._is_lattice_reference(node):
                self._mark_unresolved(LATTICE, reason)

        if definition and hasattr(statement, "name"):
            self._variables[statement.name] = _Function(statement)

    def _define_function(self, statement: ast.stmt) -> None:
        """
        Defines a top-level function without decorators: its default values
        and annotations run now, its body when it is called. A function whose
        code uses the simulation or module-level variables cannot be used
        before the first tracking call (@see check_script_functions), and a
        call of another one follows its code then.
        """
        reason = (
            f"it is a default value of '{statement.name}' in line "
            f"{statement.lineno}, which may change it"
        )
        arguments = statement.args
        defaults = arguments.defaults + [d for d in arguments.kw_defaults if d]
        for default in defaults:
            self._mark_unresolved(self._evaluate(default), reason)
        for node in _definition_time_nodes(statement):
            if not any(node is default for default in defaults):
                self._run_block_of_expression(node)  # annotations
        self._variables[statement.name] = _Function(statement)

    def _new_simulation(self, node: ast.Call) -> "_Simulation":
        """
        Creates the simulation, e.g., sim = ImpactX(). Only one is supported.
        """
        if self._simulation is not None:
            raise ValueError(
                f"Cannot import the lattice: '{self._source(node)}' in line "
                f"{node.lineno} creates a second simulation."
            )
        self._simulation = _Simulation()
        return self._simulation

    def _call_function(
        self, function: _Function, node: Optional[ast.AST] = None
    ) -> Any:
        """
        Follows a call of a function, lambda or class of the script: the lists
        and elements its code uses now may change. Calling a class creates an
        instance.
        """
        if id(function.node) not in self._running_functions:
            self._running_functions.add(id(function.node))
            try:
                self._run_block(function.node)
            finally:
                self._running_functions.discard(id(function.node))
        if isinstance(function.node, ast.ClassDef):
            return _Instance(function)
        call = self._source(node) if node is not None else "a call"
        return _Unresolved(f"'{call}' cannot be evaluated")

    def _assign(self, target: ast.expr, value: Any) -> None:
        """
        Binds a value to an assignment target, like Python.
        """
        where = f"in line {target.lineno}"

        if isinstance(target, ast.Name):
            self._variables[target.id] = value
            return

        if isinstance(target, (ast.Tuple, ast.List)):
            starred = [
                i for i, item in enumerate(target.elts) if isinstance(item, ast.Starred)
            ]
            after = len(target.elts) - 1 - starred[0] if len(starred) == 1 else 0
            if (
                isinstance(value, (list, tuple))
                and not self._unresolved_reason(value)
                and len(starred) == 1
                and len(value) >= len(target.elts) - 1
            ):
                # e.g., first, *rest = cells: rest is a new list of the same items
                stop = len(value) - after
                parts = list(value[: starred[0]])
                parts.append(self._inherit(list(value[starred[0] : stop]), value))
                parts += list(value[stop:])
                for item_target, item_value in zip(target.elts, parts):
                    if isinstance(item_target, ast.Starred):
                        item_target = item_target.value
                    self._assign(item_target, item_value)
                return
            if (
                starred
                or not (
                    isinstance(value, (list, tuple)) and len(value) == len(target.elts)
                )
                or self._unresolved_reason(value)
            ):
                # the parts may hold the lists and elements of value, e.g., a, *b = x
                self._mark_unresolved(
                    value,
                    f"it is unpacked in line {target.lineno}, which cannot be followed",
                )
                self._run_block_of_target(target)
                return
            for item_target, item_value in zip(target.elts, value):
                self._assign(item_target, item_value)
            return

        if isinstance(target, ast.Attribute):
            if target.attr == "lattice":
                self._cannot_follow(LATTICE, target, where)  # sim.lattice = [...]
            owner = self._evaluate(target.value)
            if isinstance(owner, _Element):
                self._change_element(owner, target, where)
            elif isinstance(owner, (_Lattice, _Selection)):
                self._cannot_follow(owner, target, where)
            else:
                # e.g., sim.periods = 2 does not change the lattice, but
                # cells[i].ds = 2 changes an element the dashboard cannot tell,
                # and obj.cell = cell stores a list the dashboard cannot follow
                self._run_block_of_expression(target.value)
                self._mark_unresolved(
                    value,
                    f"a list or element it uses is stored in "
                    f"'{self._source(target)}' {where}, whose changes cannot be "
                    f"followed",
                )
            return

        if isinstance(target, ast.Subscript):
            container = self._evaluate(target.value)
            if isinstance(container, (_Lattice, _Selection)):
                self._cannot_follow(container, target, where)
            elif isinstance(container, list):
                key = self._evaluate_key(target.slice)
                try:
                    if isinstance(key, slice) and isinstance(value, (list, tuple)):
                        container[key] = list(value)
                        self._inherit(container, value)
                    elif type(key) is int:
                        container[key] = value
                    else:
                        raise TypeError
                except (IndexError, TypeError, ValueError):
                    self._cannot_follow(container, target, where)
            else:
                reason = f"the change of '{self._source(target)}' {where}"
                self._mark_unresolved(container, reason)
                self._mark_unresolved(value, reason)
                self._mark_unresolved(self._evaluate_key(target.slice), reason)
            return

        self._run_block_of_target(target)

    def _run_block_of_target(self, target: ast.expr) -> None:
        """
        An assignment target whose values the dashboard does not know, e.g.,
        a, *b = cell: each part is assigned an unresolved value.
        """
        if isinstance(target, (ast.Tuple, ast.List)):
            for item in target.elts:
                self._run_block_of_target(item)
        elif isinstance(target, ast.Starred):
            self._run_block_of_target(target.value)
        else:
            reason = f"'{self._source(target)}' is assigned in line {target.lineno}"
            self._assign(target, _Unresolved(f"{reason}, which cannot be followed"))

    def _augmented_assign(self, statement: ast.AugAssign) -> None:
        """
        Runs an in-place operation like cell += [drift] or n -= 1.
        """
        target = statement.target
        value = self._evaluate(statement.value)
        where = f"in line {statement.lineno}"

        if isinstance(target, ast.Name):
            current = self._variables.get(target.id)
            result = self._binary_operation(statement.op, current, value, statement)
            if isinstance(current, list) and isinstance(result, list):
                current[:] = result  # in place, i.e., for all names of this list
                self._inherit(current, result)
            elif isinstance(current, (list, _Lattice)):
                self._cannot_follow(current, statement, where)
            else:
                self._variables[target.id] = result
            return

        container = self._evaluate(
            target.value
            if isinstance(target, (ast.Attribute, ast.Subscript))
            else target
        )
        if isinstance(container, _Element) and isinstance(target, ast.Attribute):
            self._change_element(container, target, where)
        elif isinstance(target, ast.Attribute) and target.attr == "lattice":
            self._cannot_follow(LATTICE, statement, where)
        else:
            self._cannot_follow(container, statement, where)

    def _change_element(self, element: _Element, target: ast.expr, where: str):
        """
        Marks an element whose attribute the script changes, e.g., quad.k = 2.0:
        the dashboard imports only the arguments of its constructor.
        """
        element.unresolved = element.unresolved or (
            f"'{self._source(target)}' {where} changes an element after its "
            f"construction, which the dashboard does not import"
        )

    def _take_element_attribute(self, element: _Element, node: ast.Attribute):
        """
        Marks an element whose method the script takes, e.g.,
        set_k = quad.__setattr__ or flip = quad.reverse: the method may change
        the element later. Its parameters, e.g., quad.k, other attributes that
        are not methods, e.g., quad.__class__, and the methods that only read
        it, e.g., quad.to_dict, do not.
        """
        element_class = getattr(impactx_elements, element.call.func.attr, None)
        missing = object()
        attribute = inspect.getattr_static(element_class, node.attr, missing)
        if node.attr in ELEMENT_READ_METHODS or not (
            attribute is missing or callable(attribute)
        ):
            return
        element.unresolved = element.unresolved or (
            f"'{self._source(node)}' in line {node.lineno} may change an element "
            f"after its construction, which the dashboard does not import"
        )

    def _delete(self, target: ast.expr) -> None:
        """
        Runs del on a name, an item or a slice.
        """
        where = f"in line {target.lineno}"

        if isinstance(target, ast.Name):
            self._variables.pop(target.id, None)
        elif isinstance(target, ast.Subscript):
            container = self._evaluate(target.value)
            if isinstance(container, (_Lattice, _Selection)):
                self._cannot_follow(container, target, where)
            elif isinstance(container, list):
                key = self._evaluate_key(target.slice)
                try:
                    if not (isinstance(key, slice) or type(key) is int):
                        raise TypeError
                    del container[key]
                except (IndexError, TypeError):
                    self._cannot_follow(container, target, where)
            else:
                reason = f"the change of '{self._source(target)}' {where}"
                self._mark_unresolved(container, reason)
                self._mark_unresolved(self._evaluate_key(target.slice), reason)
        elif isinstance(target, (ast.Tuple, ast.List)):
            for item in target.elts:
                self._delete(item)
        elif isinstance(target, ast.Attribute):
            owner = self._evaluate(target.value)
            if isinstance(owner, _Element):
                self._change_element(owner, target, where)
            elif target.attr == "lattice":
                self._cannot_follow(LATTICE, target, where)
            else:
                self._run_block_of_expression(target.value)

    # -----------------------------------------------------------------------------
    # Expressions
    # -----------------------------------------------------------------------------

    def _evaluate(self, node: ast.expr) -> Any:
        """
        Evaluates an expression, including its changes of lists, elements and
        sim.lattice, in the order of Python.

        :return: an _Element, a list, a tuple, an int, LATTICE, a _Selection,
                 or an _Unresolved value.
        """
        if isinstance(node, ast.Constant):
            if type(node.value) is int:
                return node.value
            return _Inert(f"'{self._source(node)}' is not an element")

        if isinstance(node, ast.Name):
            return self._variables.get(
                node.id, _Unresolved(f"'{node.id}' is not defined at the top level")
            )

        if isinstance(node, (ast.List, ast.Tuple)):
            values = []
            sources = []
            unknown_length = None
            for item in node.elts:
                if isinstance(item, ast.Starred):
                    value = self._evaluate(item.value)
                    sources.append(value)
                    if isinstance(value, (list, tuple)):
                        values.extend(value)
                    else:
                        # e.g., [*b] of an unknown b: unknown items and length
                        unknown_length = f"'{self._source(item)}' cannot be evaluated"
                        values.append(_Unresolved(unknown_length))
                else:
                    values.append(self._evaluate(item))
            result = values if isinstance(node, ast.List) else tuple(values)
            if unknown_length:
                # e.g., lines[0] of lines = [arc, *other] may be any item
                self._mark_unresolved(result, unknown_length)
            return self._inherit(result, *sources)

        if isinstance(node, ast.Attribute):
            if node.attr == "lattice":
                owner = self._evaluate(node.value)
                if (
                    isinstance(node.value, ast.Name)
                    and node.value.id not in self._script_names
                    and self._simulation is None
                ):
                    # a part of a script, without ImpactX(): the simulation is
                    # the name the script does not bind, e.g., sim.lattice
                    owner = self._variables[node.value.id] = _Simulation()
                    self._simulation = owner
                if isinstance(owner, _Simulation):
                    return LATTICE
                return _Unresolved(
                    f"'{self._source(node)}' is not the lattice of the simulation"
                )
            owner = self._evaluate(node.value)
            # e.g., append = cell.append: the method may change its owner later
            reason = f"'{self._source(node)}' in line {node.lineno} may change it"
            if isinstance(owner, list) and node.attr not in LIST_READ_METHODS:
                self._mark_unresolved(owner, reason)
            elif isinstance(owner, _Lattice) and node.attr not in LATTICE_READ_METHODS:
                self._mark_unresolved(owner, reason)
            elif isinstance(owner, _Selection) and (
                node.attr not in SELECTION_READ_METHODS
            ):
                self._mark_unresolved(owner, reason)
            elif isinstance(owner, _Element):
                self._take_element_attribute(owner, node)
            return _Unresolved(f"'{self._source(node)}' cannot be evaluated")

        if isinstance(node, ast.Subscript):
            return self._subscript(node)

        if isinstance(node, ast.UnaryOp):
            value = self._evaluate(node.operand)
            if isinstance(node.op, ast.USub) and type(value) is int:
                return -value
            if isinstance(node.op, ast.UAdd) and type(value) is int:
                return value
            if type(value) is int or isinstance(value, _Inert):
                return _Inert(f"'{self._source(node)}' is not an element")  # e.g., -1.0
            return _Unresolved(f"'{self._source(node)}' cannot be evaluated")

        if isinstance(node, ast.BinOp):
            left = self._evaluate(node.left)
            right = self._evaluate(node.right)
            return self._binary_operation(node.op, left, right, node)

        if isinstance(node, ast.NamedExpr):
            value = self._evaluate(node.value)
            self._variables[node.target.id] = value
            return value

        if isinstance(node, ast.Call):
            return self._call(node)

        if isinstance(node, ast.Lambda):
            self._run_block(node)  # like a function definition
            return _Function(node)

        if isinstance(node, ast.Compare):
            values = [self._evaluate(part) for part in [node.left] + node.comparators]
            if not all(isinstance(op, (ast.Is, ast.IsNot)) for op in node.ops):
                self._mark_operands(values, node)  # e.g., obj == cell, cell in obj
            return _Inert(f"'{self._source(node)}' is not an element")

        if isinstance(node, (ast.JoinedStr, ast.FormattedValue)):
            parts = node.values if isinstance(node, ast.JoinedStr) else [node.value]
            for part in parts:
                self._evaluate(part)
            return _Inert(f"'{self._source(node)}' is not an element")

        if isinstance(node, (ast.Dict, ast.Set)):
            parts = node.elts if isinstance(node, ast.Set) else node.keys + node.values
            reason = (
                f"a list or element it uses is put into '{self._source(node)}' in "
                f"line {node.lineno}, whose changes cannot be followed"
            )
            for part in parts:
                if part is not None:
                    self._mark_unresolved(self._evaluate(part), reason)
            return _Unresolved(f"'{self._source(node)}' is not an element")

        # other expressions, e.g., comprehensions: lists and elements they use may change
        self._run_block_of_expression(node)
        return _Unresolved(f"'{self._source(node)}' cannot be evaluated")

    def _run_block_of_expression(self, node: ast.expr) -> None:
        """
        Follows an expression that the dashboard does not evaluate.
        """
        reason = (
            f"a list or element it uses may be changed by the code in line "
            f"{node.lineno}"
        )
        self._raise_if_changes_lattice(node, f"in line {node.lineno}")
        for child in self._walk_without_lambdas(node):
            self._mark_used(child, reason, set())
            if self._is_lattice_reference(child):
                self._mark_unresolved(LATTICE, reason)
            if isinstance(child, ast.Name) and isinstance(child.ctx, ast.Store):
                # e.g., (n := 3) in a condition
                self._variables[child.id] = _Unresolved(
                    f"'{child.id}' {reason.removeprefix('a list or element it uses ')}"
                )

    def _mark_used(self, node: ast.AST, reason: str, parameters: set) -> None:
        """
        Marks the lists and elements that a name holds as unresolved, and the
        lattice, if the simulation is passed to a call, e.g., build(sim).
        Names in parameters belong to a function and are skipped.
        """
        if isinstance(node, ast.Call):
            for argument in node.args + [keyword.value for keyword in node.keywords]:
                if isinstance(argument, ast.Starred):
                    argument = argument.value
                if isinstance(argument, ast.Name) and argument.id not in parameters:
                    value = self._variables.get(argument.id)
                    if isinstance(value, _Simulation):
                        self._mark_unresolved(value, reason)
        elif isinstance(node, ast.Name) and node.id not in parameters:
            value = self._variables.get(node.id)
            if not isinstance(value, _Simulation):
                self._mark_unresolved(value, reason)

    def _binary_operation(self, op: ast.operator, left: Any, right: Any, node) -> Any:
        """
        Evaluates +, - and * of integers, and + and * of lists and tuples.
        """
        if type(left) is int and type(right) is int:
            if isinstance(op, ast.Add):
                return left + right
            if isinstance(op, ast.Sub):
                return left - right
            if isinstance(op, ast.Mult):
                return left * right
        elif isinstance(op, ast.Add):
            if isinstance(left, list) and isinstance(right, (list, tuple)):
                # list + tuple raises in Python, but list += tuple extends
                if isinstance(right, list) or not isinstance(node, ast.BinOp):
                    return self._inherit(left + list(right), left, right)
            if isinstance(left, tuple) and isinstance(right, tuple):
                return self._inherit(left + right, left, right)
        elif isinstance(op, ast.Mult):
            if type(left) is int:
                left, right = right, left
            if isinstance(left, (list, tuple)) and type(right) is int:
                items, elements = self._count_elements(left)
                reason = None
                if elements * right > MAX_REPEATED_ELEMENTS:
                    reason = f"has more than {MAX_REPEATED_ELEMENTS} elements"
                elif items * right > MAX_LATTICE_ELEMENTS:
                    reason = f"has more than {MAX_LATTICE_ELEMENTS} items"
                if reason:
                    # not kept, e.g., history = [math.nan] * n_turns, it raises only
                    # in the lattice; its lists and elements may change through it
                    reason = f"'{self._source(node)}' {reason}"
                    self._mark_unresolved(left, reason)
                    return _Unresolved(reason)
                return self._inherit(left * right, left)

        # numbers and text, e.g., 2.0 * 1e-3, run no code of the script, and
        # text formatting, e.g., "%s" % cell, only reads the list
        if isinstance(left, (int, _Inert)) and (
            isinstance(right, (int, _Inert)) or isinstance(op, ast.Mod)
        ):
            return _Inert(f"'{self._source(node)}' is not an element")

        # the result may hold the items of a list operand, e.g., [arc] * len(names),
        # or obj + cell may run code of obj with it
        reason = (
            f"'{self._source(node)}' in line {node.lineno} cannot be followed, and may "
            f"use or change it"
        )
        for value in (left, right):
            if isinstance(value, (list, tuple)):
                self._mark_unresolved(value, reason)
        return _Unresolved(f"'{self._source(node)}' cannot be evaluated")

    def _mark_operands(self, values: List[Any], node: ast.AST) -> None:
        """
        Marks the lists among the operands of an operator as unresolved, if
        another operand may run code with them, e.g., obj + cell runs the
        __radd__ method of an object the dashboard does not know.
        """
        known = KNOWN_OPERANDS + (_Element, _Inert, _Lattice, _Selection)
        if all(isinstance(value, known) for value in values):
            return
        reason = (
            f"'{self._source(node)}' in line {node.lineno} may run code that changes it"
        )
        for value in values:
            if isinstance(value, (list, tuple)):
                self._mark_unresolved(value, reason)

    def _evaluate_key(self, node: ast.expr) -> Any:
        """
        Evaluates an index, an int, or a slice of ints.
        """
        if isinstance(node, ast.Slice):
            bounds = [
                None if bound is None else self._evaluate(bound)
                for bound in (node.lower, node.upper, node.step)
            ]
            if all(bound is None or type(bound) is int for bound in bounds):
                if bounds[2] != 0:
                    return slice(*bounds)
            return None
        return self._evaluate(node)

    def _subscript(self, node: ast.Subscript) -> Any:
        """
        Evaluates an item or a slice of a list, a tuple or sim.lattice.
        """
        container = self._evaluate(node.value)
        key = self._evaluate_key(node.slice)
        if isinstance(key, (list, tuple)) and not isinstance(container, KNOWN_OPERANDS):
            self._mark_operands([container, key], node)  # e.g., obj[cell]
        items = self._lattice if isinstance(container, _Lattice) else container
        reason = self._unresolved_reason(container)
        if reason:
            if isinstance(key, slice):
                return self._inherit(list(items[key]), container)
            return _Unresolved(reason)
        if isinstance(items, (list, tuple)):
            if isinstance(key, slice):
                result = items[key]
                return list(result) if isinstance(container, _Lattice) else result
            if type(key) is int and -len(items) <= key < len(items):
                return items[key]
            # e.g., cells[len(cells) - 1] may be any item, which may change later
            self._mark_unresolved(
                container,
                f"'{self._source(node)}' in line {node.lineno} cannot be followed",
            )
        return _Unresolved(f"'{self._source(node)}' cannot be evaluated")

    def _call(self, node: ast.Call) -> Any:
        """
        Evaluates a call, with its effects on lists, elements and sim.lattice.
        """
        function = node.func

        if (isinstance(function, ast.Name) and function.id == "ImpactX") or (
            isinstance(function, ast.Attribute) and function.attr == "ImpactX"
        ):
            for argument in node.args + [keyword.value for keyword in node.keywords]:
                self._evaluate_argument(argument)
            return self._new_simulation(node)

        if self._is_element_constructor(node):
            for argument in node.args + [keyword.value for keyword in node.keywords]:
                self._evaluate(
                    argument.value if isinstance(argument, ast.Starred) else argument
                )
            return _Element(node)

        if isinstance(function, ast.Attribute):
            owner = self._evaluate(function.value)
            if (
                isinstance(function.value, ast.Attribute)
                and function.value.attr == "lattice"
                and not isinstance(owner, _Lattice)
                and function.attr not in LATTICE_READ_METHODS
            ):
                raise ValueError(
                    f"Cannot import the lattice: '{self._source(node)}' in line "
                    f"{node.lineno} changes a lattice that is not the one of the "
                    f"simulation."
                )
            arguments = [self._evaluate_argument(argument) for argument in node.args]
            keywords = [self._evaluate(keyword.value) for keyword in node.keywords]
            method = function.attr

            if isinstance(owner, _Lattice):
                return self._lattice_method(node, method, arguments, keywords)
            if isinstance(owner, _Simulation) and "lattice" in method:
                raise ValueError(
                    f"Cannot import the lattice: sim.{method}() in line "
                    f"{node.lineno} is not supported."
                )
            if isinstance(owner, (_Function, _Instance)):
                # e.g., a method of a class of the script
                self._call_function(
                    owner if isinstance(owner, _Function) else owner.function
                )
            if isinstance(owner, _Selection):
                if method in SELECTION_READ_METHODS:
                    if method == "select":
                        return _Selection()
                    return _Unresolved(f"'{self._source(node)}' is not an element")
                self._cannot_follow(LATTICE, node, f"in line {node.lineno}")
            if isinstance(owner, list):
                if keywords or any(isinstance(a, ast.Starred) for a in node.args):
                    if method not in LIST_READ_METHODS:
                        self._cannot_follow(owner, node, f"in line {node.lineno}")
                    return _Unresolved(f"'{self._source(node)}' is not an element")
                return self._list_method(node, owner, method, arguments)
            if isinstance(owner, tuple) and method in LIST_READ_METHODS:
                return _Unresolved(f"'{self._source(node)}' is not an element")
            if isinstance(owner, _Element):
                if method == "copy" and not arguments and not keywords:
                    copy = _Element(owner.call)
                    copy.unresolved = owner.unresolved
                    return copy
                if method in ELEMENT_READ_METHODS:
                    return _Unresolved(f"'{self._source(node)}' is not an element")
                owner.unresolved = (
                    f"'{self._source(node)}' in line {node.lineno} cannot be followed"
                )
                return _Unresolved(f"'{self._source(node)}' cannot be evaluated")

            if (
                isinstance(function.value, ast.Constant)
                and isinstance(function.value.value, str)
            ) or isinstance(function.value, ast.JoinedStr):
                # e.g., "{}".format(cell) formats its arguments as text
                return _Inert(f"'{self._source(node)}' is not an element")

            # a method of an object the dashboard does not know may change its arguments
            reason = (
                f"a list or element it uses is passed to '{self._source(function)}' "
                f"in line {node.lineno}, which may change it"
            )
            for value in arguments + keywords:
                self._mark_unresolved(value, reason)
            return _Unresolved(f"'{self._source(node)}' cannot be evaluated")

        arguments = [self._evaluate_argument(argument) for argument in node.args]
        keywords = [self._evaluate(keyword.value) for keyword in node.keywords]

        if isinstance(function, ast.Name) and function.id in self._variables:
            callee = self._variables[function.id]  # e.g., a function of the script
        elif isinstance(function, ast.Name):
            callee = None  # a built-in
            if (
                function.id in ("list", "tuple")
                and not keywords
                and len(arguments) <= 1
            ):
                items = arguments[0] if arguments else []
                if isinstance(items, _Lattice):
                    items = self._lattice
                if isinstance(items, (list, tuple)):
                    result = list(items) if function.id == "list" else tuple(items)
                    return self._inherit(result, items)
            if function.id in SAFE_FUNCTIONS:
                if function.id in INERT_FUNCTIONS:
                    return _Inert(f"'{self._source(node)}' is not an element")
                return _Unresolved(f"'{self._source(node)}' is not an element")
        else:
            callee = self._evaluate(function)

        if isinstance(callee, _SimulationClass):
            return self._new_simulation(node)

        result = _Unresolved(f"'{self._source(node)}' cannot be evaluated")
        if isinstance(callee, _Function):
            result = self._call_function(callee, node)
        elif isinstance(callee, _Instance):
            self._call_function(callee.function)

        # a function the dashboard does not know may change its arguments
        reason = (
            f"a list or element it uses is passed to '{self._source(function)}' "
            f"in line {node.lineno}, which may change it"
        )
        for value in arguments + keywords:
            self._mark_unresolved(value, reason)
        return result

    def _evaluate_argument(self, argument: ast.expr) -> Any:
        """
        Evaluates a positional argument; *values gives a tuple of the values.
        """
        if isinstance(argument, ast.Starred):
            value = self._evaluate(argument.value)
            if isinstance(value, _Lattice):
                return tuple(self._lattice)
            if isinstance(value, (list, tuple)):
                return self._inherit(tuple(value), value)
            return value
        return self._evaluate(argument)

    def _lattice_method(self, node: ast.Call, method: str, arguments, keywords) -> Any:
        """
        Runs a method of sim.lattice.
        """
        where = f"in line {node.lineno}"

        if method in LATTICE_READ_METHODS:
            if method == "select":
                return _Selection()
            return _Unresolved(f"'{self._source(node)}' is not an element")

        if method not in LATTICE_SUPPORTED_METHODS or keywords:
            raise ValueError(
                f"Cannot import the lattice: sim.lattice.{method}() {where} is not "
                f"supported."
            )

        expected_arguments = {"append": 1, "extend": 1}
        if len(arguments) != expected_arguments.get(method, 0) or any(
            isinstance(argument, ast.Starred) for argument in node.args
        ):
            raise ValueError(
                f"Cannot import the lattice: unsupported arguments in "
                f"'{self._source(node)}'."
            )

        if method == "clear":
            self._lattice.clear()
        elif method == "pop_back":
            if self._lattice:
                self._lattice.pop()
        elif method == "append":
            if not isinstance(arguments[0], _Element):
                self._flatten(arguments[0])  # raises with the reason
                self._cannot_follow(LATTICE, node, where)
            self._lattice.append(arguments[0])
        else:
            self._lattice.extend(self._flatten(arguments[0]))

        if len(self._lattice) > MAX_LATTICE_ELEMENTS:
            raise ValueError(
                f"Cannot import the lattice: it has more than "
                f"{MAX_LATTICE_ELEMENTS} elements."
            )
        return _Unresolved(f"'{self._source(node)}' is not an element")

    def _list_method(self, node: ast.Call, values: list, method: str, arguments) -> Any:
        """
        Runs a method of a list.
        """
        result = _Unresolved(f"'{self._source(node)}' is not an element")
        reason = self._unresolved_reason(values)
        try:
            if method in LIST_READ_METHODS:
                if method == "copy" and not arguments:
                    return self._inherit(list(values), values)
                return result
            if method == "append" and len(arguments) == 1:
                values.append(arguments[0])
            elif method == "extend" and len(arguments) == 1:
                if isinstance(arguments[0], _Lattice):
                    values.extend(self._lattice)
                elif isinstance(arguments[0], (list, tuple)):
                    values.extend(arguments[0])
                    self._inherit(values, arguments[0])
                else:
                    raise TypeError
            elif method == "insert" and len(arguments) == 2:
                if type(arguments[0]) is not int:
                    raise TypeError
                values.insert(arguments[0], arguments[1])
            elif method == "pop" and len(arguments) <= 1:
                if arguments and type(arguments[0]) is not int:
                    raise TypeError
                popped = values.pop(*arguments)
                return _Unresolved(reason) if reason else popped
            elif method == "remove" and len(arguments) == 1:
                # Python removes the first equal item, and equal elements of the
                # same kind compare equal: follow only an unambiguous removal
                for index, item in enumerate(values):
                    if item is arguments[0]:
                        del values[index]
                        break
                    if self._may_be_equal(item, arguments[0]):
                        raise ValueError
                else:
                    raise ValueError
            elif method == "reverse" and not arguments:
                values.reverse()
            elif method == "clear" and not arguments:
                values.clear()  # a mark stays: e.g., a stored append may add later
            else:
                raise TypeError
        except (IndexError, TypeError, ValueError):
            self._cannot_follow(values, node, f"in line {node.lineno}")
        return result

    # -----------------------------------------------------------------------------
    # Helpers
    # -----------------------------------------------------------------------------

    def _unresolved_reason(self, value: Any) -> Optional[str]:
        """
        Returns the reason if a list or tuple may have changed, else None.
        """
        entry = self._unresolved_sequences.get(id(value))
        return entry[1] if entry is not None and entry[0] is value else None

    def _inherit(self, result: Any, *sources: Any) -> Any:
        """
        Marks a list or tuple made from others, e.g., a slice or a copy, as
        unresolved if one of them is.
        """
        if isinstance(result, (list, tuple)):
            for source in sources:
                reason = self._unresolved_reason(source)
                if reason:
                    self._unresolved_sequences[id(result)] = (result, reason)
                    break
        return result

    @staticmethod
    def _count_elements(value: Any) -> tuple:
        """
        Counts the items of a list or tuple, with those of nested ones, and of
        these the ones that may be elements: not, e.g., numbers or None.
        """
        items = elements = 0
        pending = [(value, frozenset())]
        while pending:
            item, enclosing = pending.pop()
            items += 1  # also lists and tuples, e.g., of [[]] * n
            if isinstance(item, (list, tuple)) and id(item) not in enclosing:
                pending.extend((sub, enclosing | {id(item)}) for sub in item)
            else:
                if isinstance(item, (_Element, _Unresolved)) and not isinstance(
                    item, _Inert
                ):
                    elements += 1
            if items > MAX_LATTICE_ELEMENTS:
                break
        return items, elements

    @staticmethod
    def _may_be_equal(item: Any, value: Any) -> bool:
        """
        Checks if Python could find two values equal: elements of the same kind
        compare by their parameters, lists and tuples by their items.
        """
        if isinstance(item, _Element) and isinstance(value, _Element):
            return item.call.func.attr == value.call.func.attr
        if isinstance(item, (list, tuple)) and isinstance(value, (list, tuple)):
            return True
        if type(item) is int and type(value) is int:
            return item == value
        return isinstance(item, _Unresolved) or isinstance(value, _Unresolved)

    def _mark_unresolved(self, value: Any, reason: str) -> None:
        """
        Marks a list, with the lists and elements in it, or an element, as
        changed in a way the dashboard cannot follow. Lists change in place, so
        every name and list that holds them sees the mark.
        """
        seen = set()

        def mark(item: Any) -> None:
            if id(item) in seen:
                return
            seen.add(id(item))
            if isinstance(item, _Element):
                item.unresolved = item.unresolved or reason
            elif isinstance(item, (list, tuple)):
                for sub_item in item:
                    mark(sub_item)
                if not self._unresolved_reason(item):
                    self._unresolved_sequences[id(item)] = (item, reason)
                if isinstance(item, list):
                    item[:] = [_Unresolved(reason)]
            elif isinstance(item, (_Function, _Instance)):
                # e.g., a function passed to map() may be called
                self._call_function(
                    item if isinstance(item, _Function) else item.function
                )
            elif isinstance(item, (_Lattice, _Selection, _Simulation)):
                self._lattice_unresolved = self._lattice_unresolved or reason
                for element in self._lattice:
                    mark(element)

        mark(value)

    def _cannot_follow(self, value: Any, node: ast.AST, where: str) -> None:
        """
        Raises for sim.lattice, else marks the value as unresolved.
        """
        if isinstance(value, (_Lattice, _Selection)):
            raise ValueError(
                f"Cannot import the lattice: '{self._source(node)}' {where} is not "
                f"supported."
            )
        self._mark_unresolved(
            value, f"the change '{self._source(node)}' {where} cannot be followed"
        )

    def _is_lattice_reference(self, node: ast.AST) -> bool:
        """
        Checks for sim.lattice, or a name that holds it.
        """
        if isinstance(node, ast.Attribute):
            return node.attr == "lattice"
        if isinstance(node, ast.Name):
            return isinstance(self._variables.get(node.id), (_Lattice, _Selection))
        return False

    def _raise_if_changes_lattice(self, node: ast.AST, where: str) -> None:
        """
        Raises if node changes sim.lattice: a call of a method that changes it,
        or an assignment to sim.lattice, its items or slices.
        """
        for child in ast.walk(node):
            if isinstance(child, ast.Call) and isinstance(child.func, ast.Attribute):
                if (
                    self._is_lattice_reference(child.func.value)
                    and child.func.attr not in LATTICE_READ_METHODS
                ) or (
                    "lattice" in child.func.attr
                    and isinstance(child.func.value, ast.Name)
                    and isinstance(
                        self._variables.get(child.func.value.id), _Simulation
                    )
                ):
                    raise ValueError(
                        f"Cannot import the lattice: '{self._source(child)}' "
                        f"{where} is inside a loop, condition, function or other "
                        f"block."
                    )
            if isinstance(getattr(child, "ctx", None), (ast.Store, ast.Del)):
                if self._is_lattice_reference(child) or (
                    isinstance(child, ast.Subscript)
                    and self._is_lattice_reference(child.value)
                ):
                    raise ValueError(
                        f"Cannot import the lattice: the change of "
                        f"'{self._source(child)}' {where} is inside a loop, "
                        f"condition, function or other block."
                    )

    @staticmethod
    def _walk_without_lambdas(node: ast.AST):
        """
        Like ast.walk, but does not enter lambdas, whose code does not run where
        it is defined.
        """
        yield node
        for child in ast.iter_child_nodes(node):
            if not isinstance(child, ast.Lambda):
                yield from DashboardLatticeParser._walk_without_lambdas(child)

    def _is_element_constructor(self, call: ast.Call) -> bool:
        """
        Checks for calls like elements.Drift(...) or impactx.elements.Drift(...)
        on the module impactx.elements: not, e.g., elements.append(...) on a list
        named elements, or elements.isclose(...). A name that is no class of the
        module is an element the dashboard does not support.
        """
        function = call.func
        if not isinstance(function, ast.Attribute):
            return False
        value = getattr(impactx_elements, function.attr, None)
        if value is not None and not isinstance(value, type):
            return False
        module = function.value
        if isinstance(module, ast.Attribute) and module.attr == "elements":
            module = module.value  # e.g., impactx.elements
        elif not (isinstance(module, ast.Name) and module.id == "elements"):
            return False
        # a module that the script imports, or a name a part of a script does not bind
        return isinstance(module, ast.Name) and (
            module.id in self._module_names or module.id not in self._script_names
        )

    def _flatten(self, value: Any) -> List[_Element]:
        """
        Expands nested lists and tuples of elements into a list of elements.
        """
        flat: List[_Element] = []

        def expand(item: Any, enclosing: frozenset) -> None:
            if isinstance(item, _Element):
                if len(flat) == MAX_LATTICE_ELEMENTS:
                    raise ValueError(
                        f"Cannot import the lattice: it has more than "
                        f"{MAX_LATTICE_ELEMENTS} elements."
                    )
                flat.append(item)
            elif isinstance(item, _Lattice):
                for element in list(self._lattice):
                    expand(element, enclosing)
            elif isinstance(item, (list, tuple)):
                if id(item) in enclosing:
                    raise ValueError(
                        "Cannot import the lattice: a list contains itself."
                    )
                reason = self._unresolved_reason(item)
                if reason:
                    raise ValueError(f"Cannot import the lattice: {reason}.")
                for sub_item in item:
                    expand(sub_item, enclosing | {id(item)})
            elif isinstance(item, _Unresolved):
                raise ValueError(f"Cannot import the lattice: {item.reason}.")
            else:
                raise ValueError(
                    f"Cannot import the lattice: '{item!r}' is not an element."
                )

        expand(value, frozenset())
        return flat

    def _argument_text(self, node: ast.expr) -> str:
        """
        Returns an argument as text for the dashboard input field.
        """
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            return node.value
        if isinstance(node, (ast.List, ast.Tuple)):
            # one line, with the source code of each item, e.g., [0.0, -1.5e-3]
            items = ", ".join(self._source(item) for item in node.elts)
            if isinstance(node, ast.List):
                return f"[{items}]"
            return f"({items},)" if len(node.elts) == 1 else f"({items})"
        return self._source(node)

    def _source(self, node: ast.AST) -> str:
        """
        Returns the source code of a node, like ast.get_source_segment, but
        splits the lines of the file only once.
        """
        if self._lines is None:
            # the line breaks of Python: a form feed does not end a line
            self._lines = re.findall(r"[^\r\n]*(?:\r\n|\r|\n)|[^\r\n]+$", self._content)
        text = None
        if getattr(node, "end_lineno", None) is not None:
            first, last = node.lineno - 1, node.end_lineno - 1
            try:
                if first == last:
                    line = self._lines[first].encode()
                    text = line[node.col_offset : node.end_col_offset].decode()
                else:
                    start = self._lines[first].encode()[node.col_offset :].decode()
                    end = self._lines[last].encode()[: node.end_col_offset].decode()
                    text = start + "".join(self._lines[first + 1 : last]) + end
            except (IndexError, UnicodeDecodeError):
                text = None
        return text or ast.unparse(node)
