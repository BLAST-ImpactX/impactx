"""
A MAD-X parser that aims to be compliant with the upstream MAD-X format.

References:
    - https://mad.web.cern.ch/mad/webguide/manual.html
    - https://github.com/MethodicalAcceleratorDesign/MAD-X/blob/master/doc/latexuguide/format.tex
    - https://github.com/MethodicalAcceleratorDesign/MAD-X/blob/master/src/mad_parse.c (and related files)
"""

from __future__ import annotations

import builtins
import copy as copy
import dataclasses
import enum
import math as math
import os as os
import random as random
import typing
import warnings as warnings
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Any

__all__: list[str] = [
    "Any",
    "AttributeExpr",
    "BEAM_PARAMETERS",
    "BEAM_VALUES",
    "Beam",
    "BinaryExpr",
    "CLIGHT",
    "EMASS",
    "Element",
    "Enum",
    "EvaluationContext",
    "Expression",
    "FunctionExpr",
    "Line",
    "ListExpr",
    "MADXExpressionParser",
    "MADXInputError",
    "MADXInputWarning",
    "MADXLexer",
    "MADXParser",
    "MADXParserError",
    "NumberExpr",
    "QELECT",
    "StringExpr",
    "Token",
    "TokenType",
    "UnaryExpr",
    "Variable",
    "VariableExpr",
    "auto",
    "copy",
    "dataclass",
    "field",
    "math",
    "os",
    "random",
    "warnings",
]

class MADXParserError(Exception):
    """
    Base exception for MAD-X parser errors.
    """

class MADXInputError(MADXParserError):
    """
    Error in MAD-X input syntax or semantics.
    """
    def __init__(self, message, line_number=None, context=None, file=None): ...
    def _format_message(self): ...

class MADXInputWarning(UserWarning):
    """
    Warning for non-fatal MAD-X input issues.
    """

class _ReturnStatement(Exception):
    """
    Internal signal that a RETURN statement was encountered in a CALL'd file.
    """

class TokenType(enum.Enum):
    """
    Token types for the MAD-X lexer.
    """

    AND: typing.ClassVar[TokenType]
    ARROW: typing.ClassVar[TokenType]
    CARET: typing.ClassVar[TokenType]
    COLON: typing.ClassVar[TokenType]
    COLON_EQUALS: typing.ClassVar[TokenType]
    COMMA: typing.ClassVar[TokenType]
    EOF: typing.ClassVar[TokenType]
    EQUALS: typing.ClassVar[TokenType]
    GT: typing.ClassVar[TokenType]
    IDENTIFIER: typing.ClassVar[TokenType]
    LBRACE: typing.ClassVar[TokenType]
    LPAREN: typing.ClassVar[TokenType]
    LT: typing.ClassVar[TokenType]
    MINUS: typing.ClassVar[TokenType]
    NEWLINE: typing.ClassVar[TokenType]
    NUMBER: typing.ClassVar[TokenType]
    OR: typing.ClassVar[TokenType]
    PLUS: typing.ClassVar[TokenType]
    RBRACE: typing.ClassVar[TokenType]
    RPAREN: typing.ClassVar[TokenType]
    SEMICOLON: typing.ClassVar[TokenType]
    SLASH: typing.ClassVar[TokenType]
    STAR: typing.ClassVar[TokenType]
    STRING: typing.ClassVar[TokenType]
    _hashable_values_: typing.ClassVar[list] = [
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
    ]
    _member_map_: typing.ClassVar[dict]
    _member_names_: typing.ClassVar[list] = [
        "NUMBER",
        "STRING",
        "IDENTIFIER",
        "PLUS",
        "MINUS",
        "STAR",
        "SLASH",
        "CARET",
        "EQUALS",
        "COLON_EQUALS",
        "COLON",
        "COMMA",
        "SEMICOLON",
        "LPAREN",
        "RPAREN",
        "LBRACE",
        "RBRACE",
        "ARROW",
        "LT",
        "GT",
        "AND",
        "OR",
        "EOF",
        "NEWLINE",
    ]
    _unhashable_values_: typing.ClassVar[list] = list()
    _unhashable_values_map_: typing.ClassVar[dict] = {}
    _use_args_: typing.ClassVar[bool] = False
    _value2member_map_: typing.ClassVar[dict]
    _value_repr_ = None
    _member_type_ = builtins.object
    @staticmethod
    def _generate_next_value_(name, start, count, last_values):
        """
        Generate the next value when not given.

        name: the name of the member
        start: the initial start value or None
        count: the number of existing members
        last_values: the list of values assigned
        """
    @staticmethod
    def _new_member_(type, *args, **kwargs):
        """
        Create and return a new object.  See help(type) for accurate signature.
        """
    @classmethod
    def __new__(cls, value): ...

class Token:
    """
    A token from the MAD-X lexer.
    """

    __dataclass_fields__: typing.ClassVar[dict]
    __dataclass_params__: typing.ClassVar[dataclasses._DataclassParams]
    __hash__: typing.ClassVar[None] = None
    __match_args__: typing.ClassVar[tuple] = ("type", "value", "line", "column")
    def __eq__(self, other): ...
    def __init__(
        self, type: TokenType, value: typing.Any, line: int, column: int
    ) -> None: ...
    def __replace__(self, **changes): ...
    def __repr__(self): ...

class MADXLexer:
    """
    Lexer (tokenizer) for MAD-X input.

    Handles:
    - Single-line comments (! or //)
    - Multi-line comments (/* ... */)
    - Case insensitivity (except for quoted strings)
    - Numbers (integer, real, scientific notation)
    - Identifiers
    - Operators and delimiters
    """
    def __init__(self, text: str, file: str | None = None): ...
    def _advance(self) -> str | None:
        """
        Advance position and return the character.
        """
    def _peek(self, offset: int = 0) -> str | None:
        """
        Look at character at current position + offset.
        """
    def _read_identifier(self) -> Token:
        """
        Read an identifier (converted to lowercase).
        """
    def _read_number(self) -> Token:
        """
        Read a number (integer, real, or scientific notation).
        """
    def _read_string(self) -> Token:
        """
        Read a quoted string.
        """
    def _skip_comment(self) -> bool:
        """
        Skip comments. Returns True if a comment was skipped.
        """
    def _skip_whitespace(self):
        """
        Skip whitespace characters (but not newlines for statement tracking).
        """
    def tokenize(self) -> list[Token]:
        """
        Tokenize the entire input.
        """

class Expression:
    """
    Base class for expressions.
    """

    __dataclass_fields__: typing.ClassVar[dict] = {}
    __dataclass_params__: typing.ClassVar[dataclasses._DataclassParams]
    __hash__: typing.ClassVar[None] = None
    __match_args__: typing.ClassVar[tuple] = tuple()
    def __eq__(self, other): ...
    def __init__(self) -> None: ...
    def __replace__(self, **changes): ...
    def __repr__(self): ...

class NumberExpr(Expression):
    """
    A literal number.
    """

    __dataclass_fields__: typing.ClassVar[dict]
    __dataclass_params__: typing.ClassVar[dataclasses._DataclassParams]
    __hash__: typing.ClassVar[None] = None
    __match_args__: typing.ClassVar[tuple] = "value"
    def __eq__(self, other): ...
    def __init__(self, value: float) -> None: ...
    def __replace__(self, **changes): ...
    def __repr__(self): ...

class StringExpr(Expression):
    """
    A literal string.
    """

    __dataclass_fields__: typing.ClassVar[dict]
    __dataclass_params__: typing.ClassVar[dataclasses._DataclassParams]
    __hash__: typing.ClassVar[None] = None
    __match_args__: typing.ClassVar[tuple] = "value"
    def __eq__(self, other): ...
    def __init__(self, value: str) -> None: ...
    def __replace__(self, **changes): ...
    def __repr__(self): ...

class VariableExpr(Expression):
    """
    A variable reference.
    """

    __dataclass_fields__: typing.ClassVar[dict]
    __dataclass_params__: typing.ClassVar[dataclasses._DataclassParams]
    __hash__: typing.ClassVar[None] = None
    __match_args__: typing.ClassVar[tuple] = "name"
    def __eq__(self, other): ...
    def __init__(self, name: str) -> None: ...
    def __replace__(self, **changes): ...
    def __repr__(self): ...

class AttributeExpr(Expression):
    """
    Element/command attribute access (element->attribute).
    """

    __dataclass_fields__: typing.ClassVar[dict]
    __dataclass_params__: typing.ClassVar[dataclasses._DataclassParams]
    __hash__: typing.ClassVar[None] = None
    __match_args__: typing.ClassVar[tuple] = ("element", "attribute")
    def __eq__(self, other): ...
    def __init__(self, element: str, attribute: str) -> None: ...
    def __replace__(self, **changes): ...
    def __repr__(self): ...

class UnaryExpr(Expression):
    """
    Unary operation (+x, -x).
    """

    __dataclass_fields__: typing.ClassVar[dict]
    __dataclass_params__: typing.ClassVar[dataclasses._DataclassParams]
    __hash__: typing.ClassVar[None] = None
    __match_args__: typing.ClassVar[tuple] = ("operator", "operand")
    def __eq__(self, other): ...
    def __init__(self, operator: str, operand: Expression) -> None: ...
    def __replace__(self, **changes): ...
    def __repr__(self): ...

class BinaryExpr(Expression):
    """
    Binary operation (x op y).
    """

    __dataclass_fields__: typing.ClassVar[dict]
    __dataclass_params__: typing.ClassVar[dataclasses._DataclassParams]
    __hash__: typing.ClassVar[None] = None
    __match_args__: typing.ClassVar[tuple] = ("operator", "left", "right")
    def __eq__(self, other): ...
    def __init__(self, operator: str, left: Expression, right: Expression) -> None: ...
    def __replace__(self, **changes): ...
    def __repr__(self): ...

class FunctionExpr(Expression):
    """
    Function call.
    """

    __dataclass_fields__: typing.ClassVar[dict]
    __dataclass_params__: typing.ClassVar[dataclasses._DataclassParams]
    __hash__: typing.ClassVar[None] = None
    __match_args__: typing.ClassVar[tuple] = ("name", "arguments")
    def __eq__(self, other): ...
    def __init__(self, name: str, arguments: list[Expression]) -> None: ...
    def __replace__(self, **changes): ...
    def __repr__(self): ...

class ListExpr(Expression):
    """
    A list of elements (for LINE definitions).
    """

    __dataclass_fields__: typing.ClassVar[dict]
    __dataclass_params__: typing.ClassVar[dataclasses._DataclassParams]
    __hash__: typing.ClassVar[None] = None
    __match_args__: typing.ClassVar[tuple] = "elements"
    def __eq__(self, other): ...
    def __init__(self, elements: list) -> None: ...
    def __replace__(self, **changes): ...
    def __repr__(self): ...

class MADXExpressionParser:
    """
    Parser for MAD-X expressions.

    Implements operator precedence:
    1. ^ (exponentiation, left associative)
    2. + - (sign, applies to the whole power: -a^b^c is -((a^b)^c))
    3. * / (multiplication, division)
    4. + - (addition, subtraction)

    Beyond MAD-X, which rejects it, the right operand of a binary operator can
    carry its own sign: a*-b, a--b, a^-b. The sign of an exponent binds to that
    exponent only: a^-b^c is (a^(-b))^c.
    """

    FUNCTIONS: typing.ClassVar[dict]
    def __init__(self, tokens: list[Token], file: str | None = None): ...
    def _advance(self) -> Token:
        """
        Advance and return current token.
        """
    def _current(self) -> Token:
        """
        Get current token.
        """
    def _expect(self, *types: TokenType) -> Token:
        """
        Expect one of the given token types.
        """
    def _parse_additive(self) -> Expression:
        """
        Parse additive expression (+ -).
        """
    def _parse_exponent(self) -> Expression:
        """
        Parse the right operand of ^, which can carry a sign (a^-b).
        """
    def _parse_multiplicative(self) -> Expression:
        """
        Parse multiplicative expression (* /).
        """
    def _parse_power(self) -> Expression:
        """
        Parse power expression (^), left associative.
        """
    def _parse_primary(self) -> Expression:
        """
        Parse primary expression (numbers, variables, function calls, parenthesized).
        """
    def _parse_unary(self) -> Expression:
        """
        Parse unary expression (+x, -x), where x is a whole power.
        """
    def _peek(self, offset: int = 0) -> Token:
        """
        Look ahead at token.
        """
    def parse_expression(self) -> Expression:
        """
        Parse an expression (entry point).
        """

class Variable:
    """
    A variable in the symbol table.
    """

    __dataclass_fields__: typing.ClassVar[dict]
    __dataclass_params__: typing.ClassVar[dataclasses._DataclassParams]
    __hash__: typing.ClassVar[None] = None
    __match_args__: typing.ClassVar[tuple] = (
        "name",
        "value",
        "expression",
        "deferred",
        "constant",
    )
    constant: typing.ClassVar[bool] = False
    deferred: typing.ClassVar[bool] = False
    expression = None
    def __eq__(self, other): ...
    def __init__(
        self,
        name: str,
        value: typing.Any,
        expression: Expression | None = None,
        deferred: bool = False,
        constant: bool = False,
    ) -> None: ...
    def __replace__(self, **changes): ...
    def __repr__(self): ...

class Element:
    """
    A MAD-X element definition.
    """

    __dataclass_fields__: typing.ClassVar[dict]
    __dataclass_params__: typing.ClassVar[dataclasses._DataclassParams]
    __hash__: typing.ClassVar[None] = None
    __match_args__: typing.ClassVar[tuple] = (
        "name",
        "type",
        "attributes",
        "attribute_exprs",
    )
    def __eq__(self, other): ...
    def __init__(
        self,
        name: str,
        type: str,
        attributes: dict[str, typing.Any] = ...,
        attribute_exprs: dict[str, Expression] = ...,
    ) -> None: ...
    def __replace__(self, **changes): ...
    def __repr__(self): ...

class Line:
    """
    A MAD-X LINE definition.
    """

    __dataclass_fields__: typing.ClassVar[dict]
    __dataclass_params__: typing.ClassVar[dataclasses._DataclassParams]
    __hash__: typing.ClassVar[None] = None
    __match_args__: typing.ClassVar[tuple] = ("name", "elements")
    def __eq__(self, other): ...
    def __init__(self, name: str, elements: list) -> None: ...
    def __replace__(self, **changes): ...
    def __repr__(self): ...

class Beam:
    """
    A MAD-X beam: the result of its BEAM commands.

    The defaults are those of the MAD-X `beam` command in `src/mad_dict.c`.
    """

    __dataclass_fields__: typing.ClassVar[dict]
    __dataclass_params__: typing.ClassVar[dataclasses._DataclassParams]
    __hash__: typing.ClassVar[None] = None
    __match_args__: typing.ClassVar[tuple] = (
        "defined",
        "particle",
        "mass",
        "charge",
        "energy",
        "energy_is_explicit",
        "freq0",
        "circ",
        "bv",
        "ex",
        "ey",
        "et",
        "sigt",
        "sige",
        "kbunch",
        "npart",
        "bcurrent",
        "bunched",
        "radiate",
        "dtbyds",
        "deltap",
        "u0",
        "qs",
        "n1min",
        "pc",
        "gamma",
        "beta",
        "brho",
        "exn",
        "eyn",
        "alfa",
        "arad",
    )
    alfa: typing.ClassVar[float] = 0.0
    arad: typing.ClassVar[float] = 0.0
    bcurrent: typing.ClassVar[float] = 0.0
    beta: typing.ClassVar[float] = 0.0
    brho: typing.ClassVar[float] = 0.0
    bunched: typing.ClassVar[float] = 1.0
    bv: typing.ClassVar[float] = 1.0
    charge: typing.ClassVar[float] = 1.0
    circ: typing.ClassVar[float] = 0.0
    defined: typing.ClassVar[bool] = False
    deltap: typing.ClassVar[float] = 0.0
    dtbyds: typing.ClassVar[float] = 0.0
    energy: typing.ClassVar[float] = 1.0
    energy_is_explicit: typing.ClassVar[bool] = False
    et: typing.ClassVar[float] = 0.001
    ex: typing.ClassVar[float] = 1.0
    exn: typing.ClassVar[float] = 0.0
    ey: typing.ClassVar[float] = 1.0
    eyn: typing.ClassVar[float] = 0.0
    freq0: typing.ClassVar[float] = 0.0
    gamma: typing.ClassVar[float] = 0.0
    kbunch: typing.ClassVar[float] = 1.0
    mass: typing.ClassVar[float] = 0.00051099895
    n1min: typing.ClassVar[float] = -1.0
    npart: typing.ClassVar[float] = 1.0
    particle: typing.ClassVar[str] = "positron"
    pc: typing.ClassVar[float] = 0.0
    qs: typing.ClassVar[float] = 0.0
    radiate: typing.ClassVar[float] = 0.0
    sige: typing.ClassVar[float] = 0.001
    sigt: typing.ClassVar[float] = 1.0
    u0: typing.ClassVar[float] = 0.0
    def __eq__(self, other): ...
    def __init__(
        self,
        defined: bool = False,
        particle: str = "positron",
        mass: float = 0.00051099895,
        charge: float = 1.0,
        energy: float = 1.0,
        energy_is_explicit: bool = False,
        freq0: float = 0.0,
        circ: float = 0.0,
        bv: float = 1.0,
        ex: float = 1.0,
        ey: float = 1.0,
        et: float = 0.001,
        sigt: float = 1.0,
        sige: float = 0.001,
        kbunch: float = 1.0,
        npart: float = 1.0,
        bcurrent: float = 0.0,
        bunched: float = 1.0,
        radiate: float = 0.0,
        dtbyds: float = 0.0,
        deltap: float = 0.0,
        u0: float = 0.0,
        qs: float = 0.0,
        n1min: float = -1.0,
        pc: float = 0.0,
        gamma: float = 0.0,
        beta: float = 0.0,
        brho: float = 0.0,
        exn: float = 0.0,
        eyn: float = 0.0,
        alfa: float = 0.0,
        arad: float = 0.0,
    ) -> None: ...
    def __post_init__(self): ...
    def __replace__(self, **changes): ...
    def __repr__(self): ...
    def set_energy(self, attr: str, value: float):
        """
        Set the quantities that follow from the mass, the charge and one energy
        parameter (`energy`, `pc`, `gamma`, `beta` or `brho`), which is kept as stated.
        """
    @property
    def is_physical(self) -> bool:
        """
        Whether this can be a reference particle, unlike all that MAD-X accepts.
        """

class EvaluationContext:
    """
    Context for evaluating expressions.

    Holds variables, elements, lines, and provides expression evaluation.
    """

    PARTICLES: typing.ClassVar[dict] = {
        "positron": (0.00051099895, 1.0),
        "electron": (0.00051099895, -1.0),
        "proton": (0.93827208816, 1.0),
        "antiproton": (0.93827208816, -1.0),
        "posmuon": (0.1056583755, 1.0),
        "negmuon": (0.1056583755, -1.0),
        "ion": (0.93956542052, 1.0),
    }
    PREDEFINED_CONSTANTS: typing.ClassVar[dict] = {
        "pi": 3.141592653589793,
        "twopi": 6.283185307179586,
        "degrad": 57.29577951308232,
        "raddeg": 0.017453292519943295,
        "e": 2.718281828459045,
        "amu0": 1.2566370614359173e-06,
        "emass": 0.00051099895,
        "pmass": 0.93827208816,
        "nmass": 0.93956542052,
        "umass": 0.93149410242,
        "mumass": 0.1056583755,
        "clight": 299792458.0,
        "qelect": 1.602176634e-19,
        "hbar": 6.582119569e-25,
        "erad": 2.8179403262e-15,
        "prad": 1.5346982671888944e-18,
        "true": 1.0,
        "false": 0.0,
    }
    def __init__(self): ...
    def _error(self, message: str) -> MADXInputError:
        """
        Create a MADXInputError with current file and line context.
        """
    def _referenced_beam(self, name: str) -> Beam:
        """
        The beam of `beam` or `beam%sequence`.

        `beam` is the default beam, as is the beam of a sequence that has none.
        """
    def _warn(self, message: str):
        """
        Emit a warning with current file and line context.
        """
    def evaluate(self, expr: Expression) -> typing.Any:
        """
        Evaluate an expression.
        """
    def get_beam_attribute(self, name: str, attr_name: str) -> float:
        """
        Value of `beam->attr` or `beam%sequence->attr`.
        """
    def get_element_attribute(self, element_name: str, attr_name: str) -> typing.Any:
        """
        Get an element's attribute value.
        """
    def get_option(self, name: str, default: typing.Any = 0.0) -> typing.Any:
        """
        Get a MAD-X OPTION value.
        """
    def get_variable(self, name: str) -> typing.Any:
        """
        Get a variable value, evaluating deferred expressions.
        """
    def set_beam_attribute(self, name: str, attr_name: str, value: float):
        """
        Set `beam->attr` or `beam%sequence->attr`, as in MAD-X.

        Only this value changes. What follows from it, such as `pc` from `energy`, keeps
        its value until the next BEAM command, so that consumers cannot rely on them.
        """
    def set_option(self, name: str, value: typing.Any):
        """
        Set a MAD-X OPTION value.
        """
    def set_variable(
        self,
        name: str,
        value: typing.Any = None,
        expression: Expression = None,
        deferred: bool = False,
        constant: bool = False,
    ):
        """
        Set a variable value.
        """
    @property
    def beam(self) -> Beam:
        """
        Get the beam for the selected sequence, falling back to default.
        """

class MADXParser:
    """
    MAD-X Parser.

    Parses MAD-X input files and builds a representation of the lattice.
    """
    def __init__(self): ...
    def __str__(self) -> str:
        """
        String representation with lattice information.
        """
    def _advance(self) -> Token:
        """
        Advance and return current token.
        """
    def _apply_beam_parameters(self, beam: Beam, exprs: dict):
        """
        Apply the parameters of a BEAM command that follow from its energy.
        """
    def _check_beam_energy(self, attr: str, value: float, mass: float):
        """
        Stop on an energy parameter where MAD-X stops.
        """
    def _current(self) -> Token:
        """
        Get current token.
        """
    def _error(self, message: str) -> MADXInputError:
        """
        Create a MADXInputError with current file and line context.
        """
    def _execute_or_skip_brace_block(self, execute: bool):
        """
        Execute or skip a { ... } block, then consume trailing semicolon.
        """
    def _expect(self, *types: TokenType) -> Token:
        """
        Expect one of the given token types.
        """
    def _flatten_line(self, line_name: str) -> list[str]:
        """
        Recursively flatten a line definition to element names.
        """
    def _get_element_length(self, elem_name: str) -> float:
        """
        Get an element's length, evaluating deferred expressions if needed.
        """
    def _iter_line_references(self, elem) -> list[str]:
        """
        Yield line/element names referenced inside a LINE definition.
        """
    def _parse_add2expr_command(self):
        """
        Parse ADD2EXPR command.

        MAD-X:
            ADD2EXPR, var=name, expr=<expression>;

        Appends an expression term to an existing deferred expression.
        """
    def _parse_beam_assignment(self):
        """
        Parse `beam->name = value` or `beam%sequence->name = value`.

        As in MAD-X, also `:=` evaluates the value now.
        """
    def _parse_beam_command(self):
        """
        Parse a BEAM command, after its name, and apply it as MAD-X `update_beam`.

        MAD-X allows multiple BEAM commands for different sequences via
        BEAM, SEQUENCE=name, ...;  The beam is stored per-sequence.
        A BEAM without SEQUENCE= updates the default beam.
        The energy comes from the first of ENERGY, PC, GAMMA, BETA, BRHO.

        The values of the command are computed on a copy, which replaces the beam at the
        end: as in MAD-X, `beam->name` in the command reads the beam before the command.
        A new beam starts from the MAD-X defaults.
        """
    def _parse_beam_name(self) -> str:
        """
        Parse a name in a beam command, as in MAD-X.

        A quoted name keeps its case. Otherwise, signs and numbers continue the name,
        as in `H-` or `PB208+82`, except for a sign before another parameter, as in
        `PROTON -RADIATE`.
        """
    def _parse_beam_parameters(
        self,
        allowed: tuple = (
            "mass",
            "charge",
            "energy",
            "pc",
            "gamma",
            "beta",
            "brho",
            "ex",
            "exn",
            "ey",
            "eyn",
            "et",
            "sigt",
            "sige",
            "kbunch",
            "npart",
            "bcurrent",
            "circ",
            "dtbyds",
            "deltap",
            "alfa",
            "u0",
            "qs",
            "arad",
            "bv",
            "bunched",
            "radiate",
            "n1min",
            "freq0",
            "pdamp",
            "particle",
            "sequence",
        ),
    ) -> tuple[str | None, dict, str | None]:
        """
        Parse the parameters of a beam command through its semicolon, as in MAD-X.

        A comma is optional except after a numeric value. A parameter without a value
        takes its default in `src/mad_dict.c`: a flag is set, or cleared by `-NAME`, a
        number is zero, and PARTICLE is a positron. An unknown parameter raises, also
        where MAD-X makes it part of a preceding string value, which turns a typo in
        `PARTICLE=PROTON, ENERGIE=3` into an unknown particle with the electron mass.

        :param allowed: the parameters of the command
        :return: the particle name, the unevaluated values, and the sequence name
        """
    def _parse_beam_value(self, attr: str, expr: Expression) -> float:
        """
        Evaluate a value of a BEAM command, of which only the used ones are read.
        """
    def _parse_call_command(self):
        """
        Parse the CALL command.

        Supports both MAD-X and MAD-8 syntax:
          CALL, FILE="other_file";
          CALL, FILENAME="other_file";

        The called file is read and parsed. If the called file contains a
        RETURN statement, parsing stops at that point. Called files may
        themselves contain CALL statements to any depth.
        """
    def _parse_comparison(self) -> bool:
        """
        Parse a single comparison: arith_expr [op arith_expr].
        """
    def _parse_condition_operand(self) -> Expression:
        """
        Parse an arithmetic expression inside a condition.

        Stops at comparison operators, logical operators, or closing paren.
        """
    def _parse_const_declaration(self):
        """
        Parse a CONST declaration.
        """
    def _parse_element_attribute_update(self):
        """
        Parse top-level element attribute update statements.

        Example:
            mqxa.1r1, k1 := kqx.r1+ktqx1.r1, polarity=+1;
        """
    def _parse_element_definition(self, name: str, element_type: str):
        """
        Parse an element definition.
        """
    def _parse_exec_command(self):
        """
        Parse EXEC, MacroName($arg1, $arg2, ...);

        Substitutes $argN with the corresponding variable value,
        then parses the expanded macro body.
        """
    def _parse_expression(self, complete: bool = False) -> Expression:
        """
        Parse an expression using the expression parser.

        :param complete: raise if tokens follow the expression before the next comma
        """
    def _parse_if_condition(self) -> bool:
        """
        Parse a condition: ( logical_expr ).

        Supports compound conditions with && and ||, e.g.:
            IF ( a > 0 && b < 1 )
            IF ( x == 1 || y == 2 )

        MAD-X comparison operators (from mad_eval.c):
            ==, <>, <, >, <=, >=
        """
    def _parse_if_statement(self):
        """
        Parse a MAD-X IF statement.

        Syntax:
            IF ( condition ) { statements; };
            IF ( condition ) { statements; } ELSEIF ( condition ) { statements; } ELSE { statements; };
        """
    def _parse_labeled_definition(self):
        """
        Parse a labeled definition (element or line) or a labeled beam command.
        """
    def _parse_line_definition(self, name: str):
        """
        Parse a LINE definition.
        """
    def _parse_line_element(self):
        """
        Parse a single line element (possibly with multiplier or reversal).
        """
    def _parse_line_elements(self) -> list:
        """
        Parse the elements inside a LINE definition.
        """
    def _parse_logical_and(self) -> bool:
        """
        Parse logical AND: expr && expr && ...
        """
    def _parse_logical_or(self) -> bool:
        """
        Parse logical OR: expr || expr || ...
        """
    def _parse_macro_definition(self):
        """
        Parse a macro definition: name(params): MACRO = { body };
        """
    def _parse_option_command(self):
        """
        Parse the OPTION command and persist option values.

        Supports forms such as:
            OPTION, RBARC=true, THIN_FOC=false;
            OPTION, ECHO;
            OPTION, -WARN;
        """
    def _parse_resbeam_command(self):
        """
        Parse a RESBEAM command, after its name, which restores the beam defaults.
        """
    def _parse_sequence_definition(self, name: str):
        """
        Parse a SEQUENCE definition.

        Syntax:
            name: SEQUENCE, L=length, REFER=centre;
              elem1: type, at=pos;   ! labeled placement
              elem2, at=pos;         ! unlabeled placement (existing element)
            ENDSEQUENCE;

        The REFER attribute controls how 'at' positions are interpreted:
            centre (default) - 'at' is the element center
            entry            - 'at' is the element entry point
            exit             - 'at' is the element exit point

        Elements are collected, positions adjusted for REFER, sorted,
        and stored as a Line for compatibility with getBeamline().
        """
    def _parse_statement(self):
        """
        Parse a single statement.
        """
    def _parse_text(self, text: str):
        """
        Internal method to parse text.
        """
    def _parse_use_command(self):
        """
        Parse the USE command.
        """
    def _parse_variable_assignment(self, integer: bool = False):
        """
        Parse a variable assignment.
        """
    def _parse_while_statement(self):
        """
        Parse a MAD-X WHILE statement.

        Syntax:
            WHILE ( condition ) { statements; }
        """
    def _peek(self, offset: int = 0) -> Token:
        """
        Look ahead at token.
        """
    def _process_line_element(self, elem) -> list[str]:
        """
        Process a single line element specification.
        """
    def _resolve_beam_particle(
        self, beam: Beam, particle: str | None, exprs: dict
    ) -> tuple[str, float, float]:
        """
        Particle name, rest energy [GeV] and charge [e] after a BEAM command.
        """
    def _resolve_beamline_name(
        self, *, line_name: str | None = None, sequence_name: str | None = None
    ) -> str | None:
        """
        Resolve which LINE/SEQUENCE should be expanded into a beamline.
        """
    def _root_lines(self) -> list[str]:
        """
        Return top-level LINE/SEQUENCE definitions not referenced by other lines.
        """
    def _skip_semicolons(self):
        """
        Skip any semicolons.
        """
    def _skip_until_semicolon(self):
        """
        Skip tokens until semicolon or EOF.
        """
    def _try_parse_comparison_op(self) -> str | None:
        """
        Try to consume a comparison operator. Returns operator string or None.
        """
    def _update_context_location(self):
        """
        Update file/line on the evaluation context for warning messages.
        """
    def _warn(self, message: str):
        """
        Emit a warning with current file and line context.
        """
    def getBeamline(
        self, *, line_name: str | None = None, sequence_name: str | None = None
    ) -> list[dict]:
        """
        Get the beamline as a list of element dictionaries.

        Returns a list compatible with the old parser format.
        """
    def getElement(self, name: str) -> Element | None:
        """
        Get an element by name.
        """
    def getEtot(self) -> float:
        """
        Get total energy in GeV.
        """
    def getFreq0(self) -> float:
        """
        Get revolution frequency in Hz.
        """
    def getLine(self, name: str) -> Line | None:
        """
        Get a line by name.
        """
    def getOption(self, name: str, default: typing.Any = 0.0) -> typing.Any:
        """
        Get an OPTION value from the parsed MAD-X input.
        """
    def getOptions(self) -> dict[str, typing.Any]:
        """
        Get all parsed OPTION values.
        """
    def getParticle(self) -> str:
        """
        Get the particle type.
        """
    def getVariable(self, name: str) -> typing.Any:
        """
        Get a variable value.
        """
    def parse(self, fn: str):
        """
        Parse a MAD-X file.

        Args:
            fn: Filename to parse
        """
    def parse_string(self, text: str):
        """
        Parse MAD-X input from a string.
        """
    @property
    def beam(self) -> dict:
        """
        Get beam parameters as a dictionary (backward compatible).
        """
    @property
    def sequence(self) -> dict:
        """
        Get selected sequence (backward compatible).
        """

def _c_divide(numerator: float, denominator: float) -> float:
    """
    `numerator / denominator` in C.
    """

def _c_math_result(name: str, args: list, error: Exception) -> float:
    """
    What the C function `name` returns for the arguments a Python one raises on.
    """

def _c_pow(base: float, exponent: float) -> float:
    """
    `pow(base, exponent)` in C.
    """

def _c_sqrt(value: float) -> float:
    """
    `sqrt(value)` in C.
    """

def _madx_formatwarning(message, category, filename, lineno, line=None): ...

BEAM_PARAMETERS: tuple = (
    "mass",
    "charge",
    "energy",
    "pc",
    "gamma",
    "beta",
    "brho",
    "ex",
    "exn",
    "ey",
    "eyn",
    "et",
    "sigt",
    "sige",
    "kbunch",
    "npart",
    "bcurrent",
    "circ",
    "dtbyds",
    "deltap",
    "alfa",
    "u0",
    "qs",
    "arad",
    "bv",
    "bunched",
    "radiate",
    "n1min",
    "freq0",
    "pdamp",
    "particle",
    "sequence",
)
BEAM_VALUES: tuple = (
    "mass",
    "charge",
    "energy",
    "pc",
    "gamma",
    "beta",
    "brho",
    "ex",
    "exn",
    "ey",
    "eyn",
    "et",
    "sigt",
    "sige",
    "kbunch",
    "npart",
    "bcurrent",
    "circ",
    "dtbyds",
    "deltap",
    "alfa",
    "u0",
    "qs",
    "arad",
    "bv",
    "bunched",
    "radiate",
    "n1min",
)
CLIGHT: float = 299792458.0
EMASS: float = 0.00051099895
QELECT: float = 1.602176634e-19
