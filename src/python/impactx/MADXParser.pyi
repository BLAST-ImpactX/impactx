"""
A MAD-X parser that aims to be compliant with the upstream MAD-X format.

References:
    - https://mad.web.cern.ch/mad/webguide/manual.html
    - https://github.com/MethodicalAcceleratorDesign/MAD-X/blob/master/doc/latexuguide/format.tex
    - https://github.com/MethodicalAcceleratorDesign/MAD-X/blob/master/src/mad_parse.c (and related files)
"""

from __future__ import annotations

import builtins
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
    "Beam",
    "BinaryExpr",
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
    "StringExpr",
    "Token",
    "TokenType",
    "UnaryExpr",
    "Variable",
    "VariableExpr",
    "auto",
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
    1. ^ (exponentiation, right associative)
    2. * / (multiplication, division)
    3. + - (addition, subtraction)
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
    def _parse_multiplicative(self) -> Expression:
        """
        Parse multiplicative expression (* /).
        """
    def _parse_power(self) -> Expression:
        """
        Parse power expression (^), right associative.
        """
    def _parse_primary(self) -> Expression:
        """
        Parse primary expression (numbers, variables, function calls, parenthesized).
        """
    def _parse_unary(self) -> Expression:
        """
        Parse unary expression (+x, -x).
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
    MAD-X BEAM command parameters.
    """

    __dataclass_fields__: typing.ClassVar[dict]
    __dataclass_params__: typing.ClassVar[dataclasses._DataclassParams]
    __hash__: typing.ClassVar[None] = None
    __match_args__: typing.ClassVar[tuple] = (
        "particle",
        "energy",
        "pc",
        "mass",
        "charge",
        "freq0",
        "bv",
    )
    bv: typing.ClassVar[float] = 1.0
    charge: typing.ClassVar[float] = 0.0
    energy: typing.ClassVar[float] = 1.0
    freq0: typing.ClassVar[float] = 0.0
    mass: typing.ClassVar[float] = 0.0
    particle: typing.ClassVar[str] = ""
    pc: typing.ClassVar[float] = 0.0
    def __eq__(self, other): ...
    def __init__(
        self,
        particle: str = "",
        energy: float = 1.0,
        pc: float = 0.0,
        mass: float = 0.0,
        charge: float = 0.0,
        freq0: float = 0.0,
        bv: float = 1.0,
    ) -> None: ...
    def __replace__(self, **changes): ...
    def __repr__(self): ...

class EvaluationContext:
    """
    Context for evaluating expressions.

    Holds variables, elements, lines, and provides expression evaluation.
    """

    PREDEFINED_CONSTANTS: typing.ClassVar[dict] = {
        "pi": 3.141592653589793,
        "twopi": 6.283185307179586,
        "degrad": 57.29577951308232,
        "raddeg": 0.017453292519943295,
        "e": 2.718281828459045,
        "emass": 0.00051099895,
        "pmass": 0.93827208816,
        "nmass": 0.93956542052,
        "umass": 0.93149410242,
        "mumass": 0.1056583715,
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
    def _warn(self, message: str):
        """
        Emit a warning with current file and line context.
        """
    def evaluate(self, expr: Expression) -> typing.Any:
        """
        Evaluate an expression.
        """
    def get_beam_for_sequence(self, sequence_name: str | None = None) -> Beam:
        """
        Get or create a beam for a specific sequence.
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
    def _parse_beam_command(self):
        """
        Parse the BEAM command.

        MAD-X allows multiple BEAM commands for different sequences via
        BEAM, SEQUENCE=name, ...;  The beam is stored per-sequence.
        A BEAM without SEQUENCE= updates the default beam.
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
    def _parse_expression(self) -> Expression:
        """
        Parse an expression using the expression parser.
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
        Parse a labeled definition (element or line).
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

def _madx_formatwarning(message, category, filename, lineno, line=None): ...
