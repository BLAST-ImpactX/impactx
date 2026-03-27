#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Matthias Frey, Andreas Adelmann, Marco Garten, Axel Huebl
# License: BSD-3-Clause
#
# -*- coding: utf-8 -*-
"""
A MAD-X parser that aims to be compliant with the upstream MAD-X format.

References:
    - https://mad.web.cern.ch/mad/webguide/manual.html
    - https://github.com/MethodicalAcceleratorDesign/MAD-X/blob/master/doc/latexuguide/format.tex
    - https://github.com/MethodicalAcceleratorDesign/MAD-X/blob/master/src/mad_parse.c (and related files)
"""

import math
import os
import random
import warnings
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Any, Optional


class MADXParserError(Exception):
    """Base exception for MAD-X parser errors."""

    pass


class MADXInputError(MADXParserError):
    """Error in MAD-X input syntax or semantics."""

    def __init__(self, message, line_number=None, context=None, file=None):
        self.message = message
        self.line_number = line_number
        self.context = context
        self.file = file
        super().__init__(self._format_message())

    def _format_message(self):
        msg = self.message
        loc = ""
        if self.file:
            loc = self.file
        if self.line_number is not None:
            loc = f"{loc}:{self.line_number}" if loc else f"Line {self.line_number}"
        if loc:
            msg = f"{loc}: {msg}"
        if self.context:
            msg = f"{msg}\n  Context: {self.context}"
        return msg


class MADXInputWarning(UserWarning):
    """Warning for non-fatal MAD-X input issues."""

    pass


# Custom format for MADXInputWarning: omit the source line that Python prints by default
_original_formatwarning = warnings.formatwarning


def _madx_formatwarning(message, category, filename, lineno, line=None):
    if issubclass(category, MADXInputWarning):
        return f"[WARNING MADXParser] {str(message)}\n"
    return _original_formatwarning(message, category, filename, lineno, line)


warnings.formatwarning = _madx_formatwarning


class _ReturnStatement(Exception):
    """Internal signal that a RETURN statement was encountered in a CALL'd file."""

    pass


# =============================================================================
# Token Types
# =============================================================================


class TokenType(Enum):
    """Token types for the MAD-X lexer."""

    # Literals
    NUMBER = auto()
    STRING = auto()
    IDENTIFIER = auto()

    # Operators
    PLUS = auto()
    MINUS = auto()
    STAR = auto()
    SLASH = auto()
    CARET = auto()
    EQUALS = auto()
    COLON_EQUALS = auto()
    COLON = auto()
    COMMA = auto()
    SEMICOLON = auto()
    LPAREN = auto()
    RPAREN = auto()
    LBRACE = auto()
    RBRACE = auto()
    ARROW = auto()  # ->
    LT = auto()  # <
    GT = auto()  # >
    AND = auto()  # &&
    OR = auto()  # ||

    # Special
    EOF = auto()
    NEWLINE = auto()


@dataclass
class Token:
    """A token from the MAD-X lexer."""

    type: TokenType
    value: Any
    line: int
    column: int

    def __repr__(self):
        return f"Token({self.type.name}, {self.value!r}, line={self.line})"


# =============================================================================
# Lexer
# =============================================================================


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

    def __init__(self, text: str, file: Optional[str] = None):
        self.text = text
        self.pos = 0
        self.line = 1
        self.column = 1
        self.length = len(text)
        self.file = file

    def _peek(self, offset: int = 0) -> Optional[str]:
        """Look at character at current position + offset."""
        pos = self.pos + offset
        if pos < self.length:
            return self.text[pos]
        return None

    def _advance(self) -> Optional[str]:
        """Advance position and return the character."""
        if self.pos >= self.length:
            return None
        char = self.text[self.pos]
        self.pos += 1
        if char == "\n":
            self.line += 1
            self.column = 1
        else:
            self.column += 1
        return char

    def _skip_whitespace(self):
        """Skip whitespace characters (but not newlines for statement tracking)."""
        while self.pos < self.length and self.text[self.pos] in " \t\r\n":
            self._advance()

    def _skip_comment(self) -> bool:
        """Skip comments. Returns True if a comment was skipped."""
        # Single-line comment with !
        if self._peek() == "!":
            while self._peek() is not None and self._peek() != "\n":
                self._advance()
            return True

        # Single-line comment with //
        if self._peek() == "/" and self._peek(1) == "/":
            while self._peek() is not None and self._peek() != "\n":
                self._advance()
            return True

        # Multi-line comment /* ... */
        if self._peek() == "/" and self._peek(1) == "*":
            self._advance()  # /
            self._advance()  # *
            while self.pos < self.length:
                if self._peek() == "*" and self._peek(1) == "/":
                    self._advance()  # *
                    self._advance()  # /
                    return True
                self._advance()
            raise MADXInputError(
                "Unterminated multi-line comment", self.line, file=self.file
            )

        return False

    def _read_number(self) -> Token:
        """Read a number (integer, real, or scientific notation)."""
        start_line = self.line
        start_col = self.column
        chars = []

        # Optional leading sign is handled by the parser, not lexer
        # Read integer part
        while self._peek() is not None and self._peek().isdigit():
            chars.append(self._advance())

        # Decimal point
        if self._peek() == ".":
            chars.append(self._advance())
            while self._peek() is not None and self._peek().isdigit():
                chars.append(self._advance())

        # Exponent (E, e, D, d)
        if self._peek() is not None and self._peek().upper() in "ED":
            chars.append(self._advance())
            if self._peek() in "+-":
                chars.append(self._advance())
            while self._peek() is not None and self._peek().isdigit():
                chars.append(self._advance())

        value_str = "".join(chars).upper().replace("D", "E")
        try:
            value = float(value_str)
        except ValueError:
            raise MADXInputError(
                f"Invalid number: {value_str}", start_line, file=self.file
            )

        return Token(TokenType.NUMBER, value, start_line, start_col)

    def _read_string(self) -> Token:
        """Read a quoted string."""
        start_line = self.line
        start_col = self.column
        quote = self._advance()  # ' or "
        chars = []

        while self._peek() is not None and self._peek() != quote:
            if self._peek() == "\n":
                raise MADXInputError("Unterminated string", start_line, file=self.file)
            chars.append(self._advance())

        if self._peek() is None:
            raise MADXInputError("Unterminated string", start_line, file=self.file)

        self._advance()  # closing quote
        return Token(TokenType.STRING, "".join(chars), start_line, start_col)

    def _read_identifier(self) -> Token:
        """Read an identifier (converted to lowercase)."""
        start_line = self.line
        start_col = self.column
        chars = []

        # First character must be a letter, underscore, or $ (macro argument)
        if self._peek() is not None and (
            self._peek().isalpha() or self._peek() in "_$"
        ):
            chars.append(self._advance())

        # Subsequent characters: letters, digits, underscore, period
        while self._peek() is not None and (
            self._peek().isalnum() or self._peek() in "_."
        ):
            chars.append(self._advance())

        # MAD-X is case-insensitive for identifiers
        value = "".join(chars).lower()
        return Token(TokenType.IDENTIFIER, value, start_line, start_col)

    def tokenize(self) -> list[Token]:
        """Tokenize the entire input."""
        tokens = []

        while self.pos < self.length:
            # Skip whitespace and comments
            self._skip_whitespace()
            while self._skip_comment():
                self._skip_whitespace()

            if self.pos >= self.length:
                break

            char = self._peek()
            start_line = self.line
            start_col = self.column

            # Number
            if char.isdigit() or (
                char == "." and self._peek(1) and self._peek(1).isdigit()
            ):
                tokens.append(self._read_number())

            # String
            elif char in "\"'":
                tokens.append(self._read_string())

            # Identifier ($ prefix for macro arguments)
            elif char.isalpha() or char in "_$":
                tokens.append(self._read_identifier())

            # Two-character operators
            elif char == ":" and self._peek(1) == "=":
                self._advance()
                self._advance()
                tokens.append(
                    Token(TokenType.COLON_EQUALS, ":=", start_line, start_col)
                )

            elif char == "-" and self._peek(1) == ">":
                self._advance()
                self._advance()
                tokens.append(Token(TokenType.ARROW, "->", start_line, start_col))

            # Single-character operators
            elif char == "+":
                self._advance()
                tokens.append(Token(TokenType.PLUS, "+", start_line, start_col))

            elif char == "-":
                self._advance()
                tokens.append(Token(TokenType.MINUS, "-", start_line, start_col))

            elif char == "*":
                self._advance()
                tokens.append(Token(TokenType.STAR, "*", start_line, start_col))

            elif char == "/":
                self._advance()
                tokens.append(Token(TokenType.SLASH, "/", start_line, start_col))

            elif char == "^":
                self._advance()
                tokens.append(Token(TokenType.CARET, "^", start_line, start_col))

            elif char == "=":
                self._advance()
                tokens.append(Token(TokenType.EQUALS, "=", start_line, start_col))

            elif char == ":":
                self._advance()
                tokens.append(Token(TokenType.COLON, ":", start_line, start_col))

            elif char == ",":
                self._advance()
                tokens.append(Token(TokenType.COMMA, ",", start_line, start_col))

            elif char == ";":
                self._advance()
                tokens.append(Token(TokenType.SEMICOLON, ";", start_line, start_col))

            elif char == "(":
                self._advance()
                tokens.append(Token(TokenType.LPAREN, "(", start_line, start_col))

            elif char == ")":
                self._advance()
                tokens.append(Token(TokenType.RPAREN, ")", start_line, start_col))

            elif char == "{":
                self._advance()
                tokens.append(Token(TokenType.LBRACE, "{", start_line, start_col))

            elif char == "}":
                self._advance()
                tokens.append(Token(TokenType.RBRACE, "}", start_line, start_col))

            elif char == "<":
                self._advance()
                tokens.append(Token(TokenType.LT, "<", start_line, start_col))

            elif char == ">":
                self._advance()
                tokens.append(Token(TokenType.GT, ">", start_line, start_col))

            elif char == "&" and self._peek(1) == "&":
                self._advance()
                self._advance()
                tokens.append(Token(TokenType.AND, "&&", start_line, start_col))

            elif char == "|" and self._peek(1) == "|":
                self._advance()
                self._advance()
                tokens.append(Token(TokenType.OR, "||", start_line, start_col))

            else:
                raise MADXInputError(
                    f"Unexpected character: {char!r}", start_line, file=self.file
                )

        tokens.append(Token(TokenType.EOF, None, self.line, self.column))
        return tokens


# =============================================================================
# Expression AST and Evaluator
# =============================================================================


@dataclass
class Expression:
    """Base class for expressions."""

    pass


@dataclass
class NumberExpr(Expression):
    """A literal number."""

    value: float


@dataclass
class StringExpr(Expression):
    """A literal string."""

    value: str


@dataclass
class VariableExpr(Expression):
    """A variable reference."""

    name: str


@dataclass
class AttributeExpr(Expression):
    """Element/command attribute access (element->attribute)."""

    element: str
    attribute: str


@dataclass
class UnaryExpr(Expression):
    """Unary operation (+x, -x)."""

    operator: str
    operand: Expression


@dataclass
class BinaryExpr(Expression):
    """Binary operation (x op y)."""

    operator: str
    left: Expression
    right: Expression


@dataclass
class FunctionExpr(Expression):
    """Function call."""

    name: str
    arguments: list[Expression]


@dataclass
class ListExpr(Expression):
    """A list of elements (for LINE definitions)."""

    elements: list


# =============================================================================
# Parser
# =============================================================================


class MADXExpressionParser:
    """
    Parser for MAD-X expressions.

    Implements operator precedence:
    1. ^ (exponentiation, right associative)
    2. * / (multiplication, division)
    3. + - (addition, subtraction)
    """

    # MAD-X built-in functions
    FUNCTIONS = {
        "sqrt": (1, math.sqrt),
        "log": (1, math.log),
        "log10": (1, math.log10),
        "exp": (1, math.exp),
        "sin": (1, math.sin),
        "cos": (1, math.cos),
        "tan": (1, math.tan),
        "asin": (1, math.asin),
        "acos": (1, math.acos),
        "atan": (1, math.atan),
        "sinh": (1, math.sinh),
        "cosh": (1, math.cosh),
        "tanh": (1, math.tanh),
        "abs": (1, abs),
        "floor": (1, math.floor),
        "ceil": (1, math.ceil),
        "round": (1, round),
        "sinc": (1, lambda x: 1.0 if x == 0 else math.sin(x) / x),
        "erf": (1, math.erf),
        "erfc": (1, math.erfc),
        "frac": (1, lambda x: x - math.floor(x)),
        "atan2": (2, math.atan2),
        "max": (2, max),
        "min": (2, min),
        "mod": (2, math.fmod),
        "ranf": (0, lambda: random.random()),
        "gauss": (0, lambda: random.gauss(0, 1)),
        "tgauss": (
            1,
            lambda cut: next(
                s for s in iter(lambda: random.gauss(0, 1), None) if abs(s) < cut
            ),
        ),
    }

    def __init__(self, tokens: list[Token], file: Optional[str] = None):
        self.tokens = tokens
        self.pos = 0
        self.file = file

    def _current(self) -> Token:
        """Get current token."""
        if self.pos < len(self.tokens):
            return self.tokens[self.pos]
        return self.tokens[-1]  # EOF

    def _peek(self, offset: int = 0) -> Token:
        """Look ahead at token."""
        pos = self.pos + offset
        if pos < len(self.tokens):
            return self.tokens[pos]
        return self.tokens[-1]

    def _advance(self) -> Token:
        """Advance and return current token."""
        token = self._current()
        if self.pos < len(self.tokens) - 1:
            self.pos += 1
        return token

    def _expect(self, *types: TokenType) -> Token:
        """Expect one of the given token types."""
        token = self._current()
        if token.type not in types:
            expected = " or ".join(t.name for t in types)
            raise MADXInputError(
                f"Expected {expected}, got {token.type.name}",
                token.line,
                file=self.file,
            )
        return self._advance()

    def parse_expression(self) -> Expression:
        """Parse an expression (entry point)."""
        return self._parse_additive()

    def _parse_additive(self) -> Expression:
        """Parse additive expression (+ -)."""
        left = self._parse_multiplicative()

        while self._current().type in (TokenType.PLUS, TokenType.MINUS):
            op = self._advance().value
            right = self._parse_multiplicative()
            left = BinaryExpr(op, left, right)

        return left

    def _parse_multiplicative(self) -> Expression:
        """Parse multiplicative expression (* /)."""
        left = self._parse_power()

        while self._current().type in (TokenType.STAR, TokenType.SLASH):
            op = self._advance().value
            right = self._parse_power()
            left = BinaryExpr(op, left, right)

        return left

    def _parse_power(self) -> Expression:
        """Parse power expression (^), right associative."""
        left = self._parse_unary()

        if self._current().type == TokenType.CARET:
            self._advance()
            right = self._parse_power()  # Right associative
            left = BinaryExpr("^", left, right)

        return left

    def _parse_unary(self) -> Expression:
        """Parse unary expression (+x, -x)."""
        if self._current().type in (TokenType.PLUS, TokenType.MINUS):
            op = self._advance().value
            operand = self._parse_unary()
            return UnaryExpr(op, operand)

        return self._parse_primary()

    def _parse_primary(self) -> Expression:
        """Parse primary expression (numbers, variables, function calls, parenthesized)."""
        token = self._current()

        # Number literal
        if token.type == TokenType.NUMBER:
            self._advance()
            return NumberExpr(token.value)

        # String literal
        if token.type == TokenType.STRING:
            self._advance()
            return StringExpr(token.value)

        # Parenthesized expression or line element list
        if token.type == TokenType.LPAREN:
            self._advance()
            # Check if this is a list (for LINE definitions)
            expr = self._parse_additive()
            if self._current().type == TokenType.COMMA:
                # This is a list
                elements = [expr]
                while self._current().type == TokenType.COMMA:
                    self._advance()
                    elements.append(self._parse_additive())
                self._expect(TokenType.RPAREN)
                return ListExpr(elements)
            self._expect(TokenType.RPAREN)
            return expr

        # Identifier (variable, function call, or attribute access)
        if token.type == TokenType.IDENTIFIER:
            name = self._advance().value

            # Function call
            if self._current().type == TokenType.LPAREN:
                self._advance()
                args = []
                if self._current().type != TokenType.RPAREN:
                    args.append(self.parse_expression())
                    while self._current().type == TokenType.COMMA:
                        self._advance()
                        args.append(self.parse_expression())
                self._expect(TokenType.RPAREN)
                return FunctionExpr(name, args)

            # Attribute access (element->attribute)
            if self._current().type == TokenType.ARROW:
                self._advance()
                attr = self._expect(TokenType.IDENTIFIER).value
                return AttributeExpr(name, attr)

            # Simple variable
            return VariableExpr(name)

        # Array literal: {expr, expr, ...}
        if token.type == TokenType.LBRACE:
            self._advance()
            elements = []
            if self._current().type != TokenType.RBRACE:
                elements.append(self._parse_additive())
                while self._current().type == TokenType.COMMA:
                    self._advance()
                    elements.append(self._parse_additive())
            self._expect(TokenType.RBRACE)
            return ListExpr(elements)

        raise MADXInputError(
            f"Unexpected token in expression: {token.type.name}",
            token.line,
            file=self.file,
        )


# =============================================================================
# Symbol Table and Evaluation Context
# =============================================================================


@dataclass
class Variable:
    """A variable in the symbol table."""

    name: str
    value: Any
    expression: Optional[Expression] = None  # For deferred evaluation
    deferred: bool = False
    constant: bool = False


@dataclass
class Element:
    """A MAD-X element definition."""

    name: str
    type: str
    attributes: dict[str, Any] = field(default_factory=dict)
    attribute_exprs: dict[str, Expression] = field(default_factory=dict)


@dataclass
class Line:
    """A MAD-X LINE definition."""

    name: str
    elements: list  # List of element names or (multiplier, name) tuples


@dataclass
class Beam:
    """MAD-X BEAM command parameters."""

    particle: str = ""
    energy: float = 1.0  # GeV, MAD-X default
    pc: float = 0.0  # momentum
    mass: float = 0.0
    charge: float = 0.0
    freq0: float = 0.0  # revolution frequency [MHz]


class EvaluationContext:
    """
    Context for evaluating expressions.

    Holds variables, elements, lines, and provides expression evaluation.
    """

    # MAD-X predefined constants (from PDG)
    PREDEFINED_CONSTANTS = {
        "pi": math.pi,
        "twopi": 2 * math.pi,
        "degrad": 180.0 / math.pi,
        "raddeg": math.pi / 180.0,
        "e": math.e,
        "emass": 0.51099895000e-3,  # GeV
        "pmass": 0.93827208816,  # GeV
        "nmass": 0.93956542052,  # GeV
        "umass": 0.93149410242,  # GeV
        "mumass": 0.1056583715,  # GeV
        "clight": 2.99792458e8,  # m/s
        "qelect": 1.602176634e-19,  # A.s
        "hbar": 6.582119569e-25,  # MeV.s
        "erad": 2.8179403262e-15,  # m
        "prad": 2.8179403262e-15 * 0.51099895000e-3 / 0.93827208816,  # m
        # Compatibility aliases
        "true": 1.0,
        "false": 0.0,
    }

    def __init__(self):
        self.variables: dict[str, Variable] = {}
        self.elements: dict[str, Element] = {}
        self.lines: dict[str, Line] = {}
        # MAD-X global OPTION state (subset with defaults from user guide).
        # Stored as numeric values for compatibility with MAD-X true/false handling.
        self.options: dict[str, Any] = {
            "rbarc": 1.0,
            "thin_foc": 1.0,
        }
        self.macros: dict[
            str, tuple[list[str], str]
        ] = {}  # name -> (params, body_text)
        self.beams: dict[Optional[str], Beam] = {None: Beam()}  # key=sequence name
        self.selected_sequence: Optional[str] = None

        # Tracking for warning messages (updated by the parser)
        self.current_file: Optional[str] = None
        self.current_line: Optional[int] = None

        # Initialize predefined constants
        for name, value in self.PREDEFINED_CONSTANTS.items():
            self.variables[name] = Variable(name, value, constant=True)

    @property
    def beam(self) -> Beam:
        """Get the beam for the selected sequence, falling back to default."""
        if self.selected_sequence and self.selected_sequence in self.beams:
            return self.beams[self.selected_sequence]
        return self.beams[None]

    def get_beam_for_sequence(self, sequence_name: Optional[str] = None) -> Beam:
        """Get or create a beam for a specific sequence."""
        if sequence_name not in self.beams:
            # New sequence beam inherits from the default beam
            default = self.beams[None]
            self.beams[sequence_name] = Beam(
                particle=default.particle,
                energy=default.energy,
                pc=default.pc,
                mass=default.mass,
                charge=default.charge,
                freq0=default.freq0,
            )
        return self.beams[sequence_name]

    def _warn(self, message: str):
        """Emit a warning with current file and line context."""
        loc = ""
        if self.current_file:
            loc = self.current_file
            if self.current_line is not None:
                loc += f":{self.current_line}"
            loc += ": "
        warnings.warn(f"{loc}{message}", MADXInputWarning)

    def _error(self, message: str) -> MADXInputError:
        """Create a MADXInputError with current file and line context."""
        return MADXInputError(
            message, line_number=self.current_line, file=self.current_file
        )

    def set_variable(
        self,
        name: str,
        value: Any = None,
        expression: Expression = None,
        deferred: bool = False,
        constant: bool = False,
    ):
        """Set a variable value."""
        name = name.lower()

        # Check if trying to modify a constant
        if name in self.variables and self.variables[name].constant:
            if constant:
                # Allow redefining constants
                pass
            else:
                raise self._error(f"Cannot modify constant: {name}")

        if deferred and expression is not None:
            self.variables[name] = Variable(
                name, None, expression, deferred=True, constant=constant
            )
        else:
            if expression is not None:
                value = self.evaluate(expression)
            self.variables[name] = Variable(
                name, value, expression, deferred=False, constant=constant
            )

    def get_variable(self, name: str) -> Any:
        """Get a variable value, evaluating deferred expressions."""
        name = name.lower()
        if name not in self.variables:
            self._warn(f"Undefined variable '{name}', using 0.0")
            return 0.0

        var = self.variables[name]
        if var.deferred and var.expression is not None:
            return self.evaluate(var.expression)
        return var.value

    def get_element_attribute(self, element_name: str, attr_name: str) -> Any:
        """Get an element's attribute value."""
        element_name = element_name.lower()
        attr_name = attr_name.lower()

        if element_name not in self.elements:
            self._warn(
                f"Unknown element '{element_name}', using 0.0 for "
                f"'{element_name}->{attr_name}'"
            )
            return 0.0

        element = self.elements[element_name]

        # Check if attribute has a deferred expression
        if attr_name in element.attribute_exprs:
            return self.evaluate(element.attribute_exprs[attr_name])

        if attr_name in element.attributes:
            return element.attributes[attr_name]

        self._warn(
            f"Element '{element_name}' has no attribute '{attr_name}', using 0.0"
        )
        return 0.0

    def set_option(self, name: str, value: Any):
        """Set a MAD-X OPTION value."""
        self.options[name.lower()] = value

    def get_option(self, name: str, default: Any = 0.0) -> Any:
        """Get a MAD-X OPTION value."""
        return self.options.get(name.lower(), default)

    def evaluate(self, expr: Expression) -> Any:
        """Evaluate an expression."""
        if isinstance(expr, NumberExpr):
            return expr.value

        if isinstance(expr, StringExpr):
            return expr.value

        if isinstance(expr, VariableExpr):
            return self.get_variable(expr.name)

        if isinstance(expr, AttributeExpr):
            return self.get_element_attribute(expr.element, expr.attribute)

        if isinstance(expr, UnaryExpr):
            val = self.evaluate(expr.operand)
            if expr.operator == "-":
                return -val
            return val

        if isinstance(expr, BinaryExpr):
            left = self.evaluate(expr.left)
            right = self.evaluate(expr.right)

            if expr.operator == "+":
                return left + right
            if expr.operator == "-":
                return left - right
            if expr.operator == "*":
                return left * right
            if expr.operator == "/":
                if right == 0:
                    self._warn("Division by zero, using 0.0")
                    return 0.0
                return left / right
            if expr.operator == "^":
                return left**right

            raise self._error(f"Unknown operator: {expr.operator}")

        if isinstance(expr, FunctionExpr):
            name = expr.name.lower()

            # EXIST(var): returns 1 if variable is defined, 0 otherwise
            # Special case: checks name existence without evaluating
            if name == "exist":
                if len(expr.arguments) == 1 and isinstance(
                    expr.arguments[0], VariableExpr
                ):
                    return 1.0 if expr.arguments[0].name in self.variables else 0.0
                return 0.0

            # TABLE(): accesses MAD-X internal tables (TWISS, SURVEY, etc.)
            if name == "table":
                self._warn(
                    "TABLE() is not supported. "
                    "Use ImpactX diagnostics for beam moments."
                )
                return None

            if name not in MADXExpressionParser.FUNCTIONS:
                raise self._error(f"Unknown function: {name}")

            nargs, func = MADXExpressionParser.FUNCTIONS[name]
            if len(expr.arguments) != nargs:
                raise self._error(
                    f"Function {name} expects {nargs} arguments, "
                    f"got {len(expr.arguments)}"
                )

            args = [self.evaluate(arg) for arg in expr.arguments]
            return func(*args)

        if isinstance(expr, ListExpr):
            return [self.evaluate(e) for e in expr.elements]

        raise self._error(f"Cannot evaluate expression type: {type(expr)}")


# =============================================================================
# Main Parser
# =============================================================================


class MADXParser:
    """
    MAD-X Parser.

    Parses MAD-X input files and builds a representation of the lattice.
    """

    def __init__(self):
        self.context = EvaluationContext()
        self.tokens: list[Token] = []
        self.pos = 0
        self._current_file = None

    def _current(self) -> Token:
        """Get current token."""
        if self.pos < len(self.tokens):
            return self.tokens[self.pos]
        return self.tokens[-1]

    def _peek(self, offset: int = 0) -> Token:
        """Look ahead at token."""
        pos = self.pos + offset
        if pos < len(self.tokens):
            return self.tokens[pos]
        return self.tokens[-1]

    def _advance(self) -> Token:
        """Advance and return current token."""
        token = self._current()
        if self.pos < len(self.tokens) - 1:
            self.pos += 1
        return token

    def _expect(self, *types: TokenType) -> Token:
        """Expect one of the given token types."""
        token = self._current()
        if token.type not in types:
            expected = " or ".join(t.name for t in types)
            raise self._error(
                f"Expected {expected}, got {token.type.name} ({token.value!r})",
            )
        return self._advance()

    def _error(self, message: str) -> MADXInputError:
        """Create a MADXInputError with current file and line context."""
        line = self._current().line if self.tokens else None
        return MADXInputError(message, line_number=line, file=self._current_file)

    def _update_context_location(self):
        """Update file/line on the evaluation context for warning messages."""
        self.context.current_file = self._current_file
        self.context.current_line = self._current().line if self.tokens else None

    def _warn(self, message: str):
        """Emit a warning with current file and line context."""
        file_info = self._current_file or "<unknown>"
        line_info = self._current().line if self.tokens else "?"
        warnings.warn(f"{file_info}:{line_info}: {message}", MADXInputWarning)

    def _skip_semicolons(self):
        """Skip any semicolons."""
        while self._current().type == TokenType.SEMICOLON:
            self._advance()

    def parse(self, fn: str):
        """
        Parse a MAD-X file.

        Args:
            fn: Filename to parse
        """
        if not os.path.isfile(fn):
            raise FileNotFoundError(f"File '{fn}' not found!")

        self._current_file = fn

        with open(fn, "r") as f:
            text = f.read()

        self._parse_text(text)

    def parse_string(self, text: str):
        """Parse MAD-X input from a string."""
        self._current_file = "<string>"
        self._parse_text(text)

    def _parse_text(self, text: str):
        """Internal method to parse text."""
        lexer = MADXLexer(text, file=self._current_file)
        self.tokens = lexer.tokenize()
        self.pos = 0

        while self._current().type != TokenType.EOF:
            self._parse_statement()
            self._skip_semicolons()

    def _parse_statement(self):
        """Parse a single statement."""
        self._update_context_location()
        token = self._current()

        if token.type == TokenType.EOF:
            return

        if token.type == TokenType.SEMICOLON:
            self._advance()
            return

        if token.type != TokenType.IDENTIFIER:
            raise self._error(f"Expected identifier, got {token.type.name}")

        name = token.value.lower()

        # Check for label: definition
        if self._peek(1).type == TokenType.COLON:
            self._parse_labeled_definition()
            return

        # Check for macro definition: name(params): MACRO = { body };
        if self._peek(1).type == TokenType.LPAREN:
            # Scan ahead to see if this is name(...): MACRO
            scan = self.pos + 2
            depth = 1
            while scan < len(self.tokens) and depth > 0:
                if self.tokens[scan].type == TokenType.LPAREN:
                    depth += 1
                elif self.tokens[scan].type == TokenType.RPAREN:
                    depth -= 1
                scan += 1
            if scan < len(self.tokens) and self.tokens[scan].type == TokenType.COLON:
                self._parse_macro_definition()
                return

        # Check for variable assignment
        if self._peek(1).type in (TokenType.EQUALS, TokenType.COLON_EQUALS):
            self._parse_variable_assignment()
            return

        # Check for commands
        if name == "beam":
            self._parse_beam_command()
            return

        if name == "use":
            self._parse_use_command()
            return

        if name == "option":
            self._parse_option_command()
            return

        if name == "const":
            self._parse_const_declaration()
            return

        if name == "real":
            self._advance()  # skip 'real'
            self._parse_variable_assignment()
            return

        if name == "int":
            self._advance()  # skip 'int'
            self._parse_variable_assignment(integer=True)
            return

        if name == "call":
            self._parse_call_command()
            return

        if name == "exec":
            self._parse_exec_command()
            return

        if name == "return":
            # RETURN stops reading the current file (used in CALL'd files)
            self._skip_until_semicolon()
            raise _ReturnStatement()

        if name == "system":
            self._warn("SYSTEM command is unsafe and not supported; skipping")
            self._skip_until_semicolon()
            return

        if name in ("seqedit", "endedit", "flatten", "install", "move", "remove"):
            self._warn(
                f"'{name.upper()}' sequence editing command is not supported. "
                "Please open a GitHub issue if you need it: "
                "https://github.com/BLAST-ImpactX/impactx/issues"
            )
            self._skip_until_semicolon()
            return

        if name in ("assign", "printf", "chdir"):
            self._warn(
                f"'{name.upper()}' I/O command is not supported; "
                "use ImpactX commands for I/O. Skipping."
            )
            self._skip_until_semicolon()
            return

        if name == "if":
            self._parse_if_statement()
            return

        if name == "while":
            self._parse_while_statement()
            return

        if name in (
            "title",
            "select",
            "twiss",
            "print",
            "value",
            "show",
            "help",
            "stop",
            "quit",
            "exit",
            "endsequence",
        ):
            # Skip these commands - just consume until semicolon
            self._skip_until_semicolon()
            return

        # Unknown statement - skip it with a warning
        self._warn(f"Skipping unknown statement starting with '{name}'")
        self._skip_until_semicolon()

    def _skip_until_semicolon(self):
        """Skip tokens until semicolon or EOF."""
        while self._current().type not in (TokenType.SEMICOLON, TokenType.EOF):
            self._advance()
        if self._current().type == TokenType.SEMICOLON:
            self._advance()

    def _parse_if_condition(self) -> bool:
        """Parse a condition: ( logical_expr ).

        Supports compound conditions with && and ||, e.g.:
            IF ( a > 0 && b < 1 )
            IF ( x == 1 || y == 2 )

        MAD-X comparison operators (from mad_eval.c):
            ==, <>, <, >, <=, >=
        """
        self._expect(TokenType.LPAREN)
        result = self._parse_logical_or()
        self._expect(TokenType.RPAREN)
        return result

    def _parse_logical_or(self) -> bool:
        """Parse logical OR: expr || expr || ..."""
        result = self._parse_logical_and()
        while self._current().type == TokenType.OR:
            self._advance()
            right = self._parse_logical_and()
            result = result or right
        return result

    def _parse_logical_and(self) -> bool:
        """Parse logical AND: expr && expr && ..."""
        result = self._parse_comparison()
        while self._current().type == TokenType.AND:
            self._advance()
            right = self._parse_comparison()
            result = result and right
        return result

    def _parse_comparison(self) -> bool:
        """Parse a single comparison: arith_expr [op arith_expr]."""
        left_expr = self._parse_condition_operand()
        left_val = self.context.evaluate(left_expr)

        op = self._try_parse_comparison_op()
        if op is None:
            return bool(left_val)

        right_expr = self._parse_condition_operand()
        right_val = self.context.evaluate(right_expr)

        _cmp = {
            "==": lambda a, b: a == b,
            "<>": lambda a, b: a != b,
            "<": lambda a, b: a < b,
            ">": lambda a, b: a > b,
            "<=": lambda a, b: a <= b,
            ">=": lambda a, b: a >= b,
        }
        return _cmp[op](left_val, right_val)

    def _try_parse_comparison_op(self) -> Optional[str]:
        """Try to consume a comparison operator. Returns operator string or None."""
        tok = self._current()
        nt = self._peek(1)

        if tok.type == TokenType.EQUALS and nt.type == TokenType.EQUALS:
            self._advance()
            self._advance()
            return "=="
        if tok.type == TokenType.LT:
            if nt.type == TokenType.GT:
                self._advance()
                self._advance()
                return "<>"
            if nt.type == TokenType.EQUALS:
                self._advance()
                self._advance()
                return "<="
            self._advance()
            return "<"
        if tok.type == TokenType.GT:
            if nt.type == TokenType.EQUALS:
                self._advance()
                self._advance()
                return ">="
            self._advance()
            return ">"
        return None

    def _parse_condition_operand(self) -> Expression:
        """Parse an arithmetic expression inside a condition.

        Stops at comparison operators, logical operators, or closing paren.
        """
        start = self.pos
        paren_depth = 0

        while self._current().type != TokenType.EOF:
            tok = self._current()
            if tok.type == TokenType.LPAREN:
                paren_depth += 1
            elif tok.type == TokenType.RPAREN:
                if paren_depth == 0:
                    break
                paren_depth -= 1
            elif paren_depth == 0:
                # Stop at logical operators
                if tok.type in (TokenType.AND, TokenType.OR):
                    break
                # Stop at comparison operators: ==, <>, <, <=, >, >=
                if tok.type == TokenType.LT or tok.type == TokenType.GT:
                    break
                if (
                    tok.type == TokenType.EQUALS
                    and self._peek(1).type == TokenType.EQUALS
                ):
                    break
            self.pos += 1

        end = self.pos
        self.pos = start

        expr_tokens = self.tokens[start:end] + [Token(TokenType.EOF, None, 0, 0)]
        expr_parser = MADXExpressionParser(expr_tokens)
        expr = expr_parser.parse_expression()
        self.pos = end
        return expr

    def _parse_if_statement(self):
        """Parse a MAD-X IF statement.

        Syntax:
            IF ( condition ) { statements; };
            IF ( condition ) { statements; } ELSEIF ( condition ) { statements; } ELSE { statements; };
        """
        self._advance()  # skip 'if'

        condition = self._parse_if_condition()
        self._execute_or_skip_brace_block(condition)

        # Handle ELSEIF / ELSE
        while (
            self._current().type == TokenType.IDENTIFIER
            and self._current().value == "elseif"
        ):
            self._advance()  # skip 'elseif'
            # Always parse condition to consume tokens, even if a prior branch matched
            this_cond = self._parse_if_condition()
            execute = not condition and this_cond
            if execute:
                condition = True  # mark so subsequent branches are skipped
            self._execute_or_skip_brace_block(execute)

        # Handle ELSE
        if (
            self._current().type == TokenType.IDENTIFIER
            and self._current().value == "else"
        ):
            self._advance()  # skip 'else'
            self._execute_or_skip_brace_block(not condition)

    def _execute_or_skip_brace_block(self, execute: bool):
        """Execute or skip a { ... } block, then consume trailing semicolon."""
        self._expect(TokenType.LBRACE)
        if execute:
            while self._current().type not in (TokenType.RBRACE, TokenType.EOF):
                self._parse_statement()
                self._skip_semicolons()
        else:
            depth = 1
            while depth > 0 and self._current().type != TokenType.EOF:
                if self._current().type == TokenType.LBRACE:
                    depth += 1
                elif self._current().type == TokenType.RBRACE:
                    depth -= 1
                    if depth == 0:
                        break
                self._advance()
        self._expect(TokenType.RBRACE)
        if self._current().type == TokenType.SEMICOLON:
            self._advance()

    def _parse_while_statement(self):
        """Parse a MAD-X WHILE statement.

        Syntax:
            WHILE ( condition ) { statements; }
        """
        self._advance()  # skip 'while'

        # Save token position before the condition so we can loop back
        max_iterations = 100000  # safety limit
        loop_start = self.pos

        for _ in range(max_iterations):
            self.pos = loop_start
            condition = self._parse_if_condition()

            if not condition:
                self._execute_or_skip_brace_block(False)
                break

            self._execute_or_skip_brace_block(True)
        else:
            raise self._error(f"WHILE loop exceeded {max_iterations} iterations")

    def _parse_macro_definition(self):
        """Parse a macro definition: name(params): MACRO = { body };"""
        macro_name = self._advance().value  # name
        self._expect(TokenType.LPAREN)

        # Parse parameter names
        params = []
        if self._current().type == TokenType.IDENTIFIER:
            params.append(self._advance().value)
            while self._current().type == TokenType.COMMA:
                self._advance()
                params.append(self._expect(TokenType.IDENTIFIER).value)
        self._expect(TokenType.RPAREN)
        self._expect(TokenType.COLON)

        keyword = self._expect(TokenType.IDENTIFIER).value.lower()
        if keyword != "macro":
            raise self._error(f"Expected MACRO, got '{keyword}'")

        self._expect(TokenType.EQUALS)
        self._expect(TokenType.LBRACE)

        # Capture body text between { and }
        body_tokens = []
        depth = 1
        while depth > 0 and self._current().type != TokenType.EOF:
            if self._current().type == TokenType.LBRACE:
                depth += 1
            elif self._current().type == TokenType.RBRACE:
                depth -= 1
                if depth == 0:
                    break
            body_tokens.append(self._current())
            self._advance()

        self._expect(TokenType.RBRACE)
        self._expect(TokenType.SEMICOLON)

        # Store the macro: params and body tokens
        self.context.macros[macro_name] = (params, body_tokens)

    def _parse_exec_command(self):
        """Parse EXEC, MacroName($arg1, $arg2, ...);

        Substitutes $argN with the corresponding variable value,
        then parses the expanded macro body.
        """
        self._advance()  # skip 'exec'
        self._expect(TokenType.COMMA)

        macro_name = self._expect(TokenType.IDENTIFIER).value

        # Parse arguments
        args = []
        if self._current().type == TokenType.LPAREN:
            self._advance()
            if self._current().type != TokenType.RPAREN:
                # $var references - the $ is part of the identifier
                args.append(self._expect(TokenType.IDENTIFIER).value)
                while self._current().type == TokenType.COMMA:
                    self._advance()
                    args.append(self._expect(TokenType.IDENTIFIER).value)
            self._expect(TokenType.RPAREN)

        self._expect(TokenType.SEMICOLON)

        if macro_name not in self.context.macros:
            self._warn(f"Unknown macro '{macro_name}', skipping EXEC")
            return

        params, body_tokens = self.context.macros[macro_name]

        # Build substitution map: param_name -> arg identifier string
        # MAD-X macros do textual substitution: the argument name replaces
        # the formal parameter wherever it appears in the body tokens.
        subs = {}
        for i, param in enumerate(params):
            if i < len(args):
                # Strip leading $ from argument name
                subs[param] = args[i].lstrip("$")

        # Substitute parameter references in body tokens
        expanded = []
        for tok in body_tokens:
            if tok.type == TokenType.IDENTIFIER and tok.value in subs:
                replacement = subs[tok.value]
                # Textual substitution: replace with an identifier token
                # (the identifier will be resolved as a variable, element
                # name, or line reference during normal parsing)
                expanded.append(
                    Token(TokenType.IDENTIFIER, replacement, tok.line, tok.column)
                )
            else:
                expanded.append(tok)

        # Parse the expanded body
        expanded.append(Token(TokenType.EOF, None, 0, 0))
        saved_tokens = self.tokens
        saved_pos = self.pos
        self.tokens = expanded
        self.pos = 0

        while self._current().type != TokenType.EOF:
            self._parse_statement()
            self._skip_semicolons()

        self.tokens = saved_tokens
        self.pos = saved_pos

    def _parse_labeled_definition(self):
        """Parse a labeled definition (element or line)."""
        name = self._advance().value  # identifier
        self._expect(TokenType.COLON)  # :

        # Get the type/keyword
        keyword = self._expect(TokenType.IDENTIFIER).value.lower()

        if keyword == "line":
            self._parse_line_definition(name)
        elif keyword == "sequence":
            self._parse_sequence_definition(name)
        else:
            self._parse_element_definition(name, keyword)

    def _parse_element_definition(self, name: str, element_type: str):
        """Parse an element definition."""
        attributes = {}
        attribute_exprs = {}

        # Check if element_type is actually a reference to a previously defined element
        # (element inheritance in MAD-X, e.g., "fmagu01: fmag;")
        parent = self.context.elements.get(element_type)
        if parent is not None:
            # Inherit type and attributes from parent element
            element_type = parent.type
            attributes.update(parent.attributes)
            attribute_exprs.update(parent.attribute_exprs)

        # Parse attributes
        # MAD-X allows commas between attributes but they are optional
        while True:
            has_comma = False
            if self._current().type == TokenType.COMMA:
                self._advance()  # ,
                has_comma = True

            if self._current().type != TokenType.IDENTIFIER:
                break

            # Without a comma, only continue if the identifier is followed
            # by = or := (otherwise it's the start of a new statement)
            if not has_comma and self._peek(1).type not in (
                TokenType.EQUALS,
                TokenType.COLON_EQUALS,
            ):
                break

            attr_name = self._advance().value.lower()

            # Check for = or :=
            if self._current().type == TokenType.COLON_EQUALS:
                self._advance()
                expr = self._parse_expression()
                attribute_exprs[attr_name] = expr
                # := means deferred evaluation; don't evaluate now
                attributes[attr_name] = None
            elif self._current().type == TokenType.EQUALS:
                self._advance()
                expr = self._parse_expression()
                attributes[attr_name] = self.context.evaluate(expr)
            else:
                # Boolean flag
                attributes[attr_name] = True

        self._expect(TokenType.SEMICOLON)

        # Store the element
        element = Element(name, element_type, attributes, attribute_exprs)
        self.context.elements[name] = element

    def _parse_line_definition(self, name: str):
        """Parse a LINE definition."""
        self._expect(TokenType.EQUALS)
        self._expect(TokenType.LPAREN)

        elements = self._parse_line_elements()

        self._expect(TokenType.RPAREN)
        self._expect(TokenType.SEMICOLON)

        line = Line(name, elements)
        self.context.lines[name] = line

    def _get_element_length(self, elem_name: str) -> float:
        """Get an element's length, evaluating deferred expressions if needed."""
        elem = self.context.elements.get(elem_name)
        if elem is None:
            return 0.0
        if "l" in elem.attribute_exprs:
            return self.context.evaluate(elem.attribute_exprs["l"])
        if "l" in elem.attributes and elem.attributes["l"] is not None:
            return elem.attributes["l"]
        return 0.0

    def _parse_sequence_definition(self, name: str):
        """Parse a SEQUENCE definition.

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
        # Parse sequence-level attributes (L=..., REFER=..., etc.)
        refer = "centre"  # MAD-X default

        while self._current().type == TokenType.COMMA:
            self._advance()
            if self._current().type != TokenType.IDENTIFIER:
                break
            attr_name = self._advance().value.lower()
            if self._current().type in (TokenType.EQUALS, TokenType.COLON_EQUALS):
                self._advance()
                if attr_name == "refer":
                    # REFER value is an identifier: centre, entry, or exit
                    if self._current().type == TokenType.IDENTIFIER:
                        refer = self._advance().value.lower()
                    else:
                        self._parse_expression()  # consume unexpected form
                else:
                    self._parse_expression()  # consume but discard
            # else: boolean flag, ignore

        self._expect(TokenType.SEMICOLON)

        # Parse sequence entries until ENDSEQUENCE
        # Each entry: (at_position, element_name, from_element_or_None)
        seq_entries = []

        while self._current().type != TokenType.EOF:
            token = self._current()

            if token.type == TokenType.SEMICOLON:
                self._advance()
                continue

            if token.type != TokenType.IDENTIFIER:
                break

            entry_name = token.value.lower()

            # Check for ENDSEQUENCE
            if entry_name == "endsequence":
                self._advance()
                self._expect(TokenType.SEMICOLON)
                break

            # Two forms:
            #   1) labeled:   entry_name : type, at=pos, from=ref;
            #   2) unlabeled: entry_name, at=pos, from=ref;
            if self._peek(1).type == TokenType.COLON:
                # Labeled: parse as element definition first, then extract 'at'/'from'
                self._parse_labeled_definition()
                elem = self.context.elements.get(entry_name)
                at_pos = 0.0
                from_elem = None
                if elem:
                    if "at" in elem.attributes and elem.attributes["at"] is not None:
                        at_pos = elem.attributes["at"]
                    from_val = elem.attributes.get("from")
                    if isinstance(from_val, str):
                        from_elem = from_val
                seq_entries.append((at_pos, entry_name, from_elem))
            else:
                # Unlabeled: element_name, at = expr, from = name ;
                self._advance()  # consume element name
                at_pos = 0.0
                from_elem = None
                while self._current().type == TokenType.COMMA:
                    self._advance()
                    if self._current().type != TokenType.IDENTIFIER:
                        break
                    attr_name = self._advance().value.lower()
                    if self._current().type in (
                        TokenType.EQUALS,
                        TokenType.COLON_EQUALS,
                    ):
                        self._advance()
                        if attr_name == "from":
                            # FROM value is an element name (identifier)
                            if self._current().type == TokenType.IDENTIFIER:
                                from_elem = self._advance().value
                            else:
                                self._parse_expression()  # consume
                        else:
                            expr = self._parse_expression()
                            if attr_name == "at":
                                at_pos = self.context.evaluate(expr)
                self._expect(TokenType.SEMICOLON)
                seq_entries.append((at_pos, entry_name, from_elem))

        # Resolve FROM references: position is relative to the named element
        if any(f is not None for _, _, f in seq_entries):
            # Build map of element name -> absolute position (for entries without FROM)
            abs_positions = {}
            for at_pos, elem_name, from_elem in seq_entries:
                if from_elem is None:
                    abs_positions[elem_name] = at_pos

            resolved = []
            for at_pos, elem_name, from_elem in seq_entries:
                if from_elem is not None and from_elem in abs_positions:
                    resolved.append((at_pos + abs_positions[from_elem], elem_name))
                else:
                    resolved.append((at_pos, elem_name))
            seq_entries = [(a, n, None) for a, n in resolved]

        # Adjust positions based on REFER and sort
        if refer != "centre":
            adjusted = []
            for at_pos, elem_name, _ in seq_entries:
                elem_len = self._get_element_length(elem_name)
                if refer == "entry":
                    # 'at' marks the entry; center is at + L/2
                    adjusted.append((at_pos + elem_len / 2.0, elem_name, None))
                elif refer == "exit":
                    # 'at' marks the exit; center is at - L/2
                    adjusted.append((at_pos - elem_len / 2.0, elem_name, None))
            seq_entries = adjusted

        seq_entries.sort(key=lambda x: x[0])

        # Insert drifts between elements to fill gaps
        elements = []
        cursor = 0.0  # tracks the end of the last element

        for center_pos, elem_name, _ in seq_entries:
            elem_len = self._get_element_length(elem_name)
            entry_pos = center_pos - elem_len / 2.0

            gap = entry_pos - cursor
            if gap > 1e-15:
                drift_name = f"_drift_{name}_{len(elements)}"
                self.context.elements[drift_name] = Element(
                    drift_name, "drift", {"l": gap}
                )
                elements.append(
                    {"name": drift_name, "multiplier": 1, "reversed": False}
                )
            elif gap < -1e-10:
                self._warn(
                    f"SEQUENCE '{name}': elements overlap by {-gap:.6g} m "
                    f"at '{elem_name}' (pos={center_pos})"
                )

            elements.append({"name": elem_name, "multiplier": 1, "reversed": False})
            cursor = center_pos + elem_len / 2.0

        line = Line(name, elements)
        self.context.lines[name] = line

    def _parse_line_elements(self) -> list:
        """Parse the elements inside a LINE definition."""
        elements = []

        while True:
            elem = self._parse_line_element()
            if elem is not None:
                elements.append(elem)

            if self._current().type == TokenType.COMMA:
                self._advance()  # ,
            elif self._current().type in (
                TokenType.IDENTIFIER,
                TokenType.NUMBER,
                TokenType.MINUS,
                TokenType.LPAREN,
            ):
                # Missing comma — tolerate it (common in real MAD-X files)
                pass
            else:
                break

        return elements

    def _parse_line_element(self):
        """Parse a single line element (possibly with multiplier or reversal)."""
        # Check for reversal operator (-)
        reversed_flag = False
        if self._current().type == TokenType.MINUS:
            self._advance()
            reversed_flag = True

        # Check for multiplier: n * element or (n * subline)
        multiplier = 1

        # Check for parenthesized group with potential multiplier
        if self._current().type == TokenType.LPAREN:
            # Could be (n * element) or just (elements)
            self._advance()

            # Check if first thing is a number followed by *
            if (
                self._current().type == TokenType.NUMBER
                and self._peek(1).type == TokenType.STAR
            ):
                multiplier = int(self._advance().value)
                self._expect(TokenType.STAR)

            # Now parse the inner elements
            inner_elements = self._parse_line_elements()
            self._expect(TokenType.RPAREN)

            return {
                "multiplier": multiplier,
                "elements": inner_elements,
                "reversed": reversed_flag,
            }

        # Check for number * element or number * (group)
        has_multiplier = False
        if self._current().type == TokenType.NUMBER:
            num = self._advance().value
            if self._current().type == TokenType.STAR:
                self._advance()
                multiplier = int(num)
                has_multiplier = True
            else:
                # Shouldn't happen in valid MAD-X
                raise self._error("Expected * after number in LINE")

        # Parenthesized group after multiplier: n * (elements)
        if has_multiplier and self._current().type == TokenType.LPAREN:
            self._advance()
            inner_elements = self._parse_line_elements()
            self._expect(TokenType.RPAREN)
            return {
                "multiplier": multiplier,
                "elements": inner_elements,
                "reversed": reversed_flag,
            }

        # Now we should have an identifier
        if self._current().type == TokenType.IDENTIFIER:
            name = self._advance().value

            # Check for identifier used as multiplier: varname * element_or_group
            if self._current().type == TokenType.STAR:
                self._advance()
                multiplier = int(self.context.get_variable(name))

                # Parenthesized group: var * (elements)
                if self._current().type == TokenType.LPAREN:
                    self._advance()
                    inner_elements = self._parse_line_elements()
                    self._expect(TokenType.RPAREN)
                    return {
                        "multiplier": multiplier,
                        "elements": inner_elements,
                        "reversed": reversed_flag,
                    }

                # Single element: var * element
                if self._current().type == TokenType.IDENTIFIER:
                    elem_name = self._advance().value
                    return {
                        "name": elem_name,
                        "multiplier": multiplier,
                        "reversed": reversed_flag,
                    }

            return {
                "name": name,
                "multiplier": multiplier,
                "reversed": reversed_flag,
            }

        return None

    def _parse_expression(self) -> Expression:
        """Parse an expression using the expression parser."""
        # Find the extent of the expression (until comma, semicolon, or rparen at depth 0)
        start = self.pos
        paren_depth = 0
        brace_depth = 0

        while self._current().type != TokenType.EOF:
            if self._current().type == TokenType.LPAREN:
                paren_depth += 1
            elif self._current().type == TokenType.RPAREN:
                if paren_depth == 0:
                    break
                paren_depth -= 1
            elif self._current().type == TokenType.LBRACE:
                brace_depth += 1
            elif self._current().type == TokenType.RBRACE:
                brace_depth -= 1
            elif (
                self._current().type in (TokenType.COMMA, TokenType.SEMICOLON)
                and paren_depth == 0
                and brace_depth == 0
            ):
                break
            self.pos += 1

        end = self.pos
        self.pos = start

        # Extract tokens for expression parser
        expr_tokens = self.tokens[start:end] + [Token(TokenType.EOF, None, 0, 0)]
        expr_parser = MADXExpressionParser(expr_tokens)
        expr = expr_parser.parse_expression()

        # Advance past the expression tokens
        self.pos = end

        return expr

    def _parse_variable_assignment(self, integer: bool = False):
        """Parse a variable assignment."""
        name = self._expect(TokenType.IDENTIFIER).value

        deferred = False
        if self._current().type == TokenType.COLON_EQUALS:
            self._advance()
            deferred = True
        else:
            self._expect(TokenType.EQUALS)

        expr = self._parse_expression()
        self._expect(TokenType.SEMICOLON)

        if integer and not deferred:
            value = int(self.context.evaluate(expr))
            self.context.set_variable(
                name, value=value, expression=expr, deferred=False
            )
        else:
            self.context.set_variable(name, expression=expr, deferred=deferred)

    def _parse_const_declaration(self):
        """Parse a CONST declaration."""
        self._advance()  # skip 'const'

        # Optional 'real' or 'int'
        integer = False
        if self._current().type == TokenType.IDENTIFIER:
            if self._current().value == "int":
                integer = True
                self._advance()
            elif self._current().value == "real":
                self._advance()

        name = self._expect(TokenType.IDENTIFIER).value

        deferred = False
        if self._current().type == TokenType.COLON_EQUALS:
            self._advance()
            deferred = True
        else:
            self._expect(TokenType.EQUALS)

        expr = self._parse_expression()
        self._expect(TokenType.SEMICOLON)

        if deferred:
            self.context.set_variable(
                name, expression=expr, deferred=True, constant=True
            )
        else:
            value = self.context.evaluate(expr)
            if integer:
                value = int(value)
            self.context.set_variable(
                name, value=value, expression=expr, deferred=False, constant=True
            )

    def _parse_beam_command(self):
        """Parse the BEAM command.

        MAD-X allows multiple BEAM commands for different sequences via
        BEAM, SEQUENCE=name, ...;  The beam is stored per-sequence.
        A BEAM without SEQUENCE= updates the default beam.
        """
        self._advance()  # skip 'beam'

        # Collect all attributes first, then apply to the right beam
        beam_attrs = {}
        sequence_name = None

        while self._current().type == TokenType.COMMA:
            self._advance()

            if self._current().type != TokenType.IDENTIFIER:
                break

            attr_name = self._advance().value.lower()

            if self._current().type in (TokenType.EQUALS, TokenType.COLON_EQUALS):
                self._advance()

                if attr_name == "particle":
                    if self._current().type == TokenType.STRING:
                        beam_attrs["particle"] = self._advance().value.lower()
                    elif self._current().type == TokenType.IDENTIFIER:
                        beam_attrs["particle"] = self._advance().value.lower()
                    else:
                        raise self._error("Expected particle name")
                elif attr_name == "sequence":
                    if self._current().type in (
                        TokenType.STRING,
                        TokenType.IDENTIFIER,
                    ):
                        sequence_name = self._advance().value.lower()
                    else:
                        self._parse_expression()  # consume
                elif attr_name in ("energy", "pc", "mass", "charge", "freq0"):
                    expr = self._parse_expression()
                    beam_attrs[attr_name] = self.context.evaluate(expr)
                else:
                    # Unknown attribute, skip the value
                    self._parse_expression()

        self._expect(TokenType.SEMICOLON)

        # Apply to the appropriate beam
        beam = self.context.get_beam_for_sequence(sequence_name)
        for attr, value in beam_attrs.items():
            setattr(beam, attr, value)

    def _parse_use_command(self):
        """Parse the USE command."""
        self._advance()  # skip 'use'

        while self._current().type == TokenType.COMMA:
            self._advance()

            if self._current().type != TokenType.IDENTIFIER:
                break

            attr_name = self._advance().value.lower()

            if attr_name in ("sequence", "period"):
                self._expect(TokenType.EQUALS)
                if self._current().type == TokenType.IDENTIFIER:
                    self.context.selected_sequence = self._advance().value
                elif self._current().type == TokenType.STRING:
                    self.context.selected_sequence = self._advance().value.lower()
            else:
                if self._current().type == TokenType.EQUALS:
                    self._advance()
                    self._parse_expression()

        self._expect(TokenType.SEMICOLON)

    def _parse_option_command(self):
        """Parse the OPTION command and persist option values.

        Supports forms such as:
            OPTION, RBARC=true, THIN_FOC=false;
            OPTION, ECHO;
            OPTION, -WARN;
        """
        self._advance()  # skip 'option'

        while self._current().type == TokenType.COMMA:
            self._advance()

            negate = False
            if self._current().type == TokenType.MINUS:
                negate = True
                self._advance()

            if self._current().type != TokenType.IDENTIFIER:
                break

            opt_name = self._advance().value.lower()

            if self._current().type in (TokenType.EQUALS, TokenType.COLON_EQUALS):
                self._advance()
                expr = self._parse_expression()
                value = self.context.evaluate(expr)
            else:
                # Bare flag means "set true" unless explicitly negated.
                value = 0.0 if negate else 1.0

            # A leading "-" explicitly disables the option.
            if negate:
                value = 0.0

            self.context.set_option(opt_name, value)

        self._expect(TokenType.SEMICOLON)

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
        self._advance()  # skip 'call'

        filename = None

        while self._current().type == TokenType.COMMA:
            self._advance()  # skip comma

            if self._current().type != TokenType.IDENTIFIER:
                break

            attr_name = self._advance().value.lower()

            if attr_name in ("file", "filename"):
                self._expect(TokenType.EQUALS)
                if self._current().type == TokenType.STRING:
                    filename = self._advance().value
                elif self._current().type == TokenType.IDENTIFIER:
                    # Some MAD-X files use unquoted filenames
                    filename = self._advance().value
                else:
                    raise self._error("Expected filename after FILE=")
            else:
                # Unknown attribute, skip its value
                if self._current().type == TokenType.EQUALS:
                    self._advance()
                    self._parse_expression()

        self._expect(TokenType.SEMICOLON)

        if filename is None:
            raise self._error("CALL command requires FILE= or FILENAME= attribute")

        # Resolve the filename relative to the current file's directory
        if self._current_file and self._current_file != "<string>":
            base_dir = os.path.dirname(os.path.abspath(self._current_file))
            resolved = os.path.join(base_dir, filename)
        else:
            resolved = filename

        if not os.path.isfile(resolved):
            raise self._error(
                f"CALL: File '{filename}' not found (resolved to '{resolved}')"
            )

        # Save parser state
        saved_tokens = self.tokens
        saved_pos = self.pos
        saved_file = self._current_file

        try:
            # Parse the called file
            self._current_file = resolved
            with open(resolved, "r") as f:
                text = f.read()

            lexer = MADXLexer(text, file=self._current_file)
            self.tokens = lexer.tokenize()
            self.pos = 0

            while self._current().type != TokenType.EOF:
                self._parse_statement()
                self._skip_semicolons()
        except _ReturnStatement:
            # RETURN was encountered in the called file — stop reading it
            pass
        finally:
            # Restore parser state
            self.tokens = saved_tokens
            self.pos = saved_pos
            self._current_file = saved_file

    # =========================================================================
    # Public Interface (backward compatible with old parser)
    # =========================================================================

    @property
    def beam(self) -> dict:
        """Get beam parameters as a dictionary (backward compatible)."""
        return {
            "particle": self.context.beam.particle,
            "energy": self.context.beam.energy,
        }

    @property
    def sequence(self) -> dict:
        """Get selected sequence (backward compatible)."""
        return {"name": self.context.selected_sequence or ""}

    def _flatten_line(self, line_name: str) -> list[str]:
        """Recursively flatten a line definition to element names."""
        line_name = line_name.lower()

        if line_name not in self.context.lines:
            # It might be a direct element reference
            if line_name in self.context.elements:
                return [line_name]
            raise self._error(f"Unknown line or element: {line_name}")

        line = self.context.lines[line_name]
        result = []

        for elem in line.elements:
            if isinstance(elem, dict):
                multiplier = elem.get("multiplier", 1)
                reversed_flag = elem.get("reversed", False)

                if "elements" in elem:
                    # Nested group
                    for _ in range(multiplier):
                        inner = []
                        for inner_elem in elem["elements"]:
                            inner.extend(self._process_line_element(inner_elem))
                        if reversed_flag:
                            inner = inner[::-1]
                        result.extend(inner)
                elif "name" in elem:
                    name = elem["name"]
                    for _ in range(multiplier):
                        expanded = self._flatten_line(name)
                        if reversed_flag:
                            expanded = expanded[::-1]
                        result.extend(expanded)

        return result

    def _process_line_element(self, elem) -> list[str]:
        """Process a single line element specification."""
        if isinstance(elem, str):
            return self._flatten_line(elem)

        if isinstance(elem, dict):
            multiplier = elem.get("multiplier", 1)
            reversed_flag = elem.get("reversed", False)

            if "elements" in elem:
                inner = []
                for inner_elem in elem["elements"]:
                    inner.extend(self._process_line_element(inner_elem))
                if reversed_flag:
                    inner = inner[::-1]
                return inner * multiplier
            elif "name" in elem:
                expanded = self._flatten_line(elem["name"])
                if reversed_flag:
                    expanded = expanded[::-1]
                return expanded * multiplier

        return []

    def getBeamline(self) -> list[dict]:
        """
        Get the beamline as a list of element dictionaries.

        Returns a list compatible with the old parser format.
        """
        if not self.context.selected_sequence:
            # Try to find a line if no sequence selected
            if self.context.lines:
                self.context.selected_sequence = list(self.context.lines.keys())[0]
            else:
                return []

        # Flatten the line to get element names
        element_names = self._flatten_line(self.context.selected_sequence)

        # Build the beamline
        beamline = []
        for name in element_names:
            if name in self.context.elements:
                elem = self.context.elements[name]

                # Re-evaluate any deferred expressions
                attrs = {}
                for attr_name, value in elem.attributes.items():
                    if attr_name in elem.attribute_exprs:
                        attrs[attr_name] = self.context.evaluate(
                            elem.attribute_exprs[attr_name]
                        )
                    else:
                        attrs[attr_name] = value

                # Remove the MAD-X "type" attribute (a user annotation)
                # so it doesn't overwrite the element's actual type
                attrs.pop("type", None)

                elem_dict = {
                    "name": elem.name,
                    "type": elem.type,
                    **attrs,
                }
                beamline.append(elem_dict)
            else:
                self._warn(f"Element '{name}' referenced but not defined")

        return beamline

    def getParticle(self) -> str:
        """Get the particle type."""
        particle = self.context.beam.particle
        known_particles = [
            "positron",
            "electron",
            "proton",
            "antiproton",
            "posmuon",
            "negmuon",
            "ion",
        ]

        if particle and particle not in known_particles:
            self._warn(f"Unknown particle type '{particle}'")

        return particle

    def getEtot(self) -> float:
        """Get total energy in GeV."""
        return self.context.beam.energy

    def getFreq0(self) -> float:
        """Get revolution frequency in MHz."""
        return self.context.beam.freq0

    def getOption(self, name: str, default: Any = 0.0) -> Any:
        """Get an OPTION value from the parsed MAD-X input."""
        return self.context.get_option(name, default)

    def getOptions(self) -> dict[str, Any]:
        """Get all parsed OPTION values."""
        return dict(self.context.options)

    def getVariable(self, name: str) -> Any:
        """Get a variable value."""
        return self.context.get_variable(name)

    def getElement(self, name: str) -> Optional[Element]:
        """Get an element by name."""
        return self.context.elements.get(name.lower())

    def getLine(self, name: str) -> Optional[Line]:
        """Get a line by name."""
        return self.context.lines.get(name.lower())

    def __str__(self) -> str:
        """String representation with lattice information."""
        if not self.context.selected_sequence:
            return "No sequence selected."

        try:
            beamline = self.getBeamline()
        except Exception as e:
            return f"Error getting beamline: {e}"

        if not beamline:
            return "Empty beamline."

        length = 0.0
        type_counts = {}

        for elem in beamline:
            if "l" in elem:
                length += elem["l"]
            elem_type = elem.get("type", "unknown")
            type_counts[elem_type] = type_counts.get(elem_type, 0) + 1

        sign = "*" * 70
        lines = [
            sign,
            "MAD-X Parser Information:",
            f"         length: {length:.6f} [m]",
            f"      #elements: {len(beamline)}",
        ]

        for elem_type, count in sorted(type_counts.items()):
            lines.append(f"            * #{elem_type}: {count}")

        lines.extend(
            [
                "           beam:",
                f"            *     particle: {self.context.beam.particle}",
                f"            * total energy: {self.context.beam.energy} [GeV]",
                sign,
            ]
        )

        return "\n".join(lines)
