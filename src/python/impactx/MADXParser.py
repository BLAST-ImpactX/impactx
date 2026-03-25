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
"""

import math
import os
import warnings
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Any, Optional


class MADXParserError(Exception):
    """Base exception for MAD-X parser errors."""

    pass


class MADXInputError(MADXParserError):
    """Error in MAD-X input syntax or semantics."""

    def __init__(self, message, line_number=None, context=None):
        self.message = message
        self.line_number = line_number
        self.context = context
        super().__init__(self._format_message())

    def _format_message(self):
        msg = self.message
        if self.line_number is not None:
            msg = f"Line {self.line_number}: {msg}"
        if self.context:
            msg = f"{msg}\n  Context: {self.context}"
        return msg


class MADXInputWarning(UserWarning):
    """Warning for non-fatal MAD-X input issues."""

    pass


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

    def __init__(self, text: str):
        self.text = text
        self.pos = 0
        self.line = 1
        self.column = 1
        self.length = len(text)

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
            raise MADXInputError("Unterminated multi-line comment", self.line)

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
            raise MADXInputError(f"Invalid number: {value_str}", start_line)

        return Token(TokenType.NUMBER, value, start_line, start_col)

    def _read_string(self) -> Token:
        """Read a quoted string."""
        start_line = self.line
        start_col = self.column
        quote = self._advance()  # ' or "
        chars = []

        while self._peek() is not None and self._peek() != quote:
            if self._peek() == "\n":
                raise MADXInputError("Unterminated string", start_line)
            chars.append(self._advance())

        if self._peek() is None:
            raise MADXInputError("Unterminated string", start_line)

        self._advance()  # closing quote
        return Token(TokenType.STRING, "".join(chars), start_line, start_col)

    def _read_identifier(self) -> Token:
        """Read an identifier (converted to lowercase)."""
        start_line = self.line
        start_col = self.column
        chars = []

        # First character must be a letter
        if self._peek() is not None and (self._peek().isalpha() or self._peek() == "_"):
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

            # Identifier
            elif char.isalpha() or char == "_":
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

            else:
                raise MADXInputError(f"Unexpected character: {char!r}", start_line)

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
    }

    def __init__(self, tokens: list[Token]):
        self.tokens = tokens
        self.pos = 0

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
                f"Expected {expected}, got {token.type.name}", token.line
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
            f"Unexpected token in expression: {token.type.name}", token.line
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
    # Add other beam parameters as needed


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
        self.beam = Beam()
        self.selected_sequence: Optional[str] = None

        # Initialize predefined constants
        for name, value in self.PREDEFINED_CONSTANTS.items():
            self.variables[name] = Variable(name, value, constant=True)

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
                raise MADXInputError(f"Cannot modify constant: {name}")

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
            warnings.warn(f"Undefined variable '{name}', using 0.0", MADXInputWarning)
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
            raise MADXInputError(f"Unknown element: {element_name}")

        element = self.elements[element_name]

        # Check if attribute has a deferred expression
        if attr_name in element.attribute_exprs:
            return self.evaluate(element.attribute_exprs[attr_name])

        if attr_name in element.attributes:
            return element.attributes[attr_name]

        raise MADXInputError(f"Element '{element_name}' has no attribute '{attr_name}'")

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
                    raise MADXInputError("Division by zero")
                return left / right
            if expr.operator == "^":
                return left**right

            raise MADXInputError(f"Unknown operator: {expr.operator}")

        if isinstance(expr, FunctionExpr):
            name = expr.name.lower()
            if name not in MADXExpressionParser.FUNCTIONS:
                raise MADXInputError(f"Unknown function: {name}")

            nargs, func = MADXExpressionParser.FUNCTIONS[name]
            if len(expr.arguments) != nargs:
                raise MADXInputError(
                    f"Function {name} expects {nargs} arguments, got {len(expr.arguments)}"
                )

            args = [self.evaluate(arg) for arg in expr.arguments]
            return func(*args)

        if isinstance(expr, ListExpr):
            return [self.evaluate(e) for e in expr.elements]

        raise MADXInputError(f"Cannot evaluate expression type: {type(expr)}")


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
            raise MADXInputError(
                f"Expected {expected}, got {token.type.name} ({token.value!r})",
                token.line,
            )
        return self._advance()

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
        lexer = MADXLexer(text)
        self.tokens = lexer.tokenize()
        self.pos = 0

        while self._current().type != TokenType.EOF:
            self._parse_statement()
            self._skip_semicolons()

    def _parse_statement(self):
        """Parse a single statement."""
        token = self._current()

        if token.type == TokenType.EOF:
            return

        if token.type == TokenType.SEMICOLON:
            self._advance()
            return

        if token.type != TokenType.IDENTIFIER:
            raise MADXInputError(
                f"Expected identifier, got {token.type.name}", token.line
            )

        name = token.value.lower()

        # Check for label: definition
        if self._peek(1).type == TokenType.COLON:
            self._parse_labeled_definition()
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

        if name == "return":
            # RETURN stops reading the current file (used in CALL'd files)
            self._skip_until_semicolon()
            raise _ReturnStatement()

        if name in (
            "title",
            "option",
            "select",
            "twiss",
            "print",
            "value",
            "show",
            "help",
            "stop",
            "quit",
            "exit",
        ):
            # Skip these commands - just consume until semicolon
            self._skip_until_semicolon()
            return

        # Unknown statement - skip it with a warning
        warnings.warn(
            f"Skipping unknown statement starting with '{name}'", MADXInputWarning
        )
        self._skip_until_semicolon()

    def _skip_until_semicolon(self):
        """Skip tokens until semicolon or EOF."""
        while self._current().type not in (TokenType.SEMICOLON, TokenType.EOF):
            self._advance()
        if self._current().type == TokenType.SEMICOLON:
            self._advance()

    def _parse_labeled_definition(self):
        """Parse a labeled definition (element or line)."""
        name = self._advance().value  # identifier
        self._expect(TokenType.COLON)  # :

        # Get the type/keyword
        keyword = self._expect(TokenType.IDENTIFIER).value.lower()

        if keyword == "line":
            self._parse_line_definition(name)
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
        while self._current().type == TokenType.COMMA:
            self._advance()  # ,

            if self._current().type != TokenType.IDENTIFIER:
                break

            attr_name = self._advance().value.lower()

            # Check for = or :=
            if self._current().type == TokenType.COLON_EQUALS:
                self._advance()
                expr = self._parse_expression()
                attribute_exprs[attr_name] = expr
                attributes[attr_name] = self.context.evaluate(expr)
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

    def _parse_line_elements(self) -> list:
        """Parse the elements inside a LINE definition."""
        elements = []

        while True:
            elem = self._parse_line_element()
            if elem is not None:
                elements.append(elem)

            if self._current().type != TokenType.COMMA:
                break
            self._advance()  # ,

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

        # Check for number * element
        if self._current().type == TokenType.NUMBER:
            num = self._advance().value
            if self._current().type == TokenType.STAR:
                self._advance()
                multiplier = int(num)
            else:
                # Shouldn't happen in valid MAD-X
                raise MADXInputError(
                    "Expected * after number in LINE", self._current().line
                )

        # Now we should have an identifier
        if self._current().type == TokenType.IDENTIFIER:
            elem_name = self._advance().value
            return {
                "name": elem_name,
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

        if integer:
            value = int(self.context.evaluate(expr))
            self.context.set_variable(
                name, value=value, expression=expr, deferred=deferred
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

        value = self.context.evaluate(expr)
        if integer:
            value = int(value)

        self.context.set_variable(
            name, value=value, expression=expr, deferred=deferred, constant=True
        )

    def _parse_beam_command(self):
        """Parse the BEAM command."""
        self._advance()  # skip 'beam'

        while self._current().type == TokenType.COMMA:
            self._advance()

            if self._current().type != TokenType.IDENTIFIER:
                break

            attr_name = self._advance().value.lower()

            if self._current().type == TokenType.EQUALS:
                self._advance()

                if attr_name == "particle":
                    if self._current().type == TokenType.STRING:
                        self.context.beam.particle = self._advance().value.lower()
                    elif self._current().type == TokenType.IDENTIFIER:
                        self.context.beam.particle = self._advance().value.lower()
                    else:
                        raise MADXInputError(
                            "Expected particle name", self._current().line
                        )
                elif attr_name == "energy":
                    expr = self._parse_expression()
                    self.context.beam.energy = self.context.evaluate(expr)
                elif attr_name == "pc":
                    expr = self._parse_expression()
                    self.context.beam.pc = self.context.evaluate(expr)
                elif attr_name == "mass":
                    expr = self._parse_expression()
                    self.context.beam.mass = self.context.evaluate(expr)
                elif attr_name == "charge":
                    expr = self._parse_expression()
                    self.context.beam.charge = self.context.evaluate(expr)
                else:
                    # Unknown attribute, skip the value
                    self._parse_expression()
            elif self._current().type == TokenType.COLON_EQUALS:
                self._advance()
                # Deferred expression - just evaluate immediately for now
                expr = self._parse_expression()
                if attr_name == "energy":
                    self.context.beam.energy = self.context.evaluate(expr)
                elif attr_name == "particle":
                    val = self.context.evaluate(expr)
                    if isinstance(val, str):
                        self.context.beam.particle = val.lower()

        self._expect(TokenType.SEMICOLON)

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
                    raise MADXInputError(
                        "Expected filename after FILE=", self._current().line
                    )
            else:
                # Unknown attribute, skip its value
                if self._current().type == TokenType.EQUALS:
                    self._advance()
                    self._parse_expression()

        self._expect(TokenType.SEMICOLON)

        if filename is None:
            raise MADXInputError("CALL command requires FILE= or FILENAME= attribute")

        # Resolve the filename relative to the current file's directory
        if self._current_file and self._current_file != "<string>":
            base_dir = os.path.dirname(os.path.abspath(self._current_file))
            resolved = os.path.join(base_dir, filename)
        else:
            resolved = filename

        if not os.path.isfile(resolved):
            raise FileNotFoundError(
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

            lexer = MADXLexer(text)
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
            raise MADXInputError(f"Unknown line or element: {line_name}")

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
                warnings.warn(
                    f"Element '{name}' referenced but not defined", MADXInputWarning
                )

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
            warnings.warn(f"Unknown particle type '{particle}'", MADXInputWarning)

        return particle

    def getEtot(self) -> float:
        """Get total energy in GeV."""
        return self.context.beam.energy

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
