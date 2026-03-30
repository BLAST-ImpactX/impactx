#!/usr/bin/env python

from __future__ import division
from pyparsing import Word, alphas, ParseException, Literal, CaselessLiteral, \
Combine, Optional, nums, Or, Forward, ZeroOrMore, StringEnd, alphanums, \
restOfLine, empty, delimitedList, LineEnd, Group, QuotedString, dblQuotedString, \
removeQuotes, sglQuotedString, OneOrMore, ParserElement

import sys
import re
import math
import operator
import inspect
import types

class Printable:
    def _enum_to_string(self, enum, val):
        for member in inspect.getmembers(enum):
            name = member[0]
            if getattr(enum, name) == val:
                return name
        return str(val)

    def __repr__(self):
        retval = '<'
        first = True
        for member in inspect.getmembers(self):
            name = member[0]
            if name[0] != '_':
                if first:
                    first = False
                else:
                    retval += ','
                retval += name + ':'
                val = getattr(self, name)
                str_val = None
                if hasattr(self, '_enum_types'):
                    if name in self._enum_types:
                        str_val = self._enum_to_string(self._enum_types[name], val)
                if not str_val:
                    str_val = str(val)
                retval += str_val
        retval += '>'
        return retval

class Stack_type:
    floatnumber = 1
    ident = 2
    operator = 3
    function = 4
    unary_minus = 5
    str = 6
    subscript_ident = 7
    range = 8

class Stack_item(Printable):
    def __init__(self, type, value):
        self.type = type
        self._enum_types = {'type' : Stack_item}
        self.value = value

class Subscript_ident(Printable):
    def __init__(self, ident, subscript):
        self.ident = ident
        self.subscript = subscript

class Expression_parser:
    def __init__(self, uninitialized_warn=True,
                 uninitialized_exception=False):
        self.uninitialized_warn = uninitialized_warn
        self.uninitialized_exception = uninitialized_exception
        self.bnf = self._construct_bnf()
        self.operators = {"+" : operator.add,
                          "-" : operator.sub,
                          "*" : operator.mul,
                          "/" : operator.truediv,
                          "^" : operator.pow}
        # self.functions values are of the form [function, number of arguments]
        self.functions = {"sqrt" : [math.sqrt, 1],
                           "log" : [math.log, 1],
                           "exp" : [math.exp, 1],
                           "sin" : [math.sin, 1],
                           "cos" : [math.cos, 1],
                           "tan" : [math.tan, 1],
                           "asin" : [math.asin, 1],
                           "abs" : [abs, 1],
                           "max" : [max, 2],
                           "min" : [min, 2]}
        self.constants = {"pi" : math.pi,
                          "twopi" : 2.0 * math.pi,
                          "degrad" : 180.0 / math.pi,
                          "raddeg" : math.pi / 180.0,
                          "e" : math.e,
                          "emass" : 0.51099906e-3,
                          "pmass" : 0.93827231,
                          "clight" : 2.99792458e8}
        self.reset()

    def _construct_bnf(self):
        point = Literal(".")
        e = CaselessLiteral('e') | CaselessLiteral('d')
        plus = Literal("+")
        minus = Literal("-")
        mult = Literal("*")
        div = Literal("/")
        lpar = Literal("(").suppress()
        rpar = Literal(")").suppress()
        lbrack = Literal("[").suppress()
        rbrack = Literal("]").suppress()
        addop = plus | minus
        multop = mult | div
        expop = Literal("^")
        plusorminus = plus | minus
        number = Word(nums)
        integer = Combine(Optional(plusorminus) + number)
        self.integer = integer
        floatnumber = (Combine(((integer + \
                                 Optional(point + Optional(number))) |  \
                               (point + number)) + \
                               Optional(e + integer)) | \
                       Optional(plusorminus) + Combine(point + number) + \
                        Optional(e + integer))
        self.floatnumber = floatnumber
        ident = Word(alphas, alphanums + '_' + '.' + "'")
        self.ident = ident
        subscript_ident = Group(ident + lbrack + ident + rbrack)
        hash = Literal("#")
        hashword = hash + Word(alphanums)
        range = Combine(hashword + div + hashword)
        self.range = range

        expr = Forward()
        atom = (Optional(minus) +
                ((floatnumber)("floatnumber").setParseAction(self._push_floatnumber) | \
                 (range)("range").setParseAction(self._push_range) | \
                 (ident + lpar + delimitedList(expr) + rpar)("function").setParseAction(self._push_function) | \
                 (subscript_ident)("subscript_ident").setParseAction(self._push_subscript_ident) | \
                 (ident)("ident").setParseAction(self._push_ident) | \
                (lpar + expr.suppress() + rpar))).setParseAction(self._push_uminus)

        # by defining exponentiation as "atom [ ^ factor ]..." instead of
        # "atom [ ^ atom ]...", we get right-to-left exponents, instead of
        # left-to-right, that is, 2^3^2 = 2^(3^2), not (2^3)^2.
        factor = Forward()
        factor << atom + ZeroOrMore((expop + factor).setParseAction(self._push_operator))

        term = factor + ZeroOrMore((multop + factor).setParseAction(self._push_operator))
        expr << term + ZeroOrMore((addop + term).setParseAction(self._push_operator))

        return expr

    def _push_floatnumber(self, strg, loc, toks):
        numstr = toks[0]
        numstr = numstr.replace('d', 'e')
        numstr = numstr.replace('D', 'e')
        self.stack.append(Stack_item(Stack_type.floatnumber, float(numstr)))

    def _push_ident(self, strg, loc, toks):
        self.stack.append(Stack_item(Stack_type.ident, toks[0].lower()))

    def _push_subscript_ident(self, strg, loc, toks):
        ident = toks[0][0].lower()
        subscript = toks[0][1].lower()
        self.stack.append(Stack_item(Stack_type.subscript_ident,
                                     Subscript_ident(ident, subscript)))

    def _push_function(self, strg, loc, toks):
        self.stack.append(Stack_item(Stack_type.function, toks[0].lower()))

    def _push_uminus(self, strg, loc, toks):
        if toks and toks[0] == '-':
            self.stack.append(Stack_item(Stack_type.unary_minus, None))

    def _push_operator(self, strg, loc, toks):
        self.stack.append(Stack_item(Stack_type.operator, toks[0]))

    def _push_range(self, strg, loc, toks):
        self.stack.append(Stack_item(Stack_type.range, toks[0]))

    def evaluate_stack(self, s, variables={}, labels={}, constants=None):
        if constants == None:
            constants = self.constants
        op = s.pop()
        if op.type == Stack_type.unary_minus:
            return -self.evaluate_stack(s, variables, labels, constants)
        if op.type == Stack_type.floatnumber:
            return op.value
        if op.type == Stack_type.operator:
            op2 = self.evaluate_stack(s, variables, labels, constants)
            op1 = self.evaluate_stack(s, variables, labels, constants)
            return self.operators[op.value](op1, op2)
        elif op.type == Stack_type.function:
            ops = []
            for i in range(0, self.functions[op.value][1]):
                ops.append(self.evaluate_stack(s, variables, labels, constants))
            return self.functions[op.value][0](*ops)
        elif op.type == Stack_type.ident:
            if op.value in constants:
                return constants[op.value]
            elif op.value in variables:
                return variables[op.value]
            else:
                if self.uninitialized_exception:
                    raise ParseException('variable "%s" unitialized' % op.value)
                if self.uninitialized_warn:
                    print(
                     'warning: variable "%s" uninitialized, treating as 0.0' \
                     % op.value)
                return 0.0
        elif op.type == Stack_type.str:
            return op.value
        elif op.type == Stack_type.subscript_ident:
            retval = None
            if op.value.ident in labels:
                attributes = labels[op.value.ident].attributes
                if op.value.subscript in attributes:
                    retval = attributes[op.value.subscript]
            if retval == None:
                raise RuntimeError("Unknown subscript %s[%s]" % \
                    (op.value.ident, op.value.subscript))
            return retval
        elif op.type == Stack_type.range:
            retval = op.value
        else:
            raise RuntimeError("Unhandled expression stack item %s" % op)

    def reset(self):
        self.stack = []

    def parse(self, text):
        self.reset()
        result = self.bnf.parseString(text)
        return self.stack

class Command(Printable):
    def __init__(self, name, attributes):
        self.name = name
        self.attributes = attributes

class Mad8_parser:
    def __init__(self):
        ParserElement.setDefaultWhitespaceChars(" \t\r")
        self.expression_parser = Expression_parser()
        self.bnf = self._construct_bnf()
        self.variables = {}
        self._attributes = {}
        self.commands = []
        self.labels = {}
        self.lines = {}
        self.no_eval_attributes = ['particle', 'type', 'filename']

    def _construct_bnf(self):
        colon = Literal(':')
        equals = Literal('=')
        semicolon = Literal(';')
        bang = Literal('!')
        minus = Literal('-')
        times = Literal('*')
        ampersand = Literal('&')
        lpar = Literal('(').suppress()
        rpar = Literal(')').suppress()

        ident = self.expression_parser.ident
        expr = self.expression_parser.bnf
        integer = self.expression_parser.integer
        floatnumber = self.expression_parser.floatnumber
        range = self.expression_parser.range

        var_assign = (equals | colon + equals |
                      colon + CaselessLiteral('constant') + equals).suppress()
        attr_assign = Literal('=')
        str = dblQuotedString.setParseAction(removeQuotes) | \
            sglQuotedString.setParseAction(removeQuotes)
        str.setParseAction(self._handle_str)
        attr_value = (str | expr)
        attr = ((ident | floatnumber)  + Optional(attr_assign + attr_value) | \
                Group(Literal('-') + ident))
        attr.setParseAction(self._handle_attr)
        attr_delim = Literal(',').suppress()
        command = ident + Optional(attr_delim + \
                                   ((ident + range) | \
                                    range | \
                                    OneOrMore(delimitedList(attr))))
        command.setParseAction(self._handle_command)
#signedmodifierdef = Group(Optional(Literal('-')) + ident) + Optional(Literal("=") + expr)
#signedmodifierdef.setParseAction(handleModifier)
        var_assign = (ident + var_assign + expr)
        var_assign.setParseAction(self._handle_var_assign)
        labeled_command = (ident + colon + ident +
                          Optional(attr_delim + OneOrMore(delimitedList(attr))))
        labeled_command.setParseAction(self._handle_labeled_command)

        multiple_ident = (Optional(minus) + Optional(Group(integer + times)) +
                          (ident | lpar + Group(delimitedList(ident)) + rpar))
        line = (CaselessLiteral("line") + attr_assign + lpar +
                Group(multiple_ident +
                      ZeroOrMore(Optional(attr_delim) + multiple_ident)) +
                rpar)
        labeled_line = (ident + colon + line)
        labeled_line.setParseAction(self._handle_labeled_line)

        empty.setWhitespaceChars(" \t")
        entry = ((var_assign |
                 labeled_line |
                 labeled_command |
                 command |
                 empty) + (semicolon | LineEnd().suppress()))

        bnf = ZeroOrMore(entry) + StringEnd()
        comment = bang + restOfLine
        bnf.ignore(comment)
        continuation = ampersand + restOfLine + LineEnd()
        bnf.ignore(continuation)

        return bnf

    def _handle_var_assign(self, str, loc, toks):
        var = toks[0].lower()
        stack = self.expression_parser.stack
        value = self.expression_parser.evaluate_stack(stack, self.variables,
                                                      self.labels)
        self.variables[var] = value

    def _handle_attr(self, str, loc, toks):
        if hasattr(toks[0], 'asList'):
            attribute = "".join(toks[0]).lower()
        else:
            attribute = toks[0].lower()
#        if attribute == 'range':
#            value = toks[2] + toks[3] + toks[4]
        if len(toks) > 1:
            stack = self.expression_parser.stack
            if attribute in self.no_eval_attributes:
                value = stack.pop().value
            else:
                value = self.expression_parser.evaluate_stack(stack,
                                                              self.variables,
                                                              self.labels)
        else:
            value = None
        self._attributes[attribute] = value

    def _handle_str(self, str, loc, toks):
        self.expression_parser.stack.append(Stack_item(Stack_type.str,
                                                       toks[0]))

    def _handle_command(self, str, loc, toks):
        name = toks[0].lower()
        self.commands.append(Command(name, self._attributes))
        self._attributes = {}

    def _handle_labeled_command(self, str, loc, toks):
        label = toks[0].lower()
        name = toks[2].lower()
        self.labels[label] = Command(name, self._attributes)
        self._attributes = {}

    def _downcase_nested_list(self, lst):
        retval = []
        for elem in lst:
            if hasattr(elem, 'asList'):
                is_multiplier = False
                if len(elem) == 2:
                    if elem[1] == '*':
                        retval.append(elem[0] + '*')
                        is_multiplier = True
                if not is_multiplier:
                    retval.append(self._downcase_nested_list(elem))
            else:
                retval.append(elem.lower())
        return retval

    def _handle_labeled_line(self, str, loc, toks):
        name = toks[0].lower()
        elements = self._downcase_nested_list(toks[4])
        self.lines[name] = elements

    def parse(self, text):
        self.expression_parser.reset()
        result = self.bnf.parseString(text)

if __name__ == '__main__':
    filename = sys.argv[1]
    mp = Mad8_parser()
    mp.parse(open(filename, 'r').read())
    print("variables =", mp.variables)
    print("labels =", mp.labels)
    print("lines =", mp.lines)
    print("commands =", mp.commands)
