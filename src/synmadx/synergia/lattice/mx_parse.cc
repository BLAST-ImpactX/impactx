
#include "mx_parse.h"
#include "mx_tree.h"

#include <any>
#include <cmath>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/optional/optional.hpp>

#include <boost/phoenix/bind.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/phoenix/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>

#include <boost/fusion/include/std_pair.hpp>

#include <synergia/foundation/physical_constants.h>

namespace ascii = ::boost::spirit::ascii;
namespace phx = ::boost::phoenix;
namespace qi = ::boost::spirit::qi;

using ascii::char_;
using ascii::no_case;
using ascii::space;

// using phx::ref;
// using phx::bind;

using qi::_1;
using qi::_2;
using qi::_3;
using qi::_val;
using qi::eol;
using qi::int_;
using qi::lexeme;
using qi::lit;
// using qi::as;  -- boost 1.45+
using qi::locals;

using boost::spirit::qi::real_parser;
using boost::spirit::qi::real_policies;

using namespace qi::labels;
using namespace std;

// ------------------------------------------------------------------

namespace synergia {
    namespace detail {
        double
        pos(double v)
        {
            return v;
        }
        double
        neg(double v)
        {
            return -v;
        }

        double
        add(double v1, double v2)
        {
            return v1 + v2;
        }
        double
        sub(double v1, double v2)
        {
            return v1 - v2;
        }
        double
        mul(double v1, double v2)
        {
            return v1 * v2;
        }
        double
        div(double v1, double v2)
        {
            return v1 / v2;
        }

        bool
        op_l(double l, double r)
        {
            return l < r;
        }
        bool
        op_le(double l, double r)
        {
            return l <= r;
        }
        bool
        op_e(double l, double r)
        {
            return l == r;
        }
        bool
        op_ne(double l, double r)
        {
            return l != r;
        }
        bool
        op_g(double l, double r)
        {
            return l > r;
        }
        bool
        op_ge(double l, double r)
        {
            return l >= r;
        }
    }
}

namespace synergia {
    template <typename Iterator, typename Skip>
    struct expression;

    // parse doc into a statement tree
    template <typename Iterator, typename Skip>
    struct madx_tree_parser;
}

template <typename Iterator, typename Skip>
struct synergia::expression : qi::grammar<Iterator, mx_expr(), Skip> {
    // real policy to allow 'd' and 'D' as exponent
    struct ts_real_policies : boost::spirit::qi::real_policies<double> {
        static bool
        parse_exp(Iterator& first, Iterator const&)
        {
            if (string("eEdD").find(*first) != string::npos) {
                ++first;
                return true;
            }

            return false;
        }
    };

    // symbol table for constants like "pi"
    struct constant_
        : qi::symbols<typename std::iterator_traits<Iterator>::value_type,
                      double> {
        constant_()
        {
            this->add("pi", boost::math::constants::pi<double>())(
                "twopi", boost::math::constants::two_pi<double>())
                // ("degpi"  , boost::math::constants::radian<double>() )   --
                // boost 1.51+
                // ("raddeg" , boost::math::constants::degree<double>() )   --
                // boost 1.51+
                ("degrad", 180.0 / boost::math::constants::pi<double>())(
                    "raddeg", boost::math::constants::pi<double>() / 180.0)(
                    "e", boost::math::constants::e<double>())(
                    "emass", pconstants::me)("pmass", pconstants::mp)(
                    "mumass", pconstants::mmu)("clight", pconstants::c)(
                    "qelect", pconstants::e)("true", 1.0)("false", 0.0);
        }
    } constant;

    // symbol table for unary operator like "- 1.0"
    struct uop_
        : qi::symbols<typename std::iterator_traits<Iterator>::value_type,
                      ufunc> {
        uop_()
        {
            this->add("+", {(ufunc_t)detail::pos, op_tag::pos})(
                "-", {(ufunc_t)detail::neg, op_tag::neg});
        }
    } uop;

    // symbol table for binary operator like "1.0+2.2"
    struct bop1_
        : qi::symbols<typename std::iterator_traits<Iterator>::value_type,
                      bfunc> {
        bop1_()
        {
            this->add("+", {(bfunc_t)detail::add, op_tag::add})(
                "-", {(bfunc_t)detail::sub, op_tag::sub});
        }
    } bop1;

    // symbol table for binary operator like "1.0*2.2"
    struct bop2_
        : qi::symbols<typename std::iterator_traits<Iterator>::value_type,
                      bfunc> {
        bop2_()
        {
            this->add("*", {(bfunc_t)detail::mul, op_tag::mul})(
                "/", {(bfunc_t)detail::div, op_tag::div});
        }
    } bop2;

    // symbol table for binary operator like "1.0^2.2"
    struct bop3_
        : qi::symbols<typename std::iterator_traits<Iterator>::value_type,
                      bfunc> {
        bop3_() { this->add("^", {(bfunc_t)std::pow, op_tag::pow_op}); }
    } bop3;

    // symbol table for unary functions like "abs"
    struct ufunc_
        : qi::symbols<typename std::iterator_traits<Iterator>::value_type,
                      ufunc> {
        ufunc_()
        {
            this->add("abs", {(ufunc_t)std::abs, op_tag::abs})(
                "acos", {(ufunc_t)std::acos, op_tag::acos})(
                "asin", {(ufunc_t)std::asin, op_tag::asin})(
                "atan", {(ufunc_t)std::atan, op_tag::atan})(
                "ceil", {(ufunc_t)std::ceil, op_tag::ceil})(
                "cos", {(ufunc_t)std::cos, op_tag::cos})(
                "cosh", {(ufunc_t)std::cosh, op_tag::cosh})(
                "exp", {(ufunc_t)std::exp, op_tag::exp})(
                "floor", {(ufunc_t)std::floor, op_tag::floor})(
                "log", {(ufunc_t)std::log, op_tag::log})(
                "log10", {(ufunc_t)std::log10, op_tag::log10})(
                "sin", {(ufunc_t)std::sin, op_tag::sin})(
                "sinh", {(ufunc_t)std::sinh, op_tag::sinh})(
                "sqrt", {(ufunc_t)std::sqrt, op_tag::sqrt})(
                "tan", {(ufunc_t)std::tan, op_tag::tan})(
                "tanh", {(ufunc_t)std::tanh, op_tag::tanh});
        }
    } ufunc;

    // symbol table for binary functions like "pow"
    struct bfunc_
        : qi::symbols<typename std::iterator_traits<Iterator>::value_type,
                      bfunc> {
        bfunc_()
        {
            this->add("pow", {(bfunc_t)std::pow, op_tag::pow})(
                "atan2", {(bfunc_t)std::atan2, op_tag::atan2});
        }
    } bfunc;

    qi::rule<Iterator, mx_expr(), Skip> expr;
    qi::rule<Iterator, mx_expr(), Skip> term;
    qi::rule<Iterator, mx_expr(), Skip> factor;
    qi::rule<Iterator, mx_expr(), Skip> primary;
    qi::rule<Iterator, string(), Skip> name;
    qi::rule<Iterator, string_pair_t(), Skip> cmdref;
    qi::real_parser<double, ts_real_policies> real;

    expression() : expression::base_type(expr)
    {
        expr =
            term[_val = phx::construct<nop_t>(_1)] >>
            *(bop1 >>
              term)[_val = phx::construct<bop_t>(_1, _val, _2)]; // '+' and '-'

        term = factor[_val = phx::construct<nop_t>(_1)] >>
               *(bop2 >> factor)[_val = phx::construct<bop_t>(
                                     _1, _val, _2)]; // '*' and '/'

        factor = primary[_val = phx::construct<nop_t>(_1)] >>
                 *(bop3 >>
                   primary)[_val = phx::construct<bop_t>(_1, _val, _2)]; // '^'

        primary = real[_val = _1] |
                  ('(' >> expr >> ')')[_val = phx::construct<nop_t>(_1, true)] |
                  (uop >> primary)[_val = phx::construct<uop_t>(_1, _2)] |
                  (no_case[constant] >> !char_(".a-zA-Z_0-9"))[_val = _1] |
                  (no_case[ufunc] >> '(' >> expr >>
                   ')')[_val = phx::construct<uop_t>(_1, _2)] |
                  (no_case[bfunc] >> '(' >> expr >> ',' >> expr >>
                   ')')[_val = phx::construct<bop_t>(_1, _2, _3)] |
                  (lit("table") >> '(' >> name >> ',' >> name >>
                   ')')[_val = 1.0] // eat table()
                  | (cmdref)[_val = _1] | (name)[_val = _1];

        cmdref = name >> "->" >> name;

        name = lexeme[char_("a-zA-Z_") >> *char_(".a-zA-Z_0-9")];
    }
};

bool
synergia::parse_expression(string const& s, mx_expr& expr)
{
    if (s.empty()) return false;

    typedef string::const_iterator iter_t;
    typedef qi::rule<iter_t> ws_t;

    ws_t whitespace = space | lit('!') >> *(char_ - eol) >> eol |
                      lit("//") >> *(char_ - eol) >> eol;

    expression<iter_t, ws_t> parser;

    iter_t begin = s.begin();
    iter_t const end = s.end();

    bool b =
        qi::phrase_parse(begin, end, parser, whitespace, expr) && begin == end;

    return b;
}

bool
synergia::parse_expression(string const& s, double& result)
{
    mx_expr expr;
    bool b = parse_expression(s, expr);
    if (!b) return false;

    result = boost::apply_visitor(synergia::mx_calculator(), expr);
    return true;
}

////
// semantic actions for madx tree parser
namespace synergia {
    namespace {
        void
        set_logic(mx_logic& logic,
                  mx_expr const& lhs,
                  logic_op_t op,
                  mx_expr const& rhs)
        {
            logic.set(lhs, op, rhs);
        }

        void
        ins_if(mx_if& if_, mx_logic const& logic, mx_tree const& block)
        {
            if_.assign_if(logic, block);
        }

        void
        ins_elseif(mx_if& if_, mx_logic const& logic, mx_tree const& block)
        {
            if_.assign_elseif(logic, block);
        }

        void
        ins_else(mx_if& if_, mx_tree const& block)
        {
            if_.assign_else(block);
        }

        void
        ins_while(mx_while& while_, mx_logic const& logic, mx_tree const& block)
        {
            while_.assign(logic, block);
        }

        void
        ins_seq_member(mx_line_seq& seq,
                       boost::optional<char> m,
                       boost::optional<int> o,
                       mx_line_member const& member)
        {
            int op = m ? -1 : 1;
            if (o) op *= o.value();
            seq.insert_member(op, member);
        }

        void
        set_attr(mx_attr& attr,
                 string const& name,
                 boost::optional<char> c,
                 std::any const& v)
        {
            if (c)
                attr.set_lazy_attr(name, v);
            else
                attr.set_attr(name, v);
        }

        void
        set_flag_attr(mx_attr& attr,
                      boost::optional<char> c,
                      string const& name)
        {
            if (c)
                attr.set_attr(name, std::any(mx_expr(0.0)));
            else
                attr.set_attr(name, std::any(mx_expr(1.0)));
        }

        void
        set_cmd_label(mx_command& cmd, string const& label)
        {
            cmd.set_label(label);
        }

        void
        set_cmd_keyword(mx_command& cmd, mx_keyword const& keyword)
        {
            cmd.set_keyword(keyword);
        }

        void
        ins_cmd_attr(mx_command& cmd, mx_attr const& attr)
        {
            cmd.ins_attr(attr);
        }
    }
}

template <typename Iterator, typename Skip>
struct synergia::madx_tree_parser : qi::grammar<Iterator, mx_tree(), Skip> {
    // keywords
    struct logic_op_ : boost::spirit::qi::symbols<
                           typename std::iterator_traits<Iterator>::value_type,
                           logic_op_t> {
        logic_op_()
        {
            this->add("<", (logic_op_t)detail::op_l)("<=",
                                                     (logic_op_t)detail::op_le)(
                "==", (logic_op_t)detail::op_e)("!=",
                                                (logic_op_t)detail::op_ne)(
                ">", (logic_op_t)detail::op_g)(">=", (logic_op_t)detail::op_ge);
        }
    } logic_op;

    struct particle_keywords_
        : boost::spirit::qi::symbols<
              typename std::iterator_traits<Iterator>::value_type,
              mx_keyword> {
        particle_keywords_()
        {
            this->add("proton", mx_keyword("proton", MX_KW_PARTICLE))(
                "prot", mx_keyword("proton", MX_KW_PARTICLE))(
                "electron", mx_keyword("electron", MX_KW_PARTICLE))(
                "positron", mx_keyword("positron", MX_KW_PARTICLE))(
                "anti-proton", mx_keyword("anti-proton", MX_KW_PARTICLE))(
                "posmuon", mx_keyword("posmuon", MX_KW_PARTICLE))(
                "negmuon", mx_keyword("negmuon", MX_KW_PARTICLE));
        }
    } particle_keywords;

    struct mp_type_keywords_
        : boost::spirit::qi::symbols<
              typename std::iterator_traits<Iterator>::value_type,
              mx_keyword> {
        mp_type_keywords_()
        {
            this->add("octpn", mx_keyword("octpn", MX_KW_MP_TYPE))(
                "wgl", mx_keyword("wgl", MX_KW_MP_TYPE))(
                "sk", mx_keyword("sk", MX_KW_MP_TYPE));
        }
    } mp_type_keywords;

    struct element_keywords_
        : boost::spirit::qi::symbols<
              typename std::iterator_traits<Iterator>::value_type,
              mx_keyword> {
        element_keywords_()
        {
            this->add("drift", mx_keyword("drift", MX_KW_ELEMENT))(
                "rbend", mx_keyword("rbend", MX_KW_ELEMENT))(
                "sbend", mx_keyword("sbend", MX_KW_ELEMENT))(
                "dipedge", mx_keyword("dipedge", MX_KW_ELEMENT))(
                "quadrupole", mx_keyword("quadrupole", MX_KW_ELEMENT))(
                "sextupole", mx_keyword("sextupole", MX_KW_ELEMENT))(
                "octupole", mx_keyword("octupole", MX_KW_ELEMENT))(
                "multipole", mx_keyword("multipole", MX_KW_ELEMENT))(
                "solenoid", mx_keyword("solenoid", MX_KW_ELEMENT))(
                "nllens", mx_keyword("nllens", MX_KW_ELEMENT))(
                "hkicker", mx_keyword("hkicker", MX_KW_ELEMENT))(
                "vkicker", mx_keyword("vkicker", MX_KW_ELEMENT))(
                "tkicker", mx_keyword("tkicker", MX_KW_ELEMENT))(
                "kicker", mx_keyword("kicker", MX_KW_ELEMENT))(
                "rfcavity", mx_keyword("rfcavity", MX_KW_ELEMENT))(
                "rfmultipole", mx_keyword("rfmultipole", MX_KW_ELEMENT))(
                "crabcavity", mx_keyword("crabcavity", MX_KW_ELEMENT))(
                "elseparator", mx_keyword("elseparator", MX_KW_ELEMENT))(
                "hmonitor", mx_keyword("hmonitor", MX_KW_ELEMENT))(
                "vmonitor", mx_keyword("vmonitor", MX_KW_ELEMENT))(
                "monitor", mx_keyword("monitor", MX_KW_ELEMENT))(
                "instrument", mx_keyword("instrument", MX_KW_ELEMENT))(
                "rcollimator", mx_keyword("rcollimator", MX_KW_ELEMENT))(
                "ecollimator", mx_keyword("ecollimator", MX_KW_ELEMENT))(
                "yrotation", mx_keyword("yrotation", MX_KW_ELEMENT))(
                "srotation", mx_keyword("srotation", MX_KW_ELEMENT))(
                "beambeam", mx_keyword("beambeam", MX_KW_ELEMENT))(
                "matrix", mx_keyword("matrix", MX_KW_ELEMENT))(
                "marker", mx_keyword("marker", MX_KW_ELEMENT))(
                "constfoc", mx_keyword("constfoc", MX_KW_ELEMENT))(
                "elens", mx_keyword("elens", MX_KW_ELEMENT))(
                "foil", mx_keyword("foil", MX_KW_ELEMENT))(
                "matrix", mx_keyword("matrix", MX_KW_ELEMENT));
        }
    } element_keywords;

    struct command_keywords_
        : boost::spirit::qi::symbols<
              typename std::iterator_traits<Iterator>::value_type,
              mx_keyword> {
        command_keywords_()
        {
            this->add
                // general
                ("assign", mx_keyword("assign", MX_KW_COMMAND))(
                    "beam", mx_keyword("beam", MX_KW_COMMAND))(
                    "call", mx_keyword("call", MX_KW_COMMAND))(
                    "coguess", mx_keyword("coguess", MX_KW_COMMAND))(
                    "create", mx_keyword("create", MX_KW_COMMAND))(
                    "dumpsequ", mx_keyword("dumpsequ", MX_KW_COMMAND))(
                    "exec", mx_keyword("exec", MX_KW_COMMAND))(
                    "exit", mx_keyword("exit", MX_KW_COMMAND))(
                    "fill", mx_keyword("fill", MX_KW_COMMAND))(
                    "help", mx_keyword("help", MX_KW_COMMAND))(
                    "option", mx_keyword("option", MX_KW_COMMAND))(
                    "print", mx_keyword("print", MX_KW_COMMAND))(
                    "quit", mx_keyword("quit", MX_KW_COMMAND))(
                    "readtable", mx_keyword("readtable", MX_KW_COMMAND))(
                    "return", mx_keyword("return", MX_KW_COMMAND))(
                    "save", mx_keyword("save", MX_KW_COMMAND))(
                    "savebeta", mx_keyword("savebeta", MX_KW_COMMAND))(
                    "select", mx_keyword("select", MX_KW_COMMAND))(
                    "set", mx_keyword("set", MX_KW_COMMAND))(
                    "show", mx_keyword("show", MX_KW_COMMAND))(
                    "stop", mx_keyword("stop", MX_KW_COMMAND))(
                    "system", mx_keyword("system", MX_KW_COMMAND))(
                    "tabstring", mx_keyword("tabstring", MX_KW_COMMAND))(
                    "title", mx_keyword("title", MX_KW_COMMAND))(
                    "use", mx_keyword("use", MX_KW_COMMAND))(
                    "value", mx_keyword("value", MX_KW_COMMAND))(
                    "write", mx_keyword("write", MX_KW_COMMAND))
                // beam specification
                ("beam", mx_keyword("beam", MX_KW_COMMAND))(
                    "resbeam", mx_keyword("resbeam", MX_KW_COMMAND))
                // plot
                ("plot", mx_keyword("plot", MX_KW_COMMAND))(
                    "resplot", mx_keyword("resplot", MX_KW_COMMAND))(
                    "setplot", mx_keyword("setplot", MX_KW_COMMAND))
                // sequence editing
                ("seqedit", mx_keyword("seqedit", MX_KW_COMMAND))(
                    "flatten", mx_keyword("flatten", MX_KW_COMMAND))(
                    "install", mx_keyword("install", MX_KW_COMMAND))(
                    "move", mx_keyword("move", MX_KW_COMMAND))(
                    "remove", mx_keyword("remove", MX_KW_COMMAND))(
                    "cycle", mx_keyword("cycle", MX_KW_COMMAND))(
                    "reflect", mx_keyword("reflect", MX_KW_COMMAND))(
                    "endedit", mx_keyword("endedit", MX_KW_COMMAND))
                // other commands
                ("twiss", mx_keyword("twiss", MX_KW_COMMAND))
                // build sequence ( not present in manual )
                ("sequence", mx_keyword("sequence", MX_KW_COMMAND))(
                    "endsequence", mx_keyword("endsequence", MX_KW_COMMAND));
        }
    } command_keywords;

    // rules
    qi::rule<Iterator, mx_tree(), Skip> doc;
    qi::rule<Iterator, mx_tree(), Skip> block;
    qi::rule<Iterator, mx_statement(), Skip> statement;
    qi::rule<Iterator, mx_command(), Skip> variable;
    qi::rule<Iterator, mx_command(), locals<mx_cmd_type>, Skip> cmd;
    qi::rule<Iterator, mx_command(), Skip> command;
    qi::rule<Iterator, mx_line(), Skip> line;
    qi::rule<Iterator, mx_line_seq(), Skip> line_seq;
    qi::rule<Iterator, mx_line_member(), Skip> line_member;
    qi::rule<Iterator, mx_keyword(), Skip> ref;
    qi::rule<Iterator, mx_if(), Skip> if_flow;
    qi::rule<Iterator, mx_while(), Skip> while_flow;
    qi::rule<Iterator, mx_logic(), Skip> logic;
    qi::rule<Iterator, mx_attr(), Skip> attr;
    qi::rule<Iterator, mx_attr(), Skip> flag_attr;
    qi::rule<Iterator, mx_attr(), Skip> cmd_attr;

    qi::rule<Iterator, string(), Skip> name;
    qi::rule<Iterator, string(), Skip> dblq_str;
    qi::rule<Iterator, string(), Skip> snglq_str;

    expression<Iterator, Skip> expr;
    qi::rule<Iterator, mx_exprs(), Skip> array;
    qi::rule<Iterator, std::any(), Skip> value;

    qi::rule<Iterator, mx_statements_t(), Skip> statements;

    madx_tree_parser() : madx_tree_parser::base_type(doc), expr()
    {
        doc =
            // as<mx_statements_t>()[*statement]  -- boost 1.45+
            statements[_val = _1];

        statements = *statement;

        statement =
            // if_flow | while_flow | command  -- boost 1.45+
            if_flow[_val = _1] | while_flow[_val = _1] | command[_val = _1] |
            line[_val = _1] | ';' // empty statement
            ;

        block =
            // as<mx_statements_t>()['{' >> *statement >> '}']  -- boost 1.45+
            '{' >> statements[_val = _1] >> '}';

        logic =
            //'(' >> *(char_ - char_(')')) >> ')'
            ('(' >> expr >> logic_op >> expr >>
             ')')[phx::bind(&set_logic, _val, _1, _2, _3)];

        if_flow = (no_case["if"] >> logic >>
                   block)[phx::bind(&ins_if, _val, _1, _2)] >>
                  *(no_case["elseif"] >> logic >>
                    block)[phx::bind(&ins_elseif, _val, _1, _2)] >>
                  -(no_case["else"] >> block)[phx::bind(&ins_else, _val, _1)];

        while_flow = (no_case["while"] >> logic >>
                      block)[phx::bind(&ins_while, _val, _1, _2)];

        array =
            //'{' >> expr % ',' >> '}'  -- boost 1.45+
            '{' >> expr[phx::push_back(_val, _1)] % ',' >> '}';

        value =
            // dblq_str | snglq_str | no_case[particle_keywords] | expr | array
            // -- boost 1.45+
            dblq_str[_val = _1] | snglq_str[_val = _1] | expr[_val = _1] |
            array[_val = _1];

        attr = (no_case[qi::string("type")] // special attr 'type'
                >> -char_(':') >> '=' >>
                (name | dblq_str |
                 snglq_str))[phx::bind(&set_attr, _val, _1, _2, _3)] |
               (no_case[qi::string("particle")] // special attr 'particle'
                >> -char_(':') >> '=' >> no_case[particle_keywords])[phx::bind(
                   &set_attr, _val, _1, _2, _3)] |
               (name // generic attr
                >> -char_(':') >> '=' >>
                value)[phx::bind(&set_attr, _val, _1, _2, _3)];

        flag_attr = // flag attributes, attributes with values only, etc.
            (-char_("-") >> name)[phx::bind(&set_flag_attr, _val, _1, _2)] |
            (dblq_str)[phx::bind(&set_attr,
                                 _val,
                                 string("default"),
                                 boost::optional<char>(),
                                 _1)] |
            (snglq_str)[phx::bind(&set_attr,
                                  _val,
                                  string("default"),
                                  boost::optional<char>(),
                                  _1)] |
            (expr)[phx::bind(&set_attr,
                             _val,
                             string("default"),
                             boost::optional<char>(),
                             _1)];

        cmd_attr = attr | flag_attr;

        variable = attr[phx::bind(&ins_cmd_attr, _val, _1)];

        cmd =
            -(name >> ':')[phx::bind(&set_cmd_label, _val, _1)] // label
            >> (no_case[element_keywords] | no_case[command_keywords] |
                ref)[phx::bind(&set_cmd_keyword, _val, _1)] // keyword
            >> -(lit(',')) >>
            *(cmd_attr[phx::bind(&ins_cmd_attr, _val, _1)] % ',') // attributes
            ;

        command = (variable | cmd) >> ';';

        line = (name >> ':' >> no_case["line"] >> '=' >> line_seq >>
                ';')[_val = phx::construct<mx_line>(_1, _2)];

        line_seq =
            lit('(') >> (*lit(',')) // allows 0 or multiple ',' in front
            >> ((-char_('-') >> -(int_ >> '*') >>
                 line_member)[phx::bind(&ins_seq_member, _val, _1, _2, _3)] %
                (*lit(',')) // allows the ',' in the line_seq to be optional
                            // or multiple occurance
                ) >>
            (*lit(',')) // allows for ',' in the tail
            >> lit(')');

        line_member = name[_val = _1] | line_seq[_val = _1];

        ref = name[_val = phx::construct<mx_keyword>(_1, MX_KW_ELEMENT_REF)];

        name = lexeme[char_("a-zA-Z_") >> *char_(".a-zA-Z_0-9")];

        dblq_str = lexeme['"' >> +(char_ - '"') >> '"'];

        snglq_str = lexeme['\'' >> +(char_ - '\'') >> '\''];
    }
};

bool
synergia::parse_int_madx(string const& s, mx_tree& doc, string const& f)
{
    if (s.empty()) return true;

    typedef string::const_iterator iter_t;
    typedef qi::rule<iter_t> ws_t;

    ws_t whitespace = space | lit('!') >> *(char_ - eol) >> eol |
                      lit('&') >> *(char_ - eol) >> eol |
                      lit("//") >> *(char_ - eol) >> eol |
                      lit("/*") >> *(char_ - "*/") >> "*/";

    madx_tree_parser<iter_t, ws_t> parser;

    iter_t begin = s.begin();
    iter_t const end = s.end();

    bool b =
        qi::phrase_parse(begin, end, parser, whitespace, doc) && begin == end;

    if (!b) {
        iter_t e = distance(begin, end) < 100 ? end : begin + 100;
        string r = string(begin, e);

        throw runtime_error("MadX: parsing " + f + " failed at \n>>>\n" + r +
                            "\n<<<");
    }

    return b;
}

bool
synergia::parse_madx(string const& str, MadX& mx, string const& f)
{
    // first parse the madx doc into a statement tree
    string s = str + "\n";
    mx_tree tree;
    parse_int_madx(s, tree, f);

    // print for debug purpose
    // tree.print();

    // interpret the syntax tree into a MadX object
    tree.interpret(mx);

    return true;
}

// helper
namespace {
    void
    read_from_file(string const& fname, string& str)
    {
        ifstream file;
        file.open(fname.c_str());

        if (!file.is_open())
            throw runtime_error("Failed to open file " + fname +
                                " for parsing");

        file.seekg(0, std::ios::end);
        str.reserve(file.tellg());
        file.seekg(0, std::ios::beg);

        str.assign((istreambuf_iterator<char>(file)),
                   istreambuf_iterator<char>());

        str.append("\n");

        file.close();
    }
}

bool
synergia::parse_int_madx_file(string const& fname, mx_tree& doc)
{
    string str;
    read_from_file(fname, str);
    return parse_int_madx(str, doc, fname);
}

bool
synergia::parse_madx_file(string const& fname, MadX& mx)
{
    string str;
    read_from_file(fname, str);
    return parse_madx(str, mx, fname);
}
