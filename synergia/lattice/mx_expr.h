#ifndef MX_EXPR_H
#define MX_EXPR_H

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include <boost/variant.hpp>

#include "synergia/utils/cereal.h"

namespace synergia {
    enum class op_tag {
        // uop operators
        pos, // +
        neg, // -

        // bop operators
        add,    // +
        sub,    // -
        mul,    // *
        div,    // /
        pow_op, // ^

        // uop functions
        abs,  // abs()
        acos, // acos()
        asin,
        atan,
        ceil,
        cos,
        cosh,
        exp,
        floor,
        log,
        log10,
        sin,
        sinh,
        sqrt,
        tan,
        tanh,

        // bop functions
        pow,   // pow(lhs, rhs)
        atan2, // atan2(lhs, rhs)
    };

    typedef double (*nfunc_t)();
    typedef double (*ufunc_t)(double);
    typedef double (*bfunc_t)(double, double);

    struct nfunc {
        nfunc_t op;
        op_tag tag;
    };
    struct ufunc {
        ufunc_t op;
        op_tag tag;
    };
    struct bfunc {
        bfunc_t op;
        op_tag tag;
    };

    typedef std::pair<std::string, std::string> string_pair_t;

    struct nop_t;
    struct uop_t;
    struct bop_t;

    class MadX;
    class mx_calculator;

    // whether the expr contains a ref string
    class mx_expr_ref_checker;

    // convert the expression to a string
    class mx_expr_writer;

    typedef boost::variant<double,
                           std::string,
                           string_pair_t,
                           boost::recursive_wrapper<nop_t>,
                           boost::recursive_wrapper<uop_t>,
                           boost::recursive_wrapper<bop_t>>
        mx_expr;

    typedef std::vector<mx_expr> mx_exprs;

    // evaluate the expression
    // double mx_eval(mx_expr const & expr, double def = mx_calculator::nan);
    // double mx_eval(mx_expr const & expr, MadX const & mx, double def =
    // mx_calculator::nan);

    // whether can be evaluated to a number
    bool mx_expr_is_number(mx_expr const& expr);

    // retrieve the ref as a string
    std::string mx_expr_refstr(mx_expr const& expr);

    // retrieve the expression as an string
    std::string mx_expr_str(mx_expr const& expr);

    // parse expression
    bool parse_expr(std::string const& s, mx_expr& expr);
}

// serialization for mx_expr type
namespace cereal {
    template <class AR>
    void
    serialize(AR& ar, synergia::mx_expr& expr)
    {
        if (AR::is_saving::value) {
            std::string str = mx_expr_str(expr);
            ar(str);
        } else {
            std::string str;
            ar(str);

            parse_expr(str, expr);
        }
    }
}

struct synergia::nop_t {
    nop_t(mx_expr const& e, bool primary = false) : expr(e), primary(primary) {}

    mx_expr expr;
    bool primary;
};

struct synergia::uop_t {
    uop_t(ufunc f, mx_expr const& param) : param(param), func(f) {}

    mx_expr param;
    ufunc func;
};

struct synergia::bop_t {
    bop_t(bfunc f, mx_expr const& lhs, mx_expr const& rhs)
        : func(f), lhs(lhs), rhs(rhs)
    {}

    bfunc func;
    mx_expr lhs;
    mx_expr rhs;
};

class synergia::mx_calculator : public boost::static_visitor<double> {
  public:
    mx_calculator() : mx(NULL), def(nan) {}

    explicit mx_calculator(double def) : mx(NULL), def(def) {}

    explicit mx_calculator(MadX const& mx) : mx(&mx), def(nan) {}

    mx_calculator(MadX const& mx, double def) : mx(&mx), def(def) {}

    double operator()(double val) const;
    double operator()(std::string const& ref) const;
    double operator()(string_pair_t const& ref) const;
    double operator()(nop_t const& n) const;
    double operator()(uop_t const& u) const;
    double operator()(bop_t const& b) const;

  public:
    static double nan;

  private:
    MadX const* mx;
    double def;
};

class synergia::mx_expr_ref_checker : public boost::static_visitor<bool> {
  public:
    bool operator()(double val) const;
    bool operator()(std::string const& ref) const;
    bool operator()(string_pair_t const& ref) const;
    bool operator()(nop_t const& n) const;
    bool operator()(uop_t const& u) const;
    bool operator()(bop_t const& b) const;
};

class synergia::mx_expr_writer : public boost::static_visitor<std::string> {
  public:
    std::string operator()(double val) const;
    std::string operator()(std::string const& ref) const;
    std::string operator()(string_pair_t const& ref) const;
    std::string operator()(nop_t const& n) const;
    std::string operator()(uop_t const& u) const;
    std::string operator()(bop_t const& b) const;
};

namespace synergia {
    // evaluate the expression
    double mx_eval(mx_expr const& expr, double def = mx_calculator::nan);
    double mx_eval(mx_expr const& expr,
                   MadX const& mx,
                   double def = mx_calculator::nan);
}

#endif
