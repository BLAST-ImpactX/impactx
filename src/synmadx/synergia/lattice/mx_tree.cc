
#include "mx_tree.h"
#include "mx_parse.h"

#include <any>
#include <cassert>
#include <iostream>
#include <stdexcept>

#include <boost/algorithm/string.hpp>

#include "synergia/foundation/four_momentum.h"
#include "synergia/foundation/physical_constants.h"

using namespace synergia;
using namespace std;

using boost::get;

// helper
namespace synergia { namespace {
  // deduct an expression and replace it with a number if all through
  mx_expr
  simplify(mx_expr const& e, MadX const& mx)
  {
    try {
      double r = boost::apply_visitor(mx_calculator(mx), e);
      return mx_expr(r);
    }
    catch (...) {
      return e;
    }
  }

  template <typename T>
  void
  insert_attr(T& t, mx_attr const& attr, MadX const& mx)
  {
    if (attr.type() == MX_ATTR_STRING) {
      t.insert_attribute(attr.name(), any_cast<string>(attr.value()));
    } else if (attr.type() == MX_ATTR_PREDEFINED) {
      t.insert_attribute(attr.name(), any_cast<mx_keyword>(attr.value()).name);
    } else if (attr.type() == MX_ATTR_NUMBER) {
#if 0
        mx_expr e = simplify( any_cast<mx_expr>(attr.value()), mx );
        t.insert_attribute( attr.name(), e );
#endif
      t.insert_attribute(attr.name(), any_cast<mx_expr>(attr.value()));
    } else if (attr.type() == MX_ATTR_LAZY_NUMBER) {
      t.insert_attribute(attr.name(), any_cast<mx_expr>(attr.value()));
    } else if (attr.type() == MX_ATTR_ARRAY) {
#if 0
        mx_exprs es = any_cast<mx_exprs>(attr.value());
        for( mx_exprs::iterator it = es.begin()
           ; it != es.end(); ++it )
        {
          mx_expr e = simplify( *it, mx );
          *it = e;
        }
        t.insert_attribute( attr.name(), es );
#endif
      mx_exprs es = any_cast<mx_exprs>(attr.value());
      t.insert_attribute(attr.name(), es);
    } else if (attr.type() == MX_ATTR_LAZY_ARRAY) {
      mx_exprs es = any_cast<mx_exprs>(attr.value());
      t.insert_attribute(attr.name(), es);
    }
  }
}}

void
mx_logic::set(mx_expr const& l, logic_op_t o, mx_expr const& r)
{
  lhs = l;
  rhs = r;
  op = o;
  use_preset = false;
}

bool
mx_logic::evaluate(MadX const& mx) const
{
  if (use_preset) return pre;

  assert(op != NULL);

  double l = boost::apply_visitor(mx_calculator(mx, 0.0), lhs);
  double r = boost::apply_visitor(mx_calculator(mx, 0.0), rhs);

  return op(l, r);
}

// attributes
void
mx_attr::set_attr(std::string const& name, std::any const& val)
{
  name_ = name;
  value_ = val;

  transform(name_.begin(), name_.end(), name_.begin(), ::tolower);

  if (val.type() == typeid(std::string))
    type_ = MX_ATTR_STRING;
  else if (val.type() == typeid(mx_expr))
    type_ = MX_ATTR_NUMBER;
  else if (val.type() == typeid(mx_exprs))
    type_ = MX_ATTR_ARRAY;
  else if (val.type() == typeid(mx_keyword) &&
           (any_cast<mx_keyword>(val).tag == MX_KW_PARTICLE ||
            any_cast<mx_keyword>(val).tag == MX_KW_MP_TYPE))
    type_ = MX_ATTR_PREDEFINED;
  else
    throw std::runtime_error("Unknown attribute value type for " + name);
}

void
mx_attr::set_lazy_attr(std::string const& name, std::any const& val)
{
  name_ = name;
  value_ = val;

  transform(name_.begin(), name_.end(), name_.begin(), ::tolower);

  if (val.type() == typeid(std::string))
    type_ = MX_ATTR_STRING; // strings are immediate
  else if (val.type() == typeid(mx_expr))
    type_ = MX_ATTR_LAZY_NUMBER;
  else if (val.type() == typeid(mx_exprs))
    type_ = MX_ATTR_LAZY_ARRAY;
  else if (val.type() == typeid(mx_keyword) &&
           (any_cast<mx_keyword>(val).tag == MX_KW_PARTICLE ||
            any_cast<mx_keyword>(val).tag == MX_KW_MP_TYPE))
    type_ = MX_ATTR_PREDEFINED; // all predefines are immediate
  else
    throw std::runtime_error("Unknown lazy attribute value type for " + name);
}

// mx_line_member

void
mx_line::interpret(MadX& mx)
{
  MadX_line new_line(mx);
  seq.interpret(mx, new_line, 1);
  mx.insert_line(name, new_line);
}

void
mx_line_seq::insert_member(int op, mx_line_member const& member)
{
  members.push_back(make_pair(member, op));
}

void
mx_line_seq::interpret(MadX const& mx, MadX_line& line, int op)
{
  if (op == 0) return;

  if (op > 0) // repeat op times
  {
    for (int z = 0; z < op; ++z) {
      for (mx_line_members::const_iterator it = members.begin();
           it != members.end();
           ++it) {
        mx_line_member member = it->first;
        int op = it->second;
        member.interpret(mx, line, op);
      }
    }
    return;
  }

  if (op < 0) // reverse then repeat
  {
    for (int z = 0; z > op; --z) {
      for (mx_line_members::reverse_iterator it = members.rbegin();
           it != members.rend();
           ++it) {
        mx_line_member member = it->first;
        int op = it->second;
        member.interpret(mx, line, op);
      }
    }
    return;
  }
}

mx_line_member::mx_line_member() : member(), tag(MX_LINE_MEMBER_NAME) {}

mx_line_member::mx_line_member(string_t const& name)
  : member(name), tag(MX_LINE_MEMBER_NAME)
{}

mx_line_member::mx_line_member(mx_line_seq const& seq)
  : member(seq), tag(MX_LINE_MEMBER_SEQ)
{}

void
mx_line_member::interpret(MadX const& mx, MadX_line& line, int op)
{
  if (op == 0) return;

  if (tag == MX_LINE_MEMBER_NAME) {
    // member is a name/reference
    string name = any_cast<string>(member);
    MadX_entry_type type = mx.entry_type(name);

    // does the name refer to a madx command?
    if (type == ENTRY_COMMAND) {
      // ok it is a command
      if (mx.command(name, true).is_element()) {
        // now it is a real element
        if (op != 1)
          throw runtime_error(
            "Line op only applies to sublines, not to elements!");

        // push to the line
        line.insert_element(name);
      } else {
        throw runtime_error("Line member '" + name + "' is not an element");
      }
    }
    // or the name referes to a pre-exisitng line
    else if (type == ENTRY_LINE) {
      MadX_line const& subline = mx.line(name);
      size_t ne = subline.element_count();
      int repeat = (op > 0) ? op : -op;

      for (int z = 0; z < repeat; ++z) {
        for (size_t i = 0; i < ne; ++i)
          line.insert_element(
            subline.element_name((op > 0) ? (i) : (ne - 1 - i)));
      }
    }
    // TODO: for now, we accept the sequence name as a simple line member
    // Needs more works!
    else if (type == ENTRY_SEQUENCE) {
      line.insert_element(name);
    }
    // something we dont support
    else {
      // TODO: for now, accept all names whether it is valid or not!!
      // throw runtime_error("Line member '" + name + "' does not exist or not
      // correct type");
    }
  } else {
    // if it is not a name, it must be a seq
    any_cast<mx_line_seq>(member).interpret(mx, line, op);
  }

  return;
}

// command
void
mx_command::set_label(string const& label)
{
  label_ = label;
  labeled_ = true;
  transform(label_.begin(), label_.end(), label_.begin(), ::tolower);
}

void
mx_command::set_keyword(mx_keyword const& keyword)
{
  if (keyword.tag == MX_KW_NONE || keyword.tag == MX_KW_PARTICLE)
    throw std::runtime_error("Invalid keyword type");

  mx_cmd_type t = (keyword.tag == MX_KW_ELEMENT) ? MX_CMD_ELEMENT :
                  (keyword.tag == MX_KW_COMMAND) ? MX_CMD_EXECUTABLE :
                                                   MX_CMD_ELEMENT_REF;

  set_keyword(keyword.name, t);
}

void
mx_command::set_keyword(string const& keyword, mx_cmd_type tag)
{
  if (tag != MX_CMD_ELEMENT && tag != MX_CMD_EXECUTABLE &&
      tag != MX_CMD_ELEMENT_REF)
    throw std::runtime_error("Unknown keyword type");

  keyed_ = true;
  keyword_ = keyword;
  type_ = tag;
}

void
mx_command::ins_attr(mx_attr const& attr)
{
  attrs_.push_back(attr);
}

MadX_command
mx_command::convert(MadX const& mx) const
{
  // prepare MadX_command object
  MadX_command cmd;

  cmd.set_label(label_);
  cmd.set_name(keyword_, ELEMENT);

  for (attrs_t::const_iterator it = attrs_.begin(); it != attrs_.end(); ++it) {
    mx_attr attr = *it;

    if (attr.name() == "from") {
      string from = mx_expr_refstr(any_cast<mx_expr>(attr.value()));
      attr.set_attr("from", from);
    }

    insert_attr(cmd, attr, mx);
  }

  return cmd;
}

void
mx_command::interpret(MadX& mx)
{
  // create an MadX_command object from mx_command
  MadX_command new_cmd = convert(mx);

  // type of the command?
  if (type_ == MX_CMD_VARIABLE) {
    mx_attr attr = attrs_[0];
    insert_attr(mx, attr, mx);
  } else if (type_ == MX_CMD_ELEMENT) {
    // unlabeled class? -- it is a warning and will be skipped
    if (!labeled_) {
      cout << "statment illegal in the context (declaring the element '"
           << keyword_ << "' without a label), will be skipped\n";
      return;
    }

    // insert a label
    mx.insert_label(label_, new_cmd);

    // building a sequence?
    if (mx.building_sequence()) {
      mx_expr at = new_cmd.attribute_as_expr("at");
      string from = new_cmd.attribute_as_string("from", "");

      // label_, at, from
      mx.append_sequence_element(label_, at, from);
    }
  } else if (type_ == MX_CMD_ELEMENT_REF) {
    // type of the referenced entry
    MadX_entry_type ref_type = mx.entry_type(keyword_);

    if (labeled_) {
      // check ref type
      if (ref_type != ENTRY_COMMAND)
        throw runtime_error("unknown class type " + keyword_ + ".");

      // create a new element instance
      MadX_command ori_cmd = mx.command(keyword_);

      new_cmd.merge(ori_cmd);
      new_cmd.set_name(ori_cmd.name(), ELEMENT);

      mx.insert_label(label_, new_cmd);
    } else if (!mx.building_sequence()) {
      // check ref type
      if (ref_type != ENTRY_COMMAND) {
        cout << "warning: unknown class type " << keyword_ << "\n";
        return;
      }

      // merge with the existing instance
      mx.fuse_command(keyword_, new_cmd);
    }

    // building a sequence?
    if (mx.building_sequence()) {
      if (ref_type != ENTRY_COMMAND && ref_type != ENTRY_SEQUENCE)
        throw runtime_error("unknown class type " + keyword_ + ".");

      mx_expr at = new_cmd.attribute_as_expr("at");
      string from = new_cmd.attribute_as_string("from", "");

      // label_ or keyword_
      mx.append_sequence_element((labeled_ ? label_ : keyword_), at, from);
    }
  } else // execute the command
  {
    execute(mx);

    // mx.insert_command(new_cmd);
  }

  return;
}

void
mx_command::execute(MadX& mx)
{
  if (keyword_ == "call") {
    // pull in the sub-file
    for (attrs_t::const_iterator it = attrs_.begin(); it != attrs_.end();
         ++it) {
      if (it->name() == "file") {
        string fname = std::any_cast<string>(it->value());
        mx_tree subroutine;
        parse_int_madx_file(fname, subroutine);
        subroutine.interpret(mx);
        return;
      }
    }
    throw runtime_error("Error executing command 'call'");
  } else if (keyword_ == "sequence") {
    double length = 0.0;
    string refer = string("");
    string refpos = string("");

    // the "refer" attribute is an unquoted string, not a reference
    for (attrs_t::iterator it = attrs_.begin(); it != attrs_.end(); ++it) {
      if (it->name() == "refer") {
        if (it->value().type() == typeid(string)) {
          refer = std::any_cast<string>(it->value());
        } else if (it->value().type() == typeid(mx_expr)) {
          try {
            mx_expr ex = any_cast<mx_expr>(it->value());

            // this line causes the final extracted string to be random
            // characters if compiled with os x clang compiler ex =
            // get<nop_t>(get<nop_t>(get<nop_t>(ex).expr).expr).expr;

            // fix to the above issue
            ex = get<nop_t>(ex).expr;
            ex = get<nop_t>(ex).expr;
            ex = get<nop_t>(ex).expr;

            refer = boost::get<string>(ex);
            it->set_attr("refer", refer);
          }
          catch (...) {
            throw runtime_error("The 'refer' attribute of sequence '" + label_ +
                                "' is not a string");
          }
        } else {
          throw runtime_error("The 'refer' attribute of sequence '" + label_ +
                              "' is not a string");
        }
      }

      if (it->name() == "refpos") {
        if (it->value().type() == typeid(string)) {
          refpos = std::any_cast<string>(it->value());
        } else if (it->value().type() == typeid(mx_expr)) {
          try {
            mx_expr ex = any_cast<mx_expr>(it->value());
            ex = get<nop_t>(ex).expr;
            ex = get<nop_t>(ex).expr;
            ex = get<nop_t>(ex).expr;
            refpos = boost::get<string>(ex);
            it->set_attr("refpos", refpos);
          }
          catch (...) {
            throw runtime_error("The 'refpos' attribute of sequence '" +
                                label_ + "' is not a string");
          }
        } else {
          throw runtime_error("The 'refpos' attribute of sequence '" + label_ +
                              "' is not a string");
        }
      }

      if (it->name() == "length" || it->name() == "l") {
        try {
          mx_expr e = std::any_cast<mx_expr>(it->value());
          length = boost::apply_visitor(mx_calculator(mx), e);
        }
        catch (...) {
          throw runtime_error("The 'length' attribute of sequence '" + label_ +
                              "' is not a number");
        }
      }

    } // end of attr iter loop

    transform(refer.begin(), refer.end(), refer.begin(), ::tolower);
    transform(refpos.begin(), refpos.end(), refpos.begin(), ::tolower);

    if (length < 0.0)
      throw runtime_error("The 'length' attribute of sequence '" + label_ +
                          "' is not a valid number");

    // tells madx object to start building sequence
    mx.start_sequence(label_, length, refer, refpos);
  } else if (keyword_ == "endsequence") {
    mx.end_sequence();
  } else if (keyword_ == "beam") {
    // beam can always be referenced with 'beam'
    if (!labeled_) {
      label_ = "beam";
      labeled_ = true;
    }

    // build attributes of beam
    double mass = 0, charge = 0, energy = 0, pc = 0, gamma = 0;
    bool have_charge = false;
    for (attrs_t::const_iterator it = attrs_.begin(); it != attrs_.end();
         ++it) {
      if (it->name() == "particle") {
        string particle = any_cast<mx_keyword>(it->value()).name;
        if (particle == "proton") {
          mass = pconstants::mp;
          charge = pconstants::proton_charge;
        } else if (particle == "antiproton") {
          mass = pconstants::mp;
          charge = pconstants::antiproton_charge;
        } else if (particle == "electron") {
          mass = pconstants::me;
          charge = pconstants::electron_charge;
        } else if (particle == "positron") {
          mass = pconstants::me;
          charge = pconstants::positron_charge;
        } else if (particle == "negmuon") {
          mass = pconstants::mmu;
          charge = pconstants::muon_charge;
        } else if (particle == "posmuon") {
          mass = pconstants::mmu;
          charge = pconstants::antimuon_charge;
        } else {
          throw runtime_error("Unknown particle type " + particle);
        }
      } else if (it->name() == "mass") {
        mx_expr e = std::any_cast<mx_expr>(it->value());
        mass = boost::apply_visitor(mx_calculator(mx), e);
      } else if (it->name() == "charge") {
        mx_expr e = std::any_cast<mx_expr>(it->value());
        charge = boost::apply_visitor(mx_calculator(mx), e);
        have_charge = true;
      } else if (it->name() == "energy") {
        mx_expr e = std::any_cast<mx_expr>(it->value());
        energy = boost::apply_visitor(mx_calculator(mx), e);
      } else if (it->name() == "pc") {
        mx_expr e = std::any_cast<mx_expr>(it->value());
        pc = boost::apply_visitor(mx_calculator(mx), e);
      } else if (it->name() == "gamma") {
        mx_expr e = std::any_cast<mx_expr>(it->value());
        gamma = boost::apply_visitor(mx_calculator(mx), e);
      }
    }

    // makes no change if the particle type is absent
    if (mass != 0) {
      Four_momentum four_momentum(mass);

      if (energy > 0) four_momentum.set_total_energy(energy);
      if (pc > 0) four_momentum.set_momentum(pc);
      if (gamma > 0) four_momentum.set_gamma(gamma);

      mx_attr attr;

      if (energy == 0) {
        attr.set_attr("energy", mx_expr(four_momentum.get_total_energy()));
        ins_attr(attr);
      }

      if (pc == 0) {
        attr.set_attr("pc", mx_expr(four_momentum.get_momentum()));
        ins_attr(attr);
      }

      if (gamma == 0) {
        attr.set_attr("gamma", mx_expr(four_momentum.get_gamma()));
        ins_attr(attr);
      }

      if (!have_charge) {
        attr.set_attr("charge", mx_expr(charge));
        ins_attr(attr);
      }

      // insert a global variable brho to the madx object
      stringstream ss;
      ss.precision(18);
      ss << 1e9 / pconstants::c << "*beam->pc";

      mx_expr expr;
      parse_expression(ss.str(), expr);
      mx.insert_variable("brho", expr);
    }

    // insert beam as a label
    MadX_command new_cmd = convert(mx);
    mx.insert_label(label_, new_cmd);
  }
}

void
mx_command::print() const
{
  cout << "type(" << (int)type_ << ") ";

  if (labeled_) cout << label_ << " : ";

  if (keyed_) cout << keyword_ << ", ";

  for (attrs_t::const_iterator it = attrs_.begin(); it != attrs_.end(); ++it) {
    cout << it->name() << " = "
         << "xxx, ";
  }

  cout << "\n";
}

// if_block
bool
mx_if_block::evaluate_logic(MadX const& mx) const
{
  return logic_expr.evaluate(mx);
}

void
mx_if_block::interpret_block(MadX& mx)
{
  block.interpret(mx);
}

void
mx_if_block::print_logic() const
{
  // cout << logic_expr;
}

void
mx_if_block::print_block() const
{
  block.print();
}

// if-elseif-else
void
mx_if::assign_if(mx_logic const& logic, mx_tree const& block)
{
  if_ = mx_if_block(logic, block);
}

void
mx_if::assign_elseif(mx_logic const& logic, mx_tree const& block)
{
  elseif_.push_back(mx_if_block(logic, block));
}

void
mx_if::assign_else(mx_tree const& block)
{
  else_ = mx_if_block(true, block);
}

void
mx_if::interpret(MadX& mx)
{
  if (!if_.valid())
    throw runtime_error("mx_if::interpret() Invalid if command");

  if (if_.evaluate_logic(mx)) {
    if_.interpret_block(mx);
  } else if (!elseif_.empty()) {
    for (mx_if_block_v::iterator it = elseif_.begin(); it != elseif_.end();
         ++it) {
      if (it->evaluate_logic(mx)) { it->interpret_block(mx); }
    }
  } else if (else_.valid()) {
    else_.interpret_block(mx);
  }

  return;
}

void
mx_if::print() const
{
  if (!if_.valid()) throw runtime_error("mx_if::print() Invalid if command");

  cout << "if (";
  if_.print_logic();
  cout << ")\n{\n";
  if_.print_block();
  cout << "}\n";

  for (mx_if_block_v::const_iterator it = elseif_.begin(); it != elseif_.end();
       ++it) {
    cout << "elseif (";
    it->print_logic();
    cout << ")\n{\n";
    it->print_block();
    cout << "}\n";
  }

  if (else_.valid()) {
    cout << "else\n{\n";
    else_.print_block();
    cout << "}\n";
  }
}

// while
void
mx_while::assign(mx_logic const& logic, mx_tree const& block)
{
  while_ = mx_if_block(logic, block);
}

void
mx_while::interpret(MadX& mx)
{
  return;
}

void
mx_while::print() const
{
  if (!while_.valid())
    throw runtime_error("mx_while::print() Invalid while command");

  cout << "while (";
  while_.print_logic();
  cout << ")\n{\n";
  while_.print_block();
  cout << "}\n";
}

// statement

synergia::mx_statement::mx_statement() : value(), type(MX_NULL) {}

synergia::mx_statement::mx_statement(mx_command const& st)
  : value(st), type(MX_COMMAND)
{}

synergia::mx_statement::mx_statement(mx_if const& st) : value(st), type(MX_IF)
{}

synergia::mx_statement::mx_statement(mx_while const& st)
  : value(st), type(MX_WHILE)
{}

synergia::mx_statement::mx_statement(mx_line const& st)
  : value(st), type(MX_LINE)
{}

mx_statement_type
synergia::mx_statement::get_type() const
{
  return type;
}

void
mx_statement::assign(mx_command const& st)
{
  value = std::any(st);
  type = MX_COMMAND;
}

void
mx_statement::assign(mx_if const& st)
{
  value = std::any(st);
  type = MX_IF;
}

void
mx_statement::assign(mx_while const& st)
{
  value = std::any(st);
  type = MX_WHILE;
}

void
mx_statement::assign(mx_line const& st)
{
  value = std::any(st);
  type = MX_LINE;
}

void
mx_statement::interpret(MadX& mx)
{
  switch (type) {
    case MX_COMMAND: std::any_cast<mx_command>(value).interpret(mx); break;
    case MX_LINE: std::any_cast<mx_line>(value).interpret(mx); break;
    case MX_IF: std::any_cast<mx_if>(value).interpret(mx); break;
    case MX_WHILE: std::any_cast<mx_while>(value).interpret(mx); break;
    case MX_NULL: break;
    default:
      throw runtime_error("mx_statement::interpret()  Unknown statement type");
  }
}

void
mx_statement::print() const
{
  if (type == MX_COMMAND)
    return std::any_cast<mx_command>(value).print();
  else if (type == MX_IF)
    return std::any_cast<mx_if>(value).print();
  else if (type == MX_WHILE)
    return std::any_cast<mx_while>(value).print();
  else
    throw runtime_error("mx_statement::print()  Unknown statement type");
}

// tree
void
mx_tree::interpret(MadX& mx)
{
  for (mx_statements_t::iterator it = statements.begin();
       it != statements.end();
       ++it) {
    auto type = it->get_type();
    if (type == MX_COMMAND || type == MX_IF || type == MX_WHILE) {
      it->interpret(mx);
    }
  }

  for (mx_statements_t::iterator it = statements.begin();
       it != statements.end();
       ++it) {
    auto type = it->get_type();
    if (type != MX_COMMAND && type != MX_IF && type != MX_WHILE) {
      it->interpret(mx);
    }
  }
}

void
mx_tree::print() const
{
  for_each(statements.begin(), statements.end(), detail::print<mx_statement>);
}
