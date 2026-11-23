#include "synergia/lattice/lattice.h"
#include "synergia/lattice/lattice_element_processor.h"
#include "synergia/lattice/madx_reader.h"

#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <unordered_set>

Lattice::Lattice()
    : name("")
    , reference_particle()
    , elements()
    , updated{true, true, true}
    , tree()
{}

Lattice::Lattice(std::string const& name)
    : name(name)
    , reference_particle()
    , elements()
    , updated{true, true, true}
    , tree()
{}

Lattice::Lattice(std::string const& name, Reference_particle const& ref)
    : name(name)
    , reference_particle(ref)
    , elements()
    , updated{true, true, true}
    , tree()
{}

Lattice::Lattice(std::string const& name, Lattice_tree const& tree)
    : name(name)
    , reference_particle()
    , elements()
    , updated{true, true, true}
    , tree(tree)
{}

Lattice::Lattice(Lattice const& o)
    : name(o.name)
    , reference_particle(o.reference_particle)
    , elements(o.elements)
    , updated(o.updated)
    , tree(o.tree)
{
    for (auto& e : elements)
        e.set_lattice(*this);
}

Lattice::Lattice(Lattice&& o) noexcept
    : name(std::move(o.name))
    , reference_particle(std::move(o.reference_particle))
    , elements(std::move(o.elements))
    , updated(std::move(o.updated))
    , tree(std::move(o.tree))
{
    for (auto& e : elements)
        e.set_lattice(*this);
}

Lattice&
Lattice::operator=(Lattice const& o)
{
    name = o.name;
    reference_particle = o.reference_particle;
    elements = o.elements;
    updated = o.updated;
    tree = o.tree;

    for (auto& e : elements)
        e.set_lattice(*this);

    return *this;
}

Lattice::Lattice(Lsexpr const& lsexpr)
    : name("")
    , reference_particle()
    , elements()
    , updated{true, true, true}
    , tree()
{
    for (auto const& lse : lsexpr) {
        if (lse.is_labeled()) {
            if (lse.get_label() == "name") {
                name = lse.get_string();
            } else if (lse.get_label() == "reference_particle") {
                reference_particle = Reference_particle(lse);
            } else if (lse.get_label() == "elements") {
                for (auto const& ele : lse) {
                    append(Lattice_element(ele));
                }
            }
        }
    }
}

#if 0 // as_lsexpr() might come back some day
Lsexpr
Lattice::as_lsexpr() const
{
    Lsexpr retval;
#if 0
    retval.push_back(Lsexpr(name, "name"));
    retval.push_back(Lsexpr(element_adaptor_map_sptr->get_label(),
                            "type"));
    if (reference_particle) {
        Lsexpr ref_lsexpr(reference_particle->as_lsexpr());
        ref_lsexpr.set_label("reference_particle");
        retval.push_back(ref_lsexpr);
    }
    Lsexpr elements_lsexpr;
    elements_lsexpr.set_label("elements");
    for(Lattice_elements::const_iterator it = elements.begin();
        it != elements.end(); ++it) {
        Lsexpr element_lsexpr((*it)->as_lsexpr());
        elements_lsexpr.push_back(element_lsexpr);
    }
    retval.push_back(elements_lsexpr);
#endif
    return retval;
}
#endif

Lattice::update_flags_t
Lattice::update()
{
    update_flags_t res = updated;

    if (updated.ref || updated.structure || updated.element) {
        // TODO: reconstruct the chef lattice
        // ...

        updated = {false, false, false};
    }

    return res;
}

void
Lattice::append(Lattice_element const& element)
{
    elements.push_back(Lattice_element_processor::process(element));
    elements.back().set_lattice(*this);
    updated.structure = true;
}

void
Lattice::set_all_double_attribute(std::string const& name,
                                  double value,
                                  bool increment_revision)
{
    for (auto& e : elements)
        e.set_double_attribute(name, value, increment_revision);

    updated.element = true;
}

void
Lattice::set_all_string_attribute(std::string const& name,
                                  std::string const& value,
                                  bool increment_revision)
{
    for (auto& e : elements)
        e.set_string_attribute(name, value, increment_revision);

    updated.element = true;
}

std::list<Lattice_element> const&
Lattice::get_elements() const
{
    return elements;
}

std::list<Lattice_element>&
Lattice::get_elements()
{
    return elements;
}

void
Lattice::reset_all_markers()
{
    for (auto& e : elements)
        e.reset_markers();
}

double
Lattice::get_length() const
{
    double length = 0.0;
    for (auto const& e : elements)
        length += e.get_length();
    return length;
}

double
Lattice::get_total_angle() const
{
    double angle = 0.0;
    for (auto const& e : elements)
        angle += e.get_bend_angle();
    return angle;
}

std::string
Lattice::as_string() const
{
    std::stringstream sstream;
    sstream << name << ":\n";
    for (auto const& e : elements)
        sstream << e.as_string() << std::endl;
    return sstream.str();
}

void
Lattice::print(std::ostream& os) const
{
    os << as_string() << std::endl;
}

void
Lattice::export_madx_file(std::string const& filename, bool const sanitize) const
{
    std::ofstream mxfile(filename);

    if (!mxfile.is_open()) {
        throw std::runtime_error(
            "Lattice::as_madx_file() unable to create file");
    }

    // "beam, pc={{momentum}}, particle={{particle}}"
    if (!reference_particle.has_value()) {
        mxfile << "! This lattice has not defined a reference particle with "
                  "BEAM statement information\n\n";
    } else {
        mxfile << reference_particle.value().as_madx() << "\n\n";
    }

    // "{{element_label}} : {{element_type}}, {{attr}}={{val}}..."
    std::unordered_set<std::string> elm_names;
    for (auto const& e : elements) {
        mxfile << e.as_madx(sanitize) << "\n";
        elm_names.insert(e.get_name());
    }

    double pos = 0.0;
    int idx = 0;

    // "{{name}}: sequence, refer=entry"
    mxfile << "\n"
           << name << ": sequence, l = " << std::setprecision(16) << get_length()
           << ", refer = entry;\n";

    // "{{element_label}} : {{element}}, at={{pos}}..."
    for (auto const& e : elements) {
        std::stringstream label;
        label << name << "_" << idx;

        while (elm_names.count(label.str()))
            label << "_dup";

        elm_names.insert(label.str());

        mxfile << label.str() << ": " << e.get_name() << ", at=" << std::setprecision(16) << pos
               << ";\n";

        pos += e.get_length();
        ++idx;
    }

    // "endsequence;"
    mxfile << "endsequence;\n";

    // done
    mxfile.close();
}

Lattice
Lattice::import_madx_file(std::string const& filename, std::string const& line)
{
    return MadX_reader().get_lattice(line, filename);
}

Lattice_tree const&
Lattice::get_lattice_tree() const
{
    return tree;
}

Lattice_tree&
Lattice::get_lattice_tree()
{
    return tree;
}

void
Lattice::set_lattice_tree(Lattice_tree const& lattice_tree)
{
    tree = lattice_tree;
}

void
Lattice::set_variable(std::string const& name, double val)
{
    tree.set_variable(name, val);
}

void
Lattice::set_variable(std::string const& name, std::string const& val)
{
    tree.set_variable(name, val);
}
