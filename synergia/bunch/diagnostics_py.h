#ifndef DIAGNOSTICS_PY_H
#define DIAGNOSTICS_PY_H

#include "synergia/bunch/diagnostics.h"
#include "synergia/utils/base64.h"

#include <pybind11/pybind11.h>

class PyDiagnostics : public Diagnostics {
  public:
    explicit PyDiagnostics(std::string const& type = "PyDiagnostics",
                           std::string const& filename = "py_diag.h5",
                           bool serial = true)
        : Diagnostics(type, filename, serial), self()
    {}

    // Declared here so that we can force instantiation of typeinfo in
    // diagnostics_py.cc
    ~PyDiagnostics() noexcept;

    // hold a reference to the python instance so it wont go out of scope
    void
    reg_self()
    {
        self = pybind11::cast(this);
    }

  private:
    pybind11::object self;

    void do_update(Bunch const& bunch) override;
    void do_reduce(Commxx const& comm, int root) override;
    void do_first_write(io_device& file) override;
    void do_write(io_device& file, size_t iteration) override;

    friend class cereal::access;

    template <class Archive>
    void
    save(Archive& ar) const
    {
        ar(cereal::base_class<Diagnostics>(this));

        namespace py = pybind11;
        py::object pickle;
        try {
            pickle = py::module::import("cPickle");
        }
        catch (...) {
            pickle = py::module::import("pickle");
        }

        auto s = pickle.attr("dumps")(self, 2).cast<std::string>();
        auto base_s = base64_encode(
            reinterpret_cast<const unsigned char*>(s.c_str()), s.length());

        ar(base_s);
    }

    template <class Archive>
    void
    load(Archive& ar)
    {
        ar(cereal::base_class<Diagnostics>(this));

        std::string base_s;
        ar(base_s);

        std::string s = base64_decode(base_s);

        namespace py = pybind11;
        py::object pickle;
        try {
            pickle = py::module::import("cPickle");
        }
        catch (...) {
            pickle = py::module::import("pickle");
        }

        self = pickle.attr("loads")(py::bytes(s));
    }
};

CEREAL_REGISTER_TYPE(PyDiagnostics);
CEREAL_SPECIALIZE_FOR_ALL_ARCHIVES(PyDiagnostics,
                                   cereal::specialization::member_load_save);

#endif
