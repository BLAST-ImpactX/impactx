
#include "synergia/bunch/diagnostics_py.h"
#include "synergia/bunch/bunch.h"

PyDiagnostics::~PyDiagnostics() noexcept {}

void
PyDiagnostics::do_update(Bunch const& bunch)
{
    namespace py = pybind11;

    py::object obj = py::cast(bunch, py::return_value_policy::reference);
    self.attr("do_update")(obj);

#if 0
    PYBIND11_OVERLOAD_PURE(
            void,        // return type
            Diagnostics, // parent class
            do_update,   // name of the function in c++
            bunch        // arguments
    );
#endif
}

void
PyDiagnostics::do_reduce(Commxx const& comm, int root)
{
    self.attr("do_reduce")(comm, root);
}

void
PyDiagnostics::do_first_write(io_device& file)
{
    // self.attr("do_first_write")(file);
}

void
PyDiagnostics::do_write(io_device& file, size_t iteration)
{
    // self.attr("do_write")(file);
}
