/* Copyright 2021-2023 The ImpactX Community
 *
 * Authors: Ryan Sandberg, Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyImpactX.H"
#include <AMReX_SmallMatrix.H>

namespace py = pybind11;

namespace pybind11 {
namespace detail {

template <typename T, int NRows, int NCols, amrex::Order ORDER, int StartIndex>
struct pybind11::detail::type_caster<amrex::SmallMatrix<T, NRows, NCols, ORDER, StartIndex>> {
public:
    PYBIND11_TYPE_CASTER(amrex::SmallMatrix<T, NRows, NCols, ORDER, StartIndex>,
                         _("SmallMatrix[") + py::detail::make_caster<T>::name() + _("]"));

    // Conversion from Python to C++
    bool load(handle src, bool) {
        // Ensure we have a numpy array
        py::array_t<T> arr = py::cast<py::array_t<T>>(src);
        py::buffer_info buf = arr.request();

        // Check dimensions and shape
        if (buf.ndim != 2) {
            throw std::runtime_error("SmallMatrix requires a 2D array.");
        }
        if (buf.shape[0] != NRows || buf.shape[1] != NCols) {
            throw std::runtime_error("SmallMatrix array shape must match NRows x NCols.");
        }

        // Create a SmallMatrix and copy data
        amrex::SmallMatrix<T, NRows, NCols, ORDER, StartIndex> mat;
        T* ptr = static_cast<T*>(buf.ptr);
        for (int i = 0; i < NRows * NCols; ++i) {
            mat.m_mat[i] = ptr[i];
        }

        value = mat;
        return true;
    }

    // Conversion from C++ to Python
    static handle cast(const amrex::SmallMatrix<T, NRows, NCols, ORDER, StartIndex>& src,
                       return_value_policy /* policy */, handle /* parent */) {
        py::array_t<T> arr({NRows, NCols});
        py::buffer_info buf = arr.request();
        T* ptr = static_cast<T*>(buf.ptr);
        for (int i = 0; i < NRows * NCols; ++i) {
            ptr[i] = src.m_mat[i];
        }
        return arr.release();
    }
};

} // namespace detail
} // namespace pybind11


PYBIND11_MODULE(example, m) {
    // You can now just bind constructors and methods normally without defining conversion code:
    py::class_<amrex::SmallMatrix<double, 6, 6>>(m, "SmallMatrix6x6")
        .def(py::init<>())  // Default init
        .def("as_array", [](const amrex::SmallMatrix<double, 6, 6>& mat) {
            return mat; // Will use type_caster to return a numpy array
        });
    
    // Now Python functions expecting a SmallMatrix<double,6,6> can pass a numpy array directly:
    // def some_func(mat: SmallMatrix6x6): ...
}
