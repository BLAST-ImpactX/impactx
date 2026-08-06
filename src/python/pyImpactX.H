/* Copyright 2022-2026 The ImpactX Community
 *
 * This header is used to centrally define classes that shall not violate the
 * C++ one-definition-rule (ODR) for various Python translation units.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#ifndef IMPACTX_PYIMPACTX_H
#define IMPACTX_PYIMPACTX_H

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
// include <pybind11/numpy.h>  // not yet used
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/warnings.h>

#include <ImpactX.H>
#include <elements/All.H>

#include <functional>
#include <list>
#include <string>
#include <unordered_map>

namespace py = pybind11;
using namespace impactx;

PYBIND11_MAKE_OPAQUE(std::list<elements::KnownElements>)
PYBIND11_MAKE_OPAQUE(std::unordered_map<std::string, std::function<void(ImpactX *)> >)

/** Wrap a Python callable into a C++ CSR kick model
 *
 * The returned wrapper acquires the GIL, passes borrowed views of the beam
 * profile and the element context to the callable and converts its result
 * (a dict with the required key "pt" and the optional keys "px"/"py", or a
 * bare array meaning {"pt": array}) into per-bin kick device vectors.
 *
 * Implemented in CSRKickModel.cpp.
 *
 * @param callable a Python callable (CSRProfile, CSRElementContext) -> dict | array
 * @return the C++ kick model wrapping the callable
 */
impactx::particles::wakefields::CSRKickModel
make_csr_kick_model (pybind11::object callable);

#endif // IMPACTX_PYIMPACTX_H
