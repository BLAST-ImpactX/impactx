---
name: impactx-update-element-benchmarks
description: This skill should be used when the user asks to "update element benchmarks", "add missing element tests", "audit element test coverage", "add a benchmark for <element>", or to check that every lattice element under `src/elements/` has matching entries in the three Python element-coverage tests (benchmark, serialization, reversibility) and that spin-capable elements carry the spin/nospin parametrization. Also covers re-running those three tests after edits.
version: 0.1.0
---

# Update Element Benchmarks

Audit and extend the three Python element-coverage tests so that every lattice element exposed from `src/elements/` is represented, and that every spin-capable element is parametrized over spin tracking. Then run the three tests.

To run this skill, invoke `/impactx-update-element-benchmarks` in a new session.
The skill also auto-loads on triggers like "update element benchmarks", "add missing element tests", "audit element test coverage", or "add a benchmark for <element>".


## The three target files

1. `tests/python/test_benchmark_elements.py` — pytest-benchmark micro-benchmarks of a single `el.push(beam)` call. Uses the `sim` fixture and `pc_setup` helper.
2. `tests/python/test_element_serialization.py` — builds one `KnownElementsList` containing one instance of every element, then verifies `to_dicts()`/`from_dicts()`/`to_py()` round-trip. A new element is "covered" by being appended inside the `all_elements` fixture.
3. `tests/python/test_reversibility_elements.py` — verifies that `forward push → reverse() → forward push` returns phase space (and, when enabled, spin) to the initial state. Uses the `sim` fixture and the `roundtrip()` helper.

## Canonical element list

Source of truth: the `KnownElements` `std::variant` at the bottom of `src/elements/All.H`. Treat each entry there as one element that needs coverage (the Python binding name may differ — e.g. `ChrUniformAcc` → `elements.ChrAcc`). Use the existing test files to learn the Python-side name for each C++ type.

## Pass 1 — coverage audit

For each element in `KnownElements`, confirm an entry exists in each of the three files. Missing coverage must either be added or intentionally skipped with a one-line comment explaining why.

### Known intentional omissions (do NOT add)

Some elements are intentionally not covered. Keep them omitted and leave the existing explanatory comments in place:

- `Empty` — placeholder type in the variant; skipped across all three files (benchmark has commented stub; serialization lists it in `SKIP_ELEMENTS` with the `to_dict` bug note).
- `Source` — requires an openPMD file path, cannot be constructed in-memory; skipped across all three files (serialization lists it in `SKIP_ELEMENTS`).
- `Marker` — no push behavior to benchmark; no state change to reverse. Present in serialization only.
- `Programmable` — user-supplied Python callback, not meaningful to benchmark or roundtrip generically. Present in serialization only.
- `BeamMonitor` — I/O element, not a push. Present in serialization only (and `.finalize()` is called in fixture teardown).
- `SpinMap` — reversibility is tracked as a known-broken draft (see the commented block at the bottom of `test_reversibility_elements.py`); keep it commented out until the C++ reversal is fixed. Present in serialization only.
- `QuadEdge` benchmark — omitted pending issue #986 (see existing TODO comment). Keep omitted.

When in doubt, check whether the element already appears as a commented-out stub in the file — that is the signal that omission is intentional.

### Aperture-style elements

`Aperture` and `PolygonAperture` are both real push elements and should have benchmark + serialization coverage (`test_Aperture` is the template). They are intentionally excluded from the reversibility file — aperture filtering drops particles, which is not a reversible phase-space operation.

### Finding element parameters

Never invent physical parameters — they produce tests that compile but don't exercise the element meaningfully (or blow up at runtime). When adding a new entry:

1. First, look in `examples/` for an input deck or Python script that uses the element; grep for the element name under `examples/`. Copy realistic values from there.
2. If no example uses the element, check how the element is instantiated in the other two test files (e.g. reversibility often shares values with the benchmark file).

For `test_element_serialization.py`, every constructor keyword must be set to a non-default value so `test_element_dict_has_constructor_params` can verify `to_dict()` round-trips them — use the style already in the file (`dx=0.001`, `dy=0.002`, `rotation=0.05`, small `aperture_x/y`, `name="test_<element>"`).
For `test_reversibility_elements.py`, use the module-level `ALIGNMENT_KWARGS` (for thin kicks / rotations / edges) or `PIPE_KWARGS` (for thick elements with apertures).

### Adding entries — per-file templates

**test_benchmark_elements.py** — follow the surrounding test style:

```python
def test_<Element>(benchmark, sim):
    el = elements.<Element>(name="...", <params-from-examples>, nslice=nslice)
    benchmark.pedantic(el.push, setup=partial(pc_setup, sim), rounds=rounds)
```

Group the new test near alphabetically-adjacent existing tests. Use `mapsteps=mapsteps` for elements that take `mapsteps`, and `nslice=nslice` for thick elements.

**test_element_serialization.py** — append one `lattice.append(elements.<Element>(...))` inside the `all_elements` fixture, alphabetically ordered relative to neighbors. Include every constructor keyword (including alignment kwargs `dx`, `dy`, `rotation`, aperture kwargs if applicable, and `name=`); `test_element_dict_has_constructor_params` asserts that `to_dict()` keys match the constructor signature, so a minimal instance will fail if any constructor kwarg is omitted from `to_dict()`.

**test_reversibility_elements.py** — place the new test in the correct section:

- "Thick elements (negate ds)" for anything with `ds`.
- "Thin kick elements (negate strength)" for thin-kick elements without `ds`.
- "Rotation elements" for pure rotations.
- "RF / energy-changing elements" for RF cavities.
- "Edge elements" for edge kicks (`DipEdge`, `QuadEdge`).

Use `PIPE_KWARGS` for elements that accept `aperture_x`/`aperture_y`, otherwise `ALIGNMENT_KWARGS`. Wrap the call in `roundtrip(el, sim, ...)`. Only relax `phase_atol` when a physically motivated reason exists (multi-step map, non-symplectic split, energy change); copy the rationale comment style already present.

## Pass 2 — spin parametrization audit

An element is spin-capable iff its header under `src/elements/<Name>.H` inherits `public mixin::SpinTransport`. To find the current set in one shot:

# via Grep tool
Grep: "mixin::SpinTransport" under src/elements/ (files_with_matches)
```

For each spin-capable element that already has a test in files 1 and 3, ensure the `sim` fixture is parametrized:

```python
@pytest.mark.parametrize("sim", [True, False], indirect=True, ids=["spin", "nospin"])
def test_<Element>(benchmark, sim):
    ...
```

In the reversibility file, additionally thread `spin=sim.spin` into the `roundtrip(...)` call.
When adding alignment kwargs breaks the spin path (some elements currently use `kwargs = {} if sim.spin else PIPE_KWARGS`), follow the same conditional pattern already in the file rather than forcing a unified call and WARN the user explicitly in the summary of the change.

The serialization file (file 2) is not parametrized over spin — it only checks construction and dict round-trip. Do not add spin parametrize marks there.

Elements that are NOT spin-capable (no `SpinTransport` inheritance) must not be parametrized — the parametrization will still "work" but it just runs the same `nospin` path twice, wasting CI time.

## Pass 3 — run the tests

After editing the Python test files, install the freshly-built Python bindings into the active environment and run the three test files directly with `pytest`.

```bash
# 1. Make sure a suitable micromamba env is active. Preference order:
#    1. impactx-cpu-mpich-dev   (project-specific, preferred)
#    2. warpx-cpu-mpich-dev     (BLAST developer's global default, works as fallback)
#    Do not re-activate if one is already active.

# 2. Install/refresh the Python package from the build tree.
cmake -S . -B build_cpu -DImpactX_PYTHON=ON -DImpactX_PYTHON_IPO=OFF -DpyAMReX_IPO=OFF
cmake --build build_cpu -j 6 --target pip_install

# 3. Run the three tests directly with pytest.
python -m pytest tests/python/test_benchmark_elements.py
python -m pytest tests/python/test_element_serialization.py
python -m pytest tests/python/test_reversibility_elements.py
```

Notes:
- `pip_install` is the full install (includes dependency checks). Use `pip_install_nodeps` instead on repeated Python-only edits to skip the dependency step and save time.
- `pytest` requires the package to be importable, which is what `pip_install` ensures — running `pytest` without it will pick up a stale install or fail to import `impactx`.
- Parallel build core count: `-j 6` per user CLAUDE.md ("use at most 6 parallel cores").
- If a specific element test fails, re-run with `-k <TestName>` for isolation, e.g. `pytest test_reversibility_elements.py -k ExactCFbend -x`.
- Report per-file pass/fail counts. On failure, show the relevant assertion and the element name; do NOT loosen tolerances in the analysis/asserts to make tests pass unless the user asks for it (per repo `CLAUDE.md`).

## Deliverables checklist

- [ ] Listed the diff between `KnownElements` and each of the three test files.
- [ ] Added missing tests using parameter values sourced from `examples/` (or documented why omission is intentional).
- [ ] Added or verified the `spin/nospin` parametrize decorator on every spin-capable element in files 1 and 3.
- [ ] Ran `pip_install` and executed the three tests with `pytest`; reported results.
- [ ] Did NOT edit any `.pyi` stub files (auto-generated on `development` — see repo `CLAUDE.md`).
- [ ] WARN user about any omissions/skips/tweaks needed to pass the new tests, e.g., alignment parameter issues.
