"""
Fermilab Booster lattice translated from the MAD-X description in
``sbbooster-cooked-rfon.madx`` to ImpactX elements.

The ring has 24 super-periods (cells). The first 13 cells use a plain ``dlong``
drift in the middle; the remaining 11 cells (14-24) replace that drift with an
RF station (two ``ShortRF`` cavities separated by drifts). All cells are
otherwise identical; only the F/D bend element names carry the cell index
suffix (e.g. ``fmagu01`` ... ``fmagu24``).
"""

from impactx import elements

# ---------------------------------------------------------------------------
# shared parameters for F and D combined-function bends
# ---------------------------------------------------------------------------

_F_EDGE = dict(
    K0=1.6449340668482264,
    K1=0.0,
    K2=0.0,
    K3=0.16666666666666666,
    K4=0.0,
    K5=0.0,
    K6=0.0,
    R=1.0,
    dx=0.0,
    dy=0.0,
    g=0.0,
    model="linear",
    psi=0.0353710910981508,
    rc=40.847086,
    rotation=0.0,
)
_F_CFBEND = dict(
    aperture_x=0.0,
    aperture_y=0.0,
    ds=2.889612,
    dx=0.0,
    dy=0.0,
    int_order=6,
    k_normal=[0.024481550532148122, 0.05410921561, 0.00546850017312],
    k_skew=[0.0, 0.0, 0.0],
    mapsteps=6,
    rotation=0.0,
    unit=0,
)
_D_EDGE = dict(
    K0=1.6449340668482264,
    K1=0.0,
    K2=0.0,
    K3=0.16666666666666666,
    K4=0.0,
    K5=0.0,
    K6=0.0,
    R=1.0,
    dx=0.0,
    dy=0.0,
    g=0.0,
    model="linear",
    psi=0.03007875592383836,
    rc=48.034101,
    rotation=0.0,
)
_D_CFBEND = dict(
    aperture_x=0.0,
    aperture_y=0.0,
    ds=2.889612,
    dx=0.0,
    dy=0.0,
    int_order=6,
    k_normal=[0.020818543059648396, -0.05738855012, -0.037869227516832],
    k_skew=[0.0, 0.0, 0.0],
    mapsteps=6,
    rotation=0.0,
    unit=0,
)


# ---------------------------------------------------------------------------
# small builders for the repetitive elements inside a cell
# ---------------------------------------------------------------------------


def _drift(name, ds):
    return elements.ExactDrift(
        aperture_x=0.0,
        aperture_y=0.0,
        ds=ds,
        dx=0.0,
        dy=0.0,
        name=name,
        rotation=0.0,
    )


def _quad(name, k):
    return elements.ExactQuad(
        aperture_x=0.0,
        aperture_y=0.0,
        ds=0.024,
        dx=0.0,
        dy=0.0,
        int_order=6,
        k=k,
        mapsteps=6,
        name=name,
        rotation=0.0,
        unit=0,
    )


def _sextupole(name, k_normal, k_skew):
    return elements.ExactMultipole(
        aperture_x=0.0,
        aperture_y=0.0,
        ds=0.024,
        dx=0.0,
        dy=0.0,
        int_order=4,
        k_normal=k_normal,
        k_skew=k_skew,
        mapsteps=5,
        name=name,
        rotation=0.0,
        unit=0,
    )


def _bend_station(name, edge, cfbend):
    """[entry edge, CFbend, exit edge] for one combined-function dipole."""
    return [
        elements.DipEdge(name=f"{name}_usedge", location="entry", **edge),
        elements.ExactCFbend(name=name, **cfbend),
        elements.DipEdge(name=f"{name}_dsedge", location="exit", **edge),
    ]


def _rf_station():
    """Two ShortRF cavities with their neighboring drifts. Replaces ``dlong``
    in cells 14..24. The total length equals the length of ``dlong`` (5.581 m).
    """
    rfc_kwargs = dict(
        V=9.688990660720476e-06,
        dx=0.0,
        dy=0.0,
        freq=44704409.98524446,
        phase=-90.0,
        rotation=0.0,
    )
    return [
        _drift("drifta", 0.21),
        _drift("drrf", 1.175),
        elements.ShortRF(name="rfc", **rfc_kwargs),
        _drift("drrf", 1.175),
        _drift("driftb", 0.12),
        _drift("drrf", 1.175),
        elements.ShortRF(name="rfc", **rfc_kwargs),
        _drift("drrf", 1.175),
        _drift("dmidls", 0.551),
    ]


# ---------------------------------------------------------------------------
# one full super-period
# ---------------------------------------------------------------------------


def _cell(i, rf):
    """Build super-period ``i`` (1..24). If ``rf`` is True, substitute the
    RF station for the ``dlong`` drift (cells 14..24)."""
    suf = f"{i:02d}"
    middle = _rf_station() if rf else [_drift("dlong", 5.581)]
    return [
        # upstream short-straight section (small-aperture correctors + diagnostics)
        _drift("sa", 0.176),
        _drift("hsxx", 0.024),
        _drift("vsxx", 0.024),
        _quad("qsxx", 0.020254116087679828),
        _drift("bpms", 0.024),
        _drift("qssxx", 0.024),
        _sextupole("sxsxx", [0.0, 0.0, -0.09001070603982472], [0.0, 0.0, 0.0]),
        _drift("sssxx", 0.024),
        _drift("sb", 0.256),
        # upstream F bend
        *_bend_station(f"fmagu{suf}", _F_EDGE, _F_CFBEND),
        _drift("mins", 0.5),
        # upstream D bend
        *_bend_station(f"dmagu{suf}", _D_EDGE, _D_CFBEND),
        # middle: plain drift ``dlong`` (cells 1..13) or RF station (cells 14..24)
        *middle,
        # downstream long-straight section (large-aperture correctors + diagnostics)
        _drift("hlxx", 0.024),
        _drift("vlxx", 0.024),
        _quad("qlxx", -0.0738960996629951),
        _drift("bpms", 0.024),
        _drift("qlsxx", 0.024),
        _sextupole("sxlxx", [0.0, 0.0, 8.591792079114004], [0.0, 0.0, -0.0]),
        _drift("sslxx", 0.024),
        _drift("drifte", 0.251),
        # downstream D bend
        *_bend_station(f"dmagd{suf}", _D_EDGE, _D_CFBEND),
        _drift("mins", 0.5),
        # downstream F bend
        *_bend_station(f"fmagd{suf}", _F_EDGE, _F_CFBEND),
        _drift("sc", 0.6),
    ]


def get_lattice():
    """Return the 24-cell Fermilab Booster lattice as a list of ImpactX elements."""
    lattice = []
    for i in range(1, 25):
        lattice.extend(_cell(i, rf=(i >= 14)))
    return lattice
