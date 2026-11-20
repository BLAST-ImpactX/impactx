"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""


def get_element_color_palette(palette="cern-lhc", plot_library="mpl"):
    """Return a dictionary with colors for all elements.

    The keys are sub-strings of the element kinds: the first key that is in the kind of an
    element sets its color, ``"other"`` the color of all other elements.

    Parameters
    ----------
    palette : str
        ``"cern-lhc"``, or ``"fnal-booster"``: dipoles blue, quadrupoles red, sextupoles
        (multipoles) yellow and combined-function magnets green.
    plot_library : str
        ``"mpl"`` for matplotlib colors, else plain color names.
    """
    color_palette = {
        "cern-lhc": {
            "Quad": "tab:blue",
            "Multipole": "tab:orange",
            "Sbend": "tab:green",
            "CFbend": "tab:olive",  # TODO: improve and plot as two on top of each other
            "ConstF": "tab:red",
            "ChrPlasmaLens": "tab:red",
            "SoftSolenoid": "tab:red",
            "TaperedPL": "tab:red",
            "RFCavity": "tab:brown",
            "ShortRF": "tab:brown",
            "Buncher": "tab:purple",
            "Aperture": "black",
            "Kicker": "tab:pink",
            # 'tab:cyan'
            "other": "tab:gray",
        },
        # the Fermilab color scheme
        "fnal-booster": {
            # combined-function (gradient) magnets
            "CFbend": "tab:green",
            # dipoles
            "Sbend": "tab:blue",
            "ThinDipole": "tab:blue",
            # quadrupoles
            "Quad": "tab:red",
            # sextupoles
            "Multipole": "gold",
            # other elements: colors that differ from the magnets above
            "ConstF": "tab:cyan",
            "ChrPlasmaLens": "tab:cyan",
            "SoftSolenoid": "tab:cyan",
            "TaperedPL": "tab:cyan",
            "RFCavity": "tab:brown",
            "ShortRF": "tab:brown",
            "Buncher": "tab:purple",
            "Aperture": "black",
            "Kicker": "tab:pink",
            "other": "tab:gray",
        },
    }

    if palette not in color_palette:
        raise ValueError(
            f"Unknown palette '{palette}', choose one of: "
            + ", ".join(f"'{name}'" for name in color_palette)
            + "."
        )
    colors = color_palette[palette]

    if plot_library != "mpl":
        # plain color names: without the "tab:" prefix of matplotlib
        colors = {k: v.removeprefix("tab:") for k, v in colors.items()}

    return colors


def get_element_color(element_kind: str, palette="cern-lhc", plot_library="mpl"):
    """Get the color for a given element type string."""
    color_palette = get_element_color_palette(palette, plot_library)

    # sub-string matching of keys
    found_keys = [key for key in color_palette.keys() if key in element_kind]

    if found_keys:
        first_found = found_keys[0]
        return color_palette[first_found]
    else:
        return color_palette["other"]
