"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""


def get_element_color_palette(palette="cern-lhc", plot_library="mpl"):
    """Return a dictionary with colors for all elements.

    The key is a regex that can be matched against the element type string. TODO TODO
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
        }
    }

    colors = color_palette[palette]

    if plot_library != "mpl":
        # remove "tab:" prefix
        for k, v in colors.items():
            colors[k] = v[4:]

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
