"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

import sys

from .start import DashboardApp


def main() -> int:
    """
    Entry point for dashboard.
    """
    app = DashboardApp()
    return app.start()


if __name__ == "__main__":
    sys.exit(main())
