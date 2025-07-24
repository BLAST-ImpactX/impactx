"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

import sys

from .start import DashboardApp

if __name__ == "__main__":
    app = DashboardApp()
    sys.exit(app.start())
