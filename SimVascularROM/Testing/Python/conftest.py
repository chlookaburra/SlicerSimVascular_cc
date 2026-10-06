"""Put the packages on the path as Slicer does, and say where the optional parts come from.

`svromsetup` sits beside this module and `svmeshcomplete` beside SimVascular Mesh Prep, which
is how both import under Slicer: every scripted module's folder is on its path. Neither needs
installing for the tests.

`sv_rom_simulation` is not in this repository. The tests that need it -- the centerlines and the
solver input -- run where it is importable, from an install or from a checkout named by
`SV_ROM_SIMULATION_PATH`, and are skipped elsewhere, as are the ones that run svZeroDSolver,
which is found on the PATH or at `SVZERODSOLVER`.
"""

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE_DIR = HERE.parents[1]
sys.path.insert(0, str(MODULE_DIR))
sys.path.insert(0, str(MODULE_DIR.parent / "SimVascularMeshPrep"))
if os.environ.get("SV_ROM_SIMULATION_PATH"):
    sys.path.insert(0, os.path.expanduser(os.environ["SV_ROM_SIMULATION_PATH"]))
