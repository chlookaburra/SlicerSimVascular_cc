"""Put the packages on the path as Slicer does, and say where the optional parts come from.

`svromsetup` sits beside this module and `svmeshcomplete` beside SimVascular Mesh Prep, which
is how both import under Slicer: every scripted module's folder is on its path. Neither needs
installing for the tests.

`svromutils` is not in this repository; it is on PyPI. The tests that need it -- the centerlines
and the solver input -- run where it is installed (`pip install svromutils`, which brings pip's
VMTK with it), and are skipped elsewhere, as are the ones that run svZeroDSolver, which is found
on the PATH or at `SVZERODSOLVER`.
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE_DIR = HERE.parents[1]
sys.path.insert(0, str(MODULE_DIR))
sys.path.insert(0, str(MODULE_DIR.parent / "SimVascularMeshPrep"))
