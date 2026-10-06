"""Running svZeroDSolver on a case, and reading back what it wrote.

The solver is run as its own executable, `svzerodsolver <input.json> <output.csv>`, rather than
imported as `pysvzerod`. Two reasons. It is C++, and a model it cannot solve ends in an abort,
which in-process takes the host with it -- under Slicer, the whole application and an unsaved
scene. And `pysvzerod` is a compiled extension that has to match the host's Python exactly,
which under Slicer is Slicer's own, so it would need building for every Slicer release; the
executable needs nothing from it. A 0D model of a few dozen vessels solves in well under a
second either way, so nothing is lost by the round trip through files.
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

SOLVER_EXECUTABLE_NAME = "svzerodsolver"

# The columns svZeroDSolver writes by default: one row per vessel per time point.
RESULT_COLUMNS = ("name", "time", "flow_in", "flow_out", "pressure_in", "pressure_out")


class SolverError(RuntimeError):
    """Raised when the solver cannot be found, fails, or writes something unreadable."""


def find_solver(configured: str | None = None) -> str | None:
    """The executable to run: the one configured, or `svzerodsolver` on the PATH, or None."""
    if configured:
        configured = os.path.expanduser(configured)
        return configured if os.path.isfile(configured) and os.access(configured, os.X_OK) else None
    return shutil.which(SOLVER_EXECUTABLE_NAME)


def run_solver(config_path, results_path, executable: str, timeout: float | None = 600, env=None):
    """Solve the case and write the results, raising with the solver's own last words if it fails.

    The results file is removed first: the solver writes it only on success, so an old one left
    in place would be read as this run's.

    :param env: the environment to run it in; see `visualization.launch` for why a host that has
      changed its own passes the one it started with.
    """
    results_path = Path(results_path)
    if results_path.is_file():
        results_path.unlink()
    try:
        completed = subprocess.run(
            [executable, str(config_path), str(results_path)],
            capture_output=True, text=True, timeout=timeout,
            cwd=str(Path(config_path).parent), env=env,
        )
    except FileNotFoundError:
        raise SolverError(f"No solver at {executable}.") from None
    except subprocess.TimeoutExpired:
        raise SolverError(f"svZeroDSolver did not finish within {timeout:g} s.") from None
    if completed.returncode != 0 or not results_path.is_file():
        said = (completed.stderr or completed.stdout or "").strip().splitlines()[-8:]
        raise SolverError("svZeroDSolver failed"
                          + (f" (exit {completed.returncode})" if completed.returncode else "")
                          + (": " + " ".join(said) if said else "."))
    return completed


@dataclass(frozen=True)
class VesselResult:
    """One vessel's flows and pressures at its two ends, over the cycle the solver wrote."""

    time: np.ndarray
    flow_in: np.ndarray
    flow_out: np.ndarray
    pressure_in: np.ndarray
    pressure_out: np.ndarray


def read_results(path) -> dict:
    """`{vessel name: VesselResult}` from svZeroDSolver's CSV.

    By default the solver writes the last cardiac cycle only, which is the one past the
    transient, so a mean over what is here is a cycle average.
    """
    rows = {}
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if header is None or tuple(header[:len(RESULT_COLUMNS)]) != RESULT_COLUMNS:
            raise SolverError(f"{path} does not start with svZeroDSolver's columns "
                              f"{', '.join(RESULT_COLUMNS)}; it starts with {header}.")
        for row in reader:
            if row:
                rows.setdefault(row[0], []).append([float(value) for value in row[1:6]])
    results = {}
    for name, values in rows.items():
        array = np.asarray(values, dtype=float)
        results[name] = VesselResult(*(array[:, column] for column in range(5)))
    if not results:
        raise SolverError(f"{path} holds no results.")
    return results
