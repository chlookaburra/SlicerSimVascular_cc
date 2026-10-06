"""svZeroDVisualization: the 0D network drawn as a graph, each block's results a click away.

It is a Dash application in svZeroDSolver's `applications/svZeroDVisualization`, started as
`python visualize_simulation.py <input.json> <output dir>` and served on a local port. It is run
as a process of its own, with a Python that has its dependencies -- `pysvzerod`, `dash`,
`plotly`, `pandas`, `networkx`, `pydot` -- because none of them are Slicer's and installing a web
framework into Slicer's Python to draw one graph is the wrong way round. A conda environment
built for svZeroDSolver usually has them all.

It solves the model again itself, through `pysvzerod`, rather than reading the results this
case already has: that is how the application works, and for a model this size it costs
nothing.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

SCRIPT_NAME = "visualize_simulation.py"

# Dash's own default, which the application does not change.
URL = "http://127.0.0.1:8050/"


class VisualizationError(RuntimeError):
    """Raised when the application cannot be started."""


def find_script(configured: str) -> Path | None:
    """The application's script, given it or the folder it is in, or an svZeroDSolver checkout."""
    if not configured:
        return None
    path = Path(os.path.expanduser(configured))
    for candidate in (path, path / SCRIPT_NAME,
                      path / "applications" / "svZeroDVisualization" / SCRIPT_NAME):
        if candidate.is_file() and candidate.name == SCRIPT_NAME:
            return candidate
    return None


def launch(python: str, script, config_path, output_dir, env=None) -> subprocess.Popen:
    """Start the application on a case's solver input, returning the running process.

    Run from the script's own folder, because it imports its helpers from beside itself, and
    with its output going to a log in `output_dir`, which is where to look when the page does not
    come up: a missing dependency is reported there and nowhere else.

    :param env: the environment to run it in. A host that has changed its own -- Slicer's launcher
      points PYTHONHOME, PYTHONPATH and the library path at Slicer's Python -- has to pass the one
      it started with, or the other Python is run on Slicer's packages and finds none of its own.
    """
    python = os.path.expanduser(python or "")
    if not python or not os.path.isfile(python):
        raise VisualizationError(f"No Python at {python!r} to run svZeroDVisualization with.")
    script = Path(script)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log = open(output_dir / "svZeroDVisualization.log", "w")
    try:
        # A session of its own, so that stopping it reaches the server too: Dash's debug mode,
        # which the application runs in, serves from a child of the process started here.
        return subprocess.Popen(
            [python, str(script), str(Path(config_path).resolve()), str(output_dir.resolve())],
            cwd=str(script.parent), stdout=log, stderr=subprocess.STDOUT, env=env,
            start_new_session=(os.name == "posix"),
        )
    finally:
        log.close()


def stop(process) -> None:
    """Stop a running instance and the server its reloader started, which holds the port."""
    if process is None or process.poll() is not None:
        return
    if os.name == "posix":
        import signal
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            return
    else:
        process.terminate()
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        process.kill()
