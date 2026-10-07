# SimVascular ROM Simulation

## Summary

Sets up a 0D (reduced-order) blood flow simulation of a vascular model whose faces are named
in [SimVascular Mesh Prep](SimVascularMeshPrep.md), runs it with
[svZeroDSolver](https://github.com/SimVascular/svZeroDSolver), and shows the results in the
scene. You give one boundary condition per cap and pick the source; the module traces
centerlines from the source with [SlicerVMTK](https://github.com/vmtk/SlicerExtension-VMTK),
builds the 0D model along them with `sv_rom_simulation`
([svROMSimulation](https://github.com/ktbolt/svROMSimulation)), writes the svZeroDSolver input
file, runs the solver, and reads the CSV back.

A 0D model replaces the three-dimensional flow with a network: each vessel segment becomes a
resistance, an inductance and a capacitance worked out from its length and cross-section. It
solves in well under a second, which makes it the model to tune boundary conditions on, or to
check a 3D case against before running it.

The results come back three ways: the pressure over the last cardiac cycle of the caps you
select, in Slicer's plot view; the centerlines drawn inside the anatomy, coloured by
cycle-averaged pressure; and a CSV of every cap's flow and pressure, exported from the panel. [svZeroDVisualization](https://simvascular.github.io/documentation/rom_simulation.html#0d-solver-visualization)
can be opened on the same model from the panel.

## Tutorial

1. **Name the faces.** Open the mesh in SimVascular Mesh Prep. A mesh out of Clip Vessel
   arrives named after its clip points, and opening it there is enough to save those names on
   it. Caps are the faces named `cap_*`; the wall is `wall` or `wall_*`.
2. **Open SimVascular ROM Simulation** and select the mesh in **Labelled mesh**. This can be the
   volume mesh CFD Mesh Generator made, or a labelled surface, since a 0D model needs only the
   surface. The **Boundary conditions** table lists the caps, and selecting a row highlights
   that cap in the 3D view.
3. **Give every cap a condition.** Choose its kind in the **Type** column, then fill in
   **Values**; the hint below the table gives the order.
   - **Inflow**: double-click the values to load a waveform from a `.flow` file; the cell then
     shows that file's path. For a steady flow instead, select the cell and type `Q`. The
     waveform itself is saved with the scene, so moving or deleting the file later does not
     change the condition.
   - **RCR**: double-click and type `Rp, C, Rd, Pd`, i.e. proximal resistance, compliance,
     distal resistance and distal pressure. `Pd` can be left out, and is then 0.
   - **Resistance**: double-click and type `R`. A resistance drains to zero pressure and has no
     distal pressure, as svMultiPhysics' Resistance condition has none and the same conditions
     are written for the 3D solver; an outlet that needs a pressure to drain to is an RCR.

   Values are in cgs units: flow in mL/s, resistance in dyn·s/cm⁵, compliance in cm⁵/dyn,
   pressure in dyn/cm² (1 mmHg = 1333.22 dyn/cm²). A value that cannot be read is shown in
   red; hover over it to see why.
4. **Choose the source** in the **Centerlines** section. The centerlines are traced from this
   cap, and it needs an inflow; the largest cap is suggested. On a model with several inflows
   any of them will do: the other inflows are prescribed where their own centerlines end, and
   the answer does not depend on which one is the source. **Compute centerlines** (optional)
   shows the branch-split centerlines, so you can check that every cap was reached. Running the
   simulation computes them if they are not up to date.
5. **Set up the simulation.** Choose an **Output folder**; by default it is `rom/<mesh>` next to the
   saved scene. Check the blood density and viscosity, and the number of **Cardiac cycles** and
   **Time points per cycle**. Only the last cycle is reported, so the cycles before it give the
   compliances time to settle. The mesh is taken to be in millimetres, which is what Slicer's
   coordinates are, and is scaled to the solver's centimetres on the way.
6. **Create solver files** (optional) writes the boundary condition files and svZeroDSolver's
   `solver_0d.json` into the output folder without running anything, computing the
   centerlines first if they are not current. Look at them, or change `solver_0d.json` by hand.
7. **Run simulation.** This runs svZeroDSolver on the solver files as they are, so an edit made by
   hand since creating them is what gets solved. The files are created first if they are not
   there, or if anything in the panel has changed since they were. The status line reports the
   total flow in and out over the last cycle, which should match. Select caps in the
   **Boundary conditions** table to plot them (the source is plotted when nothing is selected);
   the centerlines in the 3D view are coloured by cycle-averaged pressure. The results model also carries the flow, which the Models module can
   colour it by instead. **Export results** writes a CSV of every cap's flow and pressure over
   the last cycle: a time column, then a flow (mL/s) and a pressure (mmHg) column per cap. Flow
   is positive the way the cap's condition drives it, into the model at an inflow and out of it
   at an outlet.
8. **Visualize results** (optional) opens svZeroDVisualization, which draws the network as a
   graph in your browser, with every block's parameters and results a click away.

The conditions and the source are saved on the mesh, so they come back when you reopen the
scene, and two anatomies in one scene keep separate conditions. The conditions are shared with
[SimVascular MultiPhysics](SimVascularMultiPhysics.md): set in either panel, they are the ones
both write, so a 0D model tuned here is the 3D case written there. An output folder that already
holds results shows them again when the mesh is selected.

### What gets written

```
<output folder>/
  mesh-complete/                the faces, as Mesh Prep's Export writes them
    mesh-complete.exterior.vtp
    mesh-surfaces/cap_*.vtp
  centerlines.vtp               branch-split centerlines
  inflow.flow                   the source's flow
  flow.dat, cap_*.flow          the other inflows
  rcrt.dat, resistance.dat      the outlets
  solver_0d.json                svZeroDSolver's input
  results.csv                   svZeroDSolver's output
```

The boundary conditions in `solver_0d.json` are named after their faces (`RCR_cap_lpa_a`,
`FLOW_cap_azygous_vein`), except the source's, which is `INFLOW`. The solver's results are per
vessel, named `branch<b>_seg<s>`; the panel works out which vessel end each cap is.

### Several inflows

A Fontan has the venae cavae, the azygous and the hepatic veins all entering, and the pulmonary
arteries all leaving. The centerlines run from the source, so every other inflow sits at a
centerline *end*. There it is prescribed as a flow entering the model. The resistances and
junctions of a 0D network do not care which way flow goes, so this is the same model as one
drawn from any other inflow.

Every inflow has to cover the same cardiac cycle. A steady inflow is spread over the period of
the waveforms next to it.

### Flow files and their sign

A `.flow` file is a time and a flow on each line, over one cycle starting at 0. Here an inflow
is positive. A file written for svMultiPhysics gives inflow as *negative*, because its flow is
measured along the outward normal. A loaded waveform whose mean is negative is therefore turned
round, and the status line says so.

## Setting up

Install **SlicerVMTK** from the Extensions Manager. Then, under **Tools** in the panel:

- **svROMSimulation**: a checkout of [svROMSimulation](https://github.com/ktbolt/svROMSimulation),
  the folder holding `sv_rom_simulation`. Not needed if it is installed into Slicer's Python.
  **Do not** pip-install its `vmtk` dependency into Slicer (see below).
- **svzerodsolver**: the svZeroDSolver executable. If left empty, `svzerodsolver` on the PATH
  is used.
- **svZeroDVisualization** and **Python environment**: an svZeroDSolver checkout (or its
  `visualize_simulation.py`), and the `python` of an environment that has `pysvzerod`, `dash`,
  `plotly`, `pandas`, `networkx` and `pydot` (a conda environment's `bin/python`, for
  instance). A conda environment built for svZeroDSolver usually has them all.

These are saved in the application's settings rather than the scene: a path belongs to one
computer, and a scene moves between them.

## Developers

### The split

`SimVascularROM.py` is an MRML adapter over `svromsetup`, the package beside it, which imports
nothing from Slicer: the conditions and SimVascular's files for them, the case folder, the calls
into `sv_rom_simulation` and the checks on what it wrote, the solver run, and the results read
back. A case set up from a terminal calls the same functions:

```python
from svromsetup import case, solver, results
from svromsetup.boundary_conditions import Inflow, Resistance

case.write_surfaces(volume_mesh, names, "case")             # {face id: name}, as Mesh Prep saves them
case.compute_centerlines("case", inlet_face_id=3)
config = case.write_solver_input("case", {"cap_RSVC": Inflow.steady(20.0), ...}, 3, "cap_RSVC",
                                 case.SimulationParameters(), "fontan")
solver.run_solver(config, "case/results.csv", "svzerodsolver")
faces = results.face_results(config, solver.read_results("case/results.csv"), "cap_RSVC")
```

`sv_rom_simulation` is imported only inside the two functions that call it, because it imports
VMTK. Everything else in `svromsetup` needs only numpy and VTK, and so do its tests.

### VMTK under Slicer

`sv_rom_simulation` imports `from vmtk import vtkvmtk`, which is how pip's VMTK is laid out.
pip's VMTK is **not** installed into Slicer, because it pins a VTK of its own (`vtk==9.6.2` for
1.5.2), and a second VTK over Slicer's breaks the application. SlicerVMTK carries the same
classes in modules of its own (`vtkvmtkComputationalGeometryPython`, …). So the panel assembles
a `vmtk.vtkvmtk` module out of them before importing the package (`ensureVmtk`). The branch
splitting the package relies on, `vtkvmtkPolyDataCenterlineBranchSplitting`, has been in
SlicerVMTK's VMTK since the update of September 2026.

### Why the solver is a process

svZeroDSolver is run as its own executable, not imported as `pysvzerod`. A model it cannot
solve ends in a C++ abort, which in-process would take Slicer and an unsaved scene with it.
`pysvzerod` is also a compiled extension that would have to be built for Slicer's own Python.
Both it and svZeroDVisualization are started with the environment Slicer was launched with
(`slicer.util.startupEnvironment()`). Slicer's launcher points `PYTHONHOME`, `PYTHONPATH` and
the library path at its own Python; inherited, another Python finds Slicer's packages and
none of its own.

### Outlets paired by position

The package matches outlet names to boundary conditions by list position. Nothing ties the
order of the centerline ends to the order of the caps, and an off-by-one list puts every
condition on its neighbour's vessel, with all the ids still valid. So the centerline ends are
paired with the cap files by where they are: one to one, closest first, and refused if an end
is further than a cap's width from its cap. `check_solver_input` then holds the written file to
what was asked: one condition per cap, on one vessel end, of the right kind, over the requested
cycles.

### Tests

The package's tests are headless (pytest, from `SimVascularROM/`):

```sh
SV_ROM_SIMULATION_PATH=~/Documents/svROMSimulation python -m pytest
```

`svromsetup.testing` builds a capped Y, with arms of different radii so that its two outlets can
be told apart by their results. The tests that need `sv_rom_simulation` (and VMTK) or the solver
skip themselves where these are missing; with only numpy, VTK and pytest, 26 run. The module's
own tests (`SimVascularROMTest`) run under Slicer: they check that conditions are kept per mesh,
that an unnamed mesh is sent to Mesh Prep, and a whole case on the Y. A `--testing` run loads
none of the installed extensions, so pass SlicerVMTK's module folders as
`--additional-module-paths`, along with `SimVascularMeshPrep` and `SimVascularROM`.
