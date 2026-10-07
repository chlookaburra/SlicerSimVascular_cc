"""Set up a 0D simulation of a labelled vascular model, run it, and look at what it says.

A reduced-order model replaces the three-dimensional flow with a network: the vessels are cut
along their centerlines into segments, each segment is a resistance, an inductance and a
capacitance worked out from its length and cross-section, and the network is closed by one
boundary condition per cap. svZeroDSolver solves it in well under a second, which is what makes
it the model to tune boundary conditions on, check a 3D case against, or run where a 3D one
cannot be afforded.

This panel is the setting up and the looking. The mesh is one whose faces SimVascular Mesh Prep
has named -- the names are read off it, never typed again here -- and the conditions are set per
cap in a table, an inflow's waveform loaded by double-clicking its values. The centerlines come
from SlicerVMTK, traced from the inlet the operator picks; the model is built and written by
`sv_rom_simulation`; svZeroDSolver is run as a process; and what it says comes back as a plot
over the cycle of the caps selected in the table, the centerlines coloured by pressure, laid
inside the anatomy they came from, and a CSV of every cap to take away.

## Where the work is

In `svromsetup`, the package beside this file, which imports nothing from Slicer: the conditions
and SimVascular's files for them, the case folder, the calls into `sv_rom_simulation` and the
checks on what it wrote, the solver run, the results read back. This file is the MRML adapter --
it reads the selected node, calls the package, and puts the answer in the scene -- so a case set
up here and one set up from a terminal are the same case.

## Several inflows

The centerlines are traced from one inlet, so a model with several inflows -- a Fontan, with the
venae cavae and the hepatic veins all entering -- has the others at centerline *ends*. They are
prescribed there as flows entering the model, which the 0D network does not mind: a vessel's
resistance and a junction's mass balance do not care which way the flow goes. Which inflow is
the inlet makes no difference to the answer, only to how the network is drawn.

## VMTK and the package under Slicer

`sv_rom_simulation` imports VMTK as `from vmtk import vtkvmtk`, which is how pip's VMTK is laid
out. pip's VMTK must not be installed into Slicer: it pins a VTK of its own, and a second VTK
over Slicer's breaks the application. SlicerVMTK carries the same classes in modules of its
own, so `vmtk.vtkvmtk` is assembled from those before the package is imported -- see
`ensureVmtk` -- and the package itself is found where it is installed or where the Tools
section says a checkout is.
"""

import importlib
import json
import logging
import os
import sys
import time
import types

import numpy as np
import qt
import slicer
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
from slicer.i18n import tr as _
from slicer.i18n import translate
from slicer.ScriptedLoadableModule import (
    ScriptedLoadableModule,
    ScriptedLoadableModuleLogic,
    ScriptedLoadableModuleTest,
    ScriptedLoadableModuleWidget,
)
from slicer.util import VTKObservationMixin

from svmeshcomplete.face_table import sanitized
from svromsetup import boundary_conditions as bcs
from svromsetup import case, results, solver, visualization

# The faces and the table of their conditions, shared with SimVascular MultiPhysics -- the names
# are Mesh Prep's and the conditions one set per mesh, whichever solver they are written for.
from SimVascularROMLib.BoundaryConditionsTable import (  # noqa: F401 (the tests use these)
    BOUNDARY_CONDITIONS_ATTRIBUTE,
    FACE_NAMES_ATTRIBUTE,
    INVALID_VALUES_COLOR,
    KIND_CHOICES,
    TYPE_COLUMN,
    VALUES_COLUMN,
    BoundaryConditionsTable,
    NamedMeshFaces,
    clearHighlight,
)

# The cap the centerlines start at, kept on the mesh beside its conditions. The 0D panel's own:
# a 3D model has no source.
INLET_FACE_ID_ATTRIBUTE = "SimVascularROM.InletFaceID"

# The nodes a run leaves in the scene, referenced from the mesh they belong to so that each
# anatomy keeps its own, and so that they are found again in a saved scene.
CENTERLINES_REFERENCE = "ROMCenterlines"
RESULTS_MODEL_REFERENCE = "ROMResultsModel"
RESULTS_TABLE_REFERENCE = "ROMResultsTable"
RESULTS_CHART_REFERENCE = "ROMResultsChart"

# The rest of the panel's state, which is the scene's.
INPUT_MESH_REFERENCE = "InputMesh"
CASE_DIRECTORY_PARAMETER = "CaseDirectory"
DENSITY_PARAMETER = "Density"
VISCOSITY_PARAMETER = "Viscosity"
CYCLES_PARAMETER = "CardiacCycles"
POINTS_PARAMETER = "PointsPerCycle"

# Where each tool is on this machine: application settings, not scene parameters, because a
# path is one computer's and a scene is copied between them.
ROM_PACKAGE_SETTING = "SimVascularROM/svROMSimulationPath"
SOLVER_SETTING = "SimVascularROM/SolverExecutable"
VISUALIZATION_SCRIPT_SETTING = "SimVascularROM/VisualizationScript"
VISUALIZATION_PYTHON_SETTING = "SimVascularROM/VisualizationPython"

# The units the mesh is in, which are Slicer's: its coordinates are millimetres whatever unit it
# displays them in. Both solvers work in centimetres, and the package scales lengths by 0.1 and
# areas by 0.01 on the way.
MESH_UNITS = "mm"

TIME_COLUMN = "time [s]"
PLOT_COLUMN_ATTRIBUTE = "SimVascularROM.Column"
# What the centerlines are coloured by and the caps are plotted in. Pressure, because it is what a
# 0D model is run to find -- the flows at the caps are mostly what the conditions set -- and the
# one quantity that varies along every branch. Flow is in the export, and on the results model
# as an array of its own for the Models module to colour by.
RESULTS_QUANTITY = "Pressure"
RESULTS_ARRAY_NAME = results.PRESSURE_ARRAY_NAME

# The modules SlicerVMTK's classes are in. pip's VMTK gathers the same classes into one module,
# `vmtk.vtkvmtk`, and that is what `sv_rom_simulation` imports.
VMTK_MODULES = (
    "vtkvmtkCommonPython",
    "vtkvmtkComputationalGeometryPython",
    "vtkvmtkDifferentialGeometryPython",
    "vtkvmtkIOPython",
    "vtkvmtkMiscPython",
    "vtkvmtkSegmentationPython",
    "vtkvmtkITKPython",
    "vtkvmtkContribPython",
)

# What a results tube is drawn in: blue low to red high, which is how pressure is read.
RESULTS_COLOR_NODES = ("vtkMRMLColorTableNodeFileColdToHotRainbow.txt",
                       "vtkMRMLColorTableNodeRainbow")

# What the mesh is faded to when results are shown inside it.
RESULTS_MESH_OPACITY = 0.25

# One colour per plotted cap, repeating after ten.
PLOT_COLORS = ((0.12, 0.47, 0.71), (1.0, 0.5, 0.05), (0.17, 0.63, 0.17), (0.84, 0.15, 0.16),
               (0.58, 0.4, 0.74), (0.55, 0.34, 0.29), (0.89, 0.47, 0.76), (0.5, 0.5, 0.5),
               (0.74, 0.74, 0.13), (0.09, 0.75, 0.81))


def ensureVmtk():
    """Make `from vmtk import vtkvmtk` work, out of SlicerVMTK's modules if need be.

    A pip VMTK, where there is one, is left as it is. Under Slicer there is none and must not be
    -- see the module docstring -- so a module of that name is assembled from SlicerVMTK's, which
    carry the same classes under the same names. Done once; the second call finds it.
    """
    try:
        from vmtk import vtkvmtk  # noqa: F401
        return
    except ImportError:
        pass
    namespace = types.ModuleType("vmtk.vtkvmtk")
    found = []
    for name in VMTK_MODULES:
        try:
            module = importlib.import_module(name)
        except ImportError:
            continue
        found.append(name)
        for attribute in dir(module):
            if attribute.startswith("vtkvmtk"):
                setattr(namespace, attribute, getattr(module, attribute))
    if not found:
        raise RuntimeError(_("VMTK is not available. Install the SlicerVMTK extension from the "
                             "Extensions Manager and restart Slicer."))
    package = types.ModuleType("vmtk")
    package.__path__ = []
    package.vtkvmtk = namespace
    sys.modules["vmtk"] = package
    sys.modules["vmtk.vtkvmtk"] = namespace


class SimVascularROM(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = _("SimVascular ROM Simulation")
        self.parent.categories = [translate("qSlicerAbstractCoreModule", "SimVascular")]
        # Mesh Prep names the faces, and its folder is where svmeshcomplete is imported from.
        self.parent.dependencies = ["SimVascularMeshPrep"]
        self.parent.contributors = ["Cardiovascular Biomechanics Computation Lab (Stanford University)"]
        self.parent.helpText = _(
            "Set up a 0D (reduced-order) simulation of a vascular model whose faces are named "
            "in SimVascular Mesh Prep: one boundary condition per cap, centerlines from the "
            "inlet, an svZeroDSolver input file written by sv_rom_simulation, the solver run, "
            "and the results shown per cap, over the cycle, and along the centerlines."
        )
        self.parent.acknowledgementText = _(
            "Developed in the Cardiovascular Biomechanics Computation Lab at Stanford "
            "University. The model is built by sv_rom_simulation and solved by svZeroDSolver; "
            "the centerlines are VMTK's, through SlicerVMTK. The setting up is the svromsetup "
            "package beside this module, which runs outside Slicer as well."
        )


class SimVascularROMWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    def __init__(self, parent=None):
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)
        self.logic = None
        self._updating = False
        # The selected mesh's faces, and the table of their conditions, which are shared with
        # SimVascular MultiPhysics: see SimVascularROMLib.BoundaryConditionsTable.
        self.namedFaces = NamedMeshFaces()
        self.conditionsTable = None
        self._inletFaceId = None
        # What the centerlines in the case folder were computed from, so that a run can tell
        # whether they are still the mesh's: changing the inlet, the folder or a name makes them
        # some other model's.
        self._centerlinesKey = None
        # What the solver files in the folder were written from, for the same reason: Run uses
        # them as they are while this still matches, which is what keeps an edit made to
        # solver_0d.json by hand between creating the files and running them.
        self._solverFilesKey = None
        self._faceResults = []
        self._resultsNote = ""
        self._visualization = None

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)
        uiWidget = slicer.util.loadUI(self.resourcePath("UI/SimVascularROM.ui"))
        self.layout.addWidget(uiWidget)
        self.ui = slicer.util.childWidgetVariables(uiWidget)
        uiWidget.setMRMLScene(slicer.mrmlScene)
        self.logic = SimVascularROMLogic()
        self.conditionsTable = BoundaryConditionsTable(
            self.ui.facesTable, self.ui.valuesHintLabel, self.setStatus,
            changed=self.updateButtons, selectionChanged=lambda _faceIds: self.plotSelected(),
            startDirectory=lambda: self.ui.caseDirectoryPathLineEdit.currentPath)
        self.conditionsTable.emphasisNote = _("The source: the centerlines start here.")

        settings = qt.QSettings()
        self.ui.romPackagePathLineEdit.currentPath = settings.value(ROM_PACKAGE_SETTING, "")
        self.ui.solverPathLineEdit.currentPath = settings.value(SOLVER_SETTING, "")
        self.ui.visualizationScriptPathLineEdit.currentPath = settings.value(
            VISUALIZATION_SCRIPT_SETTING, "")
        self.ui.visualizationPythonPathLineEdit.currentPath = settings.value(
            VISUALIZATION_PYTHON_SETTING, "")
        for lineEdit, key in ((self.ui.romPackagePathLineEdit, ROM_PACKAGE_SETTING),
                              (self.ui.solverPathLineEdit, SOLVER_SETTING),
                              (self.ui.visualizationScriptPathLineEdit, VISUALIZATION_SCRIPT_SETTING),
                              (self.ui.visualizationPythonPathLineEdit, VISUALIZATION_PYTHON_SETTING)):
            lineEdit.connect("currentPathChanged(QString)",
                             lambda path, key=key: qt.QSettings().setValue(key, path))

        self.ui.inputMeshSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onMeshChanged)
        self.ui.openMeshPrepButton.connect("clicked(bool)", self.onOpenMeshPrep)
        self.ui.inletComboBox.connect("currentIndexChanged(int)", self.onInletChanged)
        self.ui.computeCenterlinesButton.connect("clicked(bool)", self.onComputeCenterlines)
        self.ui.createSolverFilesButton.connect("clicked(bool)", self.onCreateSolverFiles)
        self.ui.runButton.connect("clicked(bool)", self.onRun)
        self.ui.exportButton.connect("clicked(bool)", self.onExport)
        self.ui.visualizationButton.connect("clicked(bool)", self.onOpenVisualization)
        self.ui.caseDirectoryPathLineEdit.connect("currentPathChanged(QString)",
                                                  self.onSimulationParameterChanged)
        for spinBox in (self.ui.densitySpinBox, self.ui.viscositySpinBox):
            spinBox.connect("valueChanged(double)", self.onSimulationParameterChanged)
        for spinBox in (self.ui.cyclesSpinBox, self.ui.pointsSpinBox):
            spinBox.connect("valueChanged(int)", self.onSimulationParameterChanged)

        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndImportEvent, self.onSceneEndImport)
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)

        self.restoreFromParameterNode()

    def enter(self):
        """Coming back, perhaps from Mesh Prep with the faces renamed, or from SimVascular
        MultiPhysics with the conditions changed: both are read again."""
        self.restoreFromParameterNode()

    def exit(self):
        clearHighlight()

    def cleanup(self):
        clearHighlight()
        visualization.stop(self._visualization)
        self.removeObservers()

    def onSceneEndImport(self, caller=None, event=None):
        self.restoreFromParameterNode()

    def onSceneEndClose(self, caller=None, event=None):
        self._centerlinesKey = None
        self._solverFilesKey = None
        self._faceResults = []
        self.onMeshChanged()

    # -- state that is saved: the conditions on the mesh, the rest on the scene ---
    def restoreFromParameterNode(self):
        node = self.logic.getParameterNode()
        self._updating = True
        try:
            directory = node.GetParameter(CASE_DIRECTORY_PARAMETER)
            if directory:
                self.ui.caseDirectoryPathLineEdit.currentPath = directory
            for parameter, spinBox, kind in (
                    (DENSITY_PARAMETER, self.ui.densitySpinBox, float),
                    (VISCOSITY_PARAMETER, self.ui.viscositySpinBox, float),
                    (CYCLES_PARAMETER, self.ui.cyclesSpinBox, int),
                    (POINTS_PARAMETER, self.ui.pointsSpinBox, int)):
                value = node.GetParameter(parameter)
                if value:
                    spinBox.value = kind(value)
            mesh = node.GetNodeReference(INPUT_MESH_REFERENCE)
            if mesh is not None:
                # Without its signal, which would measure the mesh once here and once below.
                wasBlocked = self.ui.inputMeshSelector.blockSignals(True)
                self.ui.inputMeshSelector.setCurrentNode(mesh)
                self.ui.inputMeshSelector.blockSignals(wasBlocked)
        finally:
            self._updating = False
        self.onMeshChanged()

    def saveToParameterNode(self):
        node = self.logic.getParameterNode()
        wasModifying = node.StartModify()
        node.SetParameter(CASE_DIRECTORY_PARAMETER, self.ui.caseDirectoryPathLineEdit.currentPath or "")
        node.SetParameter(DENSITY_PARAMETER, str(self.ui.densitySpinBox.value))
        node.SetParameter(VISCOSITY_PARAMETER, str(self.ui.viscositySpinBox.value))
        node.SetParameter(CYCLES_PARAMETER, str(self.ui.cyclesSpinBox.value))
        node.SetParameter(POINTS_PARAMETER, str(self.ui.pointsSpinBox.value))
        node.SetNodeReferenceID(INPUT_MESH_REFERENCE, self.ui.inputMeshSelector.currentNodeID or None)
        node.EndModify(wasModifying)

    def saveInlet(self):
        mesh = self.ui.inputMeshSelector.currentNode()
        if mesh is not None and self.conditionsTable.caps:
            mesh.SetAttribute(INLET_FACE_ID_ATTRIBUTE,
                              str(self._inletFaceId) if self._inletFaceId is not None else "")

    def onSimulationParameterChanged(self, *_args):
        if self._updating:
            return
        self.saveToParameterNode()
        self.updateButtons()

    # -- the mesh -------------------------------------------------------------------
    def onMeshChanged(self, _node=None):
        """Read the selected mesh's names, measure its caps, and bring back its conditions.

        The names are Mesh Prep's, read off the mesh every time rather than remembered, so that a
        rename there arrives here on coming back.
        """
        self._inletFaceId = None
        self._faceResults, self._resultsNote = [], ""
        self.ui.openMeshPrepButton.visible = False
        node = self.ui.inputMeshSelector.currentNode()
        if not self._updating:
            self.saveToParameterNode()
        if node is None or node.GetMesh() is None:
            self.ui.meshStatusLabel.text = ""
            self.showConditions(None)
            return
        try:
            named = self.namedFaces.read(node)
        except ValueError as error:
            self.ui.meshStatusLabel.text = str(error)
            self.showConditions(None)
            return
        self.ui.meshStatusLabel.text = self.namedFaces.describe()
        self.ui.openMeshPrepButton.visible = not named or bool(self.namedFaces.unnamed)
        if not named:
            self.showConditions(None)
            return
        self.showConditions(node)
        capIds = set(self.conditionsTable.capNames())
        try:
            saved = int(node.GetAttribute(INLET_FACE_ID_ATTRIBUTE) or "")
        except ValueError:
            saved = None
        self._inletFaceId = saved if saved in capIds else self.suggestedInlet()
        if not self.ui.caseDirectoryPathLineEdit.currentPath:
            self.ui.caseDirectoryPathLineEdit.currentPath = self.logic.suggestedCaseDirectory(node)
        self.populate()
        self.loadResultsFromCase()

    def showConditions(self, node):
        self.conditionsTable.setMesh(node, self.namedFaces if node is not None else None)
        self.populate()

    def suggestedInlet(self):
        """The largest cap given an inflow, or the largest cap: the inlet of most models.

        Any inflow would do -- the answer is the same, only the network is drawn from a different
        root -- but the largest is the trunk of the tree, and a tree drawn from its trunk is the
        one that reads as anatomy.
        """
        conditions = self.conditionsTable.conditions
        caps = self.conditionsTable.caps
        inflows = [(geometry.area, faceId) for faceId, _name, geometry in caps
                   if isinstance(conditions.get(faceId), bcs.Inflow)]
        candidates = inflows or [(geometry.area, faceId) for faceId, _name, geometry in caps]
        return max(candidates)[1] if candidates else None

    def onOpenMeshPrep(self):
        mesh = self.ui.inputMeshSelector.currentNode()
        slicer.util.selectModule("SimVascularMeshPrep")
        if mesh is not None:
            slicer.util.getModuleWidget("SimVascularMeshPrep").ui.inputMeshSelector.setCurrentNode(mesh)

    # -- the source, and the table it is emphasised in -----------------------------------
    def capName(self, faceId):
        return self.conditionsTable.capName(faceId)

    def populate(self):
        caps = self.conditionsTable.caps
        self._updating = True
        try:
            self.ui.inletComboBox.clear()
            for faceId, name, _geometry in caps:
                self.ui.inletComboBox.addItem(name, faceId)
            self.ui.inletComboBox.currentIndex = next(
                (row for row, (faceId, _name, _geometry) in enumerate(caps)
                 if faceId == self._inletFaceId), -1)
        finally:
            self._updating = False
        self.conditionsTable.emphasised = self._inletFaceId
        self.conditionsTable.populate()
        self.updateButtons()

    def onInletChanged(self, index):
        caps = self.conditionsTable.caps
        if self._updating or index < 0 or index >= len(caps):
            return
        self._inletFaceId = caps[index][0]
        self.saveInlet()
        self.populate()

    # -- centerlines and the run ---------------------------------------------------------
    def caseDirectory(self):
        directory = self.ui.caseDirectoryPathLineEdit.currentPath
        if not directory:
            raise RuntimeError(_("Choose an output folder to write the simulation into."))
        os.makedirs(directory, exist_ok=True)
        return directory

    def centerlinesKey(self):
        node = self.ui.inputMeshSelector.currentNode()
        return (node.GetID(), node.GetMesh().GetMTime(), self._inletFaceId,
                os.path.abspath(self.ui.caseDirectoryPathLineEdit.currentPath or ""),
                json.dumps(sorted(self.namedFaces.names.items())))

    def onComputeCenterlines(self):
        with slicer.util.tryWithErrorDisplay(_("Centerlines could not be computed."), waitCursor=True):
            self.computeCenterlines()

    def computeCenterlines(self):
        """Write the faces into the case, and trace and split the centerlines from the inlet."""
        node = self.ui.inputMeshSelector.currentNode()
        if node is None or not self.conditionsTable.caps:
            raise RuntimeError(_("Select a mesh whose faces are named first."))
        if self._inletFaceId is None:
            raise RuntimeError(_("Choose the source the centerlines start from."))
        directory = self.caseDirectory()
        self.logic.importRomPackage(self.ui.romPackagePathLineEdit.currentPath)
        started = time.time()
        case.write_surfaces(node.GetMesh(), self.namedFaces.names, directory, self.namedFaces.arrayName)
        centerlines = case.compute_centerlines(directory, self._inletFaceId)
        self._centerlinesKey = self.centerlinesKey()
        self.logic.showCenterlines(node, centerlines.geometry)
        branches = len({int(branch) for branch in vtk_to_numpy(
            centerlines.geometry.GetPointData().GetArray("BranchId")) if branch >= 0})
        self.ui.centerlinesStatusLabel.text = _(
            "{branches} branches from {inlet} to {outlets} other caps, each centerline end paired "
            "with the cap it ends at ({seconds:.1f} s).").format(
            branches=branches, inlet=self.capName(self._inletFaceId),
            outlets=len(centerlines.outlet_names), seconds=time.time() - started)
        return centerlines

    def simulationParameters(self):
        return case.SimulationParameters(
            cardiac_cycles=int(self.ui.cyclesSpinBox.value),
            points_per_cycle=int(self.ui.pointsSpinBox.value),
            density=float(self.ui.densitySpinBox.value),
            viscosity=float(self.ui.viscositySpinBox.value),
            units=MESH_UNITS,
        )

    def onCreateSolverFiles(self):
        with slicer.util.tryWithErrorDisplay(_("The solver files could not be created."),
                                             waitCursor=True):
            self.setStatus(self.createSolverFiles())

    def checkConditions(self):
        capNames = self.conditionsTable.capNames()
        found = bcs.problems(self.conditionsTable.conditions, capNames, self._inletFaceId)
        if found:
            raise RuntimeError(" ".join(found))
        return capNames

    def solverFilesKey(self):
        node = self.ui.inputMeshSelector.currentNode()
        return (self.centerlinesKey(), bcs.to_json(self.conditionsTable.conditions), self.simulationParameters(),
                node.GetName())

    def createSolverFiles(self):
        """Write what svZeroDSolver reads into the output folder, without running it.

        The centerlines if they are not current, SimVascular's boundary condition files and
        `solver_0d.json` -- the case as the solver will see it, to look at or change before it is
        run. The centerlines are computed again only if they are not the current mesh's from the
        current inlet into the current folder, which is what makes trying another set of
        conditions cost a second rather than a centerline run.
        """
        node = self.ui.inputMeshSelector.currentNode()
        capNames = self.checkConditions()
        directory = self.caseDirectory()
        self.logic.importRomPackage(self.ui.romPackagePathLineEdit.currentPath)
        if self._centerlinesKey != self.centerlinesKey() or case.read_centerlines(directory) is None:
            self.computeCenterlines()
        byName = self.conditionsTable.conditionsByName()
        config = case.write_solver_input(directory, byName, self._inletFaceId,
                                         capNames[self._inletFaceId], self.simulationParameters(),
                                         sanitized(node.GetName()) or "model")
        self._solverFilesKey = self.solverFilesKey()
        with open(config) as handle:
            written = json.load(handle)
        return _("Wrote {config}: {vessels} vessels, {junctions} junctions and {conditions} "
                 "boundary conditions, with SimVascular's condition files beside it.").format(
            config=config, vessels=len(written.get("vessels", [])),
            junctions=len(written.get("junctions", [])),
            conditions=len(written.get("boundary_conditions", [])))

    def onRun(self):
        with slicer.util.tryWithErrorDisplay(_("The simulation could not be run."), waitCursor=True):
            self.runSimulation()

    def runSimulation(self):
        """Solve the solver files in the output folder, and show what they say.

        The files are written first only if they are not there, or the panel has changed since
        they were: otherwise they are run as they are, so that what Create solver files wrote --
        and anything changed in it by hand since -- is what is solved.
        """
        self.checkConditions()
        directory = self.caseDirectory()
        executable = solver.find_solver(self.ui.solverPathLineEdit.currentPath or None)
        if executable is None:
            raise RuntimeError(_("No svzerodsolver: set it under Tools, or put it on the PATH."))
        config = os.path.join(directory, case.SOLVER_INPUT_NAME)
        if self._solverFilesKey != self.solverFilesKey() or not os.path.isfile(config):
            self.createSolverFiles()

        started = time.time()
        solver.run_solver(config, os.path.join(directory, case.RESULTS_NAME), executable,
                          env=self.logic.externalEnvironment())
        self.showResults(directory)
        flowIn = sum(face.mean_flow for face in self._faceResults if face.is_inflow)
        flowOut = sum(face.mean_flow for face in self._faceResults if not face.is_inflow)
        self.setStatus(_("Solved in {seconds:.1f} s: {inflow:.4g} mL/s in, {outflow:.4g} mL/s out "
                         "over the last cycle. Written to {directory}.").format(
            seconds=time.time() - started, inflow=flowIn, outflow=flowOut, directory=directory))

    def loadResultsFromCase(self):
        """Show the results a case folder already holds, if they are newer than its input.

        So that reopening a scene shows the last run again without running it. Older results than
        input are a run that failed after writing it, and are left alone.
        """
        directory = self.ui.caseDirectoryPathLineEdit.currentPath
        config = os.path.join(directory or "", case.SOLVER_INPUT_NAME)
        solved = os.path.join(directory or "", case.RESULTS_NAME)
        if not (directory and os.path.isfile(config) and os.path.isfile(solved)
                and os.path.getmtime(solved) >= os.path.getmtime(config)
                and self._inletFaceId is not None):
            self.populateResults()
            return
        try:
            self.showResults(directory)
            self._resultsNote = _("From the last run in {directory}.").format(directory=directory)
        except Exception as error:
            logging.warning("SimVascular ROM could not read the results in %s: %s", directory, error)
            self._faceResults = []
        self.populateResults()

    def showResults(self, directory):
        node = self.ui.inputMeshSelector.currentNode()
        config = os.path.join(directory, case.SOLVER_INPUT_NAME)
        solved = solver.read_results(os.path.join(directory, case.RESULTS_NAME))
        self._faceResults = results.face_results(config, solved, self.capName(self._inletFaceId))
        self._resultsNote = ""
        centerlines = case.read_centerlines(directory)
        if centerlines is not None:
            coloured = results.centerline_results(centerlines, config, solved)
            self.logic.showResults(node, coloured, RESULTS_ARRAY_NAME)
        self.logic.writeResultsTable(node, self._faceResults)
        self.populateResults()
        self.plotSelected()

    # -- results -------------------------------------------------------------------------
    def populateResults(self):
        self.ui.exportButton.enabled = bool(self._faceResults)
        self.ui.visualizationButton.enabled = bool(self._faceResults)
        if self._resultsNote:
            self.setStatus(self._resultsNote)

    def plotSelected(self):
        """Plot the caps selected in the boundary conditions table, or the inlet if none are.

        The table the conditions are set in is the one results are asked of, because it is the
        one with a row per cap already, and selecting a row there already shows that cap in the
        3D view: one selection says which cap, in the anatomy and in the plot.
        """
        node = self.ui.inputMeshSelector.currentNode()
        if node is None or not self._faceResults:
            return
        solved = {face.name for face in self._faceResults}
        names = [self.capName(faceId) for faceId in self.conditionsTable.selectedFaceIds()
                 if self.capName(faceId) in solved]
        if not names and self.capName(self._inletFaceId) in solved:
            names = [self.capName(self._inletFaceId)]
        if names:
            self.logic.plotFaces(node, names, RESULTS_QUANTITY)

    def onExport(self):
        path = self.chooseExportFile()
        if not path:
            return
        with slicer.util.tryWithErrorDisplay(_("Could not export the results."), waitCursor=True):
            self.setStatus(self.exportResults(path))

    def chooseExportFile(self):
        """Where the results are to be written, or "" if nowhere was chosen.

        Its own method, as chooseWaveformFile is, so that a test can answer it.
        """
        node = self.ui.inputMeshSelector.currentNode()
        name = sanitized(node.GetName()) if node is not None else ""
        suggested = os.path.join(self.ui.caseDirectoryPathLineEdit.currentPath or "",
                                 f"{name or 'model'}_cap_results.csv")
        return qt.QFileDialog.getSaveFileName(
            slicer.util.mainWindow(), _("Export results"), suggested,
            _("CSV files (*.csv);;All files (*)"))

    def exportResults(self, path):
        written = results.write_face_results_csv(path, self._faceResults)
        return _("Wrote {count} caps' flow and pressure over the last cycle to {path}.").format(
            count=len(self._faceResults), path=written)

    def onOpenVisualization(self):
        with slicer.util.tryWithErrorDisplay(_("svZeroDVisualization could not be started.")):
            self.openVisualization()

    def openVisualization(self):
        directory = self.caseDirectory()
        script = visualization.find_script(self.ui.visualizationScriptPathLineEdit.currentPath)
        if script is None:
            raise RuntimeError(_("Set where svZeroDVisualization is under Tools: an svZeroDSolver "
                                 "checkout, or its visualize_simulation.py."))
        visualization.stop(self._visualization)
        output = os.path.join(directory, "svZeroDVisualization")
        self._visualization = visualization.launch(
            self.ui.visualizationPythonPathLineEdit.currentPath, script,
            os.path.join(directory, case.SOLVER_INPUT_NAME), output,
            env=self.logic.externalEnvironment())
        # Dash takes a few seconds to come up, and a browser opened before it has is an error page.
        qt.QTimer.singleShot(6000, lambda: qt.QDesktopServices.openUrl(qt.QUrl(visualization.URL)))
        self.setStatus(_("Starting svZeroDVisualization at {url}; if the page does not come up, "
                         "its log is in {output}.").format(url=visualization.URL, output=output))

    # -- panel state -----------------------------------------------------------------------
    def updateButtons(self):
        if self.conditionsTable is None:
            return
        caps = self.conditionsTable.caps
        found = (bcs.problems(self.conditionsTable.conditions, self.conditionsTable.capNames(),
                              self._inletFaceId) if caps else [])
        ready = bool(caps) and not found and bool(self.ui.caseDirectoryPathLineEdit.currentPath)
        self.ui.computeCenterlinesButton.enabled = bool(caps) and self._inletFaceId is not None
        self.ui.createSolverFilesButton.enabled = ready
        self.ui.runButton.enabled = ready
        notReady = " ".join(found) or _("Choose an output folder.")
        self.ui.createSolverFilesButton.toolTip = (
            _("Write the boundary condition files and svZeroDSolver's solver_0d.json into the "
              "output folder, computing the centerlines first if they are not current, without "
              "running anything. Run simulation solves them as they are, so they can be looked "
              "at or changed by hand in between.") if ready else notReady)
        self.ui.runButton.toolTip = (
            _("Run svZeroDSolver on the solver files and show the results. The files are created "
              "first if they are not there, or the panel has changed since they were.")
            if ready else notReady)
        if caps and not self._faceResults:
            self.setStatus(" ".join(found) if found else _("Ready to run."), warning=bool(found))

    def setStatus(self, text, warning=False):
        self.ui.statusLabel.text = text
        self.ui.statusLabel.styleSheet = "QLabel { color: #d08000; }" if warning else ""


class SimVascularROMLogic(ScriptedLoadableModuleLogic):
    """What touches the scene: the names on the mesh, and the nodes a run leaves beside it."""

    def __init__(self):
        ScriptedLoadableModuleLogic.__init__(self)

    @staticmethod
    def importRomPackage(checkout=""):
        """Make `sv_rom_simulation` importable: VMTK from SlicerVMTK, the package from a checkout."""
        ensureVmtk()
        try:
            import sv_rom_simulation  # noqa: F401
            return
        except ImportError:
            pass
        if checkout and os.path.isdir(os.path.join(checkout, "sv_rom_simulation")):
            if checkout not in sys.path:
                sys.path.insert(0, checkout)
            import sv_rom_simulation  # noqa: F401
            return
        raise RuntimeError(_("sv_rom_simulation is not available. Under Tools, set where a "
                             "checkout of svROMSimulation is."))

    @staticmethod
    def externalEnvironment():
        """The environment Slicer was started in, for running anything that is not Slicer.

        Slicer's launcher points PYTHONHOME, PYTHONPATH and the library path at its own Python and
        libraries. A process started with those inherited runs on Slicer's: another Python finds
        Slicer's packages and none of its own, which is how svZeroDVisualization failed to import
        pysvzerod from a Python that has it.
        """
        return slicer.util.startupEnvironment()

    @staticmethod
    def suggestedCaseDirectory(meshNode):
        """`rom` beside the saved scene, once the scene has been saved; nowhere until then."""
        url = slicer.mrmlScene.GetURL()
        directory = os.path.dirname(url) if url else ""
        if not os.path.isdir(directory):
            return ""
        return os.path.join(directory, "rom", sanitized(meshNode.GetName()) or "model")

    # -- nodes a run leaves --------------------------------------------------------------
    @staticmethod
    def referencedNode(meshNode, role, className, name):
        """The node the mesh references in `role`, made if it has none, named after the mesh."""
        node = meshNode.GetNodeReference(role)
        if node is None:
            node = slicer.mrmlScene.AddNewNodeByClass(className, f"{meshNode.GetName()} {name}")
            meshNode.SetNodeReferenceID(role, node.GetID())
        return node

    def showCenterlines(self, meshNode, centerlines):
        node = self.referencedNode(meshNode, CENTERLINES_REFERENCE, "vtkMRMLModelNode",
                                   "ROM centerlines")
        node.SetAndObserveMesh(centerlines)
        node.SetAndObserveTransformNodeID(meshNode.GetTransformNodeID())
        if node.GetDisplayNode() is None:
            node.CreateDefaultDisplayNodes()
        display = node.GetDisplayNode()
        display.SetColor(0.1, 0.5, 1.0)
        display.SetLineWidth(3)
        display.SetScalarVisibility(False)
        display.SetVisibility(True)
        results = meshNode.GetNodeReference(RESULTS_MODEL_REFERENCE)
        if results is not None and results.GetDisplayNode() is not None:
            results.GetDisplayNode().SetVisibility(False)
        return node

    def showResults(self, meshNode, coloured, arrayName):
        """The centerlines as tubes of the inscribed sphere radius, coloured by `arrayName`.

        Tubes rather than lines because a line is one pixel wide whatever it carries, and the
        colour is the point; the inscribed sphere radius is what VMTK traced them by, so the tube
        sits inside the vessel the way the vessel does. The mesh is faded so that it can be seen
        through to them.
        """
        tube = vtk.vtkTubeFilter()
        radius = coloured.GetPointData().GetArray("MaximumInscribedSphereRadius")
        withRadius = vtk.vtkPolyData()
        withRadius.ShallowCopy(coloured)
        if radius is not None:
            withRadius.GetPointData().SetActiveScalars("MaximumInscribedSphereRadius")
            tube.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
        else:
            tube.SetRadius(0.5)
        tube.SetInputData(withRadius)
        tube.SetNumberOfSides(16)
        tube.CappingOn()
        tube.Update()

        node = self.referencedNode(meshNode, RESULTS_MODEL_REFERENCE, "vtkMRMLModelNode", "ROM results")
        node.SetAndObserveMesh(tube.GetOutput())
        node.SetAndObserveTransformNodeID(meshNode.GetTransformNodeID())
        if node.GetDisplayNode() is None:
            node.CreateDefaultDisplayNodes()
        node.GetDisplayNode().SetVisibility(True)
        self.setResultsQuantity(meshNode, arrayName)

        centerlines = meshNode.GetNodeReference(CENTERLINES_REFERENCE)
        if centerlines is not None and centerlines.GetDisplayNode() is not None:
            centerlines.GetDisplayNode().SetVisibility(False)
        meshDisplay = meshNode.GetDisplayNode()
        if meshDisplay is not None and meshDisplay.GetOpacity() > RESULTS_MESH_OPACITY:
            meshDisplay.SetOpacity(RESULTS_MESH_OPACITY)
        return node

    @staticmethod
    def setResultsQuantity(meshNode, arrayName):
        node = meshNode.GetNodeReference(RESULTS_MODEL_REFERENCE)
        if node is None or node.GetDisplayNode() is None:
            return
        display = node.GetDisplayNode()
        display.SetActiveScalarName(arrayName)
        display.SetActiveAttributeLocation(vtk.vtkAssignAttribute.POINT_DATA)
        for colorNodeId in RESULTS_COLOR_NODES:
            if slicer.mrmlScene.GetNodeByID(colorNodeId) is not None:
                display.SetAndObserveColorNodeID(colorNodeId)
                break
        display.SetScalarRangeFlag(display.UseDataScalarRange)
        display.SetScalarVisibility(True)
        # A legend, where this Slicer has them: a colour without its scale says nothing.
        try:
            colors = slicer.modules.colors.logic()
            legend = colors.GetColorLegendDisplayNode(node) or colors.AddDefaultColorLegendDisplayNode(node)
            legend.SetTitleText(arrayName)
            # Two decimals: a model's pressures often span a fraction of a mmHg, and a legend
            # reading 12.1, 12.1, 12.0 says nothing about which end is which.
            legend.SetLabelFormat("%.2f")
            legend.SetVisibility(True)
        except Exception:
            logging.debug("SimVascular ROM could not show a colour legend.", exc_info=True)

    def writeResultsTable(self, meshNode, faceResults):
        """Every cap's pressure and flow over the last cycle, one column each, for plotting."""
        node = self.referencedNode(meshNode, RESULTS_TABLE_REFERENCE, "vtkMRMLTableNode",
                                   "ROM results table")
        table = vtk.vtkTable()
        if faceResults:
            columns = [(TIME_COLUMN, faceResults[0].time)]
            for face in faceResults:
                columns.append((f"{face.name} pressure [mmHg]", face.pressure / bcs.MMHG))
                columns.append((f"{face.name} flow [mL/s]", face.flow))
            for name, values in columns:
                array = numpy_to_vtk(np.asarray(values, dtype=float), deep=True)
                array.SetName(name)
                table.AddColumn(array)
        node.SetAndObserveTable(table)
        return node

    def plotFaces(self, meshNode, names, quantity):
        """Plot `quantity` over the last cycle for the named caps, and put the plot in the layout."""
        tableNode = meshNode.GetNodeReference(RESULTS_TABLE_REFERENCE)
        if tableNode is None:
            return None
        chart = self.referencedNode(meshNode, RESULTS_CHART_REFERENCE, "vtkMRMLPlotChartNode",
                                    "ROM plot")
        unit = "mmHg" if quantity == "Pressure" else "mL/s"
        chart.SetTitle(f"{quantity} over the last cycle")
        chart.SetXAxisTitle(TIME_COLUMN)
        chart.SetYAxisTitle(f"{quantity.lower()} [{unit}]")
        chart.RemoveAllPlotSeriesNodeIDs()
        for index, name in enumerate(names):
            column = f"{name} {quantity.lower()} [{unit}]"
            series = self.plotSeries(tableNode, column)
            if series is None:
                # Named after the cap alone, which is what the plot's legend shows; found again by
                # the table and column it plots, which is what makes it this cap's.
                series = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotSeriesNode", name)
                series.SetAttribute(PLOT_COLUMN_ATTRIBUTE, column)
                series.SetAndObserveTableNodeID(tableNode.GetID())
            series.SetXColumnName(TIME_COLUMN)
            series.SetYColumnName(column)
            series.SetPlotType(slicer.vtkMRMLPlotSeriesNode.PlotTypeScatter)
            series.SetMarkerStyle(slicer.vtkMRMLPlotSeriesNode.MarkerStyleNone)
            series.SetColor(*PLOT_COLORS[index % len(PLOT_COLORS)])
            chart.AddAndObservePlotSeriesNodeID(series.GetID())
        slicer.modules.plots.logic().ShowChartInLayout(chart)
        # A maximized view -- a scene can be saved with its 3D view maximized -- covers the plot
        # the line above switched the layout to, and then asking for a plot shows nothing at all.
        layoutNode = slicer.app.layoutManager().layoutLogic().GetLayoutNode()
        if layoutNode is not None and layoutNode.GetNumberOfMaximizedViewNodes():
            layoutNode.RemoveAllMaximizedViewNodes()
        return chart

    @staticmethod
    def plotSeries(tableNode, column):
        """The series already plotting `column` of `tableNode`, or None."""
        for series in slicer.util.getNodesByClass("vtkMRMLPlotSeriesNode"):
            if (series.GetAttribute(PLOT_COLUMN_ATTRIBUTE) == column
                    and series.GetTableNodeID() == tableNode.GetID()):
                return series
        return None


class SimVascularROMTest(ScriptedLoadableModuleTest):
    """Runs under Slicer; the package's own tests run without it."""

    def setUp(self):
        slicer.mrmlScene.Clear()
        slicer.mrmlScene.SetURL("")

    def runTest(self):
        self.setUp()
        self.test_conditionsAreKeptOnTheMeshTheyAreFor()
        self.setUp()
        self.test_doubleClickingAnInflowLoadsItsWaveform()
        self.setUp()
        self.test_aMeshWithoutNamesIsSentToMeshPrep()
        self.setUp()
        self.test_aWholeCaseOnAY()

    def namedY(self, name="Y"):
        from svromsetup import testing

        node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", name)
        node.SetAndObserveMesh(testing.y_surface())
        node.CreateDefaultDisplayNodes()
        node.SetAttribute(FACE_NAMES_ATTRIBUTE, json.dumps(
            {str(faceId): faceName for faceId, faceName in testing.NAMES.items()}))
        return node

    def widget(self):
        widget = slicer.util.getModuleWidget("SimVascularROM")
        # The widget outlives the scene, so what an earlier test left on it is reset here.
        widget.ui.caseDirectoryPathLineEdit.currentPath = ""
        return widget

    def typeCondition(self, widget, faceId, kind, text):
        row = widget.conditionsTable.rowOf(faceId)
        widget.ui.facesTable.cellWidget(row, TYPE_COLUMN).currentIndex = KIND_CHOICES.index(kind)
        widget.ui.facesTable.item(row, VALUES_COLUMN).setText(text)

    def test_conditionsAreKeptOnTheMeshTheyAreFor(self):
        """Two anatomies with the same face ids keep their own conditions, and get them back."""
        from svromsetup import testing

        widget = self.widget()
        first, second = self.namedY("first"), self.namedY("second")
        widget.ui.inputMeshSelector.setCurrentNode(first)
        self.assertEqual([faceId for faceId, _name, _geometry in widget.conditionsTable.caps],
                         [testing.INLET_ID, testing.RIGHT_ID, testing.LEFT_ID])
        self.assertEqual(widget._inletFaceId, testing.INLET_ID, "the largest cap is the inlet")
        self.typeCondition(widget, testing.RIGHT_ID, bcs.RESISTANCE_KIND, "1500")
        self.assertEqual(widget.conditionsTable.conditions[testing.RIGHT_ID], bcs.Resistance(1500.0))

        widget.ui.inputMeshSelector.setCurrentNode(second)
        self.assertEqual(widget.conditionsTable.conditions, {}, "a condition set on one anatomy is not on another")
        widget.ui.inputMeshSelector.setCurrentNode(first)
        self.assertEqual(widget.conditionsTable.conditions, {testing.RIGHT_ID: bcs.Resistance(1500.0)})
        self.assertEqual(bcs.from_json(first.GetAttribute(BOUNDARY_CONDITIONS_ATTRIBUTE)),
                         {testing.RIGHT_ID: bcs.Resistance(1500.0)})

        # What cannot be read stays on screen, in red, and is not a condition.
        self.typeCondition(widget, testing.LEFT_ID, bcs.RCR_KIND, "1, 2")
        self.assertNotIn(testing.LEFT_ID, widget.conditionsTable.conditions)
        item = widget.ui.facesTable.item(widget.conditionsTable.rowOf(testing.LEFT_ID), VALUES_COLUMN)
        self.assertEqual(item.text(), "1, 2")
        self.assertEqual(item.foreground().color().red(), INVALID_VALUES_COLOR.red())
        self.assertFalse(widget.ui.runButton.enabled)
        self.delayDisplay("Conditions kept per mesh")

    def test_doubleClickingAnInflowLoadsItsWaveform(self):
        """On an inflow the values are a file to load; on anything else, values to type.

        The file dialog is answered by replacing the method that opens it, which is the one thing
        here that would wait for a person.
        """
        import tempfile

        from svromsetup import testing

        widget = self.widget()
        widget.ui.inputMeshSelector.setCurrentNode(self.namedY())
        table = widget.ui.facesTable
        # Written as svMultiPhysics writes an inflow, negative, so that turning it round is
        # checked too.
        path = os.path.join(tempfile.mkdtemp(), "inlet.flow")
        with open(path, "w") as handle:
            handle.write("0.0 -5.0\n0.4 -15.0\n0.8 -5.0\n")
        asked = []
        widget.conditionsTable.chooseWaveformFile = lambda faceId: asked.append(faceId) or path
        try:
            self.typeCondition(widget, testing.INLET_ID, bcs.INFLOW, "")
            widget.conditionsTable.onCellDoubleClicked(widget.conditionsTable.rowOf(testing.INLET_ID), VALUES_COLUMN)
            self.assertEqual(asked, [testing.INLET_ID])
            self.assertEqual(widget.conditionsTable.conditions[testing.INLET_ID],
                             bcs.Inflow((0.0, 0.4, 0.8), (5.0, 15.0, 5.0)))
            self.assertIn("turned round", widget.ui.statusLabel.text)
            self.assertEqual(table.item(widget.conditionsTable.rowOf(testing.INLET_ID), VALUES_COLUMN).text(),
                             os.path.realpath(path), "the cell names the file the waveform came from")

            # Anywhere else a double-click opens the values for typing and asks for no file.
            self.typeCondition(widget, testing.RIGHT_ID, bcs.RCR_KIND, "100, 1e-4, 1000")
            widget.conditionsTable.onCellDoubleClicked(widget.conditionsTable.rowOf(testing.RIGHT_ID), VALUES_COLUMN)
            self.assertEqual(asked, [testing.INLET_ID])
            self.assertEqual(table.state(), qt.QAbstractItemView.EditingState)
            # Closed with Escape, as a person would close it. Not with reset(), which closes
            # editors by deleting every widget in the view's cells -- the condition combo boxes
            # with them.
            for editor in table.viewport().findChildren("QLineEdit"):
                qt.QApplication.sendEvent(editor, qt.QKeyEvent(qt.QEvent.KeyPress, qt.Qt.Key_Escape,
                                                               qt.Qt.NoModifier))
            self.assertNotEqual(table.state(), qt.QAbstractItemView.EditingState)
            self.assertIsNotNone(table.cellWidget(widget.conditionsTable.rowOf(testing.RIGHT_ID), TYPE_COLUMN))
            # And a cap with no kind yet is left alone.
            widget.conditionsTable.onCellDoubleClicked(widget.conditionsTable.rowOf(testing.LEFT_ID), VALUES_COLUMN)
            self.assertEqual(asked, [testing.INLET_ID])
        finally:
            del widget.conditionsTable.chooseWaveformFile
        self.delayDisplay("Double-clicking an inflow loads its waveform")

    def test_aMeshWithoutNamesIsSentToMeshPrep(self):
        from svmeshcomplete import testing

        widget = self.widget()
        node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", "unnamed")
        node.SetAndObserveMesh(testing.cube_mesh())
        widget.ui.inputMeshSelector.setCurrentNode(node)
        self.assertEqual(widget.conditionsTable.caps, [])
        self.assertFalse(widget.ui.openMeshPrepButton.isHidden())
        self.assertIn("Mesh Prep", widget.ui.meshStatusLabel.text)
        self.delayDisplay("An unnamed mesh is sent to Mesh Prep")

    def test_aWholeCaseOnAY(self):
        """Through the package and the solver, where both are on this machine."""
        import tempfile

        from svromsetup import testing

        widget = self.widget()
        # Where the Tools section does not say -- a --testing run has settings of its own -- the
        # variables the package's headless tests read.
        if not widget.ui.romPackagePathLineEdit.currentPath:
            widget.ui.romPackagePathLineEdit.currentPath = os.environ.get("SV_ROM_SIMULATION_PATH", "")
        if not widget.ui.solverPathLineEdit.currentPath:
            widget.ui.solverPathLineEdit.currentPath = os.environ.get("SVZERODSOLVER", "")
        try:
            widget.logic.importRomPackage(widget.ui.romPackagePathLineEdit.currentPath)
        except RuntimeError as error:
            self.delayDisplay(f"Skipped: {error}")
            return
        if solver.find_solver(widget.ui.solverPathLineEdit.currentPath or None) is None:
            self.delayDisplay("Skipped: no svzerodsolver")
            return
        node = self.namedY()
        widget.ui.inputMeshSelector.setCurrentNode(node)
        widget.ui.caseDirectoryPathLineEdit.currentPath = tempfile.mkdtemp()
        self.typeCondition(widget, testing.INLET_ID, bcs.INFLOW, "10")
        self.typeCondition(widget, testing.RIGHT_ID, bcs.RESISTANCE_KIND, "2000")
        self.typeCondition(widget, testing.LEFT_ID, bcs.INFLOW, "3")
        self.assertTrue(widget.ui.runButton.enabled)
        widget.ui.cyclesSpinBox.value = 3
        widget.ui.pointsSpinBox.value = 20
        widget.runSimulation()

        faces = {face.name: face for face in widget._faceResults}
        self.assertAlmostEqual(faces["cap_right"].mean_flow, 13.0, places=4)
        self.assertIsNotNone(node.GetNodeReference(RESULTS_MODEL_REFERENCE))
        table = node.GetNodeReference(RESULTS_TABLE_REFERENCE).GetTable()
        self.assertIsNotNone(table.GetColumnByName("cap_right pressure [mmHg]"))

        # Selecting a cap in the conditions table plots it, and the export holds every cap.
        widget.ui.facesTable.selectRow(widget.conditionsTable.rowOf(testing.RIGHT_ID))
        chart = node.GetNodeReference(RESULTS_CHART_REFERENCE)
        self.assertEqual(chart.GetNthPlotSeriesNode(0).GetYColumnName(), "cap_right pressure [mmHg]")
        self.assertTrue(widget.ui.exportButton.enabled)
        exported = os.path.join(widget.ui.caseDirectoryPathLineEdit.currentPath, "caps.csv")
        widget.exportResults(exported)
        with open(exported) as handle:
            header = handle.readline().strip().split(",")
        self.assertEqual(header[0], "time [s]")
        self.assertIn("cap_left flow [mL/s]", header)

        # Created solver files are run as they are: a resistance doubled by hand in
        # solver_0d.json between creating and running is the one solved. Changing a condition in
        # the panel afterwards makes Run write them again, from the panel.
        def pressureAtRight():
            return {face.name: face for face in widget._faceResults}["cap_right"].mean_pressure_mmhg
        before = pressureAtRight()
        self.assertTrue(widget.ui.createSolverFilesButton.enabled)
        widget.createSolverFiles()
        config = os.path.join(widget.ui.caseDirectoryPathLineEdit.currentPath, case.SOLVER_INPUT_NAME)
        with open(config) as handle:
            edited = json.load(handle)
        for condition in edited["boundary_conditions"]:
            if condition["bc_name"] == "RESISTANCE_cap_right":
                condition["bc_values"]["R"] *= 2
        with open(config, "w") as handle:
            json.dump(edited, handle)
        widget.runSimulation()
        self.assertAlmostEqual(pressureAtRight(), 2 * before, places=2)
        self.typeCondition(widget, testing.RIGHT_ID, bcs.RESISTANCE_KIND, "2000")
        self.typeCondition(widget, testing.LEFT_ID, bcs.INFLOW, "3.5")
        widget.runSimulation()
        self.assertAlmostEqual(pressureAtRight(), before * 13.5 / 13.0, places=2)
        self.delayDisplay("A whole case on a Y")
