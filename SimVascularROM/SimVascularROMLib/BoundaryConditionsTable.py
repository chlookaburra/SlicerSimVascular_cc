"""What the SimVascular simulation panels share: a mesh's named faces, and the conditions on them.

Two panels set boundary conditions on the caps of a mesh SimVascular Mesh Prep has named --
SimVascular ROM Simulation for the 0D model, SimVascular MultiPhysics for the 3D one -- and they
set the *same* ones: one set per mesh, kept on the mesh node, so that a 0D model tuned in one
panel is the 3D case written in the other, with nothing copied between them to go stale. The
table they are set in is therefore here, once, rather than in each panel: two copies of it would
be two tables that drift, and the first sign would be a condition that reads one way in one
panel and another in the other.

What differs between the panels stays in them. The 0D panel has a source -- the cap its
centerlines start at -- and shows it in bold here through `emphasised`; the 3D panel has none.
Each says what a change means to it through the `changed` and `selectionChanged` callbacks.
"""

import json
import logging
import os

import ctk
import qt
import slicer
import vtk
from slicer.i18n import tr as _

from svmeshcomplete import faces
from svmeshcomplete.face_table import Face
from svromsetup import boundary_conditions as bcs

# The face names are SimVascular Mesh Prep's, kept on the mesh as `{"<face id>": "<name>"}`.
# Read here, never written: a name typed in two panels is two names, and the one in Mesh Prep is
# the one the mesh-complete folder and every case of this mesh are written under.
FACE_NAMES_ATTRIBUTE = "SimVascularMeshPrep.FaceNames"

# Where the conditions are kept: on the mesh, for the reason Mesh Prep keeps the names there. A
# scene can hold a pre-op and a post-op anatomy, each with a face 4, and conditions kept per scene
# would move between them when a selector changed. Named for the panel that first kept them, and
# kept under that name so that scenes saved before the 3D panel existed still open with theirs.
BOUNDARY_CONDITIONS_ATTRIBUTE = "SimVascularROM.BoundaryConditions"

# Where the last waveform was loaded from, shared by both panels: a case's waveforms are kept
# together, whichever solver they are being given to.
WAVEFORM_DIRECTORY_SETTING = "SimVascularROM/WaveformDirectory"

COLUMNS = ("Name", "Type", "Values")
NAME_COLUMN = COLUMNS.index("Name")
TYPE_COLUMN = COLUMNS.index("Type")
VALUES_COLUMN = COLUMNS.index("Values")

# What the type column offers, the first meaning none: a cap with no condition is shown as one,
# rather than given a default that looks like a decision.
KIND_CHOICES = (None,) + bcs.KINDS
NO_CONDITION_TEXT = "-"

INVALID_VALUES_COLOR = qt.QColor(200, 60, 40)
INVALID_COLOR_NAME = "#c83c28"
WAVEFORM_COLOR = (31, 94, 160)
WARNING_COLOR = "#d08000"


def okCancelButtons(dialog):
    """OK and Cancel, wired to the dialog. Added one by one: given Ok | Cancel, PythonQt picks
    another constructor, and the box comes up empty -- no way to answer but closing the window."""
    buttons = qt.QDialogButtonBox()
    buttons.addButton(qt.QDialogButtonBox.Ok)
    buttons.addButton(qt.QDialogButtonBox.Cancel)
    buttons.connect("accepted()", dialog.accept)
    buttons.connect("rejected()", dialog.reject)
    return buttons


class WaveformChart:
    """A waveform plotted with every point marked, and a line at zero to read its sign against.

    The tables and plots are kept on the instance: the chart holds them only weakly from Python's
    side, and one collected while shown takes the chart down with it.
    """

    def __init__(self, quantity):
        self.view = ctk.ctkVTKChartView()
        self.view.minimumSize = qt.QSize(620, 300)
        self.curve, self.zero = vtk.vtkTable(), vtk.vtkTable()
        for table in (self.curve, self.zero):
            for name in ("time (s)", quantity):
                column = vtk.vtkDoubleArray()
                column.SetName(name)
                table.AddColumn(column)
        self.zeroPlot = vtk.vtkPlotLine()
        self.zeroPlot.SetInputData(self.zero, 0, 1)
        self.zeroPlot.SetColor(150, 150, 150, 255)
        self.zeroPlot.SetWidth(1.0)
        self.linePlot = vtk.vtkPlotLine()
        self.linePlot.SetInputData(self.curve, 0, 1)
        self.linePlot.SetColor(*WAVEFORM_COLOR, 255)
        self.linePlot.SetWidth(1.5)
        self.linePlot.SetMarkerStyle(vtk.vtkPlotPoints.CIRCLE)
        chart = self.view.chart()
        chart.AddPlot(self.zeroPlot)
        chart.AddPlot(self.linePlot)
        chart.GetAxis(vtk.vtkAxis.BOTTOM).SetTitle(_("time (s)"))
        chart.GetAxis(vtk.vtkAxis.LEFT).SetTitle(_(quantity))
        chart.SetShowLegend(False)

    def plot(self, times=(), values=(), period=None):
        """The waveform, or nothing when there is none to show."""
        self.linePlot.SetMarkerSize(4.0 if len(times) <= 400 else 2.0)
        zero = (0.0, period) if period else ()
        for table, xs, ys in ((self.curve, times, values), (self.zero, zero, (0.0,) * len(zero))):
            table.SetNumberOfRows(len(xs))
            for row, (x, y) in enumerate(zip(xs, ys)):
                table.SetValue(row, 0, vtk.vtkVariant(float(x)))
                table.SetValue(row, 1, vtk.vtkVariant(float(y)))
            table.Modified()
        self.view.chart().RecalculateBounds()
        self.view.scene().SetDirty(True)
        self.view.renderWindow().Render()


class WaveformReview:
    """A waveform file's flow, plotted point by point before it is given to a cap, with the choice
    of negating it.

    Asked, not guessed. A file written for svMultiPhysics stores an inflow negative, its flow
    measured along the cap's outward normal, and one written for a 0D solver stores it positive;
    nothing in the file says which. The table used to turn round any waveform whose mean was
    negative, which suited every class file and would be wrong, without a word, for a waveform
    that is mostly backflow -- and "turned round" in a status line showed nobody what had been done
    to their numbers. So the curve is shown, every point of it, and negating is the user's choice:
    unticked until they make it, the file's numbers as they are.
    """

    def __init__(self, waveform, capName, parent=None):
        self.waveform = waveform
        self.dialog = qt.QDialog(parent or slicer.util.mainWindow())
        self.dialog.setWindowTitle(_("Inflow waveform for {name}").format(name=capName))
        layout = qt.QVBoxLayout(self.dialog)
        intro = qt.QLabel(_(
            "{file}: every point in the file, as it would be given to {name}. Here an inflow is "
            "positive, into the model; a file written for svMultiPhysics stores it negative. "
            "Negate it if the curve is upside down.").format(
            file=os.path.basename(waveform.source) or _("The waveform"), name=capName))
        intro.wordWrap = True
        layout.addWidget(intro)
        self.chart = WaveformChart("flow (mL/s)")
        self.chartView, self._curve = self.chart.view, self.chart.curve
        layout.addWidget(self.chartView)
        self.summaryLabel = qt.QLabel()
        self.summaryLabel.wordWrap = True
        layout.addWidget(self.summaryLabel)
        self.negateCheckBox = qt.QCheckBox(_("Negate flow (multiply every value by -1)"))
        self.negateCheckBox.checked = False
        layout.addWidget(self.negateCheckBox)
        layout.addWidget(okCancelButtons(self.dialog))
        self.negateCheckBox.connect("toggled(bool)", lambda _checked: self.plot())
        self.plot()
        # Sized to what it holds: left to its default, the buttons fell below its bottom edge.
        self.dialog.adjustSize()

    def given(self):
        """The waveform as it would be given to the cap: the file's, negated if that is ticked."""
        return self.waveform.flipped() if self.negateCheckBox.checked else self.waveform

    def plot(self):
        waveform = self.given()
        self.chart.plot(waveform.time, waveform.flow, waveform.period)
        text = _("{points} points over {period:g} s: mean {mean:.4g} mL/s, from {low:.4g} to "
                 "{high:.4g} mL/s.").format(points=len(waveform.time), period=waveform.period,
                                            mean=waveform.mean, low=min(waveform.flow),
                                            high=max(waveform.flow))
        negative = waveform.mean < 0
        if negative:
            text += " " + _("The mean is negative: as it is, this cap would take flow out of the model.")
        self.summaryLabel.text = text
        self.summaryLabel.styleSheet = f"QLabel {{ color: {WARNING_COLOR}; }}" if negative else ""

    def ask(self):
        """True to negate the file's flow, False to take it as it is, None if cancelled."""
        if not self.dialog.exec_():
            return None
        return bool(self.negateCheckBox.checked)


class CoronaryDialog:
    """A coronary outlet set in a window of its own, as SimVascular sets one in its Set
    Inlet/Outlet BCs dialog: Ra, Ca, Ra-micro, Cim, Rv in its order; Pim from a file, with the
    period it is stretched onto (Pim Period) and the factor its values are multiplied by (Pim
    Scaling); and, beyond SimVascular's, the venous pressure Pv that svZeroDSolver takes.

    Each value in a field of its own beside its unit, rather than SimVascular's five on one line:
    five numbers a few orders of magnitude apart, in two units, are where a value goes in the
    wrong place. And the Pim plotted as the solver will get it -- stretched and scaled -- with the
    inflow's cardiac cycle beside its period, since svZeroDSolver refuses the two differing at
    all and says only "Inconsistent cardiac cycle period".
    """

    # (attribute, label, unit, what it is), in SimVascular's order, then Pv.
    FIELDS = (("Ra", "R<sub>a</sub>", "dyn·s/cm⁵", "arterial resistance"),
              ("Ca", "C<sub>a</sub>", "cm⁵/dyn", "arterial compliance"),
              ("Ram", "R<sub>a-micro</sub>", "dyn·s/cm⁵", "microvascular resistance"),
              ("Cim", "C<sub>im</sub>", "cm⁵/dyn", "intramyocardial compliance"),
              ("Rv", "R<sub>v</sub>", "dyn·s/cm⁵", "venous resistance"),
              ("Pv", "P<sub>v</sub>", "dyn/cm²", "venous pressure, 0 if left empty"))

    def __init__(self, capName, coronary=None, cycle=None, chooseFile=None, parent=None):
        self.capName = capName
        self.cycle = cycle
        self.chooseFile = chooseFile or (lambda: "")
        self.pim = coronary.pim if coronary is not None else None
        self.dialog = qt.QDialog(parent or slicer.util.mainWindow())
        self.dialog.setWindowTitle(_("Coronary outlet: {name}").format(name=capName))
        layout = qt.QVBoxLayout(self.dialog)
        intro = qt.QLabel(_(
            "The open-loop coronary of Kim et al. (2010), svZeroDSolver's CORONARY, with its values "
            "in SimVascular's order. Pim, the intramyocardial pressure, is read from a file of time "
            "and pressure, stretched onto its period and multiplied by its scaling."))
        intro.wordWrap = True
        layout.addWidget(intro)

        form = qt.QFormLayout()
        layout.addLayout(form)
        self.fields = {}
        for attribute, label, unit, what in self.FIELDS:
            field = qt.QLineEdit()
            field.toolTip = f"{what}, {unit}"
            value = getattr(coronary, attribute) if coronary is not None else None
            field.text = f"{value:.10g}" if value is not None else ""
            if attribute == "Pv":
                field.placeholderText = "0"
            row = qt.QWidget()
            rowLayout = qt.QHBoxLayout(row)
            rowLayout.setContentsMargins(0, 0, 0, 0)
            rowLayout.addWidget(field)
            rowLayout.addWidget(qt.QLabel(unit))
            text = qt.QLabel(f"{label}:")
            text.toolTip = field.toolTip
            form.addRow(text, row)
            self.fields[attribute] = field

        fileRow = qt.QWidget()
        fileLayout = qt.QHBoxLayout(fileRow)
        fileLayout.setContentsMargins(0, 0, 0, 0)
        self.fileLabel = qt.QLabel(os.path.basename(self.pim.source) if self.pim else _("none loaded"))
        self.browseButton = qt.QPushButton(_("Browse..."))
        fileLayout.addWidget(self.fileLabel, 1)
        fileLayout.addWidget(self.browseButton)
        form.addRow(qt.QLabel(_("P<sub>im</sub> (from file):")), fileRow)

        periodRow = qt.QWidget()
        periodLayout = qt.QHBoxLayout(periodRow)
        periodLayout.setContentsMargins(0, 0, 0, 0)
        self.periodField = qt.QLineEdit()
        period = (coronary.pim_period or coronary.pim.period) if self.pim is not None else None
        self.periodField.text = f"{period:.10g}" if period else ""
        self.periodField.toolTip = _("The period Pim is stretched onto, in s: the file's own when it "
                                     "is loaded. It has to be the inflow's cardiac cycle.")
        self.cycleButton = qt.QPushButton(_("Use the inflow's {cycle:g} s").format(cycle=cycle or 0))
        self.cycleButton.visible = cycle is not None
        periodLayout.addWidget(self.periodField)
        periodLayout.addWidget(qt.QLabel("s"))
        periodLayout.addWidget(self.cycleButton)
        form.addRow(qt.QLabel(_("P<sub>im</sub> period:")), periodRow)

        self.scaleField = qt.QLineEdit()
        self.scaleField.text = f"{coronary.pim_scale:.10g}" if coronary is not None else "1"
        self.scaleField.toolTip = _("What every value in the file is multiplied by: 1333.22 turns "
                                    "a file in mmHg into dyn/cm², which the solver takes.")
        form.addRow(qt.QLabel(_("P<sub>im</sub> scaling:")), self.scaleField)

        self.chart = WaveformChart("Pim (dyn/cm²)")
        layout.addWidget(self.chart.view)
        self.summaryLabel = qt.QLabel()
        self.summaryLabel.wordWrap = True
        layout.addWidget(self.summaryLabel)
        self.messageLabel = qt.QLabel()
        self.messageLabel.wordWrap = True
        self.messageLabel.styleSheet = f"QLabel {{ color: {WARNING_COLOR}; }}"
        layout.addWidget(self.messageLabel)
        self.buttons = okCancelButtons(self.dialog)
        layout.addWidget(self.buttons)

        for field in list(self.fields.values()) + [self.periodField, self.scaleField]:
            field.connect("textChanged(QString)", lambda _text: self.update())
        self.browseButton.connect("clicked(bool)", lambda _checked: self.browse())
        self.cycleButton.connect("clicked(bool)", lambda _checked: self.useCycle())
        self.update()
        self.dialog.adjustSize()

    def browse(self):
        path = self.chooseFile()
        if not path:
            return
        try:
            self.pim = bcs.read_pressure_file(path)
        except bcs.BoundaryConditionError as error:
            self.messageLabel.text = str(error)
            return
        self.fileLabel.text = os.path.basename(path)
        self.periodField.text = f"{self.pim.period:.10g}"    # the file's own, as SimVascular does
        self.update()

    def useCycle(self):
        if self.cycle is not None:
            self.periodField.text = f"{self.cycle:.10g}"

    def _number(self, field, positive=True, empty=None):
        text = field.text.strip()
        if not text and empty is not None:
            return empty
        try:
            value = float(text)
        except ValueError:
            return None
        if value != value or value in (float("inf"), float("-inf")) or (positive and not value > 0):
            return None
        return value

    def result(self):
        """The coronary as set, or None while something is missing or does not read."""
        values = {attribute: self._number(field, positive=attribute != "Pv",
                                          empty=0.0 if attribute == "Pv" else None)
                  for attribute, field in self.fields.items()}
        period, scale = self._number(self.periodField), self._number(self.scaleField, positive=False)
        if any(value is None for value in values.values()) or self.pim is None \
                or period is None or scale is None:
            return None
        return bcs.Coronary(**values).with_pim(
            self.pim, None if period == self.pim.period else period, scale)

    def update(self):
        problems = []
        for attribute, label, _unit, _what in self.FIELDS:
            field = self.fields[attribute]
            bad = self._number(field, positive=attribute != "Pv",
                               empty=0.0 if attribute == "Pv" else None) is None
            field.styleSheet = f"color: {INVALID_COLOR_NAME};" if bad and field.text.strip() else ""
            if bad:
                name = label.replace("<sub>", "").replace("</sub>", "")
                problems.append(_("{name} has to be a number more than 0.").format(name=name)
                                if attribute != "Pv" else _("Pv has to be a number."))
        coronary = self.result()
        if self.pim is None:
            problems.append(_("Load P<sub>im</sub> from a file of time and pressure."))
        elif self._number(self.periodField) is None:
            problems.append(_("The P<sub>im</sub> period has to be a time more than 0."))
        elif self._number(self.scaleField, positive=False) is None:
            problems.append(_("The P<sub>im</sub> scaling has to be a number."))
        pim = None
        if self.pim is not None and self._number(self.periodField) and \
                self._number(self.scaleField, positive=False) is not None:
            period, scale = self._number(self.periodField), self._number(self.scaleField, positive=False)
            pim = bcs.Coronary(1, 1, 1, 1, 1).with_pim(self.pim, period, scale).applied_pim
        if pim is not None:
            self.chart.plot(pim.time, pim.pressure, pim.period)
            self.summaryLabel.text = _(
                "As the solver gets it: {points} points over {period:g} s, mean {mean:.4g} dyn/cm² "
                "({mmHg:.3g} mmHg), from {low:.4g} to {high:.4g} dyn/cm².").format(
                points=len(pim.time), period=pim.period, mean=pim.mean, mmHg=pim.mean / bcs.MMHG,
                low=min(pim.pressure), high=max(pim.pressure))
            if self.cycle is not None and not bcs.same_period(pim.period, self.cycle):
                problems.append(_("The inflow's cardiac cycle is {cycle:g} s, and P<sub>im</sub> covers "
                                  "{period:g} s: the solver refuses the two differing, so set its "
                                  "period to the inflow's.").format(cycle=self.cycle, period=pim.period))
        else:
            self.chart.plot()
            self.summaryLabel.text = ""
        self.cycleButton.enabled = self.cycle is not None and (
            pim is None or not bcs.same_period(pim.period, self.cycle))
        self.messageLabel.text = " ".join(problems)
        # A Pim off the inflow's cycle is kept, and named by the panel as what stops the run, so
        # that a coronary set before its inflow is not lost; anything else unread is not a condition.
        self.buttons.button(qt.QDialogButtonBox.Ok).enabled = coronary is not None

    def ask(self):
        """The coronary as set, or None if cancelled."""
        if not self.dialog.exec_():
            return None
        return self.result()


def faceNames(meshNode):
    """`{face id: name}` as SimVascular Mesh Prep saved them on the mesh, or {}."""
    text = meshNode.GetAttribute(FACE_NAMES_ATTRIBUTE) if meshNode is not None else None
    if not text:
        return {}
    try:
        stored = json.loads(text)
    except ValueError:
        logging.warning("Could not read the face names on %s.", meshNode.GetName())
        return {}
    names = {}
    for faceId, name in stored.items():
        try:
            names[int(faceId)] = str(name)
        except (TypeError, ValueError):
            continue
    return names


class NamedMeshFaces:
    """A mesh's faces as Mesh Prep named them, measured once for as long as the mesh is the same.

    Measured once because coming back to a panel reads the names again, which is cheap, and would
    measure the faces again with them, which on a clinical volume mesh of a million cells is not.
    The boundary is kept too: it is what the highlight is cut from.
    """

    def __init__(self):
        self._key = None
        self._measured = []
        self.boundary = None
        self.arrayName = None
        self.names = {}
        self.caps = []
        self.unnamed = []
        self.walls = 0

    def read(self, node):
        """Measure the node's faces and read their names; raises ValueError if it has no faces.

        :return: whether the faces are named at all. Unnamed faces are left out of `caps` and
          listed in `unnamed`, so a panel can say which ones Mesh Prep has still to name.
        """
        key = (node.GetID(), node.GetMesh().GetMTime())
        if self._key != key:
            self.arrayName = faces.find_face_id_array(node.GetMesh())
            self._measured = faces.measure_faces(node.GetMesh(), self.arrayName)
            self.boundary = faces.boundary_of(node.GetMesh())
            self._key = key
        self.names = faceNames(node)
        self.unnamed = [face.face_id for face in self._measured if not self.names.get(face.face_id)]
        self.caps = [(face.face_id, self.names[face.face_id], face) for face in self._measured
                     if self.names.get(face.face_id)
                     and Face(face.face_id, self.names[face.face_id]).is_cap]
        self.walls = len(self._measured) - len(self.caps) - len(self.unnamed)
        return bool(self.names)

    @property
    def count(self):
        return len(self._measured)

    def describe(self):
        """What the mesh panel says about these faces."""
        if not self.names:
            return _("{count} faces, none of them named. The names are the boundary conditions' "
                     "and the results', and they come from SimVascular Mesh Prep.").format(
                count=self.count)
        text = _("{caps} caps and {walls} wall face(s), named in SimVascular Mesh Prep.").format(
            caps=len(self.caps), walls=self.walls)
        if self.unnamed:
            text += " " + _("Face(s) {faces} have no name and are left out; name them in Mesh "
                            "Prep.").format(faces=", ".join(str(faceId) for faceId in self.unnamed))
        return text


def highlight(mesh, faceId, arrayName, boundary=None):
    """Show one face, with Mesh Prep's own highlight, so that every panel shows a face the same way."""
    from SimVascularMeshPrep import SimVascularMeshPrepLogic
    SimVascularMeshPrepLogic().highlight(mesh, faceId, arrayName or "", boundary=boundary)


def clearHighlight():
    try:
        from SimVascularMeshPrep import SimVascularMeshPrepLogic
    except ImportError:
        return
    SimVascularMeshPrepLogic.clearHighlight()


class BoundaryConditionsTable:
    """One row per cap: its name, the type of its condition, and the condition's values.

    The values of an inflow are a waveform, which nobody types: double-clicking them asks for the
    file it is in, and a steady flow is typed by selecting the cell and typing. A coronary is set
    in a window of its own, as SimVascular sets one, opened by choosing it or by double-clicking
    its values; its cell shows what was set, and is not typed in. Anything that
    cannot be read stays on screen in red, with the reason as its tooltip, rather than vanishing
    on a typo; it is not a condition until it reads.
    """

    def __init__(self, table, hintLabel, setStatus, changed=None, selectionChanged=None,
                 startDirectory=lambda: ""):
        self.table = table
        self.hintLabel = hintLabel
        self.setStatus = setStatus
        self.changed = changed or (lambda: None)
        self.selectionChanged = selectionChanged or (lambda faceIds: None)
        self.startDirectory = startDirectory
        self.mesh = None
        self.faces = None
        self.caps = []
        self.conditions = {}
        self.invalid = {}
        # The cap shown in bold, and what its tooltip says about why: the 0D panel's source.
        self.emphasised = None
        self.emphasisNote = ""
        self.updating = False

        table.setColumnCount(len(COLUMNS))
        table.setHorizontalHeaderLabels([_(name) for name in COLUMNS])
        table.verticalHeader().setVisible(False)
        # A waveform shows as its file's path, and a path cut short at its end loses the one part
        # of it that says which file it is. Without word wrap, which the table has on by default
        # and which, a path being full of places to break, elides the end whatever the mode says.
        table.setWordWrap(False)
        table.setTextElideMode(qt.Qt.ElideMiddle)
        table.horizontalHeader().setSectionResizeMode(qt.QHeaderView.ResizeToContents)
        table.horizontalHeader().setSectionResizeMode(VALUES_COLUMN, qt.QHeaderView.Stretch)
        # A double-click is this table's to answer rather than Qt's, because it means two things:
        # on an inflow or a coronary it asks for a file, anywhere else it opens the cell for typing.
        # Typing over a selected cell, which is how a steady inflow is entered, still opens it.
        table.setEditTriggers(qt.QAbstractItemView.EditKeyPressed
                              | qt.QAbstractItemView.AnyKeyPressed)
        table.connect("cellChanged(int,int)", self.onValuesEdited)
        table.connect("itemSelectionChanged()", self.onSelectionChanged)
        table.connect("cellDoubleClicked(int,int)", self.onCellDoubleClicked)

    # -- which mesh ---------------------------------------------------------------------
    def setMesh(self, mesh, namedFaces):
        """Show `mesh`'s caps, with the conditions saved on it, or nothing for None."""
        self.mesh = mesh
        self.faces = namedFaces
        self.caps = list(namedFaces.caps) if (mesh is not None and namedFaces is not None) else []
        self.invalid = {}
        capIds = {faceId for faceId, _name, _geometry in self.caps}
        self.conditions = ({faceId: condition for faceId, condition in bcs.from_json(
            mesh.GetAttribute(BOUNDARY_CONDITIONS_ATTRIBUTE)).items() if faceId in capIds}
            if self.caps else {})
        clearHighlight()
        self.populate()

    def save(self):
        """Write the conditions onto the mesh, as they are edited."""
        if self.mesh is None or not self.caps:
            return
        capIds = {faceId for faceId, _name, _geometry in self.caps}
        self.mesh.SetAttribute(BOUNDARY_CONDITIONS_ATTRIBUTE, bcs.to_json(
            {faceId: condition for faceId, condition in self.conditions.items()
             if faceId in capIds}))

    # -- what is in it ----------------------------------------------------------------------
    def capName(self, faceId):
        return next((name for capId, name, _geometry in self.caps if capId == faceId), str(faceId))

    def capNames(self):
        return {faceId: name for faceId, name, _geometry in self.caps}

    def conditionsByName(self):
        names = self.capNames()
        return {names[faceId]: condition for faceId, condition in self.conditions.items()
                if faceId in names}

    def period(self):
        return bcs.cycle_period(self.conditions.values())

    def kindOf(self, faceId):
        """The kind shown for a face: its condition's, or the one chosen for what failed to read."""
        if faceId in self.invalid:
            return self.invalid[faceId][0]
        condition = self.conditions.get(faceId)
        return condition.kind if condition is not None else None

    def rowOf(self, faceId):
        return next((row for row, (capId, _name, _geometry) in enumerate(self.caps)
                     if capId == faceId), None)

    def selectedFaceIds(self):
        rows = sorted({index.row() for index in self.table.selectedIndexes()})
        return [self.caps[row][0] for row in rows if row < len(self.caps)]

    # -- drawing it ------------------------------------------------------------------------
    def populate(self):
        self.updating = True
        try:
            table = self.table
            table.setRowCount(len(self.caps))
            for row, (faceId, name, _geometry) in enumerate(self.caps):
                nameItem = table.item(row, NAME_COLUMN)
                if nameItem is None:
                    nameItem = qt.QTableWidgetItem()
                    table.setItem(row, NAME_COLUMN, nameItem)
                nameItem.setText(name)
                nameItem.setFlags(qt.Qt.ItemIsEnabled | qt.Qt.ItemIsSelectable)
                font = nameItem.font()
                font.setBold(faceId == self.emphasised)
                nameItem.setFont(font)
                nameItem.setToolTip(self.emphasisNote if faceId == self.emphasised else "")

                combo = table.cellWidget(row, TYPE_COLUMN)
                if combo is None:
                    combo = qt.QComboBox()
                    for kind in KIND_CHOICES:
                        combo.addItem(_(kind) if kind else NO_CONDITION_TEXT)
                    combo.connect("currentIndexChanged(int)",
                                  lambda index, combo=combo: self.onKindChanged(combo, index))
                    table.setCellWidget(row, TYPE_COLUMN, combo)
                combo.setProperty("faceId", faceId)
                combo.currentIndex = KIND_CHOICES.index(self.kindOf(faceId))
                self.populateValues(row, faceId)
        finally:
            self.updating = False
        self.updateHint()

    def populateValues(self, row, faceId):
        item = self.table.item(row, VALUES_COLUMN)
        if item is None:
            item = qt.QTableWidgetItem()
            self.table.setItem(row, VALUES_COLUMN, item)
        if faceId in self.invalid:
            _kind, text, reason = self.invalid[faceId]
            item.setText(text)
            item.setForeground(qt.QBrush(INVALID_VALUES_COLOR))
            item.setToolTip(reason)
        else:
            text = bcs.to_text(self.conditions.get(faceId))
            item.setText(text)
            item.setForeground(qt.QBrush())
            # The text as well as the hint: a waveform's path is often longer than the column.
            hint = self.hintFor(self.kindOf(faceId))
            item.setToolTip(f"{text}\n{hint}" if text else hint)
        # A coronary is set in its window, not typed: the cell is what was set.
        editable = self.kindOf(faceId) not in (None, bcs.CORONARY_KIND)
        item.setFlags(qt.Qt.ItemIsEnabled | qt.Qt.ItemIsSelectable
                      | (qt.Qt.ItemIsEditable if editable else 0))

    def redrawValues(self, faceId):
        row = self.rowOf(faceId)
        if row is None:
            return
        self.updating = True
        try:
            self.populateValues(row, faceId)
        finally:
            self.updating = False

    @staticmethod
    def hintFor(kind):
        if kind == bcs.INFLOW:
            return _("Inflow: double-click the values to load a waveform from a .flow file, or "
                     "select them and type a steady flow Q in mL/s.")
        if kind == bcs.RCR_KIND:
            return _("RCR: Rp, C, Rd, Pd -- proximal resistance, compliance, distal resistance, "
                     "distal pressure (may be left out, and is then 0).")
        if kind == bcs.RESISTANCE_KIND:
            return _("Resistance: R -- the resistance, draining to zero pressure. An outlet "
                     "that needs a pressure to drain to is an RCR.")
        if kind == bcs.CORONARY_KIND:
            return _("Coronary: double-click to set it in its own window -- Ra, Ca, Ra-micro, Cim, "
                     "Rv in SimVascular's order, the venous pressure Pv, and the intramyocardial "
                     "pressure Pim from a file, with its period and scaling.")
        return _("Choose a kind of condition for the cap first.")

    def updateHint(self):
        selected = self.selectedFaceIds()
        self.hintLabel.text = self.hintFor(self.kindOf(selected[0])) if selected else ""

    # -- editing it ------------------------------------------------------------------------
    def onKindChanged(self, combo, index):
        """A cap given another kind of condition: what was typed is read again as the new kind.

        So that a value that fits the new kind carries over, and one that does not is shown in
        red rather than dropped.
        """
        if self.updating:
            return
        faceId = int(combo.property("faceId"))
        kind = KIND_CHOICES[index]
        row = self.rowOf(faceId)
        previous = self.conditions.get(faceId)
        self.invalid.pop(faceId, None)
        if kind is None:
            self.conditions.pop(faceId, None)
        elif kind == bcs.CORONARY_KIND and not isinstance(previous, bcs.Coronary):
            # Nothing typed carries over into a coronary: its window opens to set it, once the
            # combo box has finished changing -- not from inside its signal.
            self.conditions.pop(faceId, None)
            self.invalid[faceId] = (kind, "", self.hintFor(kind))
            qt.QTimer.singleShot(0, lambda: self.setCoronary(faceId))
        elif previous is None or previous.kind != kind:
            text = self.table.item(row, VALUES_COLUMN).text() if row is not None else ""
            if isinstance(previous, bcs.Inflow) and not previous.is_steady:
                text = ""
            # What the cell says about a coronary's Pim is no value of any other kind.
            text = text.split(bcs.PIM_SEPARATOR)[0].strip()
            self.setFromText(faceId, kind, text)
        self.save()
        self.redrawValues(faceId)
        self.updateHint()
        self.changed()

    def setFromText(self, faceId, kind, text):
        try:
            condition = bcs.from_text(kind, text, self.conditions.get(faceId), self.period())
        except bcs.BoundaryConditionError as error:
            self.conditions.pop(faceId, None)
            self.invalid[faceId] = (kind, text, str(error))
            return
        self.invalid.pop(faceId, None)
        if condition is None:
            self.conditions.pop(faceId, None)
            if kind is not None:
                # Chosen but not yet given values: keep the choice visible.
                self.invalid[faceId] = (kind, "", self.hintFor(kind))
        else:
            self.conditions[faceId] = condition

    def onValuesEdited(self, row, column):
        if self.updating or column != VALUES_COLUMN or row >= len(self.caps):
            return
        faceId = self.caps[row][0]
        kind = KIND_CHOICES[self.table.cellWidget(row, TYPE_COLUMN).currentIndex]
        if kind is None:
            return
        self.setFromText(faceId, kind, self.table.item(row, column).text())
        self.save()
        self.redrawValues(faceId)
        self.changed()

    def onSelectionChanged(self):
        selected = self.selectedFaceIds()
        if selected and self.mesh is not None and self.mesh.GetMesh() is not None and self.faces:
            highlight(self.mesh.GetMesh(), selected[0], self.faces.arrayName, self.faces.boundary)
        else:
            clearHighlight()
        self.updateHint()
        self.selectionChanged(selected)

    # -- a waveform, by double-clicking an inflow's values ---------------------------------
    def onCellDoubleClicked(self, row, column):
        """Ask for a waveform file on an inflow's values; open any other values for typing."""
        if column != VALUES_COLUMN or row >= len(self.caps):
            return
        faceId = self.caps[row][0]
        kind = self.kindOf(faceId)
        if kind == bcs.INFLOW:
            path = self.chooseWaveformFile(faceId)
            if not path:
                return
            waveform = None
            with slicer.util.tryWithErrorDisplay(_("Could not read the waveform."), waitCursor=True):
                waveform = bcs.read_flow_file(path)
            if waveform is None:
                return
            negate = self.reviewWaveform(waveform, faceId)
            if negate is None:
                return
            with slicer.util.tryWithErrorDisplay(_("Could not read the waveform."), waitCursor=True):
                self.setStatus(self.loadWaveform(path, faceId, negate=negate))
        elif kind == bcs.CORONARY_KIND:
            self.setCoronary(faceId)
        elif kind is not None:
            self.table.editItem(self.table.item(row, column))

    def chooseWaveformFile(self, faceId):
        """The file a cap's waveform is to be read from, or "" if none was chosen.

        Opened where the last waveform came from, since a case's waveforms are usually kept
        together, and a Fontan has seven of them to load one after another.
        """
        directory = qt.QSettings().value(WAVEFORM_DIRECTORY_SETTING, "") or self.startDirectory()
        title = (_("Intramyocardial pressure for {name}") if self.kindOf(faceId) == bcs.CORONARY_KIND
                 else _("Inflow waveform for {name}"))
        path = qt.QFileDialog.getOpenFileName(
            slicer.util.mainWindow(), title.format(name=self.capName(faceId)),
            directory, _("Waveform files (*.flow *.dat *.txt *.csv);;All files (*)"))
        if path:
            qt.QSettings().setValue(WAVEFORM_DIRECTORY_SETTING, os.path.dirname(path))
        return path

    def setCoronary(self, faceId):
        """Open a cap's coronary window, and give the cap what was set there; leave it as it was if
        the window was cancelled."""
        if self.kindOf(faceId) != bcs.CORONARY_KIND:
            return      # the kind was changed again before the window opened
        current = self.conditions.get(faceId)
        coronary = self.askCoronary(faceId, current if isinstance(current, bcs.Coronary) else None)
        if coronary is None:
            return
        self.conditions[faceId] = coronary
        self.invalid.pop(faceId, None)
        self.save()
        self.populate()
        self.changed()
        pim = coronary.applied_pim
        self.setStatus(_("{name}: coronary set, Pim of {points} points over {period:g} s, mean "
                         "{mean:.4g} dyn/cm² ({mmHg:.3g} mmHg).").format(
            name=self.capName(faceId), points=len(pim.time), period=pim.period, mean=pim.mean,
            mmHg=pim.mean / bcs.MMHG))

    def inflowCycle(self, faceId=None):
        """The cardiac cycle an inflow waveform sets, other than `faceId`'s, or None if none does."""
        return next((condition.period for capId, condition in self.conditions.items()
                     if capId != faceId and isinstance(condition, bcs.Inflow) and not condition.is_steady),
                    None)

    def askCoronary(self, faceId, current):
        """The coronary set in its window for a cap, or None if cancelled. Its own method, as
        chooseWaveformFile is, so that a test can answer it."""
        return CoronaryDialog(self.capName(faceId), current, cycle=self.inflowCycle(faceId),
                              chooseFile=lambda: self.chooseWaveformFile(faceId)).ask()

    def reviewWaveform(self, waveform, faceId):
        """Show a file's waveform, and ask whether to negate it: True, False, or None to give the
        cap nothing. Its own method, as chooseWaveformFile is, so that a test can answer it."""
        return WaveformReview(waveform, self.capName(faceId)).ask()

    def loadWaveform(self, path, faceId, negate=False):
        """Give a cap the waveform in a file, negated if asked; saying so, and saying when its mean
        is negative, which takes flow out of the model at that cap."""
        waveform = bcs.read_flow_file(path)
        if negate:
            waveform = waveform.flipped()
        self.conditions[faceId] = waveform
        self.invalid.pop(faceId, None)
        self.save()
        self.populate()
        self.changed()
        note = _("{name}: {points} points over {period:g} s, mean {mean:.4g} mL/s.").format(
            name=self.capName(faceId), points=len(waveform.time), period=waveform.period,
            mean=waveform.mean)
        if negate:
            note += " " + _("Negated: every flow in the file multiplied by -1.")
        if waveform.mean < 0:
            note += " " + _("Its mean is negative, so this cap takes flow out of the model: "
                            "double-click its values to load it again, negated.")
        return note
