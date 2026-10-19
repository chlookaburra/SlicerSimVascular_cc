"""What happens at each cap: a flow coming in, or a resistance (with or without a compliance) out.

A reduced-order model is a network of vessels, and the network is closed by one condition per
cap. Three kinds cover what a vascular model takes, and each is what svZeroDSolver calls it:

- **Inflow**, a prescribed flow entering the model -- a waveform, or a steady flow, which is a
  waveform that does not change. Positive is *into* the model. A Fontan has several: the venae
  cavae, the hepatic veins and the azygous all enter, and the pulmonary arteries all leave.
- **RCR**, a proximal resistance, a compliance and a distal resistance, draining to a distal
  pressure: the Windkessel that stands in for the vascular bed beyond an outlet.
- **Resistance**, a resistance alone, draining to zero pressure: what an outlet into a
  low-impedance bed -- the lungs, in a Fontan -- is usually given. It has no distal pressure,
  because svMultiPhysics' Resistance condition has none (it is `P = R Q`), and the same
  conditions are written for the 0D and the 3D solver: an outlet that needs a pressure to drain
  to is an RCR.
- **Coronary**, the open-loop coronary of Kim et al. (2010), svZeroDSolver's `CORONARY`: an
  arterial resistance Ra and compliance Ca, a microvascular resistance Ra-micro, an
  intramyocardial compliance Cim squeezed by the intramyocardial pressure Pim(t), and a venous
  resistance Rv draining to a venous pressure Pv -- in SimVascular's order, Ra, Ca, Ra-micro, Cim,
  Rv, so that values from a SimVascular project or a paper that follows it go in as they read. An
  RCR cannot stand in for it: coronary flow is highest in diastole, because the contracting heart
  squeezes its own vessels in systole, and that is what Pim does. Pim is a file's waveform,
  stretched onto a period and scaled, as SimVascular's Pim Period and Pim Scaling do.
  svMultiPhysics has no such condition of its own -- it takes one only coupled to svZeroDSolver
  -- so for now it is the 0D panel's alone, and the 3D panel says so.

Units are cgs throughout, as both solvers take them: flow in mL/s, resistance in dyn·s/cm⁵,
compliance in cm⁵/dyn, pressure in dyn/cm². Nothing here converts; the panel says so beside
the table, which is where a value in mmHg would be typed.

Every condition is kept per face *name*, not per face id, by the files the solvers read: those
are SimVascular's `rcrt.dat`, `resistance.dat` and `cort.dat`, all keyed by face name, and the
reader here takes the first two from a SimVascular project as they are.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Mapping

# The kinds, as the panel lists them and as they are saved.
INFLOW = "Inflow"
RCR_KIND = "RCR"
RESISTANCE_KIND = "Resistance"
CORONARY_KIND = "Coronary"
KINDS = (INFLOW, RCR_KIND, RESISTANCE_KIND, CORONARY_KIND)

# What one mmHg is in dyn/cm², for the panel to report pressures in the unit they are read in.
MMHG = 1333.22

_SEPARATORS = re.compile(r"[,;\s]+")


class BoundaryConditionError(ValueError):
    """Raised when a condition cannot be read, or the set of them cannot close the model."""


@dataclass(frozen=True)
class Inflow:
    """A flow entering the model at a cap, over one cardiac cycle.

    The cycle is the waveform's own last time: both solvers repeat a waveform over it, so a
    steady flow is two points, at 0 and at the period, carrying the same value.
    """

    time: tuple
    flow: tuple
    source: str = field(default="", compare=False)
    """The file the waveform was read from, if it was one: what the table shows for it, so that
    the condition on a cap says which file it is. Not compared -- the same numbers are the same
    inflow wherever they came from -- and not needed to run it, since the numbers are kept with
    the scene: a file moved or deleted since leaves the condition as it was."""
    kind = INFLOW

    def __post_init__(self):
        if len(self.time) != len(self.flow):
            raise BoundaryConditionError("An inflow needs as many flows as times.")
        if len(self.time) < 2:
            raise BoundaryConditionError(
                "An inflow needs two points at least: the waveform is repeated over its last "
                "time, so one point has no period to repeat over.")
        if any(later <= earlier for earlier, later in zip(self.time, self.time[1:])):
            raise BoundaryConditionError("An inflow's times have to increase.")
        if self.time[0] != 0.0:
            raise BoundaryConditionError(
                "An inflow has to start at time 0: the cycle is taken as running from 0 to "
                "its last time.")

    @classmethod
    def steady(cls, flow: float, period: float = 1.0) -> "Inflow":
        return cls((0.0, float(period)), (float(flow), float(flow)))

    @property
    def period(self) -> float:
        return float(self.time[-1])

    @property
    def is_steady(self) -> bool:
        return len(set(self.flow)) == 1

    @property
    def mean(self) -> float:
        """The cycle-averaged flow, by the trapezoidal rule over the waveform's own points."""
        total = sum((t1 - t0) * (q0 + q1) / 2.0 for t0, t1, q0, q1 in
                    zip(self.time, self.time[1:], self.flow, self.flow[1:]))
        return total / self.period

    def flipped(self) -> "Inflow":
        return Inflow(self.time, tuple(-q for q in self.flow), self.source)


@dataclass(frozen=True)
class RCR:
    Rp: float
    C: float
    Rd: float
    Pd: float = 0.0
    kind = RCR_KIND


@dataclass(frozen=True)
class Resistance:
    R: float
    kind = RESISTANCE_KIND


def check_series(time, values, what: str):
    """Refuse a waveform the solvers would read wrongly: as many values as times, two points at
    least, times increasing from 0 -- the cycle runs from 0 to the last time, and both solvers
    repeat a waveform over it."""
    if len(time) != len(values):
        raise BoundaryConditionError(f"{what} needs as many values as times.")
    if len(time) < 2:
        raise BoundaryConditionError(
            f"{what} needs two points at least: it is repeated over its last time, so one point "
            "has no period to repeat over.")
    if any(later <= earlier for earlier, later in zip(time, time[1:])):
        raise BoundaryConditionError(f"{what}'s times have to increase.")
    if time[0] != 0.0:
        raise BoundaryConditionError(
            f"{what} has to start at time 0: the cycle is taken as running from 0 to its last time.")


@dataclass(frozen=True)
class Pressure:
    """A pressure over one cardiac cycle, read from a file: a coronary's intramyocardial pressure.

    In dyn/cm², as the file has it -- nothing converts it, as nothing converts any other value
    here.
    """

    time: tuple
    pressure: tuple
    source: str = field(default="", compare=False)

    def __post_init__(self):
        check_series(self.time, self.pressure, "A pressure waveform")

    @property
    def period(self) -> float:
        return float(self.time[-1])

    @property
    def mean(self) -> float:
        total = sum((t1 - t0) * (p0 + p1) / 2.0 for t0, t1, p0, p1 in
                    zip(self.time, self.time[1:], self.pressure, self.pressure[1:]))
        return total / self.period


@dataclass(frozen=True)
class Coronary:
    """The open-loop coronary outlet, in SimVascular's order: Ra, Ca, Ram (Ra-micro), Cim, Rv --
    resistances in dyn·s/cm⁵, compliances in cm⁵/dyn -- and the venous pressure Pv in dyn/cm².

    Its Pim is kept as the file has it, with the period it is stretched onto (None for its own)
    and the scaling its values are multiplied by, as SimVascular keeps its Timed Pressure, Pim
    Period and Pim Scaling: so that opening the condition again shows what was chosen, and
    `applied_pim` is what the solver gets. A file in mmHg is scaled by `MMHG`. Kept without a Pim
    until one is loaded; `problems` says when one is still missing.
    """

    Ra: float
    Ca: float
    Ram: float
    Cim: float
    Rv: float
    Pv: float = 0.0
    pim: Pressure | None = None
    pim_period: float | None = None
    pim_scale: float = 1.0
    kind = CORONARY_KIND

    def with_pim(self, pim: Pressure, period: float | None = None, scale: float = 1.0) -> "Coronary":
        return replace(self, pim=pim, pim_period=period, pim_scale=scale)

    @property
    def applied_pim(self) -> Pressure | None:
        """The Pim the solver gets: the file's times stretched onto the period, its values scaled.

        The last time is the period itself rather than the file's last time stretched, which can
        miss it in the last bit; svZeroDSolver takes the period as the last time less the first,
        and refuses two that differ at all.
        """
        if self.pim is None:
            return None
        period = self.pim_period or self.pim.period
        factor = period / self.pim.period
        time = tuple(t * factor for t in self.pim.time[:-1]) + (float(period),)
        return Pressure(time, tuple(p * self.pim_scale for p in self.pim.pressure), self.pim.source)


# -- as typed in the table ----------------------------------------------------
# What each kind's values are, in the order they are typed. The panel shows these as the cell's
# hint, so the order is the one thing to get right and is said once, here.
VALUE_ORDER = {
    INFLOW: ("Q",),
    RCR_KIND: ("Rp", "C", "Rd", "Pd"),
    RESISTANCE_KIND: ("R",),
    CORONARY_KIND: ("Ra", "Ca", "Ram", "Cim", "Rv", "Pv"),
}
# Each kind's last value may be left out, and is then 0: an RCR's distal pressure, a coronary's
# venous pressure. Every other value has to be typed.
OPTIONAL_LAST = (RCR_KIND, CORONARY_KIND)
# What separates a coronary's typed values from what the cell says about its Pim, which is shown
# beside them and is not typed: not one of the separators the values are split on.
PIM_SEPARATOR = "|"


def to_text(condition) -> str:
    """The condition as the table shows it, which is also what parses back to it.

    A waveform is the exception: it is loaded from a file rather than typed, so it is shown as
    the file it came from, and typing a number over that replaces it with a steady flow. Only a
    waveform that came from nowhere -- made in code, or saved before its file was kept -- is
    described instead.
    """
    if condition is None:
        return ""
    if isinstance(condition, Inflow):
        if condition.source:
            return condition.source
        if condition.is_steady:
            return _number(condition.flow[0])
        return (f"waveform, {len(condition.time)} points over {_number(condition.period)} s, "
                f"mean {_number(condition.mean)}")
    if isinstance(condition, RCR):
        return ", ".join(_number(value) for value in (condition.Rp, condition.C, condition.Rd,
                                                      condition.Pd))
    if isinstance(condition, Resistance):
        return _number(condition.R)
    if isinstance(condition, Coronary):
        values = ", ".join(_number(value) for value in (condition.Ra, condition.Ca, condition.Ram,
                                                        condition.Cim, condition.Rv, condition.Pv))
        return f"{values} {PIM_SEPARATOR} Pim: {describe_pim(condition)}"
    raise BoundaryConditionError(f"Not a boundary condition: {condition!r}")


def describe_pim(coronary: Coronary) -> str:
    """A coronary's Pim in a few words: its file, and the scaling and period it is used at when
    they are not the file's own; or that it is still to be loaded."""
    pim = coronary.pim
    if pim is None:
        return "double-click to set"
    text = Path(pim.source).name if pim.source else f"{len(pim.time)} points"
    if coronary.pim_scale != 1.0:
        text += f" x {_number(coronary.pim_scale)}"
    return text + f" over {_number(coronary.applied_pim.period)} s"


def from_text(kind: str, text: str, previous=None, period: float = 1.0):
    """A condition of `kind` out of what was typed, or None for an empty cell.

    `previous` is what the cell held before. It is what keeps a loaded waveform when its
    description is left as it was -- the table hands back the cell's text whether or not it was
    edited -- and `period` is what a steady inflow is spread over, so that it can sit beside
    waveforms of the other inflows without its period disagreeing with theirs.
    """
    text = (text or "").strip()
    if not text:
        return None
    if previous is not None and getattr(previous, "kind", None) == kind and text == to_text(previous):
        return previous
    if kind == CORONARY_KIND:
        # What follows the separator says which Pim is loaded; it is not typed, and a loaded Pim
        # stays with the cap whatever its values are typed as.
        text = text.split(PIM_SEPARATOR)[0].strip()
        if not text:
            return None
    try:
        values = [float(value) for value in _SEPARATORS.split(text) if value]
    except ValueError:
        raise BoundaryConditionError(
            f"Could not read {text!r} as numbers. {kind} takes "
            f"{', '.join(VALUE_ORDER[kind])}, in that order.") from None
    names = VALUE_ORDER[kind]
    required = len(names) - 1 if kind in OPTIONAL_LAST else len(names)
    if not required <= len(values) <= len(names):
        raise BoundaryConditionError(
            f"{kind} takes {', '.join(names)}"
            + (f" ({names[-1]} may be left out, and is then 0)" if kind in OPTIONAL_LAST else "")
            + f"; {len(values)} value(s) were given.")
    if any(not math.isfinite(value) for value in values):
        raise BoundaryConditionError("Every value has to be a finite number.")
    if kind == INFLOW:
        return Inflow.steady(values[0], period)
    if kind == RCR_KIND:
        return RCR(*values)
    if kind == CORONARY_KIND:
        # Each of the five has to be more than 0: the model divides by every one of them, and a
        # zero takes a branch out of it -- a capacitance of 0 is no capacitance, a resistance of
        # 0 a short -- which svZeroDSolver solves without a word as some other condition.
        bad = [name for name, value in zip(names[:5], values[:5]) if not value > 0]
        if bad:
            raise BoundaryConditionError(f"A coronary's {', '.join(bad)} have to be more than 0.")
        coronary = Coronary(*values)
        if isinstance(previous, Coronary) and previous.pim is not None:
            coronary = coronary.with_pim(previous.pim, previous.pim_period, previous.pim_scale)
        return coronary
    return Resistance(*values)


def _number(value: float) -> str:
    return f"{value:.6g}"


# -- kept on the mesh -----------------------------------------------------------
def to_json(conditions: Mapping[int, object]) -> str:
    """`{"<face id>": {...}}`, sorted, so a scene saved twice untouched is the same bytes twice."""
    stored = {}
    for face_id, condition in sorted(conditions.items()):
        if condition is None:
            continue
        if isinstance(condition, Inflow):
            entry = {"type": INFLOW, "t": list(condition.time), "Q": list(condition.flow)}
            if condition.source:
                entry["file"] = condition.source
        elif isinstance(condition, RCR):
            entry = {"type": RCR_KIND, "Rp": condition.Rp, "C": condition.C, "Rd": condition.Rd,
                     "Pd": condition.Pd}
        elif isinstance(condition, Resistance):
            entry = {"type": RESISTANCE_KIND, "R": condition.R}
        elif isinstance(condition, Coronary):
            entry = {"type": CORONARY_KIND, "Ra": condition.Ra, "Ca": condition.Ca,
                     "Ram": condition.Ram, "Cim": condition.Cim, "Rv": condition.Rv, "Pv": condition.Pv}
            if condition.pim is not None:
                entry["t"] = list(condition.pim.time)
                entry["Pim"] = list(condition.pim.pressure)
                entry["Pim_scale"] = condition.pim_scale
                if condition.pim_period is not None:
                    entry["Pim_period"] = condition.pim_period
                if condition.pim.source:
                    entry["file"] = condition.pim.source
        else:
            raise BoundaryConditionError(f"Not a boundary condition: {condition!r}")
        stored[str(face_id)] = entry
    return json.dumps(stored, separators=(",", ":"))


def from_json(text: str) -> dict:
    """The conditions back, by face id. An entry that cannot be read is dropped, not raised:
    a face left without a condition is something the panel shows, and a scene that will not
    open over one bad value is not."""
    if not text:
        return {}
    try:
        stored = json.loads(text)
    except ValueError:
        return {}
    conditions = {}
    for face_id, entry in stored.items():
        try:
            kind = entry["type"]
            if kind == INFLOW:
                condition = Inflow(tuple(float(t) for t in entry["t"]),
                                   tuple(float(q) for q in entry["Q"]), str(entry.get("file", "")))
            elif kind == RCR_KIND:
                condition = RCR(float(entry["Rp"]), float(entry["C"]), float(entry["Rd"]),
                                float(entry.get("Pd", 0.0)))
            elif kind == RESISTANCE_KIND:
                condition = Resistance(float(entry["R"]))
            elif kind == CORONARY_KIND:
                condition = Coronary(*(float(entry[name]) for name in ("Ra", "Ca", "Ram", "Cim", "Rv")),
                                     float(entry.get("Pv", 0.0)))
                if "Pim" in entry:
                    condition = condition.with_pim(
                        Pressure(tuple(float(t) for t in entry["t"]),
                                 tuple(float(p) for p in entry["Pim"]), str(entry.get("file", ""))),
                        float(entry["Pim_period"]) if "Pim_period" in entry else None,
                        float(entry.get("Pim_scale", 1.0)))
            else:
                continue
            conditions[int(face_id)] = condition
        except (KeyError, TypeError, ValueError):     # BoundaryConditionError is a ValueError
            continue
    return conditions


# -- whether they close the model ------------------------------------------------
def problems(conditions: Mapping[int, object], cap_names: Mapping[int, str],
             inlet_face_id: int | None, *, source_required: bool = True) -> list:
    """What stops these conditions from closing the model, in words; empty when nothing does.

    Each is something the package or the solver would otherwise find later and say less
    clearly, or not say at all:

    - a cap with no condition, which the package refuses by name, after centerlines;
    - an inlet that is not an inflow -- centerlines start at the inlet, and the package
      prescribes its flow there, so anything else typed against it would be dropped;
    - no outlet at all. With every cap prescribing flow, nothing sets the pressure, and the
      solver's system is singular; it fails to converge rather than saying why;
    - inflows over different periods, which the package refuses: the cycle is the inlet's;
    - a coronary without its intramyocardial pressure, which is half the condition, or with one
      over another period than the inflows': svZeroDSolver refuses any two time-varying values
      whose periods differ at all, and says only "Inconsistent cardiac cycle period".

    `source_required` is False for a 3D model, which has no centerlines and so no source; the
    other three hold of it as they do of a 0D one.
    """
    found = []
    missing = [name for face_id, name in sorted(cap_names.items())
               if conditions.get(face_id) is None]
    if missing:
        found.append(f"{len(missing)} cap(s) have no boundary condition: {', '.join(missing)}.")
    if not source_required:
        pass
    elif inlet_face_id is None:
        found.append("Choose the source the centerlines start from.")
    elif inlet_face_id in cap_names and not isinstance(conditions.get(inlet_face_id), Inflow):
        found.append(f"The source, {cap_names[inlet_face_id]}, has to be an inflow: the centerlines "
                     "start there and its flow is what is prescribed on them.")
    given = {face_id: condition for face_id, condition in conditions.items()
             if face_id in cap_names and condition is not None}
    if given and not any(isinstance(condition, (RCR, Resistance, Coronary))
                         for condition in given.values()):
        found.append("At least one cap has to be an RCR, a resistance or a coronary: with flow "
                     "prescribed everywhere, nothing sets the pressure.")
    periods = sorted({condition.period for condition in given.values()
                      if isinstance(condition, Inflow) and not condition.is_steady})
    if len(periods) > 1:
        found.append("The inflow waveforms cover different periods ("
                     + ", ".join(f"{period:g} s" for period in periods)
                     + "); they have to cover the same cardiac cycle.")
    cycle = periods[0] if len(periods) == 1 else None
    for face_id, condition in sorted(given.items()):
        if not isinstance(condition, Coronary):
            continue
        if condition.pim is None:
            found.append(f"{cap_names[face_id]} is a coronary without its intramyocardial pressure "
                         "Pim: double-click its values to set it.")
        elif cycle is not None and not same_period(condition.applied_pim.period, cycle):
            found.append(f"The intramyocardial pressure of {cap_names[face_id]} covers "
                         f"{condition.applied_pim.period:g} s, and the inflow {cycle:g} s; they have "
                         "to cover the same cardiac cycle -- set its Pim period to the inflow's.")
    pims = sorted({condition.applied_pim.period for condition in given.values()
                   if isinstance(condition, Coronary) and condition.pim is not None})
    if cycle is None and len(pims) > 1 and not all(same_period(p, pims[0]) for p in pims):
        found.append("The intramyocardial pressures cover different periods ("
                     + ", ".join(f"{period:g} s" for period in pims)
                     + "); they have to cover the same cardiac cycle.")
    return found


def same_period(a: float, b: float) -> bool:
    """Whether two periods are one cycle, to the rounding of the files they were read from:
    0.909091 s and 0.9090909091 s are, and the writer puts a Pim on the inflow's exact period."""
    return math.isclose(a, b, rel_tol=1e-6)


def cycle_period(conditions) -> float:
    """The cardiac cycle: the period of the first waveform among them, or 1 s if all are steady.

    A steady inflow has no period of its own, so it takes this one. That is what lets a steady
    azygous flow sit beside a pulsatile hepatic one without the two disagreeing.
    """
    conditions = list(conditions)
    for condition in conditions:
        if isinstance(condition, Inflow) and not condition.is_steady:
            return condition.period
    # With every inflow steady, a coronary's Pim is the only cycle there is, and the steady
    # inflows are spread over it.
    for condition in conditions:
        if isinstance(condition, Coronary) and condition.pim is not None:
            return condition.applied_pim.period
    return 1.0


def with_period(condition, period: float):
    """A steady inflow respread over `period`, a coronary's Pim put on it exactly when it covers
    the same cycle to rounding; anything else as it is.

    Done through its Pim period, which stretches all its times: svZeroDSolver takes a waveform's
    period as its last time less its first, and refuses two that differ in the last bit, so a Pim
    read as 0.909091 s beside an inflow of 0.9090909091 s would stop the run.
    """
    if isinstance(condition, Inflow) and condition.is_steady and condition.period != period:
        return Inflow.steady(condition.flow[0], period)
    if (isinstance(condition, Coronary) and condition.pim is not None
            and condition.applied_pim.period != period
            and same_period(condition.applied_pim.period, period)):
        return replace(condition, pim_period=period)
    return condition


# -- SimVascular's files ---------------------------------------------------------
def read_flow_file(path) -> Inflow:
    """A waveform from a `.flow` file: time and flow per line, space- or comma-separated.

    As written, sign and all. A file made for svMultiPhysics carries an inflow as *negative* --
    its flow is measured along the outward normal -- and the reduced-order solvers take it as
    positive. Which one a file is cannot be read off the file, so the caller says so: the panels
    plot it and ask, and `as_inflow` guesses from the mean for a script that has to.
    """
    times, flows = _read_two_columns(path, "a flow")
    return Inflow(tuple(times), tuple(flows), str(Path(path).resolve()))


def read_pressure_file(path) -> Pressure:
    """A coronary's intramyocardial pressure from a file laid out as a `.flow` is: time and
    pressure per line, `#` lines skipped. In dyn/cm², as written."""
    times, pressures = _read_two_columns(path, "a pressure")
    return Pressure(tuple(times), tuple(pressures), str(Path(path).resolve()))


def _read_two_columns(path, what: str):
    times, values = [], []
    for number, line in enumerate(Path(path).read_text().splitlines(), start=1):
        fields = [value for value in _SEPARATORS.split(line.strip()) if value]
        if not fields or fields[0].startswith("#"):
            continue
        try:
            times.append(float(fields[0]))
            values.append(float(fields[1]))
        except (IndexError, ValueError):
            raise BoundaryConditionError(
                f"{path}, line {number}: expected a time and {what}, found {line.strip()!r}."
            ) from None
    return times, values


def as_inflow(waveform: Inflow):
    """The waveform as an inflow, positive in, and whether it had to be turned round to be one.

    A waveform whose mean is negative is taken to be one written for svMultiPhysics, where
    inflow is negative, rather than a flow out of the model: an outlet prescribing flow is
    not something the panels offer. Returned with a flag rather than flipped quietly, so that a
    caller can say it did. Not what the panels do: they plot the waveform and leave negating it to
    the user (`SimVascularROMLib.BoundaryConditionsTable.WaveformReview`), since a waveform that
    is mostly backflow has a negative mean too.
    """
    if waveform.mean < 0:
        return waveform.flipped(), True
    return waveform, False


def write_flow_file(path, inflow: Inflow) -> Path:
    path = Path(path)
    path.write_text("".join(f"{t:.10g} {q:.10g}\n" for t, q in zip(inflow.time, inflow.flow)))
    return path


def read_rcrt(path) -> dict:
    """`{face name: RCR}` from a SimVascular `rcrt.dat`.

    The file opens with the number of distal pressure points, and each face's block repeats it:
    the face name, Rp, C, Rd, then that many time / distal pressure pairs. Only a constant
    distal pressure is read, because that is what the 0D writer takes; a varying one is refused
    by name rather than averaged.
    """
    lines = [line.strip() for line in Path(path).read_text().splitlines()]
    lines = [line for line in lines if line]
    if not lines:
        return {}
    keyword = lines[0]
    found = {}
    index = 1
    while index < len(lines):
        if lines[index] != keyword:
            raise BoundaryConditionError(
                f"{path}: expected a block starting with {keyword!r}, found {lines[index]!r}.")
        try:
            name = lines[index + 1]
            rp, c, rd = (float(value) for value in lines[index + 2:index + 5])
            count = int(keyword)
            pressures = [float(lines[index + 5 + point].split()[1]) for point in range(count)]
        except (IndexError, ValueError):
            raise BoundaryConditionError(f"{path}: the block starting at line {index + 1} is "
                                         "not a face name, Rp, C, Rd and its pressures.") from None
        if len(set(pressures)) > 1:
            raise BoundaryConditionError(f"{path}: {name} has a distal pressure that varies in "
                                         "time, which the 0D writer cannot take.")
        found[name] = RCR(rp, c, rd, pressures[0] if pressures else 0.0)
        index += 5 + count
    return found


def write_rcrt(path, conditions: Mapping[str, RCR]) -> Path:
    """SimVascular's `rcrt.dat`, with the distal pressure held over the cycle 0 to 1.

    The pressure's times do not matter to the 0D writer, which takes the pressure and checks
    that both are the same; 0 and 1 are what SimVascular itself writes.
    """
    path = Path(path)
    blocks = ["2"]
    for name, rcr in conditions.items():
        blocks += ["2", name, f"{rcr.Rp:.10g}", f"{rcr.C:.10g}", f"{rcr.Rd:.10g}",
                   f"0.0 {rcr.Pd:.10g}", f"1.0 {rcr.Pd:.10g}"]
    path.write_text("\n".join(blocks) + "\n")
    return path


def read_resistance(path) -> dict:
    """`{face name: Resistance}` from a SimVascular `resistance.dat`: name, R, and a Pd of 0.

    SimVascular's file can carry a distal pressure as a third column, and a resistance here has
    none (see the module docstring). One of 0 is the same condition and is read; any other is
    refused by name rather than dropped, which would lower that outlet's pressure by it.
    """
    found = {}
    for number, line in enumerate(Path(path).read_text().splitlines(), start=1):
        values = line.split()
        if not values:
            continue
        try:
            if len(values) not in (2, 3):
                raise ValueError
            numbers = [float(value) for value in values[1:]]
        except ValueError:
            raise BoundaryConditionError(f"{path}, line {number}: expected a face name and a "
                                         f"resistance, found {line.strip()!r}.") from None
        if len(numbers) == 2 and numbers[1] != 0.0:
            raise BoundaryConditionError(
                f"{path}: {values[0]} has a distal pressure of {numbers[1]:g}, and a resistance "
                "here has none. Give that outlet an RCR instead.")
        found[values[0]] = Resistance(numbers[0])
    return found


def write_resistance(path, conditions: Mapping[str, Resistance]) -> Path:
    path = Path(path)
    path.write_text("".join(f"{name} {resistance.R:.10g}\n"
                            for name, resistance in conditions.items()))
    return path


# SimVascular's coefficients for an open-loop coronary, in the order `cort.dat` lists them.
CORT_COEFFICIENTS = ("q0", "q1", "q2", "p0", "p1", "p2", "b0", "b1", "b2", "dQinidT", "dPinidT")


def coronary_coefficients(coronary: Coronary) -> dict:
    """The coronary as SimVascular's `cort.dat` stores it: the coefficients of

        p0 P + p1 dP/dt + p2 d2P/dt2 = q0 Q + q1 dQ/dt + q2 d2Q/dt2 + b0 Pim + b1 dPim/dt + b2 d2Pim/dt2

    in terms of the five values, after Kim et al. (2010), with no venous microvascular
    resistance: SimVascular's own (`sv4guiROMSimulationUtils::CreateCORTFileContent`), and the
    relations the package checks when it works the values back out of a `cort.dat`
    (`svromutils.io_1d.coronary_sv_to_oned`). Written in them because that is the only form the
    package takes a coronary in.
    """
    Ra, Ram, Rv, Ca, Cim = coronary.Ra, coronary.Ram, coronary.Rv, coronary.Ca, coronary.Cim
    return {
        "q0": Ra + Ram + Rv,
        "q1": Ra * Ca * (Ram + Rv) + Cim * (Ra + Ram) * Rv,
        "q2": Ca * Cim * Ra * Ram * Rv,
        "p0": 1.0,
        "p1": Ram * Ca + Rv * (Ca + Cim),
        "p2": Ca * Cim * Ram * Rv,
        "b0": 0.0,
        "b1": Cim * Rv,
        "b2": 0.0,
        "dQinidT": 0.0,
        "dPinidT": 0.0,
    }


def write_cort(path, conditions: Mapping[str, Coronary]) -> Path:
    """SimVascular's `cort.dat`: the number of Pim points; then per coronary that number again,
    its face name, its eleven coefficients and its time / Pim pairs.

    One number of points for every block, the largest, because the package reads the file's first
    line as the marker each block starts with, and a block of fewer points is ended by the next
    marker all the same. The coefficients to 17 figures, so the values the package works back out
    of them are the ones set; the times to 10, as the inflow's are, so that a Pim on the inflow's
    period is read back on exactly it.
    """
    path = Path(path)
    marker = str(max(len(coronary.pim.time) for coronary in conditions.values()))
    lines = [marker]
    for name, coronary in conditions.items():
        coefficients = coronary_coefficients(coronary)
        pim = coronary.applied_pim
        lines += [marker, name] + [f"{coefficients[key]:.17g}" for key in CORT_COEFFICIENTS]
        lines += [f"{t:.10g} {p:.17g}" for t, p in zip(pim.time, pim.pressure)]
    path.write_text("\n".join(lines) + "\n")
    return path
