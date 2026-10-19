import json

import pytest

from svromsetup import boundary_conditions as bcs
from svromsetup.boundary_conditions import RCR, BoundaryConditionError, Inflow, Resistance

CAPS = {2: "cap_inlet", 3: "cap_right", 4: "cap_left"}


# -- typed in the table ----------------------------------------------------------
@pytest.mark.parametrize("condition", [Inflow.steady(12.5), RCR(100.0, 1e-4, 1500.0, 0.0),
                                       Resistance(250.0)])
def test_what_the_table_shows_reads_back_as_the_same_condition(condition):
    text = bcs.to_text(condition)
    assert bcs.from_text(condition.kind, text) == condition


def test_an_rcr_distal_pressure_may_be_left_out_and_a_resistance_has_none():
    """A resistance is R alone, as svMultiPhysics' is: a second value is refused rather than
    read as a distal pressure the 3D solver would then drop."""
    assert bcs.from_text(bcs.RCR_KIND, "100, 1e-4 1500") == RCR(100.0, 1e-4, 1500.0, 0.0)
    assert bcs.from_text(bcs.RESISTANCE_KIND, "250") == Resistance(250.0)
    with pytest.raises(BoundaryConditionError, match="takes R; 2 value"):
        bcs.from_text(bcs.RESISTANCE_KIND, "250, 6666")


def test_a_waveform_survives_its_own_description_and_is_replaced_by_a_number():
    """The table hands back the cell's text whether it was edited or not, and a waveform's text
    is a description that cannot be parsed back into one. Left alone it has to stay the waveform;
    typed over with a number, it becomes a steady flow."""
    waveform = Inflow((0.0, 0.4, 0.8), (5.0, 20.0, 5.0))
    text = bcs.to_text(waveform)
    assert "waveform" in text
    assert bcs.from_text(bcs.INFLOW, text, previous=waveform) is waveform
    assert bcs.from_text(bcs.INFLOW, "7", previous=waveform, period=0.8) == Inflow.steady(7.0, 0.8)


def test_a_waveform_read_from_a_file_is_shown_as_that_file(tmp_path):
    """The cell says which file a cap's inflow is, and leaving it as it is keeps the waveform."""
    path = tmp_path / "cap_RSVC.flow"
    path.write_text("0.0 -5.0\n0.4 -15.0\n0.8 -5.0\n")
    waveform, turned = bcs.as_inflow(bcs.read_flow_file(path))
    assert turned
    assert bcs.to_text(waveform) == str(path.resolve()), "turned round, and still that file's"
    assert bcs.from_text(bcs.INFLOW, str(path.resolve()), previous=waveform) is waveform
    # The file is where it came from, not what it is: the same numbers are the same inflow.
    assert waveform == Inflow((0.0, 0.4, 0.8), (5.0, 15.0, 5.0))
    # A steady flow read from a file is still that file's, rather than a number nobody typed.
    (tmp_path / "steady.flow").write_text("0.0 3.0\n0.8 3.0\n")
    assert bcs.to_text(bcs.read_flow_file(tmp_path / "steady.flow")).endswith("steady.flow")
    # And the file is kept with the scene, so it is shown again after reopening it.
    assert bcs.from_json(bcs.to_json({3: waveform}))[3].source == str(path.resolve())


@pytest.mark.parametrize("kind, text", [(bcs.RCR_KIND, "100, 200"), (bcs.RESISTANCE_KIND, "1 2 3"),
                                        (bcs.INFLOW, "fast"), (bcs.RCR_KIND, "1, nan, 3")])
def test_values_that_are_not_a_condition_say_what_was_expected(kind, text):
    with pytest.raises(BoundaryConditionError, match="takes|numbers|finite"):
        bcs.from_text(kind, text)


def test_an_empty_cell_is_no_condition_rather_than_an_error():
    assert bcs.from_text(bcs.RCR_KIND, "  ") is None


def test_a_waveform_has_to_be_one_cycle_from_zero():
    with pytest.raises(BoundaryConditionError, match="two points"):
        Inflow((0.0,), (1.0,))
    with pytest.raises(BoundaryConditionError, match="start at time 0"):
        Inflow((0.1, 1.0), (1.0, 1.0))
    with pytest.raises(BoundaryConditionError, match="increase"):
        Inflow((0.0, 0.5, 0.5), (1.0, 1.0, 1.0))


def test_the_mean_of_a_waveform_is_over_time_not_over_points():
    # Most of the points crowd the first tenth of the cycle; a mean over points would read 9.
    inflow = Inflow((0.0, 0.025, 0.05, 0.075, 0.1, 1.0), (10.0, 10.0, 10.0, 10.0, 10.0, 0.0))
    assert inflow.mean == pytest.approx(0.1 * 10 + 0.9 * 5.0)


# -- kept on the mesh ----------------------------------------------------------------
def test_the_conditions_survive_being_saved_and_a_bad_entry_is_dropped_not_raised():
    conditions = {2: Inflow((0.0, 0.5, 1.0), (1.0, 3.0, 1.0)), 3: RCR(1.0, 2.0, 3.0, 4.0),
                  4: Resistance(5.0)}
    text = bcs.to_json(conditions)
    assert bcs.from_json(text) == conditions
    assert bcs.to_json(conditions) == text, "saving twice untouched should be the same bytes"

    stored = json.loads(text)
    stored["3"] = {"type": "RCR", "Rp": "not a number"}
    stored["9"] = {"type": "Something else"}
    assert bcs.from_json(json.dumps(stored)) == {2: conditions[2], 4: conditions[4]}
    assert bcs.from_json("not json") == {}


# -- whether they close the model ---------------------------------------------------------
def test_a_complete_set_has_no_problems():
    conditions = {2: Inflow.steady(10.0), 3: Resistance(100.0), 4: Inflow.steady(2.0)}
    assert bcs.problems(conditions, CAPS, 2) == []


def test_every_problem_is_named():
    assert "cap_left" in " ".join(bcs.problems({2: Inflow.steady(1.0), 3: Resistance(1.0)}, CAPS, 2))
    assert "has to be an inflow" in " ".join(
        bcs.problems({2: Resistance(1.0), 3: Resistance(1.0), 4: Resistance(1.0)}, CAPS, 2))
    assert "nothing sets the pressure" in " ".join(
        bcs.problems({2: Inflow.steady(1.0), 3: Inflow.steady(1.0), 4: Inflow.steady(1.0)}, CAPS, 2))
    assert "source" in " ".join(bcs.problems({}, CAPS, None))
    periods = bcs.problems({2: Inflow((0.0, 0.4, 0.8), (1.0, 2.0, 1.0)),
                            3: Inflow((0.0, 0.5, 1.0), (1.0, 2.0, 1.0)), 4: Resistance(1.0)}, CAPS, 2)
    assert "different periods" in " ".join(periods)


def test_a_steady_inflow_takes_the_period_of_the_waveforms_beside_it():
    waveform = Inflow((0.0, 0.4, 0.8), (1.0, 2.0, 1.0))
    period = bcs.cycle_period([Inflow.steady(3.0), waveform, Resistance(1.0)])
    assert period == 0.8
    assert bcs.with_period(Inflow.steady(3.0), period) == Inflow.steady(3.0, 0.8)
    assert bcs.with_period(waveform, 2.0) is waveform
    assert bcs.cycle_period([Inflow.steady(3.0)]) == 1.0


# -- SimVascular's files -------------------------------------------------------------------
def test_rcrt_and_resistance_files_round_trip(tmp_path):
    rcrs = {"cap_right": RCR(100.0, 1e-4, 1500.0, 10.0), "cap_left": RCR(200.0, 2e-4, 2500.0, 0.0)}
    assert bcs.read_rcrt(bcs.write_rcrt(tmp_path / "rcrt.dat", rcrs)) == rcrs
    resistances = {"cap_right": Resistance(100.0), "cap_left": Resistance(200.0)}
    assert bcs.read_resistance(bcs.write_resistance(tmp_path / "resistance.dat", resistances)) \
        == resistances


def test_a_simvascular_resistance_with_a_distal_pressure_is_refused_by_name(tmp_path):
    path = tmp_path / "resistance.dat"
    path.write_text("cap_lpa 1200 0\ncap_rpa 900\n")
    assert bcs.read_resistance(path) == {"cap_lpa": Resistance(1200.0), "cap_rpa": Resistance(900.0)}
    path.write_text("cap_lpa 1200 6666\n")
    with pytest.raises(BoundaryConditionError, match="cap_lpa.*RCR"):
        bcs.read_resistance(path)


def test_an_rcrt_written_by_simvascular_is_read(tmp_path):
    path = tmp_path / "rcrt.dat"
    path.write_text("2\n2\ncap_aorta\n121\n0.000155\n1212\n0.0 0\n1.0 0\n"
                    "2\ncap_carotid\n500\n0.0001\n9000\n0.0 0\n1.0 0\n")
    assert bcs.read_rcrt(path) == {"cap_aorta": RCR(121.0, 0.000155, 1212.0, 0.0),
                                   "cap_carotid": RCR(500.0, 0.0001, 9000.0, 0.0)}


def test_an_rcrt_with_a_varying_distal_pressure_is_refused_by_name(tmp_path):
    path = tmp_path / "rcrt.dat"
    path.write_text("2\n2\ncap_aorta\n121\n0.000155\n1212\n0.0 0\n1.0 10\n")
    with pytest.raises(BoundaryConditionError, match="cap_aorta"):
        bcs.read_rcrt(path)


def test_a_flow_file_round_trips_and_a_negative_one_is_turned_round_and_said_to_be(tmp_path):
    inflow = Inflow((0.0, 0.5, 1.0), (4.0, 9.0, 4.0))
    assert bcs.read_flow_file(bcs.write_flow_file(tmp_path / "inflow.flow", inflow)) == inflow

    # svMultiPhysics' convention: inflow is negative, along the outward normal.
    (tmp_path / "3d.flow").write_text("0.0 -4.0\n0.5 -9.0\n1.0 -4.0\n")
    turned, was_turned = bcs.as_inflow(bcs.read_flow_file(tmp_path / "3d.flow"))
    assert was_turned and turned == inflow
    assert bcs.as_inflow(inflow) == (inflow, False)

    (tmp_path / "bad.flow").write_text("0.0 1.0\n0.5\n")
    with pytest.raises(BoundaryConditionError, match="line 2"):
        bcs.read_flow_file(tmp_path / "bad.flow")


def test_the_package_reads_the_files_as_they_are_written(tmp_path):
    """The files are written for `svromutils` to read, so its own reader is the judge."""
    pytest.importorskip("vmtk")
    io_1d = pytest.importorskip("svromutils.io_1d")
    from svromutils.parameters import Parameters

    bcs.write_rcrt(tmp_path / "rcrt.dat", {"cap_right": RCR(100.0, 1e-4, 1500.0, 10.0)})
    bcs.write_resistance(tmp_path / "resistance.dat", {"cap_left": Resistance(200.0)})
    parameters = Parameters()
    parameters.outflow_bc_type = ["rcrt.dat", "resistance.dat"]
    parameters.outflow_bc_file = str(tmp_path)
    values, kinds = io_1d.read_variable_outflow_bcs(parameters)
    assert values == {"cap_right": [100.0, 1e-4, 1500.0, 10.0], "cap_left": ["200", 0.0]}
    assert kinds == {"cap_right": "rcr", "cap_left": "resistance"}


# -- coronaries ---------------------------------------------------------------------
PIM = bcs.Pressure((0.0, 0.3, 0.6, 0.8), (0.0, 1.5e5, 2.0e4, 0.0), "/somewhere/pim.dat")
LCA = dict(Ra=1.8e5, Ca=1.1e-6, Ram=2.9e5, Cim=9e-6, Rv=1e5)


def test_a_coronary_is_typed_in_simvascular_order_with_a_venous_pressure_that_may_be_left_out():
    """Ra, Ca, Ra-micro, Cim, Rv, as SimVascular's dialog takes them, and Pv after them."""
    five = bcs.from_text(bcs.CORONARY_KIND, "1.8e5, 1.1e-6, 2.9e5, 9e-6, 1e5")
    assert five == bcs.Coronary(**LCA)
    assert five.pim is None, "set first, the Pim is loaded after"
    assert bcs.from_text(bcs.CORONARY_KIND, "1.8e5 1.1e-6 2.9e5 9e-6 1e5 1000").Pv == 1000.0
    with pytest.raises(BoundaryConditionError, match="Pv may be left out"):
        bcs.from_text(bcs.CORONARY_KIND, "1.8e5, 1.1e-6, 2.9e5, 9e-6")
    with pytest.raises(BoundaryConditionError, match="Ram, Cim have to be more than 0"):
        bcs.from_text(bcs.CORONARY_KIND, "1.8e5, 1.1e-6, 0, -9e-6, 1e5")


def test_a_pim_is_stretched_onto_its_period_and_scaled_as_simvascular_does():
    """SimVascular's Pim Period and Pim Scaling: a file in mmHg, over the file's own 1 s, used
    over a 0.8 s cycle in dyn/cm²."""
    mmhg = bcs.Pressure((0.0, 0.3, 1.0), (0.0, 100.0, 0.0), "/somewhere/lv.dat")
    coronary = bcs.Coronary(**LCA).with_pim(mmhg, period=0.8, scale=bcs.MMHG)
    applied = coronary.applied_pim
    assert applied.time == pytest.approx((0.0, 0.24, 0.8)) and applied.time[-1] == 0.8
    assert applied.pressure == pytest.approx((0.0, 133322.0, 0.0))
    assert coronary.pim == mmhg, "the file's own values kept, to show again"
    assert bcs.Coronary(**LCA).with_pim(PIM).applied_pim == PIM, "unstretched, unscaled by default"


def test_a_coronary_keeps_its_pim_whatever_its_values_are_typed_as():
    loaded = bcs.Coronary(**LCA).with_pim(PIM, scale=2.0)
    text = bcs.to_text(loaded)
    assert text.endswith("| Pim: pim.dat x 2 over 0.8 s")
    assert bcs.from_text(bcs.CORONARY_KIND, text, loaded) == loaded
    retyped = bcs.from_text(bcs.CORONARY_KIND, "2e5, 1.1e-6, 2.9e5, 9e-6, 1e5", loaded)
    assert retyped.Ra == 2e5 and (retyped.pim, retyped.pim_scale) == (PIM, 2.0)
    assert bcs.to_text(bcs.Coronary(1, 1, 1, 1, 1)).endswith("| Pim: double-click to set")


def test_a_coronary_survives_being_saved_with_and_without_its_pim():
    conditions = {3: bcs.Coronary(**LCA, Pv=500.0).with_pim(PIM, period=0.9, scale=1333.22),
                  4: bcs.Coronary(1.0, 2.0, 3.0, 4.0, 5.0)}
    back = bcs.from_json(bcs.to_json(conditions))
    assert back == conditions
    assert back[3].pim.source == "/somewhere/pim.dat" and back[3].pim_period == 0.9


def test_a_coronary_sets_the_pressure_but_needs_its_pim_on_the_inflows_cycle():
    inflow = Inflow((0.0, 0.4, 0.8), (5.0, 15.0, 5.0))
    coronary = bcs.Coronary(**LCA)
    found = " ".join(bcs.problems({2: inflow, 3: coronary, 4: coronary.with_pim(PIM)}, CAPS, 2))
    assert "cap_right is a coronary without its intramyocardial pressure" in found
    assert "nothing sets the pressure" not in found and "cap_left" not in found
    other = bcs.Pressure((0.0, 0.5, 1.0), (0.0, 1e5, 0.0))
    assert "intramyocardial pressure of cap_left covers 1 s, and the inflow 0.8 s" in " ".join(
        bcs.problems({2: inflow, 3: Resistance(1.0), 4: coronary.with_pim(other)}, CAPS, 2))
    assert bcs.problems({2: inflow, 3: Resistance(1.0), 4: coronary.with_pim(other, period=0.8)},
                        CAPS, 2) == [], "its Pim period set to the inflow's"
    # To rounding, it is the same cycle -- and is written on the inflow's exactly.
    rounded = bcs.Pressure((0.0, 0.4, 0.8000003), (0.0, 1e5, 0.0))
    assert bcs.problems({2: inflow, 3: Resistance(1.0), 4: coronary.with_pim(rounded)}, CAPS, 2) == []
    assert bcs.with_period(coronary.with_pim(rounded), 0.8).applied_pim.time[-1] == 0.8


def test_with_every_inflow_steady_the_pim_is_the_cycle():
    coronary = bcs.Coronary(**LCA).with_pim(PIM)
    assert bcs.cycle_period([Inflow.steady(3.0), coronary]) == 0.8
    assert bcs.cycle_period([Inflow((0.0, 0.5, 1.0), (1.0, 2.0, 1.0)), coronary]) == 1.0, \
        "an inflow waveform's cycle first; problems says when the Pim's differs"


def test_a_pim_is_read_from_a_file_as_a_flow_is(tmp_path):
    path = tmp_path / "pim.dat"
    path.write_text("# left ventricular pressure, dyn/cm2\r\n0.0\t0\r\n0.3\t150000\r\n0.8\t0\r\n")
    pim = bcs.read_pressure_file(path)
    assert pim.time == (0.0, 0.3, 0.8) and pim.pressure == (0.0, 150000.0, 0.0)
    assert pim.source == str(path.resolve())
    with pytest.raises(BoundaryConditionError, match="has to start at time 0"):
        bcs.Pressure((0.1, 0.5), (0.0, 1.0))


def test_cort_dat_holds_simvascular_coefficients():
    """The coefficients SimVascular itself writes for these five values
    (sv4guiROMSimulationUtils::CreateCORTFileContent), worked by hand for a round set."""
    c = bcs.coronary_coefficients(bcs.Coronary(Ra=1.0, Ca=2.0, Ram=3.0, Cim=4.0, Rv=5.0))
    assert (c["q0"], c["q1"], c["q2"]) == (9.0, 1 * 2 * 8 + 4 * 4 * 5, 2 * 4 * 1 * 3 * 5)
    assert (c["p0"], c["p1"], c["p2"], c["b1"]) == (1.0, 3 * 2 + 5 * 6, 2 * 4 * 3 * 5, 4 * 5)


def test_the_package_reads_a_coronary_back_as_it_was_set(tmp_path):
    """`cort.dat` holds a coronary as SimVascular's coefficients; the package works the values
    back out of them, and its own reader is the judge of whether they are the ones set."""
    pytest.importorskip("vmtk")
    io_1d = pytest.importorskip("svromutils.io_1d")
    from svromutils.parameters import Parameters

    coronaries = {"cap_right": bcs.Coronary(Ra=6.31e4, Ca=3.52e-6, Ram=1.03e5, Cim=2.86e-5,
                                            Rv=3.45e4).with_pim(PIM, scale=2.0),
                  "cap_left": bcs.Coronary(Ra=2.1e6, Ca=9e-8, Ram=3.5e6, Cim=7.5e-7, Rv=1.2e6).with_pim(
                      bcs.Pressure((0.0, 0.8), (1000.0, 1000.0)))}
    bcs.write_cort(tmp_path / "cort.dat", coronaries)
    parameters = Parameters()
    parameters.outflow_bc_type = ["cort.dat"]
    parameters.outflow_bc_file = str(tmp_path)
    values, kinds = io_1d.read_variable_outflow_bcs(parameters)
    assert kinds == {"cap_right": "coronary", "cap_left": "coronary"}
    for name, coronary in coronaries.items():
        read = values[name]["var"]
        for key, attribute in (("Ra1", "Ra"), ("Ra2", "Ram"), ("Rv1", "Rv"), ("Ca", "Ca"), ("Cc", "Cim")):
            assert read[key] == pytest.approx(getattr(coronary, attribute), rel=1e-12)
        assert values[name]["time"] == list(coronary.applied_pim.time)
        assert values[name]["pressure"] == pytest.approx(list(coronary.applied_pim.pressure))
