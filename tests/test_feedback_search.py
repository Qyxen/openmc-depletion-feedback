import csv
from types import MethodType, SimpleNamespace

import pytest

from openmc_feedbacks import FeedbackControl
from openmc_feedbacks.api import _FeedbackOperator, _keff_stats, _nominal_keff


def make_operator(tmp_path, keff_at_concentration, **control_options):
    options = {
        "material_name": "Coolant",
        "nuclides": ("B",),
        "reference_concentration": 1000.0,
        "minimum": 0.0,
        "maximum": 1000.0,
        "target_keff": 1.01,
        "tolerance": 1.0,
        "keff_tolerance": 1.0e-8,
        "max_iterations": 12,
        "history_file": "feedback.csv",
    }
    options.update(control_options)
    control = FeedbackControl(**options)
    # Bypass OpenMC initialization: these tests exercise the root-search
    # control flow with deterministic keff values, not transport physics.
    operator = object.__new__(_FeedbackOperator)
    operator.control = control
    operator.history_path = tmp_path / "feedback.csv"
    operator._step = 0
    operator._trial = 0
    operator._previous_concentration = None
    operator._previous_boc_keff = None
    operator.expected_step_keff_loss = control.expected_step_keff_loss
    operator.search_level = 0
    operator.calls = []
    with operator.history_path.open("w", newline="") as f:
        csv.writer(f).writerow(
            ["step", "trial", "phase", "concentration", "keff",
             "keff_std_dev", "target_keff", "lower_bound", "upper_bound",
             "search_level", "particles", "batches", "selected",
             "keff_error", "status", "stopping_reason"]
        )

    def transport(self, vec, source_rate, concentration, final=False, search_level=None):
        self.calls.append((concentration, final, search_level))
        value = keff_at_concentration(concentration, search_level, final)
        if isinstance(value, tuple):
            value = SimpleNamespace(nominal_value=value[0], std_dev=value[1])
        else:
            value = SimpleNamespace(nominal_value=value, std_dev=0.0)
        return SimpleNamespace(k=value)

    operator._transport = MethodType(transport, operator)
    return operator


def test_nominal_keff_supports_scalar_and_uncertain_values():
    assert _nominal_keff(SimpleNamespace(k=1.002)) == pytest.approx(1.002)
    assert _nominal_keff(
        SimpleNamespace(k=SimpleNamespace(nominal_value=1.002, std_dev=0.0001))
    ) == pytest.approx(1.002)


def test_bisection_finds_decreasing_boron_response_and_writes_history(tmp_path):
    operator = make_operator(tmp_path, lambda ppm, level, final: 1.02 - 0.00002 * ppm)

    result = operator((object(),), source_rate=1.0)

    assert _nominal_keff(result) == pytest.approx(1.01)
    assert operator.calls[-1] == (500.0, True, None)
    with (tmp_path / "feedback.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    selected = next(row for row in rows if row["selected"] == "1")
    assert selected["concentration"] == "500.00000000"
    assert selected["keff"] == "1.01000000"
    assert selected["keff_std_dev"] == "0.00000000"
    assert selected["status"] == "target_met"
    assert selected["stopping_reason"] == "keff_tolerance"


def test_bisection_handles_increasing_response(tmp_path):
    operator = make_operator(tmp_path, lambda ppm, level, final: 1.0 + 0.00002 * ppm)
    operator.control.target_keff = 1.01

    operator((object(),), source_rate=1.0)

    assert operator.calls[-1] == (500.0, True, None)


@pytest.mark.parametrize(
    "response, expected_concentration, expected_status",
    [
        (lambda ppm, level, final: 1.02 - 0.000001 * ppm, 1000.0, "upper_bound"),
        (lambda ppm, level, final: 1.00 - 0.000001 * ppm, 0.0, "lower_bound"),
    ],
)
def test_unbracketed_target_selects_closest_endpoint(
    tmp_path, response, expected_concentration, expected_status
):
    operator = make_operator(tmp_path, response)

    operator((object(),), source_rate=1.0)

    assert operator.calls[-1] == (expected_concentration, True, None)
    with (tmp_path / "feedback.csv").open(newline="") as f:
        row = next(row for row in csv.DictReader(f) if row["selected"] == "1")
    assert row["stopping_reason"] == expected_status


def test_warm_start_searches_local_interval_first(tmp_path):
    operator = make_operator(
        tmp_path,
        lambda ppm, level, final: 1.02 - (0.01 / 620.0) * ppm,
        warm_start_width=100.0,
    )
    operator._previous_concentration = 600.0

    operator((object(),), source_rate=1.0)

    assert operator.calls[:2] == [(500.0, False, None), (700.0, False, None)]
    assert operator.calls[-1][0] == pytest.approx(620.0)
    assert operator.calls[-1][1:] == (True, None)


def test_search_refines_histories_when_keff_uncertainty_is_too_large(tmp_path):
    def noisy_response(ppm, level, final):
        mean = 1.02 - 0.00002 * ppm
        sigma = 0.0004 if level in (None, 0) else 0.00005
        return mean, sigma

    operator = make_operator(
        tmp_path,
        noisy_response,
        keff_tolerance=0.0002,
        max_search_refinements=2,
        confidence_z=2.0,
    )

    result = operator((object(),), source_rate=1.0)

    assert _keff_stats(result) == pytest.approx((1.01, 0.0004))
    assert (500.0, False, 1) in operator.calls
    with (tmp_path / "feedback.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    refinements = [row for row in rows if row["phase"] == "refinement"]
    assert refinements
    assert refinements[-1]["search_level"] == "1"


def test_end_target_compensates_with_filtered_previous_step_loss(tmp_path):
    operator = make_operator(
        tmp_path,
        lambda ppm, level, final: 1.0,
        target_mode="end",
        target_keff=1.005,
        expected_step_keff_loss=0.01,
        loss_filter_alpha=0.5,
    )
    operator._previous_concentration = 300.0
    operator._previous_boc_keff = 1.02
    evaluated = {}

    def evaluate(self, vec, source_rate, concentration, target, lower, upper,
                 phase="search", level=None):
        evaluated["concentration"] = concentration
        evaluated["phase"] = phase
        return 1.012, 0.0001

    operator._evaluate = MethodType(evaluate, operator)

    target = operator._target_for_call((object(),), source_rate=1.0)

    assert evaluated == {"concentration": 300.0, "phase": "loss_measurement"}
    assert operator.expected_step_keff_loss == pytest.approx(0.009)
    assert target == pytest.approx(1.014)
