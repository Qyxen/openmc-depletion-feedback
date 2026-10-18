import csv
from types import MethodType, SimpleNamespace

import pytest

from openmc_depletion_feedback import FeedbackControl
from openmc_depletion_feedback.api import _FeedbackOperator, _nominal_keff


def make_operator(tmp_path, keff_at_concentration, **control_options):
    control = FeedbackControl(
        material_name="Coolant",
        nuclides=("B",),
        reference_concentration=1000.0,
        minimum=0.0,
        maximum=1000.0,
        target_keff=1.01,
        tolerance=1.0,
        keff_tolerance=1.0e-8,
        max_iterations=12,
        history_file="feedback.csv",
        **control_options,
    )
    # Bypass OpenMC initialization: these tests exercise the root-search
    # control flow with deterministic keff values, not transport physics.
    operator = object.__new__(_FeedbackOperator)
    operator.control = control
    operator.history_path = tmp_path / "feedback.csv"
    operator._step = 0
    operator.calls = []
    with operator.history_path.open("w", newline="") as f:
        csv.writer(f).writerow(
            ["step", "concentration", "keff", "target_keff", "status"]
        )

    def transport(self, vec, source_rate, concentration, final=False):
        self.calls.append((concentration, final))
        return SimpleNamespace(k=keff_at_concentration(concentration))

    operator._transport = MethodType(transport, operator)
    return operator


def test_nominal_keff_supports_scalar_and_uncertain_values():
    assert _nominal_keff(SimpleNamespace(k=1.002)) == pytest.approx(1.002)
    assert _nominal_keff(
        SimpleNamespace(k=SimpleNamespace(nominal_value=1.002, std_dev=0.0001))
    ) == pytest.approx(1.002)


def test_bisection_finds_decreasing_boron_response_and_writes_history(tmp_path):
    operator = make_operator(tmp_path, lambda ppm: 1.02 - 0.00002 * ppm)

    result = operator((object(),), source_rate=1.0)

    assert result.k == pytest.approx(1.01)
    assert operator.calls[-1] == (500.0, True)
    with (tmp_path / "feedback.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    assert rows == [{
        "step": "0",
        "concentration": "500.0",
        "keff": "1.01",
        "target_keff": "1.01",
        "status": "controlled",
    }]


def test_bisection_handles_increasing_response(tmp_path):
    operator = make_operator(tmp_path, lambda ppm: 1.0 + 0.00002 * ppm)
    operator.control.target_keff = 1.01

    operator((object(),), source_rate=1.0)

    assert operator.calls[-1] == (500.0, True)


@pytest.mark.parametrize(
    "response, expected_concentration, expected_status",
    [
        (lambda ppm: 1.02 - 0.000001 * ppm, 1000.0, "upper bound"),
        (lambda ppm: 1.00 - 0.000001 * ppm, 0.0, "lower bound"),
    ],
)
def test_unbracketed_target_selects_closest_endpoint(
    tmp_path, response, expected_concentration, expected_status
):
    operator = make_operator(tmp_path, response)

    operator((object(),), source_rate=1.0)

    assert operator.calls[-1] == (expected_concentration, True)
    with (tmp_path / "feedback.csv").open(newline="") as f:
        row = next(csv.DictReader(f))
    assert row["status"] == expected_status
