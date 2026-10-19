import pytest

from openmc_feedbacks import FeedbackControl


def make_control(**overrides):
    values = {
        "material_name": "Coolant",
        "nuclides": ("B",),
        "reference_concentration": 1000.0,
        "minimum": 0.0,
        "maximum": 1000.0,
        "target_keff": 1.0,
    }
    values.update(overrides)
    return FeedbackControl(**values)


def test_accepts_valid_control():
    control = make_control()

    assert control.material_name == "Coolant"
    assert control.nuclides == ("B",)
    assert control.minimum == 0.0
    assert control.maximum == 1000.0


@pytest.mark.parametrize(
    "overrides, message",
    [
        ({"material_name": ""}, "material_name"),
        ({"nuclides": ()}, "nuclide"),
        ({"reference_concentration": 0.0}, "reference_concentration"),
        ({"minimum": -1.0}, "Require 0 <= minimum"),
        ({"maximum": 1001.0}, "reference_concentration"),
        ({"tolerance": 0.0}, "tolerance and max_iterations"),
        ({"max_iterations": 0}, "tolerance and max_iterations"),
    ],
)
def test_rejects_invalid_control(overrides, message):
    with pytest.raises(ValueError, match=message):
        make_control(**overrides)
