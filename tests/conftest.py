"""Supply an import-only OpenMC stub when testing without OpenMC installed."""

import importlib.util
import sys
from types import ModuleType, SimpleNamespace


if importlib.util.find_spec("openmc") is None:
    fake_openmc = ModuleType("openmc")
    fake_openmc.__path__ = []
    fake_openmc.Model = type("Model", (), {})
    fake_openmc.lib = SimpleNamespace()

    fake_deplete = ModuleType("openmc.deplete")

    class FakeCoupledOperator:
        def __init__(self, *args, **kwargs):
            pass

    fake_deplete.CoupledOperator = FakeCoupledOperator
    fake_openmc.deplete = fake_deplete
    sys.modules["openmc"] = fake_openmc
    sys.modules["openmc.deplete"] = fake_deplete
