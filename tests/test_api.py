from pathlib import Path

import openmc
import openmc.deplete as deplete

from openmc_depletion_feedback.api import run_depletion


def test_run_depletion_passes_results_path_and_sets_operator_output_dir(
    tmp_path, monkeypatch
):
    chain_file = tmp_path / "chain.xml"
    chain_file.write_text("<chain />")
    output_dir = tmp_path / "results"
    captured = {}

    class FakeOperator:
        def __init__(self, model, chain_file):
            self.model = model
            self.chain_file = chain_file

    class FakeIntegrator:
        def __init__(self, operator, timesteps, *, power, timestep_units):
            captured["operator"] = operator

        def integrate(self, *, path):
            captured["path"] = path

    monkeypatch.setattr(deplete, "CoupledOperator", FakeOperator, raising=False)
    monkeypatch.setattr(deplete, "PredictorIntegrator", FakeIntegrator, raising=False)
    monkeypatch.setattr(deplete, "CECMIntegrator", FakeIntegrator, raising=False)
    monkeypatch.setattr(deplete, "EPCRK4Integrator", FakeIntegrator, raising=False)

    model = openmc.Model()
    model.materials = []
    run_depletion(
        model=model,
        timesteps=[1.0],
        power=1.0,
        chain_file=chain_file,
        output_dir=output_dir,
    )

    assert captured["operator"].output_dir == output_dir.resolve()
    assert Path(captured["path"]) == output_dir.resolve() / "depletion_results.h5"
