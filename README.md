# openmc-depletion-feedback

An importable wrapper around OpenMC depletion. The package accepts a prepared
OpenMC `Model` or a directory containing OpenMC XML inputs and runs depletion
with one function. Provide the reactor geometry, materials, temperatures,
volumes, cross-section data, and depletion chain for your calculation.

## Install

Install OpenMC and its nuclear data according to the [OpenMC documentation](https://docs.openmc.org/), then install this package:

```bash
pip install .
```

## Documentation

The Sphinx documentation source is in `docs/` and uses reStructuredText
(`.rst`). To build the website locally:

```bash
python -m pip install -r docs/requirements.txt
sphinx-build -b html docs docs/_build/html
```

The generated site is in `docs/_build/html`. The GitHub Actions workflow in
`.github/workflows/docs.yml` builds and publishes it with GitHub Pages when
Pages is configured to use GitHub Actions.

Run the unit tests with:

```bash
python -m pip install -e ".[test]"
pytest
```

The feedback-search tests use deterministic keff values and do not launch
OpenMC transport simulations.

## Ordinary depletion

```python
import openmc
from openmc_depletion_feedback import run_depletion

run_depletion(
    model="model",  # directory containing geometry.xml, materials.xml, settings.xml, etc.
    timesteps=[1, 7, 30, 90],
    timestep_units="d",
    power=160e6,  # watts
    chain_file="data/chain.xml",
    output_dir="depletion-results",
)
```

## Soluble-species feedback

Provide a named material containing the soluble species at the stated
reference concentration. The controller scales the listed nuclide atom
densities together, leaving every other nuclide in the material untouched.
Species concentration units are intentionally user-defined, so the same
controller can represent ppm boron, ppm gadolinium, or another soluble
absorber when the reference material has been prepared consistently.

```python
from openmc_depletion_feedback import FeedbackControl, run_depletion

control = FeedbackControl(
    material_name="Primary coolant",
    nuclides=("B10", "B11"),
    reference_concentration=1000,  # ppm represented by the model densities
    minimum=0,
    maximum=1000,
    target_keff=1.005,
    search_particles=5000,
    search_batches=30,
    search_inactive=5,
)
run_depletion(
    model="model",  # directory containing the OpenMC XML inputs
    timesteps=[1, 7, 30],
    power=160e6,
    chain_file="data/chain.xml",
    control=control,
    output_dir="controlled-results",
)
```

Set `nuclides=("Xe135",)` for a single-nuclide absorber, or use an element
symbol such as `nuclides=("B",)` to scale all isotopes of that element already
present in the control material. Multiple symbols and isotope names can also
be mixed in the tuple. The target is the keff during the transport solve at
the start of each depletion step. Feedback uses transport calculations to
search the concentration; only the selected case is used for that depletion
step. The selected case is rerun using the model's original transport
statistics. Search statistics are optional; when omitted, model settings are
used throughout.

`FeedbackControl` accepts a concentration interval, interval tolerance,
keff tolerance, iteration limit, and `history_file` path. It writes the
selected concentration, keff, target, and bound/controlled status to CSV.
Concentration feedback assumes keff changes monotonically between the
configured endpoints. Validate that assumption and the atom-density versus
concentration relationship for each reactor and species.

## API

```python
run_depletion(model, timesteps, power, chain_file, *, timestep_units="d",
              integrator="predictor", control=None, output_dir=".")
```

The supported integrators are `predictor`, `cecm`, and `epc-rk4`. All
depletable materials in the model need positive volumes. The package does not
construct reactor geometry or choose nuclear data on the user's behalf.

