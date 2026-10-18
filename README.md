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

The full user guide can be found here: [openmc-depletion-feedback documentation](https://qyxen.github.io/openmc-depletion-feedback/)

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
be mixed in the tuple. With the default `target_mode="beginning"`,
`target_keff` is the keff target at the start of each depletion step. Feedback
uses transport calculations to search the concentration; only the selected
case is used for that depletion step. The selected case is rerun using the
model's original transport statistics. Search statistics are optional; when
omitted, model settings are used for initial trials and increased for
near-target refinements.

`FeedbackControl` starts each new depletion step near the preceding step's
selected concentration, expands the bracket if needed, then uses a safeguarded
secant proposal with bisection fallback. Near-target trials can be rerun with
increasing particle and batch counts. The CSV records every trial's keff and
reported standard deviation, the selected result, remaining error, and
stopping reason. The search assumes keff changes monotonically with
concentration across the configured interval; validate that assumption and
the atom-density versus concentration relationship for each reactor and
species.

By default, `target_keff` applies at the beginning of each step. Set
`target_mode="end"` to estimate the beginning-of-step target needed to reach
the requested end-of-step keff. This mode filters measured step losses and
requires `integrator="predictor"`. Configure `expected_step_keff_loss` from a
pilot calculation and tune its clipping range and `loss_filter_alpha` for the
reactor and timestep schedule.

See the [detailed API reference](https://qyxen.github.io/openmc-depletion-feedback/api.html)
for all function parameters, feedback settings, validation rules, and outputs.

