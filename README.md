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

## Citation

Citation metadata is provided in [`CITATION.cff`](CITATION.cff).
