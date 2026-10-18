"""Public depletion API."""

from pathlib import Path
import csv

import openmc
import openmc.deplete as deplete


def _nominal_keff(result):
    value = getattr(result, "k", getattr(result, "keff", None))
    if value is None:
        raise AttributeError("OpenMC operator result does not expose k or keff")
    return float(getattr(value, "nominal_value", value))


class _FeedbackOperator(deplete.CoupledOperator):
    def __init__(self, model, chain_file, control):
        super().__init__(model, chain_file=chain_file)
        self.control = control
        self.material = None
        self.base_densities = None
        self.original_stats = None
        self.history_path = None
        self._step = 0

    def initial_condition(self):
        result = super().initial_condition()
        c = self.control
        source = next((m for m in self.model.materials if m.name == c.material_name), None)
        if source is None:
            raise ValueError(f"Control material {c.material_name!r} was not found")
        self.material = openmc.lib.materials[source.id]
        names = self.material.nuclides
        dens = self.material.densities
        resolved = []
        missing = []
        for requested in c.nuclides:
            matches = [n for n in names if n == requested or
                       (requested.isalpha() and n.startswith(requested) and
                        n[len(requested):].isdigit())]
            if not matches:
                missing.append(requested)
            else:
                resolved.extend(matches)
        if missing:
            raise ValueError(f"Control material lacks nuclides: {', '.join(missing)}")
        self.control_nuclides = tuple(dict.fromkeys(resolved))
        self.base_densities = {n: float(dens[names.index(n)]) for n in self.control_nuclides}
        settings = openmc.lib.settings
        self.original_stats = (settings.particles, settings.get_batches(), settings.inactive)
        self.history_path = Path(self.output_dir) / c.history_file
        self.history_path.parent.mkdir(parents=True, exist_ok=True)
        with self.history_path.open("w", newline="") as f:
            csv.writer(f).writerow(["step", "concentration", "keff", "target_keff", "status"])
        return result

    def _set_concentration(self, concentration):
        c = self.control
        factor = concentration / c.reference_concentration
        names, dens = self.material.nuclides, self.material.densities.copy()
        for nuclide, base in self.base_densities.items():
            dens[names.index(nuclide)] = base * factor
        self.material.set_densities(names, dens)

    def _transport(self, vec, source_rate, concentration, final=False):
        c, settings = self.control, openmc.lib.settings
        if not final:
            if c.search_particles is not None:
                settings.particles = c.search_particles
            if c.search_batches is not None:
                settings.set_batches(c.search_batches)
            if c.search_inactive is not None:
                settings.inactive = c.search_inactive
        else:
            settings.particles, batches, settings.inactive = self.original_stats
            settings.set_batches(batches)
        self._set_concentration(concentration)
        return super().__call__(vec, source_rate)

    def __call__(self, vec, source_rate):
        if source_rate == 0.0:
            return super().__call__(vec, source_rate)
        c = self.control
        lo, hi = c.minimum, c.maximum
        klo = _nominal_keff(self._transport(vec, source_rate, lo))
        khi = _nominal_keff(self._transport(vec, source_rate, hi))
        # Infer monotonic direction from the bracket endpoints.
        if abs(klo - c.target_keff) <= c.keff_tolerance:
            selected, status = lo, "controlled"
        elif abs(khi - c.target_keff) <= c.keff_tolerance:
            selected, status = hi, "controlled"
        elif (klo - c.target_keff) * (khi - c.target_keff) > 0:
            selected = lo if abs(klo - c.target_keff) < abs(khi - c.target_keff) else hi
            status = "lower bound" if selected == lo else "upper bound"
        else:
            status, selected = "controlled", (lo + hi) / 2
            increasing = khi > klo
            for _ in range(c.max_iterations):
                mid = (lo + hi) / 2
                kmid = _nominal_keff(self._transport(vec, source_rate, mid))
                selected = mid
                if abs(kmid - c.target_keff) <= c.keff_tolerance or hi - lo <= c.tolerance:
                    break
                if (kmid < c.target_keff) == increasing:
                    lo = mid
                else:
                    hi = mid
        result = self._transport(vec, source_rate, selected, final=True)
        kval = _nominal_keff(result)
        with self.history_path.open("a", newline="") as f:
            csv.writer(f).writerow([self._step, selected, kval, c.target_keff, status])
        print(f"[feedback] {status}: {selected:g}, keff={kval:.6f}")
        self._step += 1
        return result

def run_depletion(model, timesteps, power, chain_file, *, timestep_units="d",
                  integrator="predictor", control=None, output_dir="."):
    """Run an OpenMC depletion calculation.

    Parameters
    ----------
    model : openmc.Model or path-like
        A prepared OpenMC model (or directory containing OpenMC XML inputs).
    timesteps : sequence
        Step durations, interpreted using ``timestep_units``.
    power : float
        Reactor power in watts.
    chain_file : path-like
        OpenMC depletion chain XML file.
    control : FeedbackControl, optional
        Enable feedback by scaling configured nuclide densities in a named
        material. Omit for an ordinary depletion run.
    output_dir : path-like
        Directory for depletion results and feedback history.

    Returns the configured OpenMC integrator after it completes.
    """
    chain_file = Path(chain_file)
    if not chain_file.is_file():
        raise FileNotFoundError(f"Depletion chain not found: {chain_file}")
    if not isinstance(model, openmc.Model):
        model = openmc.Model.from_xml(directory=str(model))
    missing = [m.name for m in model.materials if m.depletable and (m.volume is None or m.volume <= 0)]
    if missing:
        raise ValueError("Depletable materials need positive volumes: " + ", ".join(missing))
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    operator_cls = _FeedbackOperator if control is not None else deplete.CoupledOperator
    if control is None:
        operator = operator_cls(model, chain_file=str(chain_file))
    else:
        operator = operator_cls(model, chain_file=str(chain_file), control=control)
    integrators = {"predictor": deplete.PredictorIntegrator,
                   "cecm": deplete.CECMIntegrator,
                   "epc-rk4": deplete.EPCRK4Integrator}
    try:
        integrator_cls = integrators[integrator.lower()]
    except (AttributeError, KeyError) as exc:
        raise ValueError(f"Unsupported integrator {integrator!r}; choose {', '.join(integrators)}") from exc
    calculation = integrator_cls(operator, timesteps, power=power,
                                 timestep_units=timestep_units)
    calculation.integrate(output_dir=str(output_dir))
    return calculation

