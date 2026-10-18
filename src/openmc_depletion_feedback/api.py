"""Public depletion API."""

from pathlib import Path
import csv
from types import SimpleNamespace

import openmc
import openmc.deplete as deplete


def _nominal_keff(result):
    value = getattr(result, "k", getattr(result, "keff", None))
    if value is None:
        raise AttributeError("OpenMC operator result does not expose k or keff")
    return float(getattr(value, "nominal_value", value))


def _keff_stats(result):
    """Return the keff estimate and its reported standard deviation, if any."""
    value = getattr(result, "k", getattr(result, "keff", None))
    if value is None:
        raise AttributeError("OpenMC operator result does not expose k or keff")
    mean = float(getattr(value, "nominal_value", value))
    sigma = getattr(value, "std_dev", None)
    return mean, None if sigma is None else float(sigma)


class _FeedbackOperator(deplete.CoupledOperator):
    def __init__(self, model, chain_file, control):
        super().__init__(model, chain_file=chain_file)
        self.control = control
        self.material = None
        self.base_densities = None
        self.original_stats = None
        self.history_path = None
        self._step = 0
        self._trial = 0
        self._previous_concentration = None
        self._previous_boc_keff = None
        self.expected_step_keff_loss = control.expected_step_keff_loss
        self.search_level = 0

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
            csv.writer(f).writerow([
                "step", "trial", "phase", "concentration", "keff",
                "keff_std_dev", "target_keff", "lower_bound", "upper_bound",
                "search_level", "particles", "batches", "selected",
                "keff_error", "status", "stopping_reason",
            ])
        return result

    def _set_concentration(self, concentration):
        c = self.control
        factor = concentration / c.reference_concentration
        names, dens = self.material.nuclides, self.material.densities.copy()
        for nuclide, base in self.base_densities.items():
            dens[names.index(nuclide)] = base * factor
        self.material.set_densities(names, dens)

    def _transport(self, vec, source_rate, concentration, final=False,
                   search_level=None):
        c, settings = self.control, openmc.lib.settings
        if not final:
            level = self.search_level if search_level is None else search_level
            factor = c.refinement_factor ** level
            base_particles, base_batches, base_inactive = self.original_stats
            settings.particles = round((c.search_particles or base_particles) * factor)
            settings.set_batches(round((c.search_batches or base_batches) * factor))
            if c.search_inactive is not None:
                settings.inactive = c.search_inactive
            else:
                settings.inactive = base_inactive
        else:
            settings.particles, batches, settings.inactive = self.original_stats
            settings.set_batches(batches)
        self._set_concentration(concentration)
        return super().__call__(vec, source_rate)

    def _record_trial(self, concentration, mean, sigma, target, lower, upper,
                      phase, selected=False, status="", stopping_reason=""):
        settings = getattr(openmc.lib, "settings", None)
        if settings is None:
            settings = SimpleNamespace()
        batches = settings.get_batches() if hasattr(settings, "get_batches") else ""
        with self.history_path.open("a", newline="") as f:
            csv.writer(f).writerow([
                self._step, self._trial, phase, f"{concentration:.8f}",
                f"{mean:.8f}", "" if sigma is None else f"{sigma:.8f}",
                f"{target:.8f}", f"{lower:.8f}", f"{upper:.8f}",
                self.search_level, getattr(settings, "particles", ""), batches,
                int(selected), f"{abs(mean - target):.8f}", status,
                stopping_reason,
            ])
        self._trial += 1

    def _evaluate(self, vec, source_rate, concentration, target, lower, upper,
                  phase="search", level=None):
        previous_level = self.search_level
        if level is not None:
            self.search_level = level
        result = self._transport(vec, source_rate, concentration,
                                 search_level=level)
        mean, sigma = _keff_stats(result)
        self._record_trial(concentration, mean, sigma, target, lower, upper,
                           phase)
        self.search_level = previous_level
        return mean, sigma

    def _confidence_met(self, mean, sigma, target):
        if sigma is None:
            return False
        return abs(mean - target) + self.control.confidence_z * sigma <= self.control.keff_tolerance

    def _target_for_call(self, vec, source_rate):
        c = self.control
        if c.target_mode != "end" or self._previous_concentration is None:
            return c.target_keff + (self.expected_step_keff_loss if c.target_mode == "end" else 0.0)

        # Predictor mode calls the operator once for the depleted state at the
        # start of the next step. Measure the previous step's loss at its
        # selected absorber concentration, then filter and bound the estimate.
        self.search_level = 0
        current_keff, _ = self._evaluate(
            vec, source_rate, self._previous_concentration,
            c.target_keff + self.expected_step_keff_loss,
            c.minimum, c.maximum, phase="loss_measurement", level=0,
        )
        measured_loss = self._previous_boc_keff - current_keff
        measured_loss = min(c.maximum_step_keff_loss,
                            max(c.minimum_step_keff_loss, measured_loss))
        alpha = c.loss_filter_alpha
        self.expected_step_keff_loss = (
            (1.0 - alpha) * self.expected_step_keff_loss + alpha * measured_loss
        )
        return c.target_keff + self.expected_step_keff_loss

    @staticmethod
    def _is_bracketed(klo, khi, target):
        return (klo - target) * (khi - target) <= 0.0

    def _warm_bracket(self, vec, source_rate, target):
        c = self.control
        if self._previous_concentration is None:
            lo, hi = c.minimum, c.maximum
        else:
            width = c.warm_start_width
            center = min(c.maximum, max(c.minimum, self._previous_concentration))
            lo = max(c.minimum, center - width)
            hi = min(c.maximum, center + width)

        while True:
            klo, _ = self._evaluate(vec, source_rate, lo, target, lo, hi,
                                    phase="bracket")
            khi, _ = self._evaluate(vec, source_rate, hi, target, lo, hi,
                                    phase="bracket")
            if self._is_bracketed(klo, khi, target):
                return lo, hi, klo, khi, True
            if lo == c.minimum and hi == c.maximum:
                return lo, hi, klo, khi, False
            width *= c.warm_start_expansion
            center = min(c.maximum, max(c.minimum, self._previous_concentration))
            new_lo = max(c.minimum, center - width)
            new_hi = min(c.maximum, center + width)
            if new_lo == lo and new_hi == hi:
                new_lo, new_hi = c.minimum, c.maximum
            lo, hi = new_lo, new_hi

    def _candidate(self, lo, hi, klo, khi, target):
        if self.control.use_secant and khi != klo:
            secant = lo + (target - klo) * (hi - lo) / (khi - klo)
            margin = 0.1 * (hi - lo)
            if lo + margin <= secant <= hi - margin:
                return secant
        return 0.5 * (lo + hi)

    def __call__(self, vec, source_rate):
        if source_rate == 0.0:
            return super().__call__(vec, source_rate)
        c = self.control
        self._trial = 0
        target = self._target_for_call(vec, source_rate)
        lo, hi, klo, khi, bracketed = self._warm_bracket(vec, source_rate, target)

        if not bracketed:
            selected = lo if abs(klo - target) < abs(khi - target) else hi
            selected_mean, selected_sigma = (klo, None) if selected == lo else (khi, None)
            stop_reason = "lower_bound" if selected == c.minimum else "upper_bound"
            status = "target_unmet"
        else:
            increasing = khi > klo
            selected = 0.5 * (lo + hi)
            selected_mean, selected_sigma = None, None
            status, stop_reason = "target_unmet", "max_iterations"
            for _ in range(c.max_iterations):
                candidate = self._candidate(lo, hi, klo, khi, target)
                mean, sigma = self._evaluate(vec, source_rate, candidate,
                                             target, lo, hi)
                selected, selected_mean, selected_sigma = candidate, mean, sigma

                # Near-target candidates are rerun with progressively more
                # histories until their confidence interval supports the
                # requested keff tolerance, or the refinement budget ends.
                if abs(mean - target) <= c.keff_tolerance:
                    for level in range(1, c.max_search_refinements + 1):
                        if self._confidence_met(mean, sigma, target):
                            break
                        mean, sigma = self._evaluate(
                            vec, source_rate, candidate, target, lo, hi,
                            phase="refinement", level=level,
                        )
                        selected_mean, selected_sigma = mean, sigma
                    if self._confidence_met(mean, sigma, target):
                        status, stop_reason = "target_met", "keff_tolerance"
                        break

                if (mean < target) == increasing:
                    lo, klo = candidate, mean
                else:
                    hi, khi = candidate, mean
                if hi - lo <= c.tolerance:
                    status = "target_met" if self._confidence_met(mean, sigma, target) else "target_unmet"
                    stop_reason = "concentration_tolerance"
                    if sigma is None or sigma * c.confidence_z > c.keff_tolerance:
                        stop_reason = "uncertainty_limit"
                    break
            else:
                # The loop exhausted its budget. Keep the last candidate and
                # distinguish unresolved Monte Carlo uncertainty explicitly.
                if selected_sigma is None or selected_sigma * c.confidence_z > c.keff_tolerance:
                    stop_reason = "uncertainty_limit"

        self.search_level = 0
        result = self._transport(vec, source_rate, selected, final=True)
        final_mean, final_sigma = _keff_stats(result)
        target_met = self._confidence_met(final_mean, final_sigma, target)
        if target_met:
            status = "target_met"
            if stop_reason in {"lower_bound", "upper_bound"}:
                stop_reason = "endpoint_within_tolerance"
        else:
            status = "target_unmet"
            if stop_reason not in {"lower_bound", "upper_bound"}:
                if final_sigma is None or final_sigma * c.confidence_z > c.keff_tolerance:
                    stop_reason = "uncertainty_limit"
                elif abs(final_mean - target) > c.keff_tolerance:
                    stop_reason = "final_keff_outside_tolerance"
        self._record_trial(selected, final_mean, final_sigma, target, lo, hi,
                           "selected_final", selected=True, status=status,
                           stopping_reason=stop_reason)
        print(
            f"[feedback] {status} ({stop_reason}): {selected:g}, "
            f"keff={final_mean:.6f}, target={target:.6f}"
        )
        self._previous_concentration = selected
        self._previous_boc_keff = final_mean
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
    integrators = {"predictor": deplete.PredictorIntegrator,
                   "cecm": deplete.CECMIntegrator,
                   "epc-rk4": deplete.EPCRK4Integrator}
    try:
        integrator_name = integrator.lower()
        integrator_cls = integrators[integrator_name]
    except (AttributeError, KeyError) as exc:
        raise ValueError(f"Unsupported integrator {integrator!r}; choose {', '.join(integrators)}") from exc
    if control is not None and control.target_mode == "end" and integrator_name != "predictor":
        raise ValueError("target_mode='end' currently requires integrator='predictor'")
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    operator_cls = _FeedbackOperator if control is not None else deplete.CoupledOperator
    if control is None:
        operator = operator_cls(model, chain_file=str(chain_file))
    else:
        operator = operator_cls(model, chain_file=str(chain_file), control=control)
    calculation = integrator_cls(operator, timesteps, power=power,
                                 timestep_units=timestep_units)
    calculation.integrate(output_dir=str(output_dir))
    return calculation

