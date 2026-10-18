"""Generic reactivity feedback using a soluble species in one material."""

from dataclasses import dataclass
from pathlib import Path
import csv


@dataclass
class FeedbackControl:
    """Configuration for bounded keff feedback.

    ``reference_concentration`` is the concentration represented by the
    isotope densities currently in ``material_name``. The controller scales
    only ``nuclides`` linearly from zero to that reference concentration.
    Concentration units are user-defined (ppm, g/kg, or another consistent
    unit); ``minimum`` and ``maximum`` use the same unit.
    """

    material_name: str
    nuclides: tuple[str, ...]
    reference_concentration: float
    minimum: float
    maximum: float
    target_keff: float
    target_mode: str = "beginning"
    expected_step_keff_loss: float = 0.01
    loss_filter_alpha: float = 0.3
    minimum_step_keff_loss: float = -0.01
    maximum_step_keff_loss: float = 0.02
    search_particles: int | None = None
    search_batches: int | None = None
    search_inactive: int | None = None
    max_search_refinements: int = 2
    refinement_factor: float = 2.0
    confidence_z: float = 2.0
    warm_start_width: float = 100.0
    warm_start_expansion: float = 2.0
    use_secant: bool = True
    tolerance: float = 0.001
    keff_tolerance: float = 0.0005
    max_iterations: int = 12
    history_file: str | Path = "feedback_history.csv"

    def __post_init__(self):
        if not self.material_name or not self.nuclides:
            raise ValueError("material_name and at least one nuclide are required")
        if self.reference_concentration <= 0:
            raise ValueError("reference_concentration must be positive")
        if not (0 <= self.minimum < self.maximum <= self.reference_concentration):
            raise ValueError("Require 0 <= minimum < maximum <= reference_concentration")
        if self.tolerance <= 0 or self.max_iterations < 1:
            raise ValueError("tolerance and max_iterations must be positive")
        if self.keff_tolerance <= 0:
            raise ValueError("keff_tolerance must be positive")
        for name in ("search_particles", "search_batches"):
            value = getattr(self, name)
            if value is not None and value < 1:
                raise ValueError(f"{name} must be a positive integer")
        if self.search_inactive is not None and self.search_inactive < 0:
            raise ValueError("search_inactive must be nonnegative")
        if self.target_mode not in {"beginning", "end"}:
            raise ValueError("target_mode must be 'beginning' or 'end'")
        if not 0.0 <= self.loss_filter_alpha <= 1.0:
            raise ValueError("loss_filter_alpha must be between 0 and 1")
        if self.minimum_step_keff_loss > self.maximum_step_keff_loss:
            raise ValueError("minimum_step_keff_loss must not exceed maximum_step_keff_loss")
        if self.max_search_refinements < 0:
            raise ValueError("max_search_refinements must be nonnegative")
        if self.refinement_factor <= 1.0:
            raise ValueError("refinement_factor must be greater than 1")
        if self.confidence_z <= 0:
            raise ValueError("confidence_z must be positive")
        if self.warm_start_width <= 0 or self.warm_start_expansion <= 1.0:
            raise ValueError("warm-start width must be positive and expansion must exceed 1")


