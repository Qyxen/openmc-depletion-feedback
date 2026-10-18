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
    search_particles: int | None = None
    search_batches: int | None = None
    search_inactive: int | None = None
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


