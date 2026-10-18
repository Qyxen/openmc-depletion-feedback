API reference
=============

This page documents the public objects exported by
``openmc_depletion_feedback``. The package does not create reactor geometry,
materials, or nuclear data. Supply a prepared OpenMC model, depletion chain,
cross-section configuration, power, and time schedule.

``run_depletion``
-----------------

.. py:function:: run_depletion(model, timesteps, power, chain_file, *, timestep_units="d", integrator="predictor", control=None, output_dir=".")

Run an OpenMC transport-coupled depletion calculation and return the
integrator object after the calculation completes.

Parameters
~~~~~~~~~~

``model``
   Either an ``openmc.Model`` object or a path to a directory containing the
   OpenMC XML inputs. For a path, the function calls
   ``openmc.Model.from_xml(directory=...)``. Pass the directory, not a path to
   an individual XML file. Typical inputs include ``geometry.xml``,
   ``materials.xml``, and ``settings.xml``.

``timesteps``
   Sequence of positive step durations. Values are durations, not cumulative
   times. All values use ``timestep_units``.

``power``
   Reactor power in watts. A single value is applied to all steps by the
   selected OpenMC integrator.

``chain_file``
   Path to an existing OpenMC depletion chain XML file. The chain must
   provide the decay and transmutation data needed by the depletable
   materials and reactions in the model.

``timestep_units``
   Time unit for the values in ``timesteps``. Defaults to ``"d"`` (days).
   Units supported by OpenMC include seconds, minutes, hours, and days.

``integrator``
   Depletion time-integration method. Supported values are:

   * ``"predictor"``: first-order predictor method; default.
   * ``"cecm"``: constant extrapolation/constant midpoint method.
   * ``"epc-rk4"``: extended predictor-corrector Runge-Kutta method.

   End-of-step feedback compensation (``control.target_mode="end"``) is
   currently supported only with ``"predictor"``.

``control``
   Optional ``FeedbackControl`` object. If omitted, the depletion calculation
   runs without soluble-species feedback.

``output_dir``
   Directory for OpenMC depletion output and feedback history. Created if it
   does not already exist. Defaults to the current directory.

Model requirements
~~~~~~~~~~~~~~~~~~

Every material marked depletable must have a positive volume. The
model's cross-section configuration must be available to OpenMC. ``chain_file``
must exist. The caller is responsible for choosing consistent units,
materials, volumes, power, and nuclear data.

The function raises:

* ``FileNotFoundError`` if the depletion chain path does not exist.
* ``ValueError`` if a depletable material has no positive volume, the
  integrator name is unsupported, or end-of-step feedback is requested with
  an integrator other than ``predictor``.
* OpenMC exceptions for invalid models, missing data, transport failures, or
  depletion failures.

``FeedbackControl``
-------------------

.. py:class:: FeedbackControl(material_name, nuclides, reference_concentration, minimum, maximum, target_keff, ...)

Configure bounded feedback that scales selected nuclide atom densities in
one material to move transport ``k_eff`` toward a target. The material must
exist in the model, and each requested nuclide or element must match one or
more nuclides in that material.

Required parameters
~~~~~~~~~~~~~~~~~~~

``material_name``
   Exact material name used to find the control material in the OpenMC model.

``nuclides``
   Tuple of isotope names or element symbols. An isotope name such as
   ``"B10"`` selects that isotope. An element symbol such as ``"B"`` selects
   all isotopes of that element already in the material. Multiple isotope
   names and element symbols can be combined.

``reference_concentration``
   Concentration represented by the selected nuclide densities currently
   defined in the input material. The density at a requested concentration
   is the reference density multiplied by
   ``requested_concentration / reference_concentration``. Concentration units
   are user-defined but must be consistent across all concentration fields.

``minimum`` and ``maximum``
   Lower and upper search bounds. The current implementation requires
   ``0 <= minimum < maximum <= reference_concentration``.

``target_keff``
   Desired transport eigenvalue. Its interpretation depends on
   ``target_mode`` below.

Feedback and target parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``target_mode="beginning"``
   Default. Treat ``target_keff`` as the target at the start of each
   depletion step. The controller searches for the concentration at each
   operator evaluation.

``target_mode="end"``
   Treat ``target_keff`` as the desired end-of-step value. The beginning-of-
   step target is adjusted using a filtered estimate of the expected keff
   loss. This mode requires the predictor integrator.

``expected_step_keff_loss=0.01``
   Initial estimate of the keff decrease over one timestep. Used only when
   ``target_mode="end"``. Calibrate this value for the reactor, power, and
   timestep schedule.

``loss_filter_alpha=0.3``
   Weight applied to each newly measured step loss in the exponential filter.
   Must be between zero and one. A higher value follows new measurements
   faster; a lower value smooths them more strongly.

``minimum_step_keff_loss=-0.01`` and ``maximum_step_keff_loss=0.02``
   Bounds applied to a measured step loss before it updates the filtered
   estimate. Used only in end-target mode.

Search parameters
~~~~~~~~~~~~~~~~~

``search_particles=None``
   Particle count for concentration trial solves at refinement level zero.
   ``None`` uses the particle count in the OpenMC model. Refinement levels
   multiply the current base particle count by ``refinement_factor``.

``search_batches=None``
   Batch count for trial solves at refinement level zero. ``None`` uses the
   model's batch count. Refinement levels multiply the base batch count by
   ``refinement_factor``.

``search_inactive=None``
   Inactive batches for search solves. ``None`` uses the model setting. The
   number of inactive batches is not multiplied during refinement.

``tolerance=0.001``
   Concentration interval width at which the bounded search may stop, in the
   same user-defined units as ``minimum`` and ``maximum``.

``keff_tolerance=0.0005``
   Maximum absolute keff error. A target is reported as met only when
   ``abs(keff - target_keff) + confidence_z * keff_std_dev`` is no greater
   than this tolerance.

``max_iterations=12``
   Maximum number of secant/bisection search iterations per operator call.

``max_search_refinements=2``
   Maximum number of higher-statistics reruns for a candidate whose nominal
   keff is near the target but whose reported uncertainty is too large.

``refinement_factor=2.0``
   Multiplier applied to trial particle and batch counts at each refinement
   level. Must be greater than one.

``confidence_z=2.0``
   Number of reported standard deviations included in the conservative keff
   acceptance criterion. For example, ``2.0`` uses an approximately two-
   standard-deviation margin when the estimate is approximately normal.

``warm_start_width=100.0``
   Half-width of the initial concentration interval around the previous
   selected concentration, in the configured concentration units. The first
   search uses the full configured interval.

``warm_start_expansion=2.0``
   Factor by which the local interval half-width grows when it does not
   bracket the target. Expansion continues up to the full configured range.

``use_secant=True``
   Use a secant estimate when it lies safely inside the current bracket. Set
   to ``False`` to use midpoint bisection for every iteration.

``history_file="feedback_history.csv"``
   CSV filename for feedback diagnostics. Relative paths are placed under
   ``output_dir``; absolute paths are used as supplied.

Search behavior and outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~

On the first feedback call, the full concentration interval is evaluated.
Later calls start near the previous selected concentration and expand the
bracket if necessary. Within a bracket, the search uses safeguarded secant
proposals with bisection fallback. If the target is outside the configured
range, the closer endpoint is used.

The selected concentration is evaluated again using the model's original
particle, batch, and inactive settings. Only this final transport result is
passed to the depletion integrator. Trial concentrations do not run trial
depletion calculations.

The history CSV contains one row for every endpoint, search, refinement,
loss-measurement, and final selected transport evaluation. Columns include
step and trial numbers, phase, concentration, keff estimate and reported
standard deviation, target, current bracket, search level and statistics,
selection flag, absolute keff error, status, and stopping reason. ``target_met``
means the final selected result met the uncertainty-aware keff criterion;
``target_unmet`` means it did not. Bound, concentration-width, iteration, and
uncertainty limits are distinguished in ``stopping_reason``.

The controller assumes the keff response is monotonic over the configured
concentration range. Monte Carlo uncertainties describe sampling variation;
they do not include nuclear-data uncertainty or model bias.

