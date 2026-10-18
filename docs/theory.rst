Theory and numerical method
===========================

Coupled transport and depletion
--------------------------------

During irradiation, neutron reactions change the number densities of nuclides
in each depletable material. A depletion calculation represents those changes
with a system of Bateman equations:

.. math::

   \frac{d\mathbf{N}}{dt} = \mathbf{A}(\boldsymbol{\phi})\,\mathbf{N},

where :math:`\mathbf{N}` is the vector of nuclide number densities,
:math:`\boldsymbol{\phi}` is the neutron flux, and :math:`\mathbf{A}` is the
depletion matrix. Its diagonal and off-diagonal terms represent loss and
production through radioactive decay and neutron-induced reactions. For a
reaction :math:`r` with microscopic cross section :math:`\sigma_{i,r}(E)`,
the reaction rate per target atom is determined by the spectrum-weighted
quantity

.. math::

   r_{i,r} = \int_0^\infty \sigma_{i,r}(E)\,\phi(E)\,dE.

This is the rate per target atom. Multiplying by the number of target atoms
gives the material reaction rate. The depletion matrix uses these rates to
represent the gain and loss of each nuclide. For example, neutron capture
removes a target nuclide and creates a heavier product, while radioactive
decay removes a parent and feeds one or more daughter nuclides. Fission
removes a fissile nuclide and creates fission products according to the
chain's fission yields. The rate depends on the nuclide, reaction, material,
and neutron spectrum. The chain file supplies decay and transmutation data;
cross-section data supplies transport reaction data. Both must cover the
nuclides and reactions relevant to the calculation.

The nuclide inventory over a fixed interval can be written formally as

.. math::

   \mathbf{N}(t+\Delta t) = \exp\!\left(\mathbf{A}\Delta t\right)\,
   \mathbf{N}(t),

if :math:`\mathbf{A}` is constant during the interval. In practice, OpenMC's
depletion integrator approximates the changing reaction rates and solves the
resulting stiff system numerically. The time step is therefore part of the
physical approximation: a long interval assumes the chosen integrator's rate
approximation remains adequate as composition and spectrum change.

``run_depletion`` uses OpenMC's ``CoupledOperator`` and one of the selected
integrators: ``predictor``, ``cecm``, or ``epc-rk4``. In predictor mode,
transport reaction rates calculated from the beginning-of-step composition
are used to advance the inventory through that step. Higher-order integrators
may request reaction rates at additional stages in the step. The feedback
controller is called whenever the integrator requests rates, so for feedback
calculations the predictor integrator gives the clearest interpretation of
one concentration selection per depletion step.

The power argument is in watts. In power-normalized depletion, OpenMC scales
the neutron source/reaction rates so that the model represents the requested
power. This requires consistent model volumes, material compositions, fission
energies, and power assumptions. Only materials marked depletable in the
input model are included in depletion, and they must have positive volumes.
Other materials can still affect neutron transport but their inventories are
not advanced by depletion.

The available integration choices differ in how they estimate reaction rates
over each interval:

* ``predictor`` is a first-order method. It uses the depletion matrix evaluated
  from the beginning-of-step state for the interval.
* ``cecm`` is the second-order constant-extrapolation/constant-midpoint
  predictor-corrector method. It evaluates transport at the beginning and at
  an estimated midpoint state.
* ``epc-rk4`` is an extended predictor-corrector method based on the classical
  fourth-order Runge-Kutta stages, with transport evaluations at intermediate
  estimated states.

Higher-order methods can improve the approximation when reaction rates vary
over a step, at the cost of additional transport solves. The actual accuracy
also depends on timestep size, Monte Carlo uncertainty, nuclear data, and
model fidelity; the method name alone does not guarantee a particular
accuracy for a reactor calculation.

The transport calculation and depletion solve have different roles. Transport
estimates neutron reaction rates for the current composition; the depletion
solver uses these rates and the chain data to update nuclide inventories.
The updated inventory changes subsequent transport calculations. This
sequence is the transport-depletion coupling. The method is not a continuous
feedback simulation within each transport history or within an infinitesimal
time interval; it advances the state at discrete integration stages.

Effective multiplication factor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The transport eigenvalue :math:`k_{\mathrm{eff}}` is the ratio of neutron
production in one generation to neutron losses by absorption and leakage,
under the eigenvalue transport model. A value near one indicates a critical
steady neutron population in that model. A value above or below one indicates
that neutron production is greater or less than losses, respectively. This
eigenvalue describes a transport calculation for a particular material state
and set of transport conditions; it does not directly measure the depletion
rate or how many operating days remain.

Adding a soluble absorber generally increases neutron absorption and lowers
:math:`k_{\mathrm{eff}}`; removing absorber generally has the opposite
effect. This provides the feedback mechanism: as fuel composition changes,
the controller changes absorber atom densities to keep the calculated
eigenvalue near the requested target. The relation depends on the reactor
state and neutron spectrum because the absorber competes with other
absorbers, fuel, and structural materials for neutrons.

Choosing a time-step schedule
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use shorter steps where composition or reactivity changes rapidly, such as
early irradiation when short-lived fission products build up, after a power
change, or when a soluble absorber changes substantially. Longer steps may
be appropriate after rates of change become slower. A timestep refinement
study compares results from a schedule with smaller intervals against the
intended schedule, using quantities such as nuclide inventories, keff, and
absorber concentration. Agreement within the accuracy needed for the study is
evidence that the schedule is sufficiently resolved; a single calculation
does not establish timestep adequacy by itself.

Soluble-species feedback
------------------------

Feedback is optional. When ``FeedbackControl`` is supplied, the package
adjusts the specified absorber concentration when the coupled operator is
asked to calculate transport reaction rates. With the predictor integrator,
this occurs at the beginning of each depletion step. The concentration is
chosen to bring the calculated transport eigenvalue :math:`k_{\mathrm{eff}}`
near ``target_keff``. Without feedback, depletion runs with the input model's
material compositions.

The controller changes the atom densities of the selected nuclides in one
named material. If the material's selected nuclides represent a reference
concentration :math:`C_{\mathrm{ref}}`, the density applied at concentration
:math:`C` is

.. math::

   N_i(C) = N_i(C_{\mathrm{ref}})\,\frac{C}{C_{\mathrm{ref}}},

for each selected nuclide :math:`i`. All unselected nuclide densities in that
material remain as they were. An element symbol, such as ``"B"``, selects all
isotopes of that element present in the material; explicit nuclide names can
be used to select only certain isotopes.

This scaling means the reference concentration must match the composition
already in the model. The controller does not convert mass concentration to
atom density, perform chemical speciation, or add/remove solvent to preserve
total material density. Prepare the reference material and choose the
concentration units consistently before running.

For example, if the input material contains B10 and B11 at atom densities
corresponding to 1000 ppm total boron, set ``reference_concentration=1000``.
At a trial concentration of 500 ppm, the controller halves the B10 and B11
densities individually, preserving their isotope ratio. It does not infer an
isotopic ratio from natural abundance or recalculate the coolant mass density.
The material's other nuclides, such as hydrogen and oxygen, are not rescaled.

Bounded search
--------------

At each feedback call, the controller calculates :math:`k_{\mathrm{eff}}` at
On the first call, it searches the full concentration range. On later calls,
it starts with an interval around the previous selected concentration. If
that interval does not bracket the target, it expands geometrically until it
does or reaches the configured full range. This warm start often reduces
transport solves when adjacent depletion steps need similar concentrations;
it still falls back to the full range if the target has moved farther.

Each trial is a new OpenMC transport solve at the same nuclide inventory,
with only the selected absorber densities changed. No trial depletion solve
is performed. Once a bracket is available, the search proposes a safeguarded
secant estimate when it lies safely inside the bracket, and otherwise uses
the midpoint. The secant estimate uses a straight-line approximation between
the endpoint results; the bracket and midpoint fallback preserve the
robustness of bisection when that approximation is poor.

The search stops when its uncertainty-aware keff criterion is met, the
concentration interval is within ``tolerance``, or ``max_iterations`` is
reached. It then runs the selected concentration again using the model's
original transport statistics; those reaction rates are returned to the
depletion integrator.

For a monotonically decreasing relation between concentration and keff, the
lower concentration endpoint has the higher keff. The implementation also
checks the endpoint direction from the results, so a monotonically increasing
relationship can be bracketed as well. It does not check for non-monotonic
behavior inside the interval. Monte Carlo uncertainty can also make measured
trial values appear locally non-monotonic, especially when search runs use
few histories.

If both full-range bounds are on the same side of the target, the controller
selects the endpoint whose calculated keff is closer. It records a
``lower_bound`` or ``upper_bound`` stopping reason, and the target is not
reachable within that range for the current state. For a bracketed search,
history distinguishes ``target_met`` from ``target_unmet`` and records the
stopping reason, including iteration limit and uncertainty limit. The
selected result's remaining absolute keff error is included in the history.
The search assumes keff varies monotonically with concentration over the
configured interval. Check this assumption for each model and absorber.

Trial searches can use reduced particle, batch, and inactive-generation
settings through ``search_particles``, ``search_batches``, and
``search_inactive``. Any setting left as ``None`` retains the model's current
value. Reduced histories lower the cost of each search transport solve, but
also increase Monte Carlo uncertainty in the estimated keff and selected
concentration. The final selected transport solve restores the particle,
batch, and inactive settings captured from the model when depletion starts.

When a trial's nominal keff is within ``keff_tolerance`` but its uncertainty
is too large to establish the target, the controller repeats that candidate
with progressively higher particle and batch counts. The counts are
multiplied by ``refinement_factor`` for each level, up to
``max_search_refinements``. A target is reported as met only when

.. math::

   |\hat{k}_{\mathrm{eff}} - k_{\mathrm{target}}|
   + z\,\sigma_k \leq \varepsilon_k,

where :math:`\hat{k}_{\mathrm{eff}}` is the estimated keff,
:math:`\sigma_k` is its reported standard deviation, :math:`z` is
``confidence_z``, and :math:`\varepsilon_k` is ``keff_tolerance``. This is a
conservative uncertainty check under the usual approximately normal
interpretation of Monte Carlo estimates. Increasing histories reduces
sampling uncertainty but does not address nuclear-data or model uncertainty.

``keff_tolerance`` and ``tolerance`` measure different things. The former is
an absolute difference in keff; the latter is a width in the user-defined
concentration units. ``warm_start_width`` controls the initial local bracket
around the prior concentration, and ``warm_start_expansion`` controls its
geometric growth. ``use_secant`` enables the safeguarded secant proposal;
setting it to ``False`` uses bisection throughout. Increasing
``max_iterations`` or ``max_search_refinements`` permits additional transport
calculations. Search tolerances smaller than the statistical uncertainty in
trial keff do not ensure a correspondingly accurate physical concentration.

The history CSV contains one row for each bracket, search, refinement, and
final selected transport solve. It records concentration, estimated keff,
reported standard deviation, current bracket, search settings, target error,
whether that solve supplied the selected result, and final stopping reason.

Timestep target interpretation
------------------------------

With ``target_mode="beginning"`` (the default), ``target_keff`` applies to
the transport calculation at the beginning of each depletion step. It is a
beginning-of-step target; the controller does not predict end-of-step keff.

With ``target_mode="end"``, ``target_keff`` is interpreted as a desired
end-of-step keff. The controller adds an estimate of the keff loss over the
upcoming interval to form its beginning-of-step target. After each predictor
step, it evaluates the new composition using the prior step's absorber
concentration, measures the keff loss, clips that measurement to the
configured range, and updates the estimate with an exponential filter:

.. math::

   L_{n+1} = (1-\alpha)L_n + \alpha\,\Delta k_n.

Here :math:`L_n` is the filtered expected step loss, :math:`\Delta k_n` is
the measured loss, and :math:`\alpha` is ``loss_filter_alpha``. The initial
estimate is ``expected_step_keff_loss``; the clipping interval is set by
``minimum_step_keff_loss`` and ``maximum_step_keff_loss``. This compensation
currently requires the predictor integrator because it assumes one operator
feedback call per depletion step. It remains an estimate: power changes,
timestep changes, and strong spectrum or composition effects can make past
loss a poor predictor of the next step.

Output and reproducibility
--------------------------

OpenMC writes depletion results to ``output_dir``. With feedback enabled,
``history_file`` records every evaluated concentration and keff estimate,
its reported uncertainty, the final selection and stopping reason. Archive
the input model, cross-section and chain data identifiers,
power, timestep schedule, controller configuration, and transport settings
alongside the results when comparing calculations.

