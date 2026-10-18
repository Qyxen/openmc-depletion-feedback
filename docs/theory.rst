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
the lower and upper concentration bounds. If the target is bracketed, it
repeatedly tests the midpoint and narrows the interval. Each trial is a new
OpenMC transport solve at the same nuclide inventory, with only the selected
absorber densities changed. No trial depletion solve is performed. The
search stops when the calculated keff error is within ``keff_tolerance``, the
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

If both bounds are on the same side of the target, the controller selects the
endpoint whose calculated keff is closer and records ``lower bound`` or
``upper bound``. The target is then outside the reachable range for that
state and those bounds. For a bracketed search, the history status is
``controlled`` even if ``max_iterations`` is reached before either numerical
tolerance is met; inspect the selected keff and concentration to judge the
remaining error. The search assumes keff varies monotonically with
concentration over the configured interval. Check this assumption for each
model and absorber.

Trial searches can use reduced particle, batch, and inactive-generation
settings through ``search_particles``, ``search_batches``, and
``search_inactive``. Any setting left as ``None`` retains the model's current
value. Reduced histories lower the cost of each search transport solve, but
also increase Monte Carlo uncertainty in the estimated keff and selected
concentration. The final selected transport solve restores the particle,
batch, and inactive settings captured from the model when depletion starts.

``keff_tolerance`` and ``tolerance`` measure different things. The former is
an absolute difference in keff; the latter is a width in the user-defined
concentration units. Reaching either condition stops the search. Increasing
``max_iterations`` permits additional bisection steps but increases the
number of transport calculations. Search tolerances smaller than the
statistical uncertainty in trial keff do not ensure a correspondingly
accurate physical concentration.

Timestep target interpretation
------------------------------

With the predictor integrator, the feedback target applies to the transport
calculation at the beginning of each depletion step. It is therefore a
beginning-of-step target; the controller does not predict the end-of-step
keff. During the interval, fuel depletion, fission-product buildup, and
decay change reactivity. If an end-of-step condition matters, use shorter
steps or calibrate a beginning-of-step target against the expected loss for
the reactor, power, and schedule. Such calibration is external to this
package.

Output and reproducibility
--------------------------

OpenMC writes depletion results to ``output_dir``. With feedback enabled,
``history_file`` records the step number, selected concentration, resulting
keff, target keff, and whether the target was controlled or a bound was
selected. Archive the input model, cross-section and chain data identifiers,
power, timestep schedule, controller configuration, and transport settings
alongside the results when comparing calculations.

