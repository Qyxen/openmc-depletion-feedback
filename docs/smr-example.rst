SMR depletion example
=====================

Use this example to run depletion for an SMR model with soluble-boron
feedback. You provide the OpenMC model, nuclear data paths, depletion chain,
reactor power, and timestep schedule. The same procedure applies to other
reactors and soluble absorbers.

1. Prepare the SMR OpenMC model
-------------------------------

Prepare the SMR as an ``openmc.Model`` object or export its inputs as OpenMC
XML files. Before depletion, confirm that:

#. Fuel materials intended for depletion have ``depletable = True``.
#. Every depletable material has a positive volume.
#. The coolant material that will be controlled has a stable, unique name.
#. The controlled coolant contains the isotopes of the soluble species at a
   known reference concentration.
#. OpenMC's cross-section library is configured and includes the nuclides in
   the model and depletion chain.

In this example, the borated coolant material is named
``Coolant with 1000 ppm boron``. Its B10 and B11 atom densities represent
1000 ppm boron. Set the cross-section and depletion-chain paths to files
available in your OpenMC environment.

Generate the OpenMC XML inputs for your SMR, then load the directory
containing those files with ``openmc.Model.from_xml``. Pass the directory, not
an individual XML filename. For example, if the files are in ``smr-xml/``:

.. code-block:: python

   import openmc

   model = openmc.Model.from_xml(directory="path/to/smr-xml")

That directory should contain the OpenMC inputs, typically ``geometry.xml``,
``materials.xml``, ``settings.xml``, and (if used) ``tallies.xml``. You can
pass either the loaded ``model`` object or the directory itself to
``run_depletion``. A direct filename such as
``path/to/smr-xml/model.xml`` is not accepted by the current API.

If you build the SMR model in Python, pass the resulting ``openmc.Model``
object directly to ``run_depletion`` instead.

2. Configure boron feedback
---------------------------

Install ``openmc-feedbacks`` in your OpenMC environment, then create a
``FeedbackControl``. The element symbol ``B`` selects all boron isotopes found
in the named coolant material, so B10 and B11 are adjusted together while
hydrogen, oxygen, and other coolant nuclides are left unchanged.

.. code-block:: python

   from openmc_feedbacks import FeedbackControl, run_depletion

   boron_feedback = FeedbackControl(
       material_name="Coolant with 1000 ppm boron",
       nuclides=("B",),
       reference_concentration=1000.0,  # ppm represented by model's B densities
       minimum=0.0,
       maximum=1000.0,
       target_keff=1.005,
       target_mode="end",  # omit to target beginning-of-step keff
       expected_step_keff_loss=0.01,  # calibrate with a pilot run
       search_particles=35000,
       search_batches=50,
       search_inactive=7,
       warm_start_width=100.0,
       max_search_refinements=2,
       tolerance=50.0,          # ppm interval width
       keff_tolerance=0.0005,
       max_iterations=10,
       history_file="boron_history.csv",
   )

The model's original particle, batch, and inactive-generation settings are
used for the final transport calculation at the selected concentration. The
lower search statistics above are used only for concentration trials.

3. Run depletion
----------------

Set ``power`` to the SMR's thermal power in watts and provide the timestep
durations you want to simulate. The values below show a short example at
160 MW thermal power; replace them with your operating power and full
irradiation schedule for a production calculation.

.. code-block:: python

   run_depletion(
       model=model,  # alternatively: model="path/to/smr-xml"
       timesteps=[1, 7, 30],
       timestep_units="d",
       power=160.0e6,  # watts
       chain_file="/path/to/chain_endfb80_pwr.xml",
       integrator="predictor",
       control=boron_feedback,
       output_dir="smr-depletion-results",
   )

At each depletion step, the controller searches for a boron concentration
using transport calculations, then uses the selected transport result for
that step's depletion calculation. It does not run trial depletion steps.
OpenMC results and the ``boron_history.csv`` log are written to
``smr-depletion-results``.

4. Adapt the example
--------------------

To use another reactor, provide its input model, chain file, thermal power,
timestep schedule, and feedback material. For another soluble absorber, set
``nuclides`` to an element symbol present in that material (for example
``"Gd"``) or to explicit isotope names (for example ``("Gd155", "Gd157")``).
Set ``reference_concentration`` to the concentration represented by the input
material densities. The feedback controller scales those listed nuclide
densities linearly; it does not convert between concentration units or
rebuild chemical compounds.

This example sets ``target_mode="end"``. The configured ``target_keff`` is the
desired end-of-step keff; the controller estimates the beginning-of-step
target by adding the filtered expected step loss. If you omit ``target_mode``
or set it to ``"beginning"``, ``target_keff`` applies directly to the
beginning-of-step transport calculation instead. End-of-step compensation is
an estimate, so calibrate it for your timestep schedule and operating
conditions.

