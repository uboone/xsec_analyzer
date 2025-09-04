Writing Your Own Selection
==========================

.. attention::

   This documentation page is incomplete. Please view tutorial materiel for
   this framework for now.

Organising Your Selection
-------------------------

Here we describe a selection called ``YourSelection``. This looks like a class that

.. cpp:class:: YourSelection : SelectionBase

   YourSelection is a class that implements all of the required methods
   (described below) that SelectionBase and the framework requires to operate a
   cross-section analysis.

   This class should also have member variables added that correspond to any
   desired additional branches you'd like available in the output tree.

   .. cpp:function:: void define_constants()

      Here you should define any particular constants that might be helpful,
      e.g.

      - Threshold values that your selection cuts depends upon.
      - PDG codes of interest.

   .. cpp:function:: bool define_signal(AnalysisEvent* event)

      This method returns a boolean indicating if :cpp:var:`event` passed a
      signal definition.

      Truth information should be used to define the signal.

   .. cpp:function:: int categorize_event(AnalysisEvent* event)

      Here you can add logic that classifies events via their
      truth-information.

      A classification scheme can be as coarse or finely defined as you want,
      using MC interaction codes, final-state particles or a combination of
      both.

   .. cpp:function:: bool selection(AnalysisEvent* event)

      This method should contain logic to operate your selection cut on an
      event.

      This method can be run for files with and without truth-information.

   .. cpp:function:: void compute_reco_observables(AnalysisEvent* event)
   .. cpp:function:: void compute_true_observables(AnalysisEvent* event)
   .. cpp:function:: void define_output_branches()
   .. cpp:function:: void reset()

      Run after each event is processed, this should reset any class variables
      to default values so the values for the previous event are not persisted
      into the new event.


Event Variables in ``AnalysisEvent``
------------------------------------

.. rubric:: ``include/XSecAnalyzer/AnalysisEvent.hh``

.. cpp:class:: AnalysisEvent

   Represents a single interaction event with all the associated event
   variables as class instance members.

   Initially read from the PeLEE ntuples.

   .. cpp:member:: float topological_score_
   .. cpp:member:: float cosmic_impact_parameter_
   .. cpp:member:: float contained_fraction_

   .. cpp:member:: float nu_completeness_from_pfp_
   .. cpp:member:: float nu_purity_from_pfp_
   .. cpp:member:: int nu_pdg_

   .. cpp:member:: int nslice_

      Number of neutrino slices identified by the SliceID. Allowed values
      are zero or one.

   .. cpp:member:: float nu_vx_

      Reco neutrino vertex coordinates (cm). Space charge corrections have
      been applied for these.
   .. cpp:member:: float nu_vy_
   .. cpp:member:: float nu_vz_

   .. cpp:member:: int num_pf_particles_
   .. cpp:member:: int num_tracks_
   .. cpp:member:: int num_showers_

NuMI-Specific Configuration
---------------------------

Several changes need to be made to change from the default BNB-specific
operating mode to a NuMI-specific mode.

Adding Mock Weights to NuMI Dirt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

NuMI Dirt files must have additional mock systematic weights added.
``xsec_analyzer`` expects a full set of systematic weights for all input PeLEE
ntuples. While BNB dirt files have a full set of weights included, NuMI dirt
files do not (instead a large uncertainty is applied).

In order to apply these additional mock weights to NuMI dirt files, the
:ref:`exec-addfakeweights` program can be invoked. The new dirt file should be
used in the remaining workflow.

Compile Time Switch
^^^^^^^^^^^^^^^^^^^

At compile-time a variable should be changed in the framework:

.. rubric:: ``include/XSecAnalyzer/Constants.hh``

.. cpp:var:: bool useNuMI = true
   :no-contents-entry:

   Enable NuMI processing mode for the framework.
