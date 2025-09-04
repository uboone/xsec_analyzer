.. xsec_analyzer documentation master file, created by
   sphinx-quickstart on Thu Aug 28 16:49:36 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

xsec_analyzer
=============

A framework for operating MicroBooNE [#ub]_ cross-section analyses.

This framework uses the PeLEE-format ROOT ntuples produced by the relevant
groups in MicroBooNE and provides scaffolding for physics analyzers.

This allows analyzers to write code to define their targeted signal using true
variables and select their events of interest from reconstructed quantities.
The framework also handles propagating systematic uncertainties from different
sources and carries out the process of cross-section unfolding across
user-defined kinematic variables.

The general approach taken by the framework towards extracting neutrino-argon
cross-sections is described mathematically in Section B. of [Gardiner]_.

Getting Started
---------------

- Find out about the :doc:`executables` that the framework provides and how to
  use them.
- View the :ref:`PeLEE branch to AnalysisEvent member mapping
  <pelee-to-analysis-event>` to see how to access event variables.
- See existing :doc:`examples`.

.. toctree::
   :maxdepth: 2
   :hidden:

   Introduction <self>
   executables
   configuration-files
   adding-selection
   examples
   glossary

.. rubric:: Errata

.. [#ub] The `Micro Booster Neutrino Experiment <https://microboone.fnal.gov>`_.

.. rubric:: References

.. [Gardiner]
      Steven Gardiner -- `Mathematical Methods for Cross-Section Extraction` 2025,
      `arXiv:2401.04065 <https://arxiv.org/abs/2401.04065>`_ [hep-ex]
