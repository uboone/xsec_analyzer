Configuration Files
===================

Several text-configuration files are used to control the operation of the
framework.

.. _config-file-properties:

``file_properties.txt``
-----------------------

This file declares a series of input files produced by the ``ProcessNTuples``
program.

Example
^^^^^^^

.. code-block::

   # Full beam-on data
   /home/me/data/CrossSections/XSec_Analyzer/xsec-ntuples/xsec-ana-numi_beam_on_pion_ntuples_run1_fhc.root 1 onBNB 5748692 2.192e+20
   # Beam-off data
   /home/me/data/CrossSections/XSec_Analyzer/xsec-ntuples/xsec-ana-numi_beam_off_pion_ntuples_run1_fhc.root 1 extBNB 4582248 0

Format
^^^^^^

Here variables to replace are labelled as ``{such}``. This can be repeated on
each line for each declared file.

.. code-block::

   # a comment (ignored by the framework)
   {file-path} {run-number} {file-type} [{triggers} {POT}]
   ...

.. option:: file-path

   An absolute path to an ``stv-tree`` ntuple file as produced by
   :ref:`exec-processntuples`.

.. option:: run-number

   An integer corresponding to a MicroBooNE running period. Between 1 and 5
   inclusive.

.. option:: file-type

   A string corresponding to a category describing the input file. One of
   :ref:`sample-type-table`. Only one of a single file type can be declared in
   the configuration file.

   .. _sample-type-table:

   .. table:: ``file-type`` options

      +--------------------+-------------------------------------+
      | Identifier         | Description                         |
      +====================+=====================================+
      | ``onBNB``          | Beam-on data (BNB or NuMI)          |
      +--------------------+-------------------------------------+
      | ``extBNB``         | Beam-off data (BNB or NuMI)         |
      +--------------------+-------------------------------------+
      | ``numuMC``         | MC files containing :math:`\nu_\mu` |
      |                    | interactions                        |
      +--------------------+-------------------------------------+
      | ``nueMC``          | MC files containing :math:`\nu_e`   |
      |                    | interactions                        |
      +--------------------+-------------------------------------+
      | ``dirtMC``         | MC files with neutrino interactions |
      |                    | occuring outside of the fiducial    |
      |                    | volume                              |
      +--------------------+-------------------------------------+
      | ``detVarCV``       | Detector variation file with the    |
      |                    | central value prediction            |
      +--------------------+-------------------------------------+
      | ``detVarLYatten``  | Detector variation file with the    |
      |                    | attenuated light yield              |
      +--------------------+-------------------------------------+
      | ``detVarLYdown``   | NtupleFileType::kDetVarMCLYdown     |
      +--------------------+-------------------------------------+
      | ``detVarLYrayl``   | NtupleFileType::kDetVarMCLYrayl     |
      +--------------------+-------------------------------------+
      | ``detVarRecomb2``  | Detector variation with alternate   |
      |                    | recombination model                 |
      +--------------------+-------------------------------------+
      | ``detVarSCE``      | NtupleFileType::kDetVarMCSCE        |
      +--------------------+-------------------------------------+
      | ``detVarWMAngleXZ``| NtupleFileType::kDetVarMCWMAngleXZ  |
      +--------------------+-------------------------------------+
      | ``detVarWMAngleYZ``| NtupleFileType::kDetVarMCWMAngleYZ  |
      +--------------------+-------------------------------------+
      | ``detVarWMdEdx``   | NtupleFileType::kDetVarMCWMdEdx     |
      +--------------------+-------------------------------------+
      | ``detVarWMX``      | NtupleFileType::kDetVarMCWMX        |
      +--------------------+-------------------------------------+
      | ``detVarWMYZ``     | NtupleFileType::kDetVarMCWMYZ       |
      +--------------------+-------------------------------------+
      | ``detVarCVExtra``  | NtupleFileType::kDetVarMCCVExtra    |
      +--------------------+-------------------------------------+
      | ``altCVMC``        | NtupleFileType::kAltCVMC            |
      +--------------------+-------------------------------------+

The remaining options are only valid for data (beam-on or -off) or fake data
files:

.. option:: triggers

   The number of recorded triggers for the data file.

.. option:: POT

   The protons-on-target for the data file, this is the POT value all other
   MC-simulated samples are scaled to.

   - For data files this is normally extracted via a script.
   - For MC files this should be empty.
   - For beam-off data files this is 0.

.. _config-bin-config:

``bin-config.txt``
------------------

This file defines a binning scheme, marking the start and end values for each
histogram bin, these are given for each kinematic variable.

.. tip::

   Instead of writing this file by hand, you can have it generated for you via
   the :ref:`exec-binscheme` command. See its documentation for more
   information.

Example
^^^^^^^

.. code-block::

   nuecc_bin_config
   stv_tree
   NuMICC1e
   15
   0 0 "MC_Signal && mc_electron_energy >= 0.03 && mc_electron_energy < 0.3"
   0 0 "MC_Signal && mc_electron_energy >= 0.3 && mc_electron_energy < 0.47"
   0 0 "MC_Signal && mc_electron_energy >= 0.47 && mc_electron_energy < 0.70"
   0 0 "MC_Signal && mc_electron_energy >= 0.70 && mc_electron_energy < 0.99"
   0 0 "MC_Signal && mc_electron_energy >= 0.99 && mc_electron_energy < 1.43"
   0 0 "MC_Signal && mc_electron_energy >= 1.43 && mc_electron_energy < 3.0"
   0 0 "MC_Signal && mc_electron_energy >= 3.0"
   0 1 "MC_Signal"
   1 -1 "EventCategory == 1"
   1 -1 "EventCategory == 2"
   1 -1 "EventCategory == 3"
   1 -1 "EventCategory == 4"
   1 -1 "EventCategory == 5"
   1 -1 "EventCategory == 6"
   1 -1 "EventCategory == 7"
   8
   0 0 "sel_nu_e_cc && reco_electron_energy >= 0.03 && reco_electron_energy < 0.3"
   0 0 "sel_nu_e_cc && reco_electron_energy >= 0.3 && reco_electron_energy < 0.47"
   0 0 "sel_nu_e_cc && reco_electron_energy >= 0.47 && reco_electron_energy < 0.70"
   0 0 "sel_nu_e_cc && reco_electron_energy >= 0.70 && reco_electron_energy < 0.99"
   0 0 "sel_nu_e_cc && reco_electron_energy >= 0.99 && reco_electron_energy < 1.43"
   0 0 "sel_nu_e_cc && reco_electron_energy >= 1.43 && reco_electron_energy < 3.0"
   0 0 "sel_nu_e_cc && reco_electron_energy >= 3.0"
   0 1 "sel_nu_e_cc"


Format
^^^^^^

.. code-block::

   {???}
   {ntuple-name}
   {selection-name}
   {number-of-bins}
   {?} {?} {cut-condition}
   ...


.. _config-syst-config:

``syst-config.txt``
-------------------

Example
^^^^^^^

Format
^^^^^^
