Liang's Analysis Note
=====================

.. sectionauthor:: Liang Liu

.. _installation:

Installation
------------

The `xsec-analyzer <https://github.com/uboone/xsec_analyzer>`_
framework is designed to be installed and run on the uboonegpvms. In principle
it can be used elsewhere if you have access to the 'PeLEE' ntuples. The only
external dependency is ROOT. To install the code, navigate to a writeable area
in your local area or (if one the :abbr:`GPVMs (general purpose virtual
machines)`) :file:`/exp/uboone/app` area.

.. code-block:: console

    $ cd your_workarea
    $ git clone https://github.com/uboone/xsec_analyzer.git
    $ git checkout –t origin/tutorial-umn
    $ source setup_xsec_analyzer.sh
    $ make

The setup script will attempt to set up ROOT. It will auto-detect whether
you're running :abbr:`SL7 (Scientific Linux 7)` or :abbr:`AL9 (AlmaLinux 9)`.
If you're running SL7, the setup script will set up ROOT via uboonecode. If
you're running AL9, ROOT will be set up by spack.

The setup script will set the path to source code in the environment variable
``XSEC_ANALYZER_DIR``. This is used by the framework to find certain files, add
``${XSEC_ANALYZER_DIR}/bin`` into the environment variable ``PATH``, add
``${XSEC_ANALYZER_DIR}/lib`` into the environment variable ``LD_LIBRARY_PATH``.

On a fresh login, just navigate to the same folder and source the
``setup_xsec_analyzer.sh`` script.

Workflow
--------

The `xsec-analyzer <https://github.com/uboone/xsec_analyzer>`_
provides six major analysis tools. They are executables in
``${XSEC_ANALYZER_DIR}/bin``:

.. code-block:: console

   $ ls ${XSEC_ANALYZER_DIR}/bin
   BinScheme  ProcessNTuples   SlicePlots  StandaloneUnfold  Unfolder  univmake

Here, I use ``CC1muXp0pi`` channel as an example to show the workflow of
xsec-analyzer.

``ProcessNTuples``
^^^^^^^^^^^^^^^^^^

The basic usage of ``ProcessNtuples``:

.. code-block:: console

  $ ProcessNTuples INPUT_PELEE_NTUPLE_FILE CC1mu1p0pi OUTPUT_FILE

Here, I use a script to handle the all the different type of ntuples for
different runs.

.. code-block:: console

  $ ./ReprocessNTuples.sh -h
      Usage: ./ReprocessNTuples.sh [-o <out directory>] [-v version] [-r runnumbers] [-s samples]
        -c configure   A text file includes all the PeLEE samples.
        -s selection   Name of your selection algorithm.
        -o directory   Specify the output directory. It should be in your data area e.g. /exp/uboone/data/users/liangliu
        -r run #       Specify runs to process (comma-separated, e.g., 1,2,3).
        -v version     Specify the version in teck note (e.g. v00_00_01)
        -t type        Specify the type of ntuples (e.g., numuMC,nueMC,dirtMC,extBNB,onBNB,openBNB,detVarCV ... )
        -h help        Print help info
  # Post-process the nu overlay of run 1
  $ ./ReprocessNTuples.sh -s CC1muXp0pi -r 1 -t numuMC
      numuMC,File Name: /exp/uboone/data/users/liangliu/hadd/run1/run1_bnb_nu_overlay_generator_resc_ntuple_ntuple_ana.root
      Index: 1
      Sample Type: numuMC
      Number of Events:
      Scaling Factor:
      Output file name: /exp/uboone/data/users/liangliu/xsec/xsec_CC1muXp0pi_v00_00_01_run1_bnb_nu_overlay_generator_resc_ntuple_ntuple_ana.root
      time ProcessNTuples /exp/uboone/data/users/liangliu/hadd/run1/run1_bnb_nu_overlay_generator_resc_ntuple_ntuple_ana.root CC1muXp0pi /exp/uboone/data/users/liangliu/xsec/xsec_CC1muXp0pi_v00_00_01_run1_bnb_nu_overlay_generator_resc_ntuple_ntuple_ana.root
      --------------------------------
  #  Post-process all the ntuples
  $ ./ReprocessNTuples.sh

The configuration file which collects all the PeLEE ntuples for my
``CC1muXp0pi`` analysis is
``${XSEC_ANALYZER_DIR}/configs/files_to_process_liangliu.txt``

- ``INPUT_PELEE_NTUPLE_FILE`` should be a PeLEE ntuples, e.g. `PeLEE Samples
  2023
  <https://docs.google.com/spreadsheets/d/1dX-W4DGTHeZbJLt2HvwXS4QDNeEwYKveHHSCkVrJcSU/edit?gid=0#gid=0>`_

- The second input argument is the name of selection algorithm. The
  ``CC1mu1p0pi`` is implemented in:

  - :file:`$XSEC_ANALYZER_DIR/include/XSecAnalyzer/Selections/CC1mu1p0pi.hh`
  - :file:`$XSEC_ANALYZER_DIR/src/selections/CC1mu1p0pi.cxx`

- ``OUTPUT_FILE`` is the post-processed ntuple. It includes a TParameter called
  ``summed_pot`` and a TTree called ``stv_tree``

  For MC samples, the ``summed_pot`` gives the simulated POT needed to
  scale to data. For real data, it is always zero. The POT and trigger
  counts for real data files are stored elsewhere (more on that later).
  This scaling is handled automatically by later parts of the framework,
  so no need to worry much about this item.

  Many branches in ``stv_tree`` are copied over directly from the PeLEE
  ntuples, some are new and analysis-specific. Name is a hold-over from a
  much older incarnation of the code

``BinScheme``
^^^^^^^^^^^^^

Plot the smearing matrix

.. code-block:: console

   $ BinScheme -c TutorialBinScheme

Save binning configuration into text files

.. code-block:: console

   $ BinScheme -s TutorialBinScheme

     ------------------------------------------------------------------
    | Welcome to ROOT 6.28/12                        https://root.cern |
    | (c) 1995-2024, The ROOT Team; conception: R. Brun, F. Rademakers |
    | Built for linuxx8664gcc on Jan 30 2024, 08:17:35                 |
    | From tags/v6-28-12@v6-28-12                                      |
    | With g++ (Spack GCC) 12.2.0                                      |
    | Try '.help'/'.?', '.demo', '.license', '.credits', '.quit'/'.q'  |
     ------------------------------------------------------------------

   non-option ARGV-elements: CCXp0piBinScheme
   muon_2d_bin
   stv_tree
   CC1muXp0pi
   174
   0 0 "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range >= 0.000 && CC1muXp0pi_sig_mc_num_proton_in_momentum_range < 1.000"
   0 0 "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range >= 1.000 && CC1muXp0pi_sig_mc_num_proton_in_momentum_range < 2.000"
   0 0 "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range >= 2.000 && CC1muXp0pi_sig_mc_num_proton_in_momentum_range < 3.000"
   0 0 "CC1muXp0pi_MC_Signal && CC1muXp0pi_sig_mc_num_proton_in_momentum_range >= 3.000 && CC1muXp0pi_sig_mc_num_proton_in_momentum_range < 10.000"
   ......
   163 1 164
   164 1 165
   165 1 166
   166 1 167
   Save universes bin configuration into => /exp/uboone/app/users/liangliu/analysis-code/tutorial/xsec_analyzer_eaf/configs/ccxp0pi_TKI_2D_bin_config.txt
   Save slice configuration into         => /exp/uboone/app/users/liangliu/analysis-code/tutorial/xsec_analyzer_eaf/configs/ccxp0pi_TKI_2D_slice_config.txt
   root [0]


``univmake``
^^^^^^^^^^^^

Using the output from step 1 and 2 to run univmake

.. code-block:: console

     # Usage:
     #  univmake LIST_FILE UNIVMAKE_CONFIG_FILE OUTPUT_ROOT_FILE
     $ univmake  $XSEC_ANALYZER_DIR/configs/file_properties_CC1muXp0pi_v00_00_01.txt $XSEC_ANALYZER_DIR/configs/ccxp0pi_TKI_2D_bin_config.txt /exp/uboone/data/users/liangliu/ntuple/

``SlicePlots``
^^^^^^^^^^^^^^

Once the unimake finished, we can plot the distributions that we configured
in Bin Scheme.

.. code-block:: console

   # Usage:
   # SlicePlots FILE_PROPERTIES SYS_CALC_CONF SLICE_CONF UNIV_FILE SLICE_OUTPUT_DIR
   $ SlicePlots ${XSEC_ANALYZER_DIR}/configs/file_properties_fsi_current_run3.txt ${XSEC_ANALYZER_DIR}/configs/systcalc.conf ${XSEC_ANALYZER_DIR}/configs/ccxp0pi_TKI_2D_slice_config.txt /exp/uboone/data/users/liangliu/workarea/fsi/univmake_tki_2d/univmake_tki_2d.root `pwd`/output

.. note::
   ``FILE_PROPERTIES`` is similar to ``LIST_FILE`` but they are different.
   ``LIST_FILE`` is just tell the analyzer framework the available samples and
   make universe for each of them. In ``FILE_PROPERTIES``, you need to
   configure the universe, to be precise, the universe of detvars to plot the
   distributions.

   - Runs 3, 4 and 5 have 9 different detvars
   - Run 2 have no detvar ntuples, we use run 1 and run3 to estimate run 2 detvars
   - Run 1 doesn't need LY Attenuation
   - MC generated for Run 4a with a special flux that models the misalignment of
     the beam -- so it's important to use the specific MC for that period weighted
     to the Run 4a POT (from Patrick)


``Unfolder``
^^^^^^^^^^^^

Subtract backgrounds, correct for inefficiency and bin-to-bin-smearing,
convert to cross-section units.

.. code-block:: console

   # Usage:
   #   Unfolder XSEC_CONF SLICE_CONF XSEC_OUTPUT_ROOT_FILE
   $ Unfolder xsec_config_fakedata_dagostini.txt ${XSEC_ANALYZER_DIR}/configs/ccxp0pi_TKI_2D_slice_config.txt xsec_muon_proton_fakedata_dagostini.root

.. code-block:: console

   UnivFile /exp/uboone/data/users/liangliu/workarea/fsi/univmake_tki_2d/univmake_tki_2d.root
   SystFile /exp/uboone/app/users/liangliu/analysis-code/tutorial/xsec_analyzer_eaf/configs/systcalc.conf
   FPFile /exp/uboone/app/users/liangliu/analysis-code/tutorial/xsec_analyzer_eaf/configs/file_properties_fsi_current_run12345_fakedata.txt
   Unfold DAgostini fm 0.025
   #Unfold WienerSVD 1 second-deriv
   Prediction uBTune "MicroBooNE Tune" univ CV
   Prediction FakeData "Fake data" univ FakeData
   #Prediction gv2 "GENIE 2.12.10" file /exp/uboone/app/users/gardiner/temp-gen/BuildEventGenerators/ubmc/comp-all/comp-gv2.root MicroBooNE_CC1MuNp_XSec_2D_PpCosp_nu_MC
   #Prediction gv3 "GENIE 3.0.6" file /exp/uboone/app/users/gardiner/temp-gen/BuildEventGenerators/ubmc/comp-all/comp-gv3.root MicroBooNE_CC1MuNp_XSec_2D_PpCosp_nu_MC
   #Prediction g1802a "GENIE 3.2.0 G18_02a" file /exp/uboone/app/users/gardiner/temp-gen/BuildEventGenerators/ubmc/comp-all/comp-gv3-g1802a.root MicroBooNE_CC1MuNp_XSec_2D_PpCosp_nu_MC
   #Prediction g2111a "GENIE 3.2.0 G21_11a" file /exp/uboone/app/users/gardiner/temp-gen/BuildEventGenerators/ubmc/comp-all/comp-gv3-g2111a.root MicroBooNE_CC1MuNp_XSec_2D_PpCosp_nu_MC
   #Prediction g2111b "GENIE 3.2.0 G21_11b" file /exp/uboone/app/users/gardiner/temp-gen/BuildEventGenerators/ubmc/comp-all/comp-gv3-g2111b.root MicroBooNE_CC1MuNp_XSec_2D_PpCosp_nu_MC
   #Prediction neut "NEUT 5.4.0.1" file /exp/uboone/app/users/gardiner/temp-gen/BuildEventGenerators/ubmc/comp-all/comp-neut.root MicroBooNE_CC1MuNp_XSec_2D_PpCosp_nu_MC
   #Prediction nuwro "NuWro 19.02.1" file /exp/uboone/app/users/gardiner/temp-gen/BuildEventGenerators/ubmc/comp-all/comp-nuwro.root MicroBooNE_CC1MuNp_XSec_2D_PpCosp_nu_MC
   #Prediction gibuu "GiBUU 2021.1" file /exp/uboone/app/users/gardiner/temp-gen/BuildEventGenerators/ubmc/mygibuu3.root MicroBooNE_CC1MuNp_XSec_2D_PpCosp_nu_MC


StandaloneUnfold
""""""""""""""""

Steven has presented `Standalone unfolding tutorial slides
<https://microboone-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=42842&filename=Unfolding-Tutorial-uB-Retreat-UMN.pdf&version=8>`_
at a MicroBooNE analysis retreat.



