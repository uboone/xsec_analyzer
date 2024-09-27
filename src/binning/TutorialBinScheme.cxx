// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/TutorialBinScheme.hh"

TutorialBinScheme::TutorialBinScheme() : BinSchemeBase( "TutorialBinScheme" ) {}

void TutorialBinScheme::DefineBlocks() {

  /////// Set some standard variables before managing the blocks

  // TTree name for the post-processed ntuples
  ntuple_ttree_name_ = "stv_tree";

  // Run numbers to use when plotting migration matrices
  runs_to_use_ = { 1 };

  // Prefix for the output bin and slice configuration text files
  out_config_prefix_ = "tutorial_";

  // Selection to use with this binning scheme
  selection_name_ = "TutorialCC1mu";

  // TDirectory file name to use when producing the univmake output histograms
  out_tdir_name_ = "tutorial_muon_1D_bin";

  /////// Define the blocks of bins in both true and reco space

  // First block: cos_theta_mu in 1D
  std::vector< double > cos_theta_mu_1D_edges = { -1., -0.85, -0.775, -0.7,
    -0.625, -0.55, -0.475, -0.4, -0.325, -0.25, -0.175, -0.1, -0.025, 0.05,
    0.125, 0.2, 0.275, 0.35, 0.425, 0.5, 0.575, 0.65, 0.725, 0.8, 0.85,
    0.875, 0.9, 0.925, 0.950, 0.975, 1. };

  Block1D* b1t = new Block1D( "TutorialCC1mu_true_p3_mu.CosTheta()",
    "muon cos#theta", "\\cos\\theta_{\\mu}", cos_theta_mu_1D_edges,
    "TutorialCC1mu_MC_Signal", kSignalTrueBin );

  Block1D* b1r = new Block1D( "TutorialCC1mu_reco_p3_mu.CosTheta()",
    "muon cos#theta", "\\cos\\theta_{\\mu}", cos_theta_mu_1D_edges,
    "TutorialCC1mu_Selected", kOrdinaryRecoBin );

  vect_block.emplace_back( b1t, b1r );

  // Second block: p_mu in 1D
  std::vector< double > pmu_1D_edges = { 0.1, 0.175, 0.2, 0.225, 0.25,
    0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.55, 0.6,
    0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2 };

  Block1D* b2t = new Block1D( "TutorialCC1mu_true_p3_mu.Mag()",
    "p_{#mu}; (GeV)", "p_{\\mu}; (GeV)", pmu_1D_edges,
    "TutorialCC1mu_MC_Signal", kSignalTrueBin );

  Block1D* b2r = new Block1D( "TutorialCC1mu_reco_p3_mu.Mag()",
    "p_{#mu}; (GeV)", "p_{\\mu}; (GeV)", pmu_1D_edges,
    "TutorialCC1mu_Selected", kOrdinaryRecoBin );

  vect_block.emplace_back( b2t, b2r );
}
