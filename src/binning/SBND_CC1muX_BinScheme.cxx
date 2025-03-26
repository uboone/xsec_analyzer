// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/SBND_CC1muX_BinScheme.hh"

#include <limits>

SBND_CC1muX_BinScheme::SBND_CC1muX_BinScheme() : BinSchemeBase( "SBND_CC1muX_BinScheme" ) {}

void SBND_CC1muX_BinScheme::DefineBlocks() {

  /////// Set some standard variables before managing the blocks

  // TTree name for the post-processed ntuples
  ntuple_ttree_name_ = "stv_tree";

  // Run numbers to use when plotting migration matrices
  runs_to_use_ = { 1 };

  // Prefix for the output bin and slice configuration text files
  out_config_prefix_ = "SBND_CC1muX_";

  // Selection to use with this binning scheme
  selection_name_ = "SBND_CC1muX";

  // TDirectory file name to use when producing the univmake output histograms
  out_tdir_name_ = "SBND_CC1muX_1D";

  /////// Define the blocks of bins in both true and reco space

  // First block: cos_theta_mu in 1D
  std::vector< double > cos_theta_mu_1D_edges = { -1., -0.5, 0., 0.27,
  0.45, 0.62, 0.76, 0.86, 0.94, 1.}; 

  Block1D* b1t = new Block1D( "true_leading_muon_costheta",
    "muon cos#theta", "\\cos\\theta_{\\mu}", cos_theta_mu_1D_edges,
    "MC_Signal", kSignalTrueBin );

  Block1D* b1r = new Block1D( "leading_muon_costheta",
    "muon cos#theta", "\\cos\\theta_{\\mu}", cos_theta_mu_1D_edges,
    "Selected", kOrdinaryRecoBin );

  vect_block.emplace_back( b1t, b1r );

  // Second block: p_mu in 1D
  std::vector< double > pmu_1D_edges = { 0., 0.3, 0.5, 0.7, 0.9, 1.1, 1.3,
  1.5, 2.0, 3.0, 10000.};

  Block1D* b2t = new Block1D( "true_leading_muon_momentum",
    "p_{#mu}; (GeV)", "p_{\\mu}; (GeV)", pmu_1D_edges,
    "MC_Signal", kSignalTrueBin );

  Block1D* b2r = new Block1D( "leading_muon_momentum",
    "p_{#mu}; (GeV)", "p_{\\mu}; (GeV)", pmu_1D_edges,
    "Selected", kOrdinaryRecoBin );

  vect_block.emplace_back( b2t, b2r );
}
