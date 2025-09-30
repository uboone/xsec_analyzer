// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/CC1muNp1piBinScheme.hh"

CC1muNp1piBinScheme::CC1muNp1piBinScheme() : BinSchemeBase( "CC1muNp1piBinScheme" ) {}


void CC1muNp1piBinScheme::DefineBlocks() {

  /////// Set some standard variables before managing the blocks

  // TTree name for the post-processed ntuples
  ntuple_ttree_name_ = "stv_tree";

  // Run numbers to use when plotting migration matrices
  runs_to_use_ = { 1 };

  // Prefix for the output bin and slice configuration text files
  out_config_prefix_ = "CC1muNp1pi_xsec";

  // Selection to use with this binning scheme
  selection_name_ = "CC1muNp1pi";

  // TDirectory file name to use when producing the univmake output histograms
  out_tdir_name_ = "CC1muNp1pi";

  /////// Define the blocks of bins in both true and reco space

  // First block: cos_theta_mu in 1D
  /*
  std::vector< double > phi3D_edges = {0., 62., 98., 122., 140., 160., 180.};

  std::vector< double > p3_mu_edges = {};
  for (double i = 0.15; i <= 1.2; i += (1.2 - 0.15) / 40) {
    p3_mu_edges.push_back(i);
  }

  std::vector< double > p3_cpi_edges = {};
  for (double i = 0.05; i <= 0.6; i += (0.6 - 0.05) / 40) {
    p3_cpi_edges.push_back(i);
  }

  std::vector< double > p3_lead_p_edges = {};
  for (double i = 0.3; i <= 1.0; i += (1.0 - 0.3) / 40) {
    p3_lead_p_edges.push_back(i);
  }
*/
std::vector< double > pn_edges = {0., 0.12, 0.24, 0.39, 0.54, 0.69, 1.2};
std::vector< double > alpha3D_edges = {0., 62., 98., 122., 140., 160., 180.};
  // pn block
  
  Block1D* b1t = new Block1D( "true_gki_Pn", 
    "p_{n} [GeV]", "p_{n} [GeV]", pn_edges,
    "MC_Signal", kSignalTrueBin );

  Block1D* b1r = new Block1D( "reco_gki_Pn",
    "p_{n} [GeV]", "p_{n} [GeV]", pn_edges,
    "Selected", kOrdinaryRecoBin );

  vect_block.emplace_back( b1t, b1r );

  Block1D *b2t = new Block1D( "true_gki_DeltaAlpha3D",
    "#alpha_{3D}", "\\alpha_{3D}", alpha3D_edges,
    "MC_Signal", kSignalTrueBin );

  Block1D *b2r = new Block1D( "reco_gki_DeltaAlpha3D",
    "#alpha_{3D}", "\\alpha_{3D}", alpha3D_edges,
    "Selected", kOrdinaryRecoBin );

  vect_block.emplace_back( b2t, b2r );

 /*
  Block1D *b3t = new Block1D( "true_gki_DeltaPhi3D",
    "#phi_{3D}", "\\phi_{3D}", phi3D_edges,
    "MC_Signal", kSignalTrueBin );

  Block1D *b3r = new Block1D( "reco_gki_DeltaPhi3D",
    "#phi_{3D}", "\\phi_{3D}", phi3D_edges,
    "Selected", kOrdinaryRecoBin );

  vect_block.emplace_back( b3t, b3r );


  // p3_mu block
  Block1D* b4t = new Block1D( "true_p3_mu.Mag()",
    "p_{#mu} [GeV]", "p_{\\mu} [GeV]", p3_mu_edges,
    "MC_Signal", kSignalTrueBin );
  
  Block1D* b4r = new Block1D( "reco_p3_mu.Mag()",
    "p_{#mu} [GeV]", "p_{\\mu} [GeV]", p3_mu_edges,
    "Selected", kOrdinaryRecoBin );

  vect_block.emplace_back( b4t, b4r );

  // p3_cpi block
  Block1D* b5t = new Block1D( "true_p3_cpi.Mag()",
    "p_{#pi} [GeV]", "p_{\\pi} [GeV]", p3_cpi_edges,
    "MC_Signal", kSignalTrueBin );

  Block1D* b5r = new Block1D( "reco_p3_cpi.Mag()",
    "p_{#pi} [GeV]", "p_{\\pi} [GeV]", p3_cpi_edges,
    "Selected", kOrdinaryRecoBin );

  vect_block.emplace_back( b5t, b5r );

  // p3_lead_p block
  Block1D* b6t = new Block1D( "true_p3_lead_p.Mag()",
    "p_{p} [GeV]", "p_{p} [GeV]", p3_lead_p_edges,
    "MC_Signal", kSignalTrueBin );

  Block1D* b6r = new Block1D( "reco_p3_lead_p.Mag()",
    "p_{p} [GeV]", "p_{p} [GeV]", p3_lead_p_edges,
    "Selected", kOrdinaryRecoBin );

  vect_block.emplace_back( b6t, b6r );
*/
}
