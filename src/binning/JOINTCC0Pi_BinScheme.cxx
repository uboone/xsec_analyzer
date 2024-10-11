// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/JOINTCC0Pi_BinScheme.hh"

JOINTCC0Pi_BinScheme::JOINTCC0Pi_BinScheme() : BinSchemeBase( "JOINTCC0Pi_BinScheme" ) {}

void JOINTCC0Pi_BinScheme::DefineBlocks() {

  /////// Set some standard variables before managing the blocks

  // TTree name for the post-processed ntuples
  ntuple_ttree_name_ = "stv_tree";

  // Run numbers to use when plotting migration matrices
  runs_to_use_ = { 1 };

  // Prefix for the output bin and slice configuration text files
  out_config_prefix_ = "JOINTCC0pi_";

  // Selection to use with this binning scheme
  selection_name_ = "JOINTCC0pi";

  // TDirectory file name to use when producing the univmake output histograms
  out_tdir_name_ = "JOINTCC0pi";

  /////// Define the blocks of bins in both true and reco space
  std::string branchexpr_reco_scheme1,
  branchexpr_true_scheme2, 
  branchexpr_true_scheme1,
  branchexpr_reco_scheme2, 
  branchexpr_sideband, 
  title, textitle,
  selection_true, selection_reco;

  selection_true = "mc_is_cc0pi_signal";
  selection_reco = "sel_CC0pi && sel_muon_contained";

 branchexpr_reco_scheme1 = "p3_mu.Mag();GeV/c;p3_mu.CosTheta();";
 branchexpr_reco_scheme2 = "p3_mu.CosTheta();;p3_mu.Mag();GeV/c";


 branchexpr_true_scheme1 = "mc_p3_mu.Mag();GeV/c;mc_p3_mu.CosTheta();";
 branchexpr_true_scheme2 = "mc_p3_mu.CosTheta();;mc_p3_mu.Mag();GeV/c";

  // First block: cos_theta_mu in 2D
  
  std::map< double, std::vector<double> > MUON_2D_BIN_EDGES = {
    { 0.1, {  -1, -0.55,    0.,  0.45, 1.00 } },
    { 0.24, { -1, -0.55, -0.25,    0., 0.25,  0.45,  0.7, 1.00 } },
    { 0.3,  { -1,  -0.4,  -0.1,   0.1, 0.35,  0.5,   0.7, 0.85, 1. } },
    { 0.38, { -1,     0,   0.5,  0.65, 0.8,   0.92,  1.00 } },
    { 0.48, { -1,   0.2,   0.5,  0.65, 0.8,   0.875, 0.950, 1.00 } },
    { 0.7,  { -1,  0.65,   0.8, 0.875, 0.950, 1.00 } },
    { 0.85, { -1,  0.85,   0.9, 0.950, 1.00} },
    { 1.28, { -1,  0.85,  1.00} },
    { 1.58, { -1,  0.85,  1.00} },
    { 2.0, {} }

  };

  Block2D *b2dmuon_mom_costheta_true_scheme1 = new Block2D(branchexpr_true_scheme1, 
      "muon mom; GeV/c; cos#theta; ",
      "P_{\\mu}; GeV/c; \\cos\\theta; ",
      MUON_2D_BIN_EDGES, selection_true, kSignalTrueBin);
      
  Block2D *b2dmuon_mom_costheta_reco_scheme1 = new Block2D(branchexpr_reco_scheme1,
      "muon mom; GeV/c; cos#theta; ", 
      "P_{\\mu}; GeV/c; \\cos\\theta; ",
      MUON_2D_BIN_EDGES, selection_reco, kOrdinaryRecoBin);
      
      
  vect_block.emplace_back(b2dmuon_mom_costheta_true_scheme1, b2dmuon_mom_costheta_reco_scheme1);


  std::map< double, std::vector<double> > MUON_2D_BIN_EDGES_inclusive = {
// the 2D binning of inclusive but took the max limit of 2 GeV and not 2.5 GeV which its given  
{ -1.00, 
     {0.1, 0.24, 0.38, 2.0} },
{ -0.50,
    {0.1, 0.24, 0.3, 0.38, 0.48, 2.0} },
 { 0.00, 
   {0.1, 0.24, 0.3, 0.38, 0.48, 2.0} },
 { 0.27, 
     {0.1, 0.3, 0.48, 0.85, 2.0} },
 { 0.45, 
     {0.1, 0.38, 0.48, .68, 0.85, 2.0} },
 { 0.62, 
     {0.1, 0.38, 0.48, 0.68, 0.85, 1.28, 2.0} },
 { 0.76, 
     {0.1, 0.3, 0.48, 0.68, 0.85, 2.0} },
 { 0.86, 
     {0.1, 0.48, 0.68, 0.85, 1.28,  2.0} },
 { 0.94, 
     {0.1, 0.38, 0.48, 0.68, 0.85, 1.28, 2.0} },
 { 1.0, {} }
};



  // the title "title; unit" is used in plot in root style
  title ="cos#theta; p_{#mu}; GeV/c;";
  // the tex title "tex title; units" is used in latex format
  textitle = "P_{\\mu}; GeV/c; \\cos\\theta; ";
  // selection


Block2D *b2dmuon_mom_costheta_true_scheme2 = new Block2D(branchexpr_true_scheme2,
      title,
      textitle,
      MUON_2D_BIN_EDGES,
      selection_true, kSignalTrueBin);
      
  Block2D *b2dmuon_mom_costheta_reco_scheme2 = new Block2D(branchexpr_reco_scheme2,
       title,
      textitle,
      MUON_2D_BIN_EDGES,
      selection_reco, kOrdinaryRecoBin);
      

  vect_block.emplace_back(b2dmuon_mom_costheta_true_scheme2, b2dmuon_mom_costheta_reco_scheme2);

  // CATEGORY is the branchexpr
  // background_index is vector of background categories.
  CATEGORY = "JOINTCC0Pi_EventCategory";
  background_index = {0, 5, 7, 8, 9, 10, 11};


/*

 will leave out for now
 std::vector< double > Mult= {-0.5, 0.5, 1.5, 2.5, 3.5};

  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_true_lead_proton_delta_pT ; GeV/c";
  // the title "title; unit" is used in plot in root style
  title = "leading proton #Delta p_{T}; GeV/c";
  // the tex title "tex title; units" is used in latex format
  textitle = "leading proton \\Delta \\p_{T}; GeV/c";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_nProtons_in_Momentum_range > 0";

  Block1D *b1dt_delta_pT = new Block1D(branchexpr, title, textitle, delta_pT_1D_edges, selection, kSignalTrueBin);

*/


const std::string NC_SIDEBAND_SELECTION =
  "sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
  " && !sel_has_muon_candidate && sel_topo_cut_passed"
  " && sel_no_reco_showers && sel_has_p_candidate"
  " && sel_protons_contained"
  " && sel_n_bdt_other == 0" 
  " && sel_n_bdt_invalid == 0"
  " && sel_lead_p_passed_mom_cuts"
  " && sel_num_proton_candidates > 0";


const std::string CCNPI_SIDEBAND_SELECTION =
  "sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
  " && sel_has_muon_candidate && sel_topo_cut_passed"
  " && sel_no_reco_showers && sel_muon_passed_mom_cuts"
  " && sel_has_pion_candidate"
  " && sel_num_pion_candidates > 0"
  " && sel_n_bdt_other == 0" 
  " && sel_n_bdt_invalid == 0";
  
    
   std::vector<double> Pmu_binningSideband = { 0.1, 0.24, 0.3, 0.38, 0.48, 0.7, 0.85, 1.28, 1.58, 2.0};
   std::vector<double> Costheta_sideBand = {-1.00, -0.50, 0.00,  0.27,  0.45, 0.62, 0.76, 0.86, 0.94, 1.0};
   std::vector<double> proton_binningSideband = {0.25, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.8, 0.87, 0.93, 1.2};


 branchexpr_sideband = "p3_lead_p.Mag(); GeV/c";
   // the title "title; unit" is used in plot in root style
  title ="p_{p}; GeV/c";
  // the tex title "tex title; units" is used in latex format
  textitle = "p_{p}; GeV/c";

  Block1D *b1r_sideband_NC = new Block1D(branchexpr_sideband, title, textitle,proton_binningSideband, NC_SIDEBAND_SELECTION, kSidebandRecoBin);
  vect_sideband.emplace_back(b1r_sideband_NC);



 branchexpr_sideband = "p3_mu.Mag(); GeV/c";
   // the title "title; unit" is used in plot in root style
  title ="p_{#mu}; GeV/c";
  // the tex title "tex title; units" is used in latex format
  textitle = "p_{\\mu}; GeV/c";


  Block1D *b1r_sideband_Npi_pMu = new Block1D(branchexpr_sideband, title, textitle, Pmu_binningSideband, CCNPI_SIDEBAND_SELECTION, kSidebandRecoBin);
  vect_sideband.emplace_back(b1r_sideband_Npi_pMu);


 branchexpr_sideband = "p3_mu.CosTheta(); GeV/c";
   // the title "title; unit" is used in plot in root style
  title ="cos#theta;";
  // the tex title "tex title; units" is used in latex format
  textitle = "\\cos\\theta;";


  Block1D *b1r_sideband_Npi_costheta = new Block1D(branchexpr_sideband, title, textitle, Costheta_sideBand, CCNPI_SIDEBAND_SELECTION, kSidebandRecoBin);
  vect_sideband.emplace_back(b1r_sideband_Npi_costheta);

}
