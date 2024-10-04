// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/CCXp0piBinScheme.hh"

CCXp0piBinScheme::CCXp0piBinScheme() : BinSchemeBase( "CCXp0piBinScheme" ) {}

void CCXp0piBinScheme::DefineBlocks() {

  /////// Set some standard variables before managing the blocks

  // TTree name for the post-processed ntuples
  ntuple_ttree_name_ = "stv_tree";

  // Run numbers to use when plotting migration matrices
  runs_to_use_ = { 1 };

  // Prefix for the output bin and slice configuration text files
  out_config_prefix_ = "ccxp0pi_";

  // Selection to use with this binning scheme
  selection_name_ = "CC1muXp0pi";

  // TDirectory file name to use when producing the univmake output histograms
  out_tdir_name_ = "muon_2d_bin";

  /////// Define the blocks of bins in both true and reco space

// First block: cos_theta_mu in 1D
  std::vector< double > cos_theta_mu_1D_edges = { -1., -0.85, -0.775, -0.7,
    -0.625, -0.55, -0.475, -0.4, -0.325, -0.25, -0.175, -0.1, -0.025, 0.05,
    0.125, 0.2, 0.275, 0.35, 0.425, 0.5, 0.575, 0.65, 0.725, 0.8, 0.85,
    0.875, 0.9, 0.925, 0.950, 0.975, 1. };

  Block1D* b1t = new Block1D( "CC1muXp0pi_true_muon_costh",
    "muon cos#theta", "\\cos\\theta_{\\mu}", cos_theta_mu_1D_edges,
    "CC1muXp0pi_MC_Signal", kSignalTrueBin );

  Block1D* b1r = new Block1D( "CC1muXp0pi_reco_muon_costh",
    "muon cos#theta", "\\cos\\theta_{\\mu}", cos_theta_mu_1D_edges,
    "CC1muXp0pi_Selected", kOrdinaryRecoBin );

  vect_block.emplace_back( b1t, b1r );

  // Second block: p_mu in 1D
  std::vector< double > pmu_1D_edges = { 0.1, 0.175, 0.2, 0.225, 0.25,
    0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.55, 0.6,
    0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2 };

  Block1D* b2t = new Block1D( "CC1muXp0pi_true_muon_p",
    "p_{#mu}; (GeV)", "p_{\\mu}; (GeV)", pmu_1D_edges,
    "CC1muXp0pi_MC_Signal", kSignalTrueBin );

  Block1D* b2r = new Block1D( "CC1muXp0pi_reco_muon_p",
    "p_{#mu}; (GeV)", "p_{\\mu}; (GeV)", pmu_1D_edges,
    "CC1muXp0pi_Selected", kOrdinaryRecoBin );

  vect_block.emplace_back( b2t, b2r );


  // Using floating-point numbers as std::map keys is admittedly evil, but
  // it's safe in this case: all we'll do with this map is iterate over the
  // elements. Keys are muon momentum bin edges, values are muon scattering
  // cosine bin edges.
  std::map< double, std::vector<double> > MUON_2D_BIN_EDGES = {
//
    // No need for an underflow bin: due to the signal definition, all muons
    // with reco momentum below 0.1 GeV/c will be lost
    { 0.1, { -1, -0.55, -0.25, 0., 0.25, 0.45, 0.7, 1.00 }, },
    { 0.24, { -1, -0.55, -0.25, 0., 0.25, 0.45, 0.7, 1.00 } },
    { 0.3,  { -1, -0.4, -0.1, 0.1, 0.35, 0.5, 0.7, 0.85, 1. } },
    { 0.38, { -1, 0, 0.5, 0.65, 0.8, 0.92, 1.00 } },
    { 0.48, { -1, 0.2, 0.5, 0.65, 0.8, 0.875, 0.950, 1.00 } },
    { 0.7, { -1, 0.65, 0.8, 0.875, 0.950, 1.00 } },
    { 0.85, { -1, 0.85, 0.9, 0.950, 1.00 } },

    // Upper edge of the last bin. Due to the signal definition, no overflow
    // bin is needed for muons above 1.2 GeV/c
    { 1.2, {} }

  };

  Block2D *b2dmuon_mom_costheta_true = new Block2D("CC1muXp0pi_true_muon_p; GeV/c; CC1muXp0pi_true_muon_costh; ", 
      "muon mom; GeV/c; cos#theta; ",
      "P_{\\mu}; GeV/c; \\cos\\theta; ",
      MUON_2D_BIN_EDGES, "CC1muXp0pi_MC_Signal", kSignalTrueBin);
  Block2D *b2dmuon_mom_costheta_reco = new Block2D("CC1muXp0pi_reco_muon_p; GeV/c; CC1muXp0pi_reco_muon_costh; ",
      "muon mom; GeV/c; cos#theta; ", 
      "P_{\\mu}; GeV/c; \\cos\\theta; ",
      MUON_2D_BIN_EDGES, "CC1muXp0pi_Selected", kOrdinaryRecoBin);
  vect_block.emplace_back(b2dmuon_mom_costheta_true, b2dmuon_mom_costheta_reco);

  std::string branchexpr, title, textitle, selection;

  std::vector< double > delta_pT_1D_edges = { 0., 0.06, 0.12, 0.18, 0.24, 0.32,
      0.4, 0.48, 0.55, 0.68, 0.75, 0.9 };
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_true_lead_proton_delta_pT ; GeV/c";
  // the title "title; unit" is used in plot in root style
  title = "leading proton #Delta p_{T}; GeV/c";
  // the tex title "tex title; units" is used in latex format
  textitle = "leading proton \\Delta \\p_{T}; GeV/c";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_nProtons_in_Momentum_range > 0";

  Block1D *b1dt_delta_pT = new Block1D(branchexpr, title, textitle, delta_pT_1D_edges, selection, kSignalTrueBin);

  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_reco_lead_proton_delta_pT ; GeV/c";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_num_proton_candidates > 0";

  // only the name of branch and the selection is different from true.
  Block1D *b1dr_delta_pT = new Block1D(branchexpr, title, textitle, delta_pT_1D_edges, selection, kOrdinaryRecoBin);

  vect_block.emplace_back(b1dt_delta_pT, b1dr_delta_pT);


  std::vector< double > delta_pTx_1D_edges = { -0.8, -0.6, -0.45, -0.35, -0.25,
      -0.15, -0.075, 0, 0.075, 0.15, 0.25, 0.35, 0.45, 0.6, 0.8 };
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_true_lead_proton_delta_pTx ; GeV/c";
  // the title "title; unit" is used in plot in root style
  title = "leading proton #Delta p_{Tx}; GeV/c";
  // the tex title "tex title; units" is used in latex format
  textitle = "leading proton \\Delta \\p_{Tx}; GeV/c";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_nProtons_in_Momentum_range > 0";

  Block1D *b1dt_delta_pTx = new Block1D(branchexpr, title, textitle, delta_pTx_1D_edges, selection, kSignalTrueBin);

  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_reco_lead_proton_delta_pTx ; GeV/c";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_num_proton_candidates > 0";

  // only the name of branch and the selection is different from true.
  Block1D *b1dr_delta_pTx = new Block1D(branchexpr, title, textitle, delta_pTx_1D_edges, selection, kOrdinaryRecoBin);

  vect_block.emplace_back(b1dt_delta_pTx, b1dr_delta_pTx);


  std::vector< double > delta_pTy_1D_edges = { -1.0, -0.6, -0.45, -0.35, -0.25,
      -0.15, -0.075, 0, 0.075, 0.15, 0.25, 0.35, 0.6 };
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_true_lead_proton_delta_pTy ; GeV/c";
  // the title "title; unit" is used in plot in root style
  title = "leading proton #Delta p_{Ty}; GeV/c";
  // the tex title "tex title; units" is used in latex format
  textitle = "leading proton \\Delta \\p_{Ty}; GeV/c";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_nProtons_in_Momentum_range > 0";

  Block1D *b1dt_delta_pTy = new Block1D(branchexpr, title, textitle, delta_pTy_1D_edges, selection, kSignalTrueBin);

  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_reco_lead_proton_delta_pTy ; GeV/c";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_num_proton_candidates > 0";

  // only the name of branch and the selection is different from true.
  Block1D *b1dr_delta_pTy = new Block1D(branchexpr, title, textitle, delta_pTy_1D_edges, selection, kOrdinaryRecoBin);

  vect_block.emplace_back(b1dt_delta_pTy, b1dr_delta_pTy);


  std::vector< double > delta_phiT_1D_edges = { 0., 0.1, 0.25, 0.45, 0.7, 1.0,
      1.3, 1.6, 1.9, 2.2, 2.5, 2.8, TMath::Pi()};
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_true_lead_proton_delta_phiT ;";
  // the title "title; unit" is used in plot in root style
  title = "leading proton #Delta #phi_{T};";
  // the tex title "tex title; units" is used in latex format
  textitle = "leading proton \\Delta \\phi T;";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_nProtons_in_Momentum_range > 0";

  Block1D *b1dt_delta_phiT = new Block1D(branchexpr, title, textitle, delta_phiT_1D_edges, selection, kSignalTrueBin);

  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_reco_lead_proton_delta_phiT ;";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_num_proton_candidates > 0";

  // only the name of branch and the selection is different from true.
  Block1D *b1dr_delta_phiT = new Block1D(branchexpr, title, textitle, delta_phiT_1D_edges, selection, kOrdinaryRecoBin);

  vect_block.emplace_back(b1dt_delta_phiT, b1dr_delta_phiT);


  std::vector< double > delta_alphaT_1D_edges = { 0, 0.261799, 0.523599, 0.785398, 1.0472, 1.309, 1.5708, 1.8326, 2.0944, 2.35619, 2.61799, 2.87979, TMath::Pi() };
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_true_lead_proton_delta_alphaT ;";
  // the title "title; unit" is used in plot in root style
  title = "leading proton #Delta #alpha_{T};";
  // the tex title "tex title; units" is used in latex format
  textitle = "leading proton \\Delta \\alpha T;";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_nProtons_in_Momentum_range > 0";

  Block1D *b1dt_delta_alphaT = new Block1D(branchexpr, title, textitle, delta_alphaT_1D_edges, selection, kSignalTrueBin);

  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_reco_lead_proton_delta_alphaT ;";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_num_proton_candidates > 0";

  // only the name of branch and the selection is different from true.
  Block1D *b1dr_delta_alphaT = new Block1D(branchexpr, title, textitle, delta_alphaT_1D_edges, selection, kOrdinaryRecoBin);

  vect_block.emplace_back(b1dt_delta_alphaT, b1dr_delta_alphaT);


  std::vector< double > sum_delta_phiT_1D_edges = { 0., 0.1, 0.25, 0.45, 0.7, 1.0,
      1.3, 1.6, 1.9, 2.2, 2.5, 2.8, TMath::Pi()};
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_true_hadron_delta_phiT ; ";
  // the title "title; unit" is used in plot in root style
  title = "sum of protons #Delta #phi_{T}; ";
  // the tex title "tex title; units" is used in latex format
  textitle = "sum of protons \\Delta \\phi T; ";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_nProtons_in_Momentum_range > 0";

  Block1D *b1dt_sum_delta_phiT = new Block1D(branchexpr, title, textitle, sum_delta_phiT_1D_edges, selection, kSignalTrueBin);

  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_reco_hadron_delta_phiT ; ";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_num_proton_candidates > 0";

  // only the name of branch and the selection is different from true.
  Block1D *b1dr_sum_delta_phiT = new Block1D(branchexpr, title, textitle, sum_delta_phiT_1D_edges, selection, kOrdinaryRecoBin);

  vect_block.emplace_back(b1dt_sum_delta_phiT, b1dr_sum_delta_phiT);


  std::vector< double > sum_delta_alphaT_1D_edges = { 0, 0.261799, 0.523599, 0.785398, 1.0472, 1.309, 1.5708, 1.8326, 2.0944, 2.35619, 2.61799, 2.87979, TMath::Pi() };
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_true_hadron_delta_alphaT ; ";
  // the title "title; unit" is used in plot in root style
  title = "sum of protons #Delta #alpha_{T};";
  // the tex title "tex title; units" is used in latex format
  textitle = "sum of protons \\Delta \\alpha T; ";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_nProtons_in_Momentum_range > 0";

  Block1D *b1dt_sum_delta_alphaT = new Block1D(branchexpr, title, textitle, sum_delta_alphaT_1D_edges, selection, kSignalTrueBin);

  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_reco_hadron_delta_alphaT ; ";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_num_proton_candidates > 0";

  // only the name of branch and the selection is different from true.
  Block1D *b1dr_sum_delta_alphaT = new Block1D(branchexpr, title, textitle, sum_delta_alphaT_1D_edges, selection, kOrdinaryRecoBin);

  vect_block.emplace_back(b1dt_sum_delta_alphaT, b1dr_sum_delta_alphaT);

  std::vector< double > theta_mu_p_1D_edges = { 0, 0.261799, 0.523599, 0.785398, 1.0472, 1.309, 1.5708, 1.8326, 2.0944, 2.35619, 2.61799, 2.87979, TMath::Pi() };
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_true_lead_proton_theta_mu_p ; ";
  // the title "title; unit" is used in plot in root style
  title = "leading proton #theta_{#mu p};";
  // the tex title "tex title; units" is used in latex format
  textitle = "leading proton \\theta_{\\mu p}; ";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_nProtons_in_Momentum_range > 0";

  Block1D *b1dt_theta_mu_p = new Block1D(branchexpr, title, textitle, theta_mu_p_1D_edges, selection, kSignalTrueBin);

  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_reco_lead_proton_theta_mu_p ; ";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_num_proton_candidates > 0";

  // only the name of branch and the selection is different from true.
  Block1D *b1dr_theta_mu_p = new Block1D(branchexpr, title, textitle, theta_mu_p_1D_edges, selection, kOrdinaryRecoBin);

  vect_block.emplace_back(b1dt_theta_mu_p, b1dr_theta_mu_p);

  std::vector< double > sum_theta_mu_p_1D_edges = { 0, 0.261799, 0.523599, 0.785398, 1.0472, 1.309, 1.5708, 1.8326, 2.0944, 2.35619, 2.61799, 2.87979, TMath::Pi() };
  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_true_hadron_theta_mu_p ; ";
  // the title "title; unit" is used in plot in root style
  title = "sum of protons #theta_{#mu p};";
  // the tex title "tex title; units" is used in latex format
  textitle = "sum of protons \\theta_{\\mu p}; ";
  // selection
  selection = "CC1muXp0pi_MC_Signal && CC1muXp0pi_nProtons_in_Momentum_range > 0";

  Block1D *b1dt_sum_theta_mu_p = new Block1D(branchexpr, title, textitle, sum_theta_mu_p_1D_edges, selection, kSignalTrueBin);

  // the branch name of the truth of proton momentum in stv root files; unit is optional with format "branch name; unit"
  branchexpr = "CC1muXp0pi_reco_hadron_theta_mu_p ; ";
  // the title "title; unit" is used in plot in root style
  // 
  // the tex title "tex title; units" is used in latex format
  // selection
  selection = "CC1muXp0pi_Selected && CC1muXp0pi_num_proton_candidates > 0";

  // only the name of branch and the selection is different from true.
  Block1D *b1dr_sum_theta_mu_p = new Block1D(branchexpr, title, textitle, sum_theta_mu_p_1D_edges, selection, kOrdinaryRecoBin);

  vect_block.emplace_back(b1dt_sum_theta_mu_p, b1dr_sum_theta_mu_p);

  // a example to use sideband bins; now we only use sideband bins for reco. 
  // you can definie many blocks of sideband with different selections
  std::vector< double > test_sideband = {1,2,3,4};
  
  Block1D *b1r_sideband = new Block1D("sideband", "sideband", "sideband",test_sideband , "selection", kSidebandRecoBin);

  vect_sideband.emplace_back(b1r_sideband);

  // CATEGORY is the branchexpr
  // background_index is vector of background categories.
  CATEGORY = "CC1muXp0pi_EventCategory";
  background_index = {17, 18, 19, 20, 21, 22};

        // 
}
