#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class CC1muNp1pi : public SelectionBase {

public:

  CC1muNp1pi();

  virtual int categorize_event( AnalysisEvent* Event ) override final;
  virtual bool selection( AnalysisEvent* Event ) override final;
  virtual bool define_signal( AnalysisEvent* Event ) override final;
  virtual void compute_reco_observables( AnalysisEvent* Event ) override final;
  virtual void compute_true_observables( AnalysisEvent* Event ) override final;
  virtual void define_output_branches() override final;
  virtual void define_constants() override final;
  virtual void define_category_map() override final;
  virtual void reset() override final;

private:

  /* SIGNAL DEFINITION FLAGS */
  //================================================================================================================

  bool sig_isNuMu_;
  bool sig_inFV_;
  bool sig_leadProtonMomInRange_;
  bool sig_muonInMomRange_;
  bool sig_noFSMesons_;
  bool sig_mc_no_fs_pi0_;
  int sig_nProtons_in_Momentum_range_;
  int sig_nCpi_in_Momentum_range_;


  /* SELECTION DEF FLAGS */
  //================================================================================================================

  bool sel_reco_vertex_in_FV_;
  bool sel_pfp_starts_in_PCV_;
  bool sel_cosmic_ip_cut_passed_;
  bool sel_topo_cut_passed_;
  bool sel_nu_mu_cc_;

  bool sel_muon_contained_;
  bool sel_pion_contained_;
  bool sel_protons_contained_;

  bool sel_no_reco_showers_;
  bool sel_min_3_tracks_;
  bool sel_2_non_proton_;

  bool sel_all_pfp_in_vtx_proximity_;
  bool sel_all_pfp_contained_;

  bool sel_has_muon_candidate_;
  bool sel_has_pion_candidate_;
  bool sel_has_p_candidate_;
  
  //bool sel_passed_proton_pid_cut_;
  bool sel_muon_passed_mom_cuts_;
  bool sel_pion_passed_mom_cuts_;
  bool sel_lead_p_passed_mom_cuts_;
  bool sel_muon_quality_ok_;
  

  /* RECO VARIABLES */
  //================================================================================================================

  int n_reco_tracks_;
  int n_non_proton_like_;

  int lead_p_candidate_idx_;
  int muon_candidate_pid_idx_;
  int muon_candidate_length_idx_;
  int muon_candidate_bdt_idx_;
  int pion_candidate_idx_;

  double delta_pT_;
  double delta_phiT_;
  double delta_alphaT_;
  double delta_pL_;
  double pn_;
  double delta_pTx_;
  double delta_pTy_;
  double theta_mu_p_;
  double theta_mu_cpi_;

  MyPointer< TVector3 > p3mu_;
  MyPointer< TVector3 > p3cpi_;
  MyPointer< TVector3 > p3p_;
  MyPointer< std::vector< TVector3 > > p3_p_vec_;
  
  /* MC VARIABLES */
  //================================================================================================================
 
  bool mc_golden_;
  bool mc_pi_stopping_;
  double mc_delta_pT_;
  double mc_delta_phiT_;
  double mc_delta_alphaT_;
  double mc_delta_pL_;
  double mc_pn_;
  double mc_delta_pTx_;
  double mc_delta_pTy_;
  double mc_theta_mu_p_;
  double mc_theta_mu_cpi_;

  int mc_n_protons_;

  MyPointer< TVector3 > mc_p3mu_;
  MyPointer< TVector3 > mc_p3p_;
  MyPointer< TVector3 > mc_p3cpi_;

  MyPointer< std::vector< TVector3 > > mc_p3_p_vec_;
  MyPointer< std::vector< TVector3 > > mc_p3_cpi_vec_;

  STVCalcType calc_type;
};
