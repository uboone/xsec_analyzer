#ifndef __CC1muXp0pi_h__
#define __CC1muXp0pi_h__


#include "XSecAnalyzer/Selections/SelectionBase.hh"

#define DEFAULT_NUM 20

class CC1muXp0pi : public SelectionBase {
 public:
  CC1muXp0pi();

  int CategorizeEvent(AnalysisEvent* Event);
  bool Selection(AnalysisEvent* Event);
  bool DefineSignal(AnalysisEvent* Event);
  void ComputeRecoObservables(AnalysisEvent* Event);
  void ComputeTrueObservables(AnalysisEvent* Event);
  void DefineOutputBranches();
  void DefineConstants();
  void DefineCategoryMap();
  
private:
  bool sig_isNuMu_;
  bool sig_inFV_;
  bool sig_leadProtonMomInRange_;
  bool sig_muonInMomRange_;
  bool sig_noFSMesons_;
  bool sig_mc_no_fs_pi0_;
  bool sig_mc_no_fs_eta_;
  bool sig_mc_no_charged_pi_above_threshold_;
  int sig_nProtons_in_Momentum_range;

  int sig_num_gamma;
  int sig_num_muplus;
  int sig_num_muminus;
  int sig_num_eplus;
  int sig_num_eminus;
  int sig_num_nu;
  int sig_num_antinu;
  int sig_num_proton;
  int sig_num_neutron;
  
  bool sel_reco_vertex_in_FV_;
  bool sel_pfp_starts_in_PCV_;
  bool sel_has_muon_candidate_;
  bool sel_topo_cut_passed_;
  bool sel_nu_mu_cc_;
  bool sel_muon_contained_;
  bool sel_muon_passed_mom_cuts_;
  bool sel_no_reco_showers_;
  bool sel_has_p_candidate_;
  bool sel_muon_quality_ok_;
  bool sel_protons_contained_;
  bool sel_passed_proton_pid_cut_;
  bool sel_lead_p_passed_mom_cuts_;
  bool sel_cosmic_ip_cut_passed_;
  bool sel_has_non_proton_particles_;
  
  int sel_num_muon_candidates;
  int sel_muon_candidate_index[DEFAULT_NUM];
  double sel_muon_candidate_range_mom[DEFAULT_NUM];
  double sel_muon_candidate_mcs_mom[DEFAULT_NUM];
  double sel_muon_candidate_px[DEFAULT_NUM];
  double sel_muon_candidate_py[DEFAULT_NUM];
  double sel_muon_candidate_pz[DEFAULT_NUM];
  int sel_muon_bt_pdg[DEFAULT_NUM];

  double sel_muon_quality[DEFAULT_NUM];
  int sel_num_proton_candidates;
  int sel_proton_candidate_index[DEFAULT_NUM];
  double sel_proton_candidate_e[DEFAULT_NUM];
  double sel_proton_candidate_px[DEFAULT_NUM];
  double sel_proton_candidate_py[DEFAULT_NUM];
  double sel_proton_candidate_pz[DEFAULT_NUM];
  int lead_p_candidate_idx_;
  int muon_candidate_idx_;
  int muon_candidate_filter_idx_;


  int sel_num_rest_particles;
  double sel_rest_llr_pid_score[DEFAULT_NUM];
  double sel_rest_track_score[DEFAULT_NUM];
  double sel_rest_track_length[DEFAULT_NUM];
  double sel_rest_proton_mom[DEFAULT_NUM];
  int sel_rest_bt_pdg[DEFAULT_NUM];

  double hadron_delta_pT_;
  double hadron_delta_phiT_;
  double hadron_delta_alphaT_;
  double hadron_delta_alpha3Dq_;
  double hadron_delta_alpha3DMu_;
  double hadron_delta_phi3D_;
  double hadron_delta_pL_;
  double hadron_pn_;
  double hadron_delta_pTx_;
  double hadron_delta_pTy_;
  double hadron_theta_mu_p_;

  double lead_proton_delta_pT_;
  double lead_proton_delta_phiT_;
  double lead_proton_delta_alphaT_;
  double lead_proton_delta_alpha3Dq_;
  double lead_proton_delta_alpha3DMu_;
  double lead_proton_delta_phi3D_;
  double lead_proton_delta_pL_;
  double lead_proton_pn_;
  double lead_proton_delta_pTx_;
  double lead_proton_delta_pTy_;
  double lead_proton_theta_mu_p_;






  MyPointer<TVector3> p3mu;
  MyPointer<TVector3> p3p;
  MyPointer<std::vector<TVector3>> p3_p_vec_;


  double muon_energy_;
  double muon_costh_;
  double muon_p_;

  double mc_muon_energy_;
  double mc_muon_costh_;
  double mc_muon_p_;

  double mc_lead_proton_energy_;
  double mc_lead_proton_costh_;
  double mc_lead_proton_p_;

  double mc_hadron_energy_;
  double mc_hadron_costh_;
  double mc_hadron_p_;
  
  double mc_hadron_delta_pT_;
  double mc_hadron_delta_phiT_;
  double mc_hadron_delta_alphaT_;
  double mc_hadron_delta_alpha3Dq_;
  double mc_hadron_delta_alpha3DMu_;
  double mc_hadron_delta_phi3D_;
  double mc_hadron_delta_pL_;
  double mc_hadron_pn_;
  double mc_hadron_delta_pTx_;
  double mc_hadron_delta_pTy_;
  double mc_hadron_theta_mu_p_;

  double mc_lead_proton_delta_pT_;
  double mc_lead_proton_delta_phiT_;
  double mc_lead_proton_delta_alphaT_;
  double mc_lead_proton_delta_alpha3Dq_;
  double mc_lead_proton_delta_alpha3DMu_;
  double mc_lead_proton_delta_phi3D_;
  double mc_lead_proton_delta_pL_;
  double mc_lead_proton_pn_;
  double mc_lead_proton_delta_pTx_;
  double mc_lead_proton_delta_pTy_;
  double mc_lead_proton_theta_mu_p_;

  double mc_lead_proton2_delta_pT_;
  double mc_lead_proton2_delta_phiT_;
  double mc_lead_proton2_delta_alphaT_;
  double mc_lead_proton2_delta_pL_;
  double mc_lead_proton2_pn_;
  double mc_lead_proton2_delta_pTx_;
  double mc_lead_proton2_delta_pTy_;
  double mc_lead_proton2_theta_mu_p_;

  MyPointer<TVector3> mc_p3mu;
  MyPointer<TVector3> mc_p3p;
  MyPointer<std::vector<TVector3>> mc_p3_p_vec_;


  std::map<double, int> proton_index;
  std::map<double, double> proton_pid_score, proton_e, proton_px, proton_py, proton_pz;

  STVCalcType CalcType;
};

#endif
