#ifndef __CC1muNnXp0pi_h__
#define __CC1muNnXp0pi_h__


#include "XSecAnalyzer/Selections/SelectionBase.hh"

#define DEFAULT_NUM 20

class CC1muNnXp0pi : virtual SelectionBase {
 public:
  CC1muNnXp0pi();

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

  //Begin Burke Additions
  double mc_lead_neutron_energy_;
  double mc_lead_neutron_costh_;
  double mc_lead_neutron_px_;
  double mc_lead_neutron_py_;
  double mc_lead_neutron_p_;
  double mc_lead_neutron_phi_;
  double mc_lead_neutron_theta_;

  bool sig_is_signal_;
  bool sig_mc_FS_neutron_;
  int sig_num_neutron_above_100MeV;
  bool sel_has_secondary_proton_candidate_;
  int sel_num_secondary_proton_cands;
  std::map<int, double> pfp_proximity_map;
  std::map<int, double> pfp_direction_map;
  std::map<int, double> pfp_vtx_disp_map;
  std::map<int, double> pfp_len_map;
  std::map<int, double> pfp_cos_theta_map;
  std::map<int, double> pfp_theta_map;
  std::map<int, double> pfp_phi_map;
  std::map<int, double> pfp_start_x_map;
  std::map<int, double> pfp_start_y_map;
  std::map<int, double> pfp_start_z_map;
  std::map<int, double> pfp_end_x_map;
  std::map<int, double> pfp_end_y_map;
  std::map<int, double> pfp_end_z_map;
  //MyPointer<std::vector<int>> secondary_proton_index;
  int secondary_proton_trk_candidate_indices[DEFAULT_NUM];
  float secondary_proton_proximity[DEFAULT_NUM];
  float secondary_proton_direction[DEFAULT_NUM];
  float secondary_proton_vtx_disp[DEFAULT_NUM];
  float secondary_proton_len[DEFAULT_NUM];
  float secondary_proton_costheta[DEFAULT_NUM];
  float secondary_proton_theta[DEFAULT_NUM];
  float secondary_proton_phi[DEFAULT_NUM];
  float secondary_proton_startx[DEFAULT_NUM];
  float secondary_proton_starty[DEFAULT_NUM];
  float secondary_proton_startz[DEFAULT_NUM];
  float secondary_proton_endx[DEFAULT_NUM];
  float secondary_proton_endy[DEFAULT_NUM];
  float secondary_proton_endz[DEFAULT_NUM];

  double event_missing_trans_px;
  double event_missing_trans_py;
  double event_missing_trans_p_mag;
  double event_missing_p_L;
  double event_missing_p_mag;
  double sel_cos_alpha_trans;
  double sel_alpha_trans;
  double sel_cos_alpha_3D;
  double sel_alpha_3D;
  double event_missing_trans_px_leading_p;
  double event_missing_trans_py_leading_p;
  double event_missing_trans_p_mag_leading_p;
  double event_missing_p_L_leading_p;
  double event_missing_p_mag_leading_p;
  double sel_cos_alpha_trans_leading_p;
  double sel_alpha_trans_leading_p;
  double sel_cos_alpha_3D_leading_p;
  double sel_alpha_3D_leading_p;

  STVCalcType CalcType;
};

#endif
