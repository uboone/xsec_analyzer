#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

#define DEFAULT_NUM 20

class CC1muNp0piNn : virtual SelectionBase {
 public:
  CC1muNp0piNn();

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
  bool sig_mc_no_charged_pi_above_threshold_;
  int sig_nProtons_in_Momentum_range;

  //Begin Burke Edit
  int sig_num_neutron_above_100MeV;
  bool sig_mc_FS_neutron_;
  int sel_num_secondary_proton_cands;
  bool sel_has_secondary_proton_candidate_;
  //End Burke Edit

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

  int lead_p_candidate_idx_;
  int muon_candidate_idx_;

  double delta_pT_;
  double delta_phiT_;
  double delta_alphaT_;
  double delta_pL_;
  double pn_;
  double delta_pTx_;
  double delta_pTy_;
  double theta_mu_p_;

  MyPointer<TVector3> p3mu;
  MyPointer<TVector3> p3p;
  MyPointer<std::vector<TVector3>> p3_p_vec_;

  double mc_delta_pT_;
  double mc_delta_phiT_;
  double mc_delta_alphaT_;
  double mc_delta_pL_;
  double mc_pn_;
  double mc_delta_pTx_;
  double mc_delta_pTy_;
  double mc_theta_mu_p_;

  MyPointer<TVector3> mc_p3mu;
  MyPointer<TVector3> mc_p3p;
  MyPointer<std::vector<TVector3>> mc_p3_p_vec_;

  //Begin Burke Additions

  MyPointer<std::vector<double>> weight_neutron_argon_xsec;
  MyPointer<std::vector<double>> weight_neutron_reint;
  MyPointer<std::vector<double>> neutron_track_weights;
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
