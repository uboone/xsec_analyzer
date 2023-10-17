#ifndef __CC1mu1p0pi_h__
#define __CC1mu1p0pi_h__

#include "SelectionBase.h"

class CC1mu1p0pi : virtual SelectionBase {
 public:
  CC1mu1p0pi();

  bool Selection(AnalysisEvent* Event);
  EventCategory CategorizeEvent(AnalysisEvent* Event);
  void ComputeRecoObservables(AnalysisEvent* Event);
  void ComputeTrueObservables(AnalysisEvent* Event);
  void DefineOutputBranches();
  bool DefineSignal(AnalysisEvent* Event);
  void DefineConstants();
  
private:
  bool sel_nslice_eq_1_;
  bool sel_nshower_eq_0_;
  bool sel_ntrack_eq_2_;
  bool sel_muoncandidate_tracklike_;
  bool sel_protoncandidate_tracklike_;
  bool sel_nuvertex_contained_;
  bool sel_muoncandidate_above_p_thresh;
  bool sel_protoncandidate_above_p_thresh;
  bool sel_muoncandidate_contained;
  bool sel_protoncandidate_contained;
  bool sel_muon_momentum_quality;
  bool sel_no_flipped_tracks_;
  bool sel_proton_cand_passed_LLRCut;
  bool sel_muon_momentum_in_range;
  bool sel_muon_costheta_in_range;
  bool sel_muon_phi_in_range;
  bool sel_proton_momentum_in_range;
  bool sel_proton_costheta_in_range;
  bool sel_proton_phi_in_range;

  bool sig_truevertex_in_fv_;
  bool sig_ccnc_;
  bool sig_is_numu_;
  bool sig_one_muon_above_thresh_;
  bool sig_one_proton_above_thresh_;
  bool sig_no_pions_;
  bool sig_no_heavy_mesons_;
  int sig_mc_n_threshold_muon;
  int sig_mc_n_threshold_proton;
  int sig_mc_n_threshold_pion0;
  int sig_mc_n_threshold_pionpm;
  int sig_mc_n_heaviermeson;
  
  int CandidateMuonIndex;
  int CandidateProtonIndex;
};

#endif
