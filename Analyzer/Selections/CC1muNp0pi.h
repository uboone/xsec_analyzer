#ifndef __CC1muNp0pi_h__
#define __CC1muNp0pi_h__

#include "SelectionBase.h"

class CC1muNp0pi : virtual SelectionBase {
 public:
  CC1muNp0pi();
  
  bool Selection(AnalysisEvent* Event);
  void ComputeObservables(AnalysisEvent* Event);
  void DefineBranches();
  bool DefineSignal(AnalysisEvent* Event);
  void DefineConstants();
  
private:
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
  
  int lead_p_candidate_idx_;
  int muon_candidate_idx_;
};

#endif
