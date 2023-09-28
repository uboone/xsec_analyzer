#ifndef __CC1mu2p0pi_h__
#define __CC1mu2p0pi_h__

#include "SelectionBase.h"

class CC1mu2p0pi : virtual SelectionBase {
 public:
  CC1mu2p0pi();
  
  bool Selection(AnalysisEvent* Event);
  EventCategory CategorizeEvent(AnalysisEvent* Event);
  void ComputeRecoObservables(AnalysisEvent* Event);
  void ComputeTrueObservables(AnalysisEvent* Event);
  void DefineOutputBranches();
  bool DefineSignal(AnalysisEvent* Event);
  void DefineConstants();
  
private:
  bool sel_reco_vertex_in_FV_;
  bool sel_has_muon_candidate_;
  bool sel_nu_mu_cc_;
  bool sel_npfps_eq_3;
  bool sel_ntracks_eq_3;
  bool sel_containedparticles;
  bool sel_correctparticles;
  bool sel_momentum_threshold_passed_;
  
  int muon_candidate_idx_;
};

#endif
