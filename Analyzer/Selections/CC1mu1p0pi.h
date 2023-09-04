#ifndef __CC1mu1p0pi_h__
#define __CC1mu1p0pi_h__

#include "SelectionBase.h"

class CC1mu1p0pi : virtual SelectionBase {
 public:
  CC1mu1p0pi();
  
  bool Selection(AnalysisEvent* Event);
  void ComputeObservables(AnalysisEvent* Event);
  void DefineBranches();
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
};

#endif
