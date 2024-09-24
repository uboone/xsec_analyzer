#ifndef __DUMMYSELECTION_h__
#define __DUMMYSELECTION_h__

#include "SelectionBase.h"

class DummySelection : virtual SelectionBase {
 public:
  DummySelection();

  EventCategory CategorizeEvent(AnalysisEvent* Event);
  bool Selection(AnalysisEvent* Event);
  bool DefineSignal(AnalysisEvent* Event);
  void ComputeRecoObservables(AnalysisEvent* Event);
  void ComputeTrueObservables(AnalysisEvent* Event);
  void DefineOutputBranches();
  void DefineConstants();
  
private:
};

#endif
