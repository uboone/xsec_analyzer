#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class DummySelection : public SelectionBase {
 public:
  DummySelection();

  int CategorizeEvent(AnalysisEvent* Event);
  bool Selection(AnalysisEvent* Event);
  bool DefineSignal(AnalysisEvent* Event);
  void ComputeRecoObservables(AnalysisEvent* Event);
  void ComputeTrueObservables(AnalysisEvent* Event);
  void DefineOutputBranches();
  void DefineConstants();
  void DefineCategoryMap();

private:
};
