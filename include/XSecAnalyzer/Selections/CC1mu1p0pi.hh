#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/KICalculator.hh"
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class CC1mu1p0pi : public SelectionBase {

public:

  CC1mu1p0pi();

  virtual std::string categorize_event( AnalysisEvent& event ) override final;
  virtual bool is_selected( AnalysisEvent& event ) override final;

private:

  KICalculator::CalcType ki_calc_type_;

};
