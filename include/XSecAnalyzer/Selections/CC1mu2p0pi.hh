#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"
#include "XSecAnalyzer/KICalculator.hh"

class CC1mu2p0pi : public SelectionBase {

public:

  CC1mu2p0pi();

  virtual std::string categorize_event( AnalysisEvent& ev ) override final;
  virtual bool is_selected( AnalysisEvent& ev ) override final;

private:

  KICalculator::CalcType ki_calc_type_;
};
