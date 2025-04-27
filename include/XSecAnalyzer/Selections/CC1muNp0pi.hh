#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"
#include "XSecAnalyzer/KICalculator.hh"

class CC1muNp0pi : public SelectionBase {

public:

  CC1muNp0pi();

  virtual std::string categorize_event( AnalysisEvent& ev ) override final;
  virtual bool is_selected( AnalysisEvent& ev ) override final;

private:

  KICalculator::CalcType stv_calc_type_;
};
