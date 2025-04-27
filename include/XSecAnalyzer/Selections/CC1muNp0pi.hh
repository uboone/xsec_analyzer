#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class CC1muNp0pi : public SelectionBase {

public:

  CC1muNp0pi();

  virtual std::string categorize_event( AnalysisEvent& ev ) override final;
  virtual bool is_selected( AnalysisEvent& ev ) override final;

private:

  STVCalcType stv_calc_type_;
};
