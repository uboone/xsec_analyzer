#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class CC1muXpWC : public SelectionBase {

public:

  CC1muXpWC();

  virtual std::string categorize_event( AnalysisEvent& ev ) override final;
  virtual bool is_selected( AnalysisEvent& ev ) override final;

};
