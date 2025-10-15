#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class NC0pi : public SelectionBase {

public:

  NC0pi();

  virtual bool is_selected( AnalysisEvent& ev ) override final;
  virtual std::string categorize_event( AnalysisEvent& ev ) override final;

};
