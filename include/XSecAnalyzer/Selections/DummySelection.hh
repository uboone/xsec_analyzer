#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class DummySelection : public SelectionBase {

public:

  DummySelection();

  virtual bool is_selected( AnalysisEvent& ev ) override final;
  virtual std::string categorize_event( AnalysisEvent& ev ) override final;

};
