#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class DummySelection : public SelectionBase {

public:

  DummySelection();

  virtual bool is_selected( TreeHandler& th ) override final;
  virtual bool is_signal( TreeHandler& th ) override final;
  virtual const std::string categorize_event( TreeHandler& th ) override final;
  virtual void compute_reco_observables( TreeHandler& th ) override final;
  virtual void compute_true_observables( TreeHandler& th ) override final;

};
