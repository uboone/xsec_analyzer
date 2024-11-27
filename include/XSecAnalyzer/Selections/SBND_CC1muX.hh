#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"
#include "XSecAnalyzer/Selections/SBND_EventCategories.hh"

class SBND_CC1muX : public SelectionBase {

public:

  SBND_CC1muX();

  virtual int categorize_event( AnalysisEvent* event ) override final;
  virtual bool selection( AnalysisEvent* event ) override final;
  virtual bool define_signal( AnalysisEvent* event ) override final;
  virtual void compute_reco_observables( AnalysisEvent* event ) override final;
  virtual void compute_true_observables( AnalysisEvent* event ) override final;
  virtual void define_output_branches() override final;
  virtual void define_constants() override final;
  virtual void define_category_map() override final;
  virtual void reset() override final;

private:

  double leading_muon_momentum;
  double leading_muon_costheta;
  double leading_muon_momentum_truth;
  double leading_muon_costheta_truth;

};
