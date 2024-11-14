// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SBND_CC1muX.hh"

SBND_CC1muX::SBND_CC1muX() : SelectionBase( "SBND_CC1muX" ) {
}

void SBND_CC1muX::define_constants() {
  // Define reco and true fiducial volumes alongside any other constants used
  // within selection cuts
}

void SBND_CC1muX::compute_reco_observables( AnalysisEvent* event ) {
  // Calculate reconstructed kinematic variables to be saved in the output
  leading_muon_momentum = event->leading_muon_momentum_;
  leading_muon_costheta = event->leading_muon_costheta_;
}

void SBND_CC1muX::compute_true_observables( AnalysisEvent* event ) {
  // Calculate true kinematic variables to be saved in the output
  leading_muon_momentum_truth = event->leading_muon_momentum_truth_;
  leading_muon_costheta_truth = event->leading_muon_costheta_truth_;
}

int SBND_CC1muX::categorize_event( AnalysisEvent* event ) {
  // Assign the event category of the selected event
  return event->event_type_;
}

bool SBND_CC1muX::define_signal(AnalysisEvent* event) {
  // Determine whether an input MC event fulfills the signal definition.
  // Only truth information should be used to determine the answer.
  return event->event_type_ == 0;
}

bool SBND_CC1muX::selection( AnalysisEvent* event ) {
  // Determine whether an input event satisfies the selection criteria.
  // Only reco information should be used to determine the answer.
  //return event->is_signal_;
  return static_cast<int>(event->is_signal_);
}

void SBND_CC1muX::define_output_branches() {
  // Call set_branch() for every new variable to be saved to the output TTree
  set_branch( &leading_muon_momentum, "leading_muon_momentum" );
  set_branch( &leading_muon_costheta, "leading_muon_costheta" );
  set_branch( &leading_muon_momentum_truth, "leading_muon_momentum_truth" );
  set_branch( &leading_muon_costheta_truth, "leading_muon_costheta_truth" );
}

void SBND_CC1muX::reset() {
  // Set variables managed by this class to their default values. This function
  // is called in preparation for analyzing each new input event.
  leading_muon_momentum = BOGUS;
  leading_muon_costheta = BOGUS;
  leading_muon_momentum_truth = BOGUS;
  leading_muon_costheta_truth = BOGUS;
}

void SBND_CC1muX::define_category_map() {
  // Define the mapping between each integer event category and
  // a string label / color integer code pair
  categ_map_ = CC1muX_MAP;
}
