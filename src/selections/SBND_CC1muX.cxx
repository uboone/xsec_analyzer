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
  pmu = event->leading_muon_momentum_;
  cos_theta_mu = event->leading_muon_costheta_;
}

void SBND_CC1muX::compute_true_observables( AnalysisEvent* event ) {
  // Calculate true kinematic variables to be saved in the output
  mc_pmu = event->leading_muon_momentum_truth_;
  mc_cos_theta_mu = event->leading_muon_costheta_truth_;
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
  return event->is_signal_;
}

void SBND_CC1muX::define_output_branches() {
  // Call set_branch() for every new variable to be saved to the output TTree
  set_branch( &pmu, "pmu" );
  set_branch( &cos_theta_mu, "cos_theta_mu" );
  set_branch( &mc_pmu, "mc_pmu" );
  set_branch( &mc_cos_theta_mu, "mc_cos_theta_mu" );
}

void SBND_CC1muX::reset() {
  // Set variables managed by this class to their default values. This function
  // is called in preparation for analyzing each new input event.
  pmu = BOGUS;
  cos_theta_mu = BOGUS;
  mc_pmu = BOGUS;
  mc_cos_theta_mu = BOGUS;
}

void SBND_CC1muX::define_category_map() {
  // Define the mapping between each integer event category and
  // a string label / color integer code pair
  categ_map_ = CC1muX_MAP;
}
