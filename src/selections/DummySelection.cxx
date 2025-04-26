// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/DummySelection.hh"

DummySelection::DummySelection() : SelectionBase( "Dummy" ) {
}

void DummySelection::compute_reco_observables( AnalysisEvent& ev ) {
  // Calculate reconstructed kinematic variables to be saved in the output
}

void DummySelection::compute_true_observables( AnalysisEvent& ev ) {
  // Calculate true kinematic variables to be saved in the output
}

std::string DummySelection::categorize_event( AnalysisEvent& ev ) {
  // Assign the event category of the selected event
  return "Unknown";
}

bool DummySelection::is_signal( AnalysisEvent& ev ) {
  // Determine whether an input MC event fulfills the signal definition.
  // Only truth information should be used to determine the answer.
  return false;
}

bool DummySelection::is_selected( AnalysisEvent& ev ) {
  // Determine whether an input event satisfies the selection criteria.
  // Only reco information should be used to determine the answer.
  return false;
}
