// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/DummySelection.hh"

DummySelection::DummySelection() : SelectionBase( "Dummy" ) {
}

void DummySelection::compute_reco_observables( TreeHandler& th ) {
  // Calculate reconstructed kinematic variables to be saved in the output
}

void DummySelection::compute_true_observables( TreeHandler& th ) {
  // Calculate true kinematic variables to be saved in the output
}

const std::string DummySelection::categorize_event( TreeHandler& th ) {
  // Assign the event category of the selected event
  return "Unknown";
}

bool DummySelection::is_signal( TreeHandler& th ) {
  // Determine whether an input MC event fulfills the signal definition.
  // Only truth information should be used to determine the answer.
  return false;
}

bool DummySelection::is_selected( TreeHandler& th ) {
  // Determine whether an input event satisfies the selection criteria.
  // Only reco information should be used to determine the answer.
  return false;
}
