// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/DummySelection.hh"

DummySelection::DummySelection() : SelectionBase( "DummySel" ) {
}

void DummySelection::define_constants() {
  // Define reco and true fiducial volumes alongside any other constants used
  // within selection cuts
}

void DummySelection::compute_reco_observables( AnalysisEvent* event ) {
  // Calculate reconstructed kinematic variables to be saved in the output
}

void DummySelection::compute_true_observables( AnalysisEvent* event ) {
  // Calculate true kinematic variables to be saved in the output
}

int DummySelection::categorize_event( AnalysisEvent* event ) {
  // Assign the event category of the selected event
  return 0;
}

bool DummySelection::define_signal(AnalysisEvent* event) {
  // Determine whether an input MC event fulfills the signal definition.
  // Only truth information should be used to determine the answer.
  return false;
}

bool DummySelection::selection( AnalysisEvent* event ) {
  // Determine whether an input event satisfies the selection criteria.
  // Only reco information should be used to determine the answer.
  return false;
}

void DummySelection::define_output_branches() {
  // Call set_branch() for every new variable to be saved to the output TTree
}

void DummySelection::reset() {
  // Set variables managed by this class to their default values. This function
  // is called in preparation for analyzing each new input event.
}

void DummySelection::define_category_map() {
  // Define the mapping between each integer event category and
  // a string label / color integer code pair
  // The color codes are documented at
  // https://root.cern.ch/doc/master/classTColor.html
  std::map< int, std::pair< std::string, int > >
    temp_map = { { 0, { "Unknown", 0 } } };
}
