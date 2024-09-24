// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/DummySelection.hh"

DummySelection::DummySelection() : SelectionBase("DummySel") {
}

void DummySelection::DefineConstants() {
  //Define Reco&True FV, alongside any other constants used within selection cuts
}

void DummySelection::ComputeRecoObservables(AnalysisEvent* Event) {
  //Define the reconstructed kinematic variables the xsec measurement will be provided in/of interest
}

void DummySelection::ComputeTrueObservables(AnalysisEvent* Event) {
  //Define the true kinematic variables of interest
}

int DummySelection::CategorizeEvent(AnalysisEvent* Event) {
  //Define the event category of the selected event
  return 0;
}

bool DummySelection::DefineSignal(AnalysisEvent* Event) {
  //Define the MC True signal definition
  return false;
}

bool DummySelection::Selection(AnalysisEvent* Event) {
  //Apply the selection cuts on Reco variables
  return false;
}

void DummySelection::DefineOutputBranches() {
  //Save any additional variables to output TTree
}

void DummySelection::DefineCategoryMap() {
  // Define the mapping between each integer event category and
  // a string label / color integer code pair
  std::map< int, std::pair< std::string, int > >
    temp_map = { { 0, { "Unknown", 0 } } };
}
