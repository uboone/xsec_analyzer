#include "DummySelection.h"

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

EventCategory DummySelection::CategorizeEvent(AnalysisEvent* Event) {
  //Define the event category of the selected event
  return kUnknown;
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

