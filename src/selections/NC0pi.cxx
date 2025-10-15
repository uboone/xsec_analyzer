// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/NC0pi.hh"

NC0pi::NC0pi() : SelectionBase( "NC0pi" ) {
}

std::string NC0pi::categorize_event( AnalysisEvent& ev ) {
  // Assign the event category of the selected event. Only MC truth
  // information should be used to determine the answer. This function
  // is not called by the framework when processing real data.
  return "Unknown";
}

bool NC0pi::is_selected( AnalysisEvent& ev ) {
  // Determine whether an input event satisfies the selection criteria.
  // Only reco information should be used to determine the answer so that
  // this function can be applied to both real data and MC events.
  return false;
}
