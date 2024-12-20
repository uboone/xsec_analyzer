/*
Description: Factory function to create appropriate event based on experiment
Date: 2024-10-16
Authors: Brinden Carlson (bcarlson1@ufl.edu)
*/

#include "AnalysisEvent.hh"
#include "SBNDEvent.hh"
#include "uBooNEEvent.hh"

// Factory function to create appropriate event based on experiment type
std::unique_ptr<AnalysisEvent> create_event(const std::string& experiment) {
  if (experiment == "sbnd") {
    return std::make_unique<SBNDEvent>();
  }
  else if (experiment == "uboone") {
    return std::make_unique<uBooNEEvent>();
  }
  // Add more experiments here with else if
  std::cerr << "Invalid experiment name: " << experiment << std::endl;
  return nullptr; // or throw an error
}