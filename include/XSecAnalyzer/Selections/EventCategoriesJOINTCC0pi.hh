#pragma once

// Standard library includes
#include <map>
#include <string>

// Enum used to label event categories of interest for analysis plots in
// the CC0pi 1p/2p/Np/Xp analyses
enum EventCategoryJOINTCC0Pi {

  // Unable to categorize (e.g., because the event is real data and thus
  // has no MC truth information)
  kUnknown = 0,

  // Signal events broken down by underlying reaction mode
  kSignalCCQE = 1,
  kSignalCCMEC = 2,
  kSignalCCRES = 3,
  kSignalOther = 4,

  // True numu CC event with at least one final-state pion above threshold
  kNuMuCCNpi = 5,

  // True numu CC event with zero final-state pions above threshold and
  // zero final-state protons above threshold
  //  Will keep  kNuMuCC0pi0p enum but removed from categories and now this type is funnled into signal region
  
  // Any true numu CC event which does not satisfy the criteria for inclusion
  // in one of the other categories above
  kNuMuCCOther = 7,

  // True nue CC event
  kNuECC = 8,

  // True neutral current event for any neutrino flavor
  kNC = 9,

  // True neutrino vertex (any reaction mode and flavor combination) is outside
  // of the fiducial volume
  kOOFV = 10,

  // Seperate numubar CC from kOther
  //kNuMuBarCC = 13,
  // All events that do not fall within any of the other categories (e.g.,)
  kOther = 11,
  
  
};

extern std::map< int, std::pair< std::string, int > > JOINTCC0Pi_MAP;



