#pragma once

// Standard library includes
#include <map>
#include <string>

// Enum used to label event categories of interest for analysis plots in
// the CC0pi 1p/2p/Np/Xp analyses
enum EventCategoryNp1pi {

  // Unable to categorize (e.g., because the event is real data and thus
  // has no MC truth information)
  kUnknown = 0, //

  // Signal events broken down by underlying reaction mode
  kSignalCCQE = 1, //
  kSignalCCMEC = 2, //
  kSignalCCRES = 3, //
  kSignalCCDIS = 4, //
  kSignalCCCOH = 5, //
  kSignalOther = 6, //

  // True numu CC event with pi0s
  kNuMuCCpi0 = 7,

  // True numu CC event with no final-state pions above threshold and
  // at least one proton above threshold
  kNuMuCC0piXp = 8,


  // Any true numu CC event which does not satisfy the criteria for inclusion
  // in one of the other categories above
  kNuMuCCOther = 9,

  // True nue CC event
  kNuECC = 10, //

  // True neutral current event for any neutrino flavor
  kNC = 11, //

  // True neutrino vertex (any reaction mode and flavor combination) is outside
  // of the fiducial volume
  kOOFV = 12, //

  // All events that do not fall within any of the other categories (e.g.,
  // numubar CC)
  kOther = 13 //

};

extern std::map< int, std::pair< std::string, int > > CC1muNp1pi_MAP;
