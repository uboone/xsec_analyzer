#pragma once

// Standard library includes
#include <map>
#include <string>

// Enum used to label event categories of interest for analysis plots in
// the transverse kin imbalance analysis with neutron information
enum EventCategoryXpXn {

  // Unable to categorize (e.g., because the event is real data and thus
  // has no MC truth information)
  kUnknown = 0,

  // Signal events broken down by underlying reaction mode
  kNuMuCC0p0pi0n = 1,
  kNuMuCC1p0pi0n = 2,
  kNuMuCCNp0pi0n = 3,

  //1 final state neutron above 100 MeV
  kNuMuCC0p0pi1n = 4,
  kNuMuCC1p0pi1n = 5,
  kNuMuCCNp0pi1n = 6,

  //>1 final state neutron above 100 MeV 
  kNuMuCC0p0piNn = 7,
  kNuMuCC1p0piNn = 8,
  kNuMuCCNp0piNn = 9,

  // True numu CC event with at least one final-state pion above threshold
  kNuMuCCNpi = 10,

  // Any true numu CC event which does not satisfy the criteria for inclusion
  // in one of the other categories above
  kNuMuCCOther = 11,

  // True nue CC event
  kNuECC = 12,

  // True neutral current event for any neutrino flavor
  kNC = 13,

  // True neutrino vertex (any reaction mode and flavor combination) is outside
  // of the fiducial volume
  kOOFV = 14,

  // All events that do not fall within any of the other categories (e.g.,
  // numubar CC)
  kOther = 15,
};

extern std::map< int, std::pair< std::string, int > > CC1muXpXn_MAP;
