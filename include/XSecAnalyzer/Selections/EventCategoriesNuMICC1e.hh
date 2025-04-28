#pragma once

// Standard library includes
#include <map>
#include <string>

// Enum used to label event categories of interest for analysis plots
enum EventCategoryNuMICC1e {

  // True electron (anti)neutrino CC event
  // signal
  kNuECC = 0,
  // other, non-signal
  kNuECCOther = 1,

  // True muon (anti)neutrino CC event
  kNuMuCCPi0 = 2,
  kNuMuCCOther = 3,

  // True neutral current event for any neutrino flavor
  kNCPi0 = 4,
  kNCOther = 5,

  // True neutrino vertex (any reaction mode and flavor combination) is outside
  // of the fiducial volume
  kOOFV = 6,

  // EXT
  kEXT = 7,

  // all events that do not fall within any of the other categories
  // should not be populated
  kUnknown = 8,
};

extern std::map< int, std::pair< std::string, int > > NuMICC1e_MAP;