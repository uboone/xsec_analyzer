#pragma once

// Standard library includes
#include <map>
#include <string>

// Enum used to label event categories of interest for analysis plots in
// the CC0pi 1p/2p/Np/Xp analyses
enum SBND_EventCategory {

  // Unable to categorize (e.g., because the event is real data and thus
  // has no MC truth information)
  kUnknown = -1,

  // Inclusive numu CC event
  kNuMuCC = 0,
  kNC = 1,
  kNuECC = 2,
  kCosmic = 3,
  kOOFV = 4,

};

extern std::map< int, std::pair< std::string, int > > CC1muX_MAP;
