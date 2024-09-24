#pragma once

// Standard library includes
#include <map>
#include <string>

// Enum used to label event categories of interest for analysis plots in
// the CC0pi 1p/2p/Np/Xp analyses
enum EventCategoryXp {

  // Unable to categorize (e.g., because the event is real data and thus
  // has no MC truth information)
  kUnknown = 0,

  // Signal events broken down by underlying reaction mode
  kNuMuCC0p0pi_CCQE = 1,
  kNuMuCC0p0pi_CCMEC = 2,
  kNuMuCC0p0pi_CCRES = 3,
  kNuMuCC0p0pi_Other = 4,

  kNuMuCC1p0pi_CCQE = 5,
  kNuMuCC1p0pi_CCMEC = 6,
  kNuMuCC1p0pi_CCRES = 7,
  kNuMuCC1p0pi_Other = 8,

  kNuMuCC2p0pi_CCQE = 9,
  kNuMuCC2p0pi_CCMEC = 10,
  kNuMuCC2p0pi_CCRES = 11,
  kNuMuCC2p0pi_Other = 12,

  // M = >2
  kNuMuCCMp0pi_CCQE = 13,
  kNuMuCCMp0pi_CCMEC = 14,
  kNuMuCCMp0pi_CCRES = 15,
  kNuMuCCMp0pi_Other = 16,

  // True numu CC event with at least one final-state pion above threshold
  kNuMuCCNpi = 17,

  // Any true numu CC event which does not satisfy the criteria for inclusion
  // in one of the other categories above
  kNuMuCCOther = 18,

  // True nue CC event
  kNuECC = 19,

  // True neutral current event for any neutrino flavor
  kNC = 20,

  // True neutrino vertex (any reaction mode and flavor combination) is outside
  // of the fiducial volume
  kOOFV = 21,

  // All events that do not fall within any of the other categories (e.g.,
  // numubar CC)
  kOther = 22,
};

extern std::map< int, std::pair< std::string, int > > CC1muXp_MAP;
