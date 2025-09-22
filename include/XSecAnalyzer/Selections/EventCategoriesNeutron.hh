#pragma once

// Standard library includes
#include <map>
#include <string>

// Enum used to label event categories of interest for analysis plots in
// the CC0pi 1p/2p/Np/Xp analyses
enum EventCategoryNeutron {

  // Unable to categorize (e.g., because the event is real data and thus
  // has no MC truth information)
  kUnknown = 0,

  // Signal events broken down by underlying reaction mode
  kNuMuCC1n0pi0p_CCQE = 1,
  kNuMuCC1n0pi0p_CCMEC = 2,
  kNuMuCC1n0pi0p_CCRES = 3,
  kNuMuCC1n0pi0p_Other = 4,

  kNuMuCCNn0pi0p_CCQE = 5,
  kNuMuCCNn0pi0p_CCMEC = 6,
  kNuMuCCNn0pi0p_CCRES = 7,
  kNuMuCCNn0pi0p_Other = 8,

  kNuMuCC1n0piNp_CCQE = 9,
  kNuMuCC1n0piNp_CCMEC = 10,
  kNuMuCC1n0piNp_CCRES = 11,
  kNuMuCC1n0piNp_Other = 12,

  kNuMuCCNn0piNp_CCQE = 13,
  kNuMuCCNn0piNp_CCMEC = 14,
  kNuMuCCNn0piNp_CCRES = 15,
  kNuMuCCNn0piNp_Other = 16,

  //BEGIN BACKGROUNDS
  // True numu CC with no neutrons, no pions, and at least one proton above threshold
  kNuMuCC0n0piNp = 17,

  // True numu CC event with at least one final-state pion above threshold
  kNuMuCCNpi = 18,

  // Any true numu CC event which does not satisfy the criteria for inclusion
  // in one of the other categories above
  kNuMuCCOther = 19,

  // True nue CC event
  kNuECC = 20,

  // True neutral current event for any neutrino flavor
  kNC = 21,

  // True neutrino vertex (any reaction mode and flavor combination) is outside
  // of the fiducial volume
  kOOFV = 22,

  // All events that do not fall within any of the other categories (e.g.,
  // numubar CC)
  kOther = 23,
};

extern std::map< int, std::pair< std::string, int > > CC1muNnXp_MAP;

