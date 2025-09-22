// ROOT includes
#include "TH1.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/EventCategoriesXp.hh"

std::map< int, std::pair< std::string, int > > CC1muXp_MAP = {
  { kUnknown, { "Unknown", kGray } },
  { kNuMuCC0p0pi_CCQE, { "CCmu0p0pi (CCQE)", kBlue - 2 } },
  { kNuMuCC0p0pi_CCMEC, { "CCmu0p0pi (CCMEC)", kBlue - 6 } },
  { kNuMuCC0p0pi_CCRES, { "CCmu0p0pi (CCRES)", kBlue - 9 } },
  { kNuMuCC0p0pi_Other, { "CCmu0p0pi (Other)", kBlue - 10 } },
  { kNuMuCC1p0pi_CCQE, { "CCmu1p0pi (CCQE)", kOrange + 4 } },
  { kNuMuCC1p0pi_CCMEC, { "CCmu1p0pi (CCMEC)", kOrange + 5 } },
  { kNuMuCC1p0pi_CCRES, { "CCmu1p0pi (CCRES)", kOrange + 6 } },
  { kNuMuCC1p0pi_Other, { "CCmu1p0pi (Other)", kOrange + 7 } },
  { kNuMuCC2p0pi_CCQE, { "CCmu2p0pi (CCQE)", kCyan - 3 } },
  { kNuMuCC2p0pi_CCMEC, { "CCmu2p0pi (CCMEC)", kCyan - 6 } },
  { kNuMuCC2p0pi_CCRES, { "CCmu2p0pi (CCRES)", kCyan - 4 } },
  { kNuMuCC2p0pi_Other, { "CCmu2p0pi (Other)", kCyan - 9 } },
  { kNuMuCCMp0pi_CCQE, { "CCmuMp0pi (CCQE)", kGreen } },
  { kNuMuCCMp0pi_CCMEC, { "CCmuMp0pi (CCMEC)", kGreen + 1 } },
  { kNuMuCCMp0pi_CCRES, { "CCmuMp0pi (CCRES)", kGreen + 2 } },
  { kNuMuCCMp0pi_Other, { "CCmuMp0pi (Other)", kGreen + 3 } },
  { kNuMuCCNpi, { "#nu_{#mu} CCN#pi", kAzure - 2 } },
  { kNuMuCCOther, { "Other #nu_{#mu} CC", kAzure } },
  { kNuECC, { "#nu_{e} CC", kViolet } },
  { kNC, { "NC", kOrange } },
  { kOOFV, {"Out FV", kRed + 3 } },
  { kOther, { "Other", kRed + 1 } }
};
