// ROOT includes
#include "TH1.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/EventCategoriesXpXn.hh"

std::map< int, std::pair< std::string, int > > CC1muXpXn_MAP = {
  { kUnknown, { "Unknown", kGray } },
  { kNuMuCC0p0pi0n, { "CCmu0p0pi0n", kMagenta} },
  { kNuMuCC1p0pi0n, { "CCmu1p0pi0n", kMagenta + 2 } },
  { kNuMuCCNp0pi0n, { "CCmuNp0pi0n", kMagenta + 4 } },
  { kNuMuCC0p0pi1n, { "CCmu0p0pi1n", kGreen} },
  { kNuMuCC1p0pi1n, { "CCmu1p0pi1n", kGreen + 2 } },
  { kNuMuCCNp0pi1n, { "CCmuNp0pi1n", kGreen + 4 } },
  { kNuMuCC0p0piNn, { "CCmu0p0piNn", kBlue - 4} },
  { kNuMuCC1p0piNn, { "CCmu1p0piNn", kBlue} },
  { kNuMuCCNp0piNn, { "CCmuNp0piNn", kBlue + 2} },
  { kNuMuCCNpi, { "#nu_{#mu} CCN#pi", kOrange -4 } },
  { kNuMuCCOther, { "Other #nu_{#mu} CC", kOrange -3 } },
  { kNuECC, { "#nu_{e} CC", kRed - 7 } },
  { kNC, { "NC", kOrange } },
  { kOOFV, {"Out FV", kRed + 3 } },
  { kOther, { "Other", kRed + 1 } }
};
