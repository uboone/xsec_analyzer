// ROOT includes
#include "TH1.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/EventCategoriesNp1pi.hh"

std::map< int, std::pair< std::string, int > > CC1muNp1pi_MAP = {
  { kUnknown, { "Unknown", kGray } },
  { kSignalCCQE, { "Signal (CCQE)", kGreen } },
  { kSignalCCMEC, { "Signal (CCMEC)", kGreen + 1 } },
  { kSignalCCRES, { "Signal (CCRES)", kGreen + 2 } },
  { kSignalCCDIS, { "Signal (CCDIS)", kGreen + 3 } },
  { kSignalCCCOH, { "Signal (CCCOH)", kGreen + 4 } },
  { kSignalOther, { "Signal (Other)", kGreen + 5 } },
  { kNuMuCCpi0, { "#nu_{#mu} CCN#pi^{0}", kAzure - 2 } },
  { kNuMuCC0piXp, { "#nu_{#mu} CC0#piXp", kAzure - 1 } },
  { kNuMuCCOther, { "Other #nu_{#mu} CC", kAzure + 3 } },
  { kNuECC, { "#nu_{e} CC", kViolet } },
  { kNC, { "NC", kOrange } },
  { kOOFV, { "OOFV", kRed + 3 } },
  { kOther, { "Other", kRed + 1 } }
};
