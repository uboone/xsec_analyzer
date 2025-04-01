// ROOT includes
#include "TH1.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/EventCategoriesNuMICC1e.hh"

std::map< int, std::pair< std::string, int > > NuMICC1e_MAP = {
  
  { kNuECC, {"#nu_{e} CC (Signal)", kRed + 1 } },
  { kNuECCOther, {"#nu_{e} CC Other", kRed - 6 } },
  { kNuMuCCPi0, { "#nu_{#mu} CC #pi^{0}", kBlue+3 } },
  { kNuMuCCOther, { "#nu_{#mu} CC Other", kBlue-3 } },
  { kNCPi0, { "NC", kMagenta+2 } },
  { kNCOther, { "NC", kMagenta-2 } },
  { kOOFV, {"Out FV", kOrange +1 } },
  { kEXT, { "Other", kGray + 2 } },
  { kUnknown, { "Unknown", kBlack } }
};