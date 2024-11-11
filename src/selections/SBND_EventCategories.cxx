// ROOT includes
#include "TH1.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SBND_EventCategories.hh"

std::map< int, std::pair< std::string, int > > CC1muX_MAP = {
  { kUnknown, { "Unknown", kGray } },
  { kNuMuCC, { "NuMuCC", kBlue } },
  { kNC, { "NC", kOrange } },
  { kNuECC, { "NuECC", kViolet } },
  { kCosmic, { "Cosmic", kRed } },
  { kOOFV, { "OOFV", kRed } }
};
