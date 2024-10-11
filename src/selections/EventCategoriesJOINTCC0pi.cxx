// ROOT includes
#include "TH1.h"



// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/EventCategoriesJOINTCC0pi.hh"

std::map< int, std::pair< std::string, int > > JOINTCC0Pi_MAP = {

      { kUnknown, {"Unknown", kGray }},
      { kSignalCCQE, {"Signal (CCQE)" ,kGreen }},
      { kSignalCCMEC, {"Signal (CCMEC)",kGreen + 1  }},
      { kSignalCCRES, {"Signal (CCRES)",kGreen + 2  }},
      { kSignalOther, {"Signal (Other)",kGreen + 3  }},
      { kNuMuCCNpi, {"#nu_{#mu} CCN#pi",kAzure - 2 }},
      //{kNuMuBarCC ,{ "#bar{#nu_{#mu}}-CC"}},
      { kNuMuCCOther, {"Other #nu_{#mu} CC",kAzure }},
      { kNuECC, {"#nu_{e} CC",kViolet }},
      { kNC, {"NC",kOrange }},
      { kOOFV, {"Out FV",kRed + 3 }},
      { kOther, {"Other",kRed + 1 }}
  
  
};
