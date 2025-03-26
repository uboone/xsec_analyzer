/**
 * @file    syst.h
 * @brief   Definitions of parameters to do systematic weights.
 * @details This file contains definitions of parameters to store systematic weights.
 * @author brindenc@fnal.gov
 */

#ifndef SYST_H
#define SYST_H

#include "sbnana/CAFAna/Core/ISyst.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/Cut.h"

#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Systs/UniverseOracle.h"

#include "slice_cuts.h"
#include "slice_variables.h"

using namespace ana;
namespace syst{
  //const Var kTrueLeadingMuonCostheta = SIMPLEVAR(slicevars::kTrueLeadingMuonCostheta);
  //const Var kTrueE = SIMPLEVAR(truth.E);

  // Get systematic weights and names
  const std::vector<std::string> genie_names = GetSBNGenieWeightNames();
  const std::vector<const ISyst*> genie_systs = GetSBNGenieWeightSysts();
  //const std::vector<const ISyst*> detsysts = GetDetectorSysts();
  
  // Define systematic names
  const std::vector<std::string> flux_names{ "expskin_Flux", "horncurrent_Flux", "nucleoninexsec_Flux", "nucleonqexsec_Flux", "nucleontotxsec_Flux", "pioninexsec_Flux", "pionqexsec_Flux", "piontotxsec_Flux", "piplus_Flux", "piminus_Flux", "kplus_Flux", "kminus_Flux", "kzero_Flux" };
  const std::vector<std::string> xsec_multisim_names{
    "GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_RPA_CCQE",
    "GENIEReWeight_SBN_v1_multisim_CoulombCCQE",
    "GENIEReWeight_SBN_v1_multisim_NormCCMEC",
    "GENIEReWeight_SBN_v1_multisim_NormNCMEC",
    "GENIEReWeight_SBN_v1_multisim_NCELVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi",
    "GENIEReWeight_SBN_v1_multisim_COHVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse",
    "GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse",
  };

}
#endif // SYST_H
