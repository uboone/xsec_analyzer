#pragma once

// Needed for the weight range limits used in DEFAULT_MC_EVENT_WEIGHT. We pull
// them from the UniverseMaker header file to ensure consistency
// with the safe_weight() function.
#include "UniverseMaker.hh"

// **** Helper code to facilitate making histograms ****

// By default, weight the MC events using the MicroBooNE CV tune. The
// TTree::Draw() expression below includes some workarounds for problematic
// weights. They are equivalent to the safe_weight() function defined in
// UniverseMaker.hh.
const std::string DEFAULT_MC_EVENT_WEIGHT = "(std::isfinite("
  "weightSplineTimesTune) && weightSplineTimesTune >= "
  + std::to_string( MIN_WEIGHT ) + " && weightSplineTimesTune"
  " <= " + std::to_string( MAX_WEIGHT )
  + " ? weightSplineTimesTune : 1)";

// Generates a vector of bin low edges equivalent to the approach used by the
// TH1 constructor that takes xmin and xmax in addition to the number of bins
std::vector<double> get_bin_low_edges( double xmin, double xmax, int Nbins );
