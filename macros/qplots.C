// XSecAnalyzer includes
#include "XSecAnalyzer/FilePropertiesManager.hh"
#include "XSecAnalyzer/QuickPlotter.hh"

// Only load the library in this way if we're using this code from inside the
// ROOT C++ interpreter. We could check for __CINT__ as well, but the specific
// R__LOAD_LIBRARY approach used here only works with ROOT 6.
#ifdef __CLING__
R__LOAD_LIBRARY(../lib/libXSecAnalyzer.so)
#endif

// Initialize the FilePropertiesManager with the samples of interest
// (currently uses hadd-ed ntuples from the April 2025 mini-retreat)
void initialize_file_properties( const std::string& sample_name ) {
  auto& fpm = FilePropertiesManager::Instance();
  std::string xsec_analyzer_dir = gSystem->Getenv( "XSEC_ANALYZER_DIR" );
  if ( sample_name == "mcc9" ) {
    fpm.load_file_properties( xsec_analyzer_dir
      + "/configs/mcc9_file_properties.txt" );
    return;
  }
  else if ( sample_name == "mcc9.10 Pandora" ) {
    fpm.load_file_properties( xsec_analyzer_dir
      + "/configs/file_properties.txt" );
    return;
  }
  else if ( sample_name == "mcc9.10 Super Unified" ) {
    fpm.load_file_properties( xsec_analyzer_dir
      + "/configs/super_unified_file_properties.txt" );
    return;
  }
  throw std::runtime_error( "Unrecognized sample name \""
    + sample_name + "\"" );
}

void qplots() {

  // **** Choose the sample you want to look at
  //initialize_file_properties( "mcc9" );
  //initialize_file_properties( "mcc9.10 Pandora" );
  initialize_file_properties( "mcc9.10 Super Unified" );

  // **** Initializes a "QuickPlotter" object to look at Run 4 results.
  // The N-proton selection is used to define the event categories.
  QuickPlotter qp( 4, "CC1muNp0pi" );

  // **** Make a plot using the data and MC ntuples
  //
  // 1st argument: Branch expression to plot (reco muon momentum in this case)
  // 2nd argument: Selection cuts to use (full CC 0-pion N-proton selection)
  // 3rd argument: Lower bound for the x-axis
  // 4th argument: Upper bound for the x-axis
  // 5th argument: Number of bins to use
  // 6th argument: x-axis caption
  // 7th argument: y-axis caption
  // 8th argument: Plot title
  qp.plot( "reco_p3_mu.Mag()", "CC1muNp0pi_Selected",
    0.1, 1.2, 40, "p_{#mu} [GeV/c]", "events", "Run 4b" );
}
