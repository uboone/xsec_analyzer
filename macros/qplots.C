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

  // Choose the sample you want to look at
  initialize_file_properties( "mcc9" );
  //initialize_file_properties( "mcc9.10 Pandora" );
  //initialize_file_properties( "mcc9.10 Super Unified" );

  QuickPlotter qp( 4, "CC1muNp0pi" );

  qp.plot( "topological_score",
    "CC1muNp0pi_reco_vertex_in_FV && CC1muNp0pi_pfp_starts_in_PCV",
    0., 1., 40, "topological score", "events", "Run 4b" );
}
