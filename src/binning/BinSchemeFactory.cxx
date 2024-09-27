// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/BinSchemeFactory.hh"
#include "XSecAnalyzer/Binning/TutorialBinScheme.hh"

BinSchemeFactory::BinSchemeFactory() {
}

BinSchemeBase* BinSchemeFactory::CreateBinScheme(
  const std::string& bin_scheme_name )
{
  BinSchemeBase* bs = nullptr;

  if ( bin_scheme_name == "TutorialBinScheme" ) {
    bs = new TutorialBinScheme;
    bs->Init();
  }
  else {
    std::cerr << "Binning scheme name requested: " << bin_scheme_name
      << " is not implemented in " << __FILE__ << '\n';
    throw;
  }

  return bs;
}
