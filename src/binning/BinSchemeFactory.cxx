// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/BinSchemeFactory.hh"
#include "XSecAnalyzer/Binning/TutorialBinScheme.hh"
#include "XSecAnalyzer/Binning/CCXp0piBinScheme.hh"
#include "XSecAnalyzer/Binning/JOINTCC0Pi_BinScheme.hh"

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
  else if(bin_scheme_name == "CCXp0piBinScheme"){
    bs = new CCXp0piBinScheme;
    bs->Init();
  }
    else if(bin_scheme_name == "JOINTCC0Pi_BinScheme"){
    bs = new JOINTCC0Pi_BinScheme;
    bs->Init();
  }
  else {
    std::cerr << "Binning scheme name requested: " << bin_scheme_name
      << " is not implemented in " << __FILE__ << '\n';
    throw;
  }

  return bs;
}
