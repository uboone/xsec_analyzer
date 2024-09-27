#pragma once

// Standard library includes
#include <string>

// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/BinSchemeBase.hh"

class BinSchemeFactory {
  public:
    BinSchemeFactory();
    BinSchemeBase* CreateBinScheme( const std::string& BinSchemeName );
};
