#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/BinSchemeBase.hh"

class CCXp0piBinScheme : public BinSchemeBase {

  public:

    CCXp0piBinScheme();
    virtual void DefineBlocks() override;
};
