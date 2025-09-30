#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/BinSchemeBase.hh"

class CC1muNp1piBinScheme : public BinSchemeBase {

  public:

    CC1muNp1piBinScheme();
    virtual void DefineBlocks() override;
};
