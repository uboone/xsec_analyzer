#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/BinSchemeBase.hh"

class SBND_CC1muX_BinScheme : public BinSchemeBase {

  public:

    SBND_CC1muX_BinScheme();
    virtual void DefineBlocks() override;
};
