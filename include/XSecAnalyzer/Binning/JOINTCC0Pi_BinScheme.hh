#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/BinSchemeBase.hh"

class JOINTCC0Pi_BinScheme : public BinSchemeBase {

  public:

    JOINTCC0Pi_BinScheme();
    virtual void DefineBlocks() override;
};
