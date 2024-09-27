#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/BinSchemeBase.hh"

class TutorialBinScheme : public BinSchemeBase {

  public:

    TutorialBinScheme();
    virtual void DefineBlocks() override;
};
