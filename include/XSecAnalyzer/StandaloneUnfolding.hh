#pragma once

// Standard library includes
#include <memory>
#include <set>
#include <stdexcept>
#include <string>

// ROOT includes
#include "TFile.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Unfolder.hh"

class StandaloneUnfolding {

  public:

    StandaloneUnfolding( const std::string& config_file_name );

    void run_unfolding() const;

  protected:

    std::unique_ptr< TFile > input_root_file_;
    std::unique_ptr< TFile > output_root_file_;
    std::unique_ptr< Unfolder > unfolder_;
    std::string blocks_file_;
};
