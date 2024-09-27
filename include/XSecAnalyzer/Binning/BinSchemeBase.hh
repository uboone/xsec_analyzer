#pragma once

// Standard library includes
#include <vector>

// XSecAnalyzerIncludes
#include "XSecAnalyzer/Binning/Block.hh"
#include "XSecAnalyzer/UniverseMaker.hh"

class BinSchemeBase {

  public:

    inline BinSchemeBase( const std::string& bin_scheme_name )
      : bin_scheme_name_( bin_scheme_name ) {}

    inline void Init() { this->DefineBlocks(); }

    virtual ~BinSchemeBase() = default;

    virtual void DefineBlocks() = 0;

    /// Name that uniquely identifies this binning scheme
    std::string bin_scheme_name_;

    /// Blocks of true + reco bin definitions
    std::vector< BlockTrueReco > vect_block;

    /// Name of the input TTree that appears in the post-processed ntuple files
    std::string ntuple_ttree_name_ = "stv_tree";

    /// The run numbers to use when plotting migration matrices
    std::set< int > runs_to_use_ = { 1 };

    /// Prefix for the output bin/slice configuration text files
    std::string out_config_prefix_;

    /// Name of the selection to use with this binning scheme
    std::string selection_name_;

    /// Name of the TDirectoryFile that will store the output histograms
    /// when univmake is run for this binning scheme
    std::string out_tdir_name_;
};
