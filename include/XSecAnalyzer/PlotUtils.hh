#pragma once

// Standard library includes
#include <map>
#include <string>
#include <vector>

// ROOT includes
#include "TGraph.h"
#include "TH1D.h"

/// Helper function that produces the standard MicroBooNE plot legend title
/// with the BNB POT displayed
std::string get_legend_title( double bnb_pot );

/// Helper function that dumps 1D histogram contents to a map of pgfplotstable
/// columns
void dump_1d_histogram( const std::string& hist_col_prefix,
  const TH1D& hist,
  std::map< std::string, std::vector<std::string> >& pgf_plots_hist_table,
  bool include_yerror = true, bool include_x_coords = false );

/// Helper function that dumps TGraph contents to a map of pgfplotstable columns
// TODO: add support for TGraphErrors
void dump_tgraph( const std::string& col_prefix, const TGraph& graph,
  std::map< std::string, std::vector<std::string> >& pgf_plots_table );

void write_pgfplots_file( const std::string& out_filename,
  std::map< std::string, std::vector<std::string> >& pgfplots_table );
