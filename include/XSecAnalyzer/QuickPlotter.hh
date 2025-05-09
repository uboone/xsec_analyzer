#pragma once

// Standard library includes
#include <memory>
#include <set>
#include <string>

// XSecAnalyzer includes
#include "XSecAnalyzer/HistUtils.hh"

// Forward-declare the SelectionBase class
class SelectionBase;

/// Implements a simple interface for quick data/MC comparisons using the
/// post-processed ntuples. Lacks a treatment of systematic uncertainties
/// since evaluating them requires processing event weights via univmake.
class QuickPlotter {

  public:

    /// Standard constructor that specifies a set of runs and a selection
    /// name used to assign MC event categories
    QuickPlotter( const std::set< int >& runs, const std::string& sel_name );

    /// Overloaded version that uses a single run
    inline QuickPlotter( int run, const std::string& sel_name )
      : QuickPlotter( std::set< int >( {run} ), sel_name ) {}

    /// Generate a plot using an arbitrary set of bin edges
    void plot( const std::string& branchexpr, const std::string& cutexpr,
      const std::vector< double >& bin_low_edges,
      const std::vector< std::string >& extra_tree_names = {},
      const std::string& x_axis_label = "",
      const std::string& y_axis_label = "",
      const std::string& title = "" ) const;

    /// Generate a plot with uniform binning over a specified x-axis range
    void plot( const std::string& branchexpr, const std::string& cutexpr,
      double xmin, double xmax, int Nbins,
      const std::vector< std::string >& extra_tree_names = {},
      const std::string& x_axis_label = "",
      const std::string& y_axis_label = "",
      const std::string& title = "" ) const;

    /// Set a custom string expression to use for assigning event weights
    /// in calls to TTree::Draw()
    inline void set_mc_event_weight_string( const std::string& weight_expr )
      { mc_event_weight_ = weight_expr; }

  protected:

    /// Selection used to define MC categories for making the stacked plots
    std::unique_ptr< SelectionBase > sel_for_categories_;

    /// Run numbers to use when generating plots
    std::set< int > runs_;

    /// String expression used with TTree::Draw() to apply event weights to the
    /// central-value MC prediction
    std::string mc_event_weight_ = DEFAULT_MC_EVENT_WEIGHT;
};
