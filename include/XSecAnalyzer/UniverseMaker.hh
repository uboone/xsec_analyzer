#pragma once

// Standard library includes
#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

// ROOT includes

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TTreeFormula.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/WeightHandler.hh"

#include "Selections/SelectionBase.hh"
#include "Selections/SelectionFactory.hh"

// Only load the library in this way if we're using this code from inside the
// ROOT C++ interpreter. We could check for __CINT__ as well, but the specific
// R__LOAD_LIBRARY approach used here only works with ROOT 6.
#ifdef __CLING__
// Pre-load the definition of the TTreeFormula::EvalInstance function from the
// TreePlayer shared library. This approach is based on the trick mentioned
// here: https://tinyurl.com/2s4yuzxm
R__LOAD_LIBRARY(libTreePlayer.so)
#endif

// Keys used to identify true and reco bin configurations in a universe output
// ROOT file
const std::string TRUE_BIN_SPEC_NAME = "true_bin_spec";
const std::string RECO_BIN_SPEC_NAME = "reco_bin_spec";

// Converts the name of an analysis ntuple file (typically with the full path)
// into a TDirectoryFile name to use as a subfolder of the main output
// TDirectoryFile used for saving universes. Since the forward slash
// character '/' cannot be used in a TDirectoryFile name, this function
// replaces all instances of this character by a '+' instead. The technique
// used here is based on https://stackoverflow.com/a/2896627/4081973
inline std::string ntuple_subfolder_from_file_name(
  const std::string& file_name )
{
  constexpr char FORWARD_SLASH = '/';
  constexpr char PLUS = '+';

  std::string result = file_name;
  std::replace( result.begin(), result.end(), FORWARD_SLASH, PLUS );
  return result;
}

// Branch names for special event weights
const std::string SPLINE_WEIGHT_NAME = "weight_splines_general_Spline";
const std::string TUNE_WEIGHT_NAME = "weight_TunedCentralValue_UBGenie";

// Special weight name to store the unweighted event counts
const std::string UNWEIGHTED_NAME = "unweighted";

constexpr double MIN_WEIGHT = 0.;
constexpr double MAX_WEIGHT = 30.;

// Event weights that are below MIN_WEIGHT, above MAX_WEIGHT, infinite, or NaN
// are reset to unity by this function. Other weights are returned unaltered.
inline double safe_weight( double w ) {
  if ( std::isfinite(w) && w >= MIN_WEIGHT && w <= MAX_WEIGHT ) return w;
  else return 1.0;
}

// Utility function used to check endings of (trimmed) weight labels based on
// branch names in the weights TTree
inline bool string_has_end( const std::string& str, const std::string& end ) {
  if ( str.length() >= end.length() ) {
    int comp = str.compare( str.length() - end.length(),
      end.length(), end );
    bool test_result = ( comp == 0 );
    return test_result;
  }
  return false;
}

// Multiplies a given event weight by extra correction factors as appropriate.
// TODO: include the rootino_fix weight as a correction to the central value
inline void apply_cv_correction_weights( const std::string& wgt_name,
  double& wgt, double spline_weight, double tune_weight )
{
  if ( string_has_end(wgt_name, "UBGenie") ) {
    wgt *= spline_weight;
  }
  else if ( wgt_name == "weight_flux_all"
    || wgt_name == "weight_reint_all"
    || wgt_name == "weight_xsr_scc_Fa3_SCC"
    || wgt_name == "weight_xsr_scc_Fv3_SCC" )
  {
    wgt *= spline_weight * tune_weight;
  }
  else if ( wgt_name == SPLINE_WEIGHT_NAME ) {
    // No extra weight factors needed
    return;
  }
  else throw std::runtime_error( "Unrecognized weight name" );
}

// Enum used to label bin types in true space
enum TrueBinType {

  // Bin contains signal events
  kSignalTrueBin = 0,

  // Bin contains background events
  kBackgroundTrueBin = 1,

  // Placeholder for undefined cases
  kUnknownTrueBin = 2

};

// Enum used to label bin types in reco space
enum RecoBinType {

  // Bin intended for use in the primary measurement of interest
  kOrdinaryRecoBin = 0,

  // Bin intended for use in a sideband constraint
  kSidebandRecoBin = 1,

  // Placeholder for undefined cases
  kUnknownRecoBin = 2

};

// Defines a histogram bin in true space
struct TrueBin {

  public:

    inline TrueBin( const std::string& cuts = "",
      TrueBinType bin_type = kSignalTrueBin,
      int block = -1 )
      : signal_cuts_( cuts ), type_( bin_type ), block_index_( block ) {}

    // Cuts to use in TTree::Draw for filling the bin. Any overall event weight
    // included here will be ignored. It is up to the user to ensure that only
    // true quantities are used for true bin cuts.
    std::string signal_cuts_;

    // Indicates the nature of the events contained in this bin
    TrueBinType type_;

    // Index used to group true signal bins together for purposes of
    // unfolding. The numerical value is otherwise ignored.
    int block_index_;
};

inline std::ostream& operator<<( std::ostream& out, const TrueBin& tb ) {
  out << tb.type_ << ' ' << tb.block_index_ << " \""
    << tb.signal_cuts_ << '\"';
  return out;
}

inline std::istream& operator>>( std::istream& in, TrueBin& tb ) {

  int temp_type;
  in >> temp_type;

  tb.type_ = static_cast< TrueBinType >( temp_type );

  int block_idx;
  in >> block_idx;

  tb.block_index_ = block_idx;

  // Use two calls to std::getline using a double quote delimiter
  // in order to get the contents of the next double-quoted string
  std::string temp_line;
  std::getline( in, temp_line, '\"' );
  std::getline( in, temp_line, '\"' );
  tb.signal_cuts_ = temp_line;

  return in;
}

// Defines a histogram bin in reco space
struct RecoBin {

  public:

    inline RecoBin( const std::string& cuts = "",
      RecoBinType bin_type = kOrdinaryRecoBin,
      int block_idx = -1 )
      : type_( bin_type ), selection_cuts_( cuts ),
      block_index_( block_idx ) {}

    // Cuts to use in TTree::Draw for filling the bin. Any overall event weight
    // included here will be ignored. It is up to the user to ensure that only
    // reco quantities are used for reco bin cuts.
    std::string selection_cuts_;

    // The kind of reco bin represented by this object
    RecoBinType type_;

    // Index used to group ordinary reco bins together for purposes of
    // unfolding. The numerical value is otherwise ignored.
    int block_index_;

};

inline std::ostream& operator<<( std::ostream& out, const RecoBin& rb ) {
  out << rb.type_ << ' ' << rb.block_index_ << " \""
    << rb.selection_cuts_ << '\"';
  return out;
}

inline std::istream& operator>>( std::istream& in, RecoBin& rb ) {

  int temp_type;
  in >> temp_type;

  rb.type_ = static_cast< RecoBinType >( temp_type );

  int block_idx;
  in >> block_idx;

  rb.block_index_ = block_idx;

  // Use two calls to std::getline using a double quote delimiter
  // in order to get the contents of the next double-quoted string
  std::string temp_line;
  std::getline( in, temp_line, '\"' );
  std::getline( in, temp_line, '\"' );
  rb.selection_cuts_ = temp_line;

  return in;
}

// Provides a set of histograms used to store summed bin counts (with
// associated MC statistical uncertainties) in a given systematic variation
// universe
class Universe {

  public:

    inline Universe( const std::string& universe_name,
      size_t universe_index, int num_true_bins, int num_reco_bins )
      : universe_name_( universe_name ), index_( universe_index )
    {
      std::string hist_name_prefix = universe_name + '_'
        + std::to_string( universe_index );

      hist_true_ = std::make_unique< TH1D >(
        (hist_name_prefix + "_true").c_str(), "; true bin number; events",
        num_true_bins, 0., num_true_bins );

      hist_reco_ = std::make_unique< TH1D >(
        (hist_name_prefix + "_reco").c_str(), "; reco bin number; events",
        num_reco_bins, 0., num_reco_bins );

      hist_2d_ = std::make_unique< TH2D >( (hist_name_prefix + "_2d").c_str(),
        "; true bin number; reco bin number; counts", num_true_bins, 0.,
        num_true_bins, num_reco_bins, 0., num_reco_bins );

      hist_reco2d_ = std::make_unique< TH2D >(
        (hist_name_prefix + "_reco2d").c_str(),
        "; reco bin number; reco bin number; counts", num_reco_bins, 0.,
        num_reco_bins, num_reco_bins, 0., num_reco_bins );

      hist_true2d_ = std::make_unique< TH2D >(
        (hist_name_prefix + "_true2d").c_str(),
        "; true bin number; true bin number; counts", num_true_bins, 0.,
        num_true_bins, num_true_bins, 0., num_true_bins );

      hist_categ_ = std::make_unique< TH2D >(
        (hist_name_prefix + "_categ").c_str(),
        "; true event category; reco bin number; counts", num_categories_, 0.,
        num_categories_, num_reco_bins, 0., num_reco_bins );

      // Store summed squares of event weights (for calculations of the MC
      // statistical uncertainty on bin contents)
      hist_true_->Sumw2();
      hist_reco_->Sumw2();
      hist_2d_->Sumw2();
      hist_categ_->Sumw2();
      hist_reco2d_->Sumw2();
      hist_true2d_->Sumw2();
    }

    // Note: the new Universe object takes ownership of the histogram
    // pointers passed to this constructor
    inline Universe( const std::string& universe_name,
      size_t universe_index, TH1D* hist_true, TH1D* hist_reco, TH2D* hist_2d,
      TH2D* hist_categ, TH2D* hist_reco2d, TH2D* hist_true2d )
      : universe_name_( universe_name ), index_( universe_index ),
      hist_true_( hist_true ), hist_reco_( hist_reco ), hist_2d_( hist_2d ),
      hist_categ_( hist_categ ), hist_reco2d_( hist_reco2d ),
      hist_true2d_( hist_true2d )
    {
      hist_true_->SetDirectory( nullptr );
      hist_reco_->SetDirectory( nullptr );
      hist_2d_->SetDirectory( nullptr );
      hist_categ_->SetDirectory( nullptr );
      hist_reco2d_->SetDirectory( nullptr );
      hist_true2d_->SetDirectory( nullptr );
    }

    inline std::unique_ptr< Universe > clone() const {
      int num_true_bins = hist_2d_->GetXaxis()->GetNbins();
      int num_reco_bins = hist_2d_->GetYaxis()->GetNbins();
      auto result = std::make_unique< Universe >( universe_name_,
        index_, num_true_bins, num_reco_bins );

      result->hist_true_->Add( this->hist_true_.get() );
      result->hist_reco_->Add( this->hist_reco_.get() );
      result->hist_2d_->Add( this->hist_2d_.get() );
      result->hist_categ_->Add( this->hist_categ_.get() );
      result->hist_reco2d_->Add( this->hist_reco2d_.get() );
      result->hist_true2d_->Add( this->hist_true2d_.get() );

      return result;
    }

    inline static void set_num_categories( const int count )
      { num_categories_ = count; }

    std::string universe_name_;
    size_t index_;
    std::unique_ptr< TH1D > hist_true_;
    std::unique_ptr< TH1D > hist_reco_;
    std::unique_ptr< TH2D > hist_2d_;
    std::unique_ptr< TH2D > hist_categ_;
    std::unique_ptr< TH2D > hist_reco2d_;
    std::unique_ptr< TH2D > hist_true2d_;

    static size_t num_categories_;
};

class UniverseMaker {

  public:

    // Initialize the UniverseMaker with true and reco bin definitions
    // stored in a configuration file
    UniverseMaker( const std::string& config_file_name );

    // Overloaded constructor that reads the configuration settings using
    // an existing input stream
    UniverseMaker( std::istream& config_stream );

    // Add an ntuple input file to the owned TChain
    void add_input_file( const std::string& input_file_name );

    // Access the bin definitions
    inline const auto& true_bins() const { return true_bins_; }
    inline const auto& reco_bins() const { return reco_bins_; }

    // Access the owned TChain
    inline auto& input_chain() { return input_chain_; }

    // Does the actual calculation of the event histograms across the
    // various systematic universes. The optional argument points to a vector
    // of branch names that will be used to retrieve systematic universe
    // weights. If it is omitted, all available ones will be auto-detected and
    // used.
    void build_universes(
      const std::vector<std::string>* universe_branch_names = nullptr );

    // Overloaded version of the function that takes a reference the
    // vector of universe branch names (for convenience). The behavior
    // is the same as the original, but in this case the explicit vector of
    // branch names definitely exists.
    void build_universes(
      const std::vector<std::string>& universe_branch_names );

    // Writes the universe histograms to an output ROOT file
    void save_histograms( const std::string& output_file_name,
      const std::string& subdirectory_name, bool update_file = true );

    // Provides read-only access to the map of Universe objects
    const auto& universe_map() const { return universes_; }

    // Returns the name of the TDirectoryFile that will be used to hold the
    // universe histograms when they are written to the output ROOT file
    const std::string& dir_name() const { return output_directory_name_; }

  protected:

    // Helper function used by the constructors
    void init( std::istream& in_file );

    // Helper struct that keeps track of bin indices and TTreeFormula weights
    // when filling universe histograms
    struct FormulaMatch {
      // The weight given here is from an evaluation of a TTreeFormula for
      // filling an individual bin (as opposed to an overall event weight which
      // may be given separately)
      inline FormulaMatch( size_t bin_idx, double wgt )
        : bin_index_( bin_idx ), weight_( wgt ) {}

      size_t bin_index_;
      double weight_;
    };

    // Prepares the TTreeFormula objects needed to test each entry for
    // membership in each bin
    void prepare_formulas();

    // Prepares the Universe objects needed to store summed event weights for
    // each bin in each systematic variation universe
    void prepare_universes( const WeightHandler& wh );

    // Bin definitions in true space
    std::vector< TrueBin > true_bins_;

    // Bin definitions in reco space
    std::vector< RecoBin > reco_bins_;

    // A TChain containing MC event ntuples that will be used to compute the
    // universe histograms
    TChain input_chain_;

    // TTreeFormula objects used to test whether the current TChain entry falls
    // into each true bin
    std::vector< std::unique_ptr<TTreeFormula> > true_bin_formulas_;

    // TTreeFormula objects used to test whether the current TChain entry falls
    // into each reco bin
    std::vector< std::unique_ptr<TTreeFormula> > reco_bin_formulas_;

    // TTreeFormula objects used to test whether the current TChain entry falls
    // into each true EventCategory
    std::vector< std::unique_ptr<TTreeFormula> > category_formulas_;

    // Stores Universe objects used to accumulate event weights
    std::map< std::string, std::vector<Universe> > universes_;

    // Root TDirectoryFile name to use when writing the universes to an output
    // ROOT file
    std::string output_directory_name_;

    // Selection whose event category definitions will be used to
    // populate the category histograms in Universes
    // std::unique_ptr< SelectionBase > sel_for_categories_;
    //FIXME: using normal pointer to avoid invalid pointer error
    SelectionBase *sel_for_categories_;
};
