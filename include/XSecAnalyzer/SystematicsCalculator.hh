#pragma once

// Standard library includes
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>

// ROOT includes
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TParameter.h"

// XSecAnalyzer includes
#include "FilePropertiesManager.hh"
#include "UniverseMaker.hh"

#include "Selections/SelectionBase.hh"
#include "Selections/SelectionFactory.hh"

// Helper function template that retrieves an object from a TDirectoryFile
// and loads a pointer to it into a std::unique_ptr of the correct type
template< typename T > std::unique_ptr< T > get_object_unique_ptr(
  const std::string& namecycle, TDirectory& df )
{
  T* temp_ptr = nullptr;
  df.GetObject( namecycle.c_str(), temp_ptr );

  // Set the directory to nullptr in the case of ROOT histograms. This will
  // avoid extra deletion attempts
  TH1* temp_hist_ptr = dynamic_cast< TH1* >( temp_ptr );
  if ( temp_hist_ptr ) temp_hist_ptr->SetDirectory( nullptr );

  return std::unique_ptr< T >( temp_ptr );
}

// Helper function to use on a new Universe object's histograms. Prevents
// auto-deletion problems by disassociating the Universe object's histograms
// from any TDirectory. Also turns off the stats box when plotting for
// convenience.
void set_stats_and_dir( Universe& univ );

// Tests whether a string ends with another string. Taken from
// https://stackoverflow.com/a/874160/4081973. In C++20, you could use
// std::string::ends_with() instead.
bool has_ending( const std::string& fullString, const std::string& ending );

// Simple container for a TH2D that represents a covariance matrix
struct CovMatrix {

  inline CovMatrix() {}

  inline CovMatrix( TH2D* cov_mat ) : cov_matrix_( cov_mat ) {}

  CovMatrix( const TMatrixD& matrix );

  std::unique_ptr< TH2D > cov_matrix_;

  // Helper function for operator+=
  void add_or_clone( std::unique_ptr<TH2D>& mine, TH2D* other );

  inline CovMatrix& operator+=( const CovMatrix& other ) {

    add_or_clone( cov_matrix_, other.cov_matrix_.get() );

    return *this;
  }

  inline std::unique_ptr< TMatrixD > get_matrix() const {
    // Note that ROOT histogram bin indices are one-based to allow for
    // underflow. The TMatrixDSym element indices, on the other hand,
    // are zero-based.
    int num_cm_bins = cov_matrix_->GetNbinsX();
    auto result = std::make_unique< TMatrixD >( num_cm_bins, num_cm_bins );
    // TODO: consider doing something more efficient than setting each
    // element manually
    for ( int a = 0; a < num_cm_bins; ++a ) {
      for ( int b = 0; b < num_cm_bins; ++b ) {
        result->operator()( a, b ) = cov_matrix_->GetBinContent( a + 1, b + 1 );
      }
    }
    return result;
  }

};

// Container that holds the results of subtracting the EXT+MC background from
// the measured data event counts in ordinary reco bins
struct MeasuredEvents {

  inline MeasuredEvents() {}

  inline MeasuredEvents( TMatrixD* bkgd_subtracted_data, TMatrixD* bkgd,
    TMatrixD* mc_plus_ext, TMatrixD* cov_mat )
    : reco_signal_( bkgd_subtracted_data ), reco_bkgd_( bkgd ),
    reco_mc_plus_ext_( mc_plus_ext ), cov_matrix_( cov_mat )
  {
    if ( bkgd_subtracted_data->GetNcols() != 1 ) throw std::runtime_error(
      "Non-row-vector background-subtracted signal passed to MeasuredEvents" );

    if ( bkgd->GetNcols() != 1 ) throw std::runtime_error( "Non-row-vector"
      " background prediction passed to MeasuredEvents" );

    if ( mc_plus_ext->GetNcols() != 1 ) throw std::runtime_error(
      "Non-row-vector MC+EXT prediction passed to MeasuredEvents" );

    int num_ordinary_reco_bins = bkgd_subtracted_data->GetNrows();
    if ( cov_mat->GetNcols() != num_ordinary_reco_bins
      || cov_mat->GetNrows() != num_ordinary_reco_bins )
    {
      throw std::runtime_error( "Bad covariance matrix dimensions passed to"
        " MeasuredEvents" );
    }
  }

  // Background-subtracted data event counts in the ordinary reco bins
  std::unique_ptr< TMatrixD > reco_signal_;

  // Background that was subtracted from each reco bin to form the signal
  // measurement
  std::unique_ptr< TMatrixD > reco_bkgd_;

  // Total MC+EXT prediction in each reco bin (used to estimate the covariance
  // matrix on the data)
  std::unique_ptr< TMatrixD > reco_mc_plus_ext_;

  // Covariance matrix for the background-subtracted data
  std::unique_ptr< TMatrixD > cov_matrix_;

};

using CovMatrixMap = std::map< std::string, CovMatrix >;
using NFT = NtupleFileType;

class SystematicsCalculator {

  public:

    SystematicsCalculator( const std::string& input_respmat_file_name,
      const std::string& syst_cfg_file_name = "",
      const std::string& respmat_tdirectoryfile_name = "" );

    void load_universes( TDirectoryFile& total_subdir );

    void build_universes( TDirectoryFile& root_tdir );

    void save_universes( TDirectoryFile& out_tdf );

    const Universe& cv_universe() const {
      return *rw_universes_.at( CV_UNIV_NAME ).front();
    }

    const std::unique_ptr< Universe >& fake_data_universe() const {
      return fake_data_universe_;
    }

    std::unique_ptr< CovMatrixMap > get_covariances() const;

    // Returns a background-subtracted measurement in all ordinary reco bins
    // with the total covariance matrix and the background event counts that
    // were subtracted.
    // NOTE: this function assumes that the ordinary reco bins are all listed
    // before any sideband reco bins
    virtual MeasuredEvents get_measured_events() const;

    // Utility functions to help with unfolding
    inline std::unique_ptr< TMatrixD > get_cv_smearceptance_matrix() const {
      const auto& cv_univ = this->cv_universe();
      return this->get_smearceptance_matrix( cv_univ );
    }

    std::unique_ptr< TMatrixD > get_smearceptance_matrix(
      const Universe& univ ) const;

    std::unique_ptr< TMatrixD > get_cv_true_signal() const;

    // Returns the expected background in each ordinary reco bin (including
    // both EXT and the central-value MC prediction for beam-correlated
    // backgrounds)
    // NOTE: This function assumes that all "ordinary" reco bins are listed
    // before the sideband ones.
    inline std::unique_ptr< TMatrixD > get_cv_ordinary_reco_bkgd() const
      { return this->get_cv_ordinary_reco_helper( true ); }

    // Returns the expected signal event counts in each ordinary reco bin
    // NOTE: This function assumes that all "ordinary" reco bins are listed
    // before the sideband ones.
    inline std::unique_ptr< TMatrixD > get_cv_ordinary_reco_signal() const
      { return this->get_cv_ordinary_reco_helper( false ); }

    inline size_t get_num_signal_true_bins() const
      { return num_signal_true_bins_; }

    // Dumps vectors of the observables evaluated in each of the systematic
    // universes to a text file for easy inspection / retrieval
    void dump_universe_observables( const std::string& out_file_name ) const;

    const SelectionBase& get_selection_for_categories() const
      { return *sel_for_categ_; }

  //protected:

    // Implements both get_cv_ordinary_reco_bkgd() and
    // get_cv_ordinary_reco_signal() in order to reduce code duplication. If
    // return_bkgd is false (true), then the background (signal) event counts
    // in each ordinary reco bin will be returned as a column vector.
    std::unique_ptr< TMatrixD > get_cv_ordinary_reco_helper(
      bool return_bkgd ) const;

    // Returns true if a given Universe represents a detector variation or
    // false otherwise
    bool is_detvar_universe( const Universe& univ ) const;

    // Overload for special cases in which the N*N covariance matrix does not
    // have dimension parameter N equal to the number of reco bins
    inline virtual size_t get_covariance_matrix_size() const
      { return reco_bins_.size(); }

    CovMatrix make_covariance_matrix( const std::string& hist_name ) const;

    // Evaluate the observable described by the covariance matrices in
    // a given universe and reco-space bin. NOTE: the reco bin index given
    // as an argument to this function is zero-based.
    virtual double evaluate_observable( const Universe& univ, int reco_bin,
      int flux_universe_index = -1 ) const = 0;

    // Evaluate a covariance matrix element for the data statistical
    // uncertainty on the observable of interest for a given pair of reco bins.
    // In cases where every event falls into a unique reco bin, only the
    // diagonal covariance matrix elements are non-vanishing. Do the
    // calculation either for BNB data (use_ext = false) or for EXT data
    // (use_ext = true). NOTE: the reco bin indices consumed by this function
    // are zero-based.
    virtual double evaluate_data_stat_covariance( int reco_bin_a,
      int reco_bin_b, bool use_ext ) const = 0;

    // Evaluate a covariance matrix element for the MC statistical uncertainty
    // on the observable of interest (including contributions from both signal
    // and background events) for a given pair of reco bins within a particular
    // universe. Typically the CV universe should be used, but MC statistical
    // uncertainties are tracked for all Universe objects in case they are
    // needed. In cases where every event falls into a unique reco bin, only
    // the diagonal covariance matrix elements are non-vanishing. NOTE: the
    // reco bin indices consumed by this function are zero-based.
    virtual double evaluate_mc_stat_covariance( const Universe& univ,
      int reco_bin_a, int reco_bin_b ) const = 0;

    // Utility function used by dump_universe_observables() to prepare the
    // output
    void dump_universe_helper( std::ostream& out, const Universe& univ,
      int flux_u_index = -1 ) const;

    // Central value universe name
    const std::string CV_UNIV_NAME = "weight_TunedCentralValue_UBGenie";

    // Beginning of the subdirectory name for the TDirectoryFile containing the
    // POT-summed histograms for the various universes across all analysis
    // ntuples. The full name is formed from this prefix and the name of the
    // FilePropertiesManager configuration file that is currently active.
    const std::string TOTAL_SUBFOLDER_NAME_PREFIX = "total_";

    // Holds reco-space histograms for data (BNB and EXT) bin counts
    std::map< NFT, std::unique_ptr<TH1D> > data_hists_;
    std::map< NFT, std::unique_ptr<TH2D> > data_hists2d_;

    // Holds universe objects for reweightable systematics
    std::map< std::string, std::vector< std::unique_ptr<Universe> > >
      rw_universes_;

    // Detector systematic universes (and the detVar CV) will be indexed using
    // ntuple file type values. We're currently scaling one set of ntuples
    // (from Run 3b) to the full dataset.
    // TODO: revisit this procedure if new detVar samples become available
    std::map< NFT, std::unique_ptr<Universe> > detvar_universes_;

    // If we are working with fake data, then this will point to a Universe
    // object containing full reco and truth information for the MC portion
    // (as opposed to the EXT contribution which is added to the reco MC
    // counts in the BNB "data" histogram). If we are working with real data,
    // then this will be a null pointer.
    std::unique_ptr< Universe > fake_data_universe_ = nullptr;

    // "Alternate CV" universes for assessing unisim systematics related to
    // interaction modeling
    std::map< NFT, std::unique_ptr<Universe> > alt_cv_universes_;

    // True bin configuration that was used to compute the universes
    std::vector< TrueBin > true_bins_;

    // Reco bin configuration that was used to compute the universes
    std::vector< RecoBin > reco_bins_;

    // Total POT exposure for the analyzed BNB data
    double total_bnb_data_pot_ = 0.;

    // Name of the systematics configuration file that should be used
    // when computing covariance matrices
    std::string syst_config_file_name_;

    // Number of entries in the reco_bins_ vector that are "ordinary" bins
    size_t num_ordinary_reco_bins_ = 0u;

    // Number of entries in the true_bins_ vector that are "signal" bins
    size_t num_signal_true_bins_ = 0u;

    // Selection used to assign event categories in the universes
    // std::unique_ptr< SelectionBase > sel_for_categ_;
    // FIXME: using normal pointer to avoid invalid pointer error
    SelectionBase *sel_for_categ_;
};
