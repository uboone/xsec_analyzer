#pragma once

// Standard library includes
#include <memory>
#include <set>
#include <stdexcept>

// ROOT includes
#include "TMatrixD.h"

// Forward-delcare some needed objects
struct TrueBin;
struct RecoBin;
class SystematicsCalculator;

// Simple container for the output of Unfolder::unfold()
struct UnfoldedMeasurement {
  UnfoldedMeasurement( TMatrixD* unfolded_signal, TMatrixD* cov_matrix,
    TMatrixD* unfolding_matrix, TMatrixD* err_prop_matrix,
    TMatrixD* add_smear_matrix, TMatrixD* smearcept )
    : unfolded_signal_( unfolded_signal ), cov_matrix_( cov_matrix ),
    unfolding_matrix_( unfolding_matrix ), err_prop_matrix_( err_prop_matrix ),
    add_smear_matrix_( add_smear_matrix ), response_matrix_( smearcept )
    {}

  std::unique_ptr< TMatrixD > unfolded_signal_;
  std::unique_ptr< TMatrixD > cov_matrix_;
  std::unique_ptr< TMatrixD > unfolding_matrix_;
  std::unique_ptr< TMatrixD > err_prop_matrix_;
  std::unique_ptr< TMatrixD > add_smear_matrix_;
  std::unique_ptr< TMatrixD > response_matrix_;
};

// Container for mapping block indices to bin indices in
// Unfolder::blockwise_unfold()
struct BlockBins {
  BlockBins() {}
  std::vector< size_t > true_bin_indices_;
  std::vector< size_t > reco_bin_indices_;
};

// Abstract base class for objects that implement an algorithm for unfolding
// measured background-subtracted event counts from reco space to
// true space, possibly with regularization.
class Unfolder {

  public:

    Unfolder() {}

    // Function that actually implements a specific unfolding algorithm
    virtual UnfoldedMeasurement unfold( const TMatrixD& data_signal,
      const TMatrixD& data_covmat, const TMatrixD& smearcept,
      const TMatrixD& prior_true_signal ) const = 0;

    virtual UnfoldedMeasurement unfold( const TMatrixD& data_signal,
      const TMatrixD& data_covmat, const TMatrixD& smearcept,
      const TMatrixD& prior_true_signal,
      const std::string& block_defs_file ) const final;

    virtual UnfoldedMeasurement unfold(
      const SystematicsCalculator& syst_calc ) const final;

    // Helper function that sets up unfolding for multiple blocks of bins,
    // then combines the results
    virtual UnfoldedMeasurement blockwise_unfold(
      const TMatrixD& data_signal, const TMatrixD& data_covmat,
      const TMatrixD& smearcept, const TMatrixD& prior_true_signal,
      const std::vector< TrueBin >& true_bins,
      const std::vector< RecoBin >& reco_bins ) const final;

    // Modified version that gets the block definitions from a text input file.
    // Used for standalone blockwise unfolding without using
    // SystematicsCalculator
    virtual UnfoldedMeasurement blockwise_unfold(
      const TMatrixD& data_signal, const TMatrixD& data_covmat,
      const TMatrixD& smearcept, const TMatrixD& prior_true_signal,
      std::ifstream& in_block_file ) const final;

  protected:

    // Helper function that does some sanity checks on the dimensions of the
    // input matrices passed to unfold()
    static void check_matrices( const TMatrixD& data_signal,
      const TMatrixD& data_covmat, const TMatrixD& smearcept,
      const TMatrixD& prior_true_signal );
};


