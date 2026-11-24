#pragma once

// Combination matrices for a global goodness of fit across blocks that share
// events.
//
// A measurement reported as several blocks (1D and sliced 2D distributions)
// built from the same selected events fills exactly one bin per block for
// every event. The indicator vector of an event (one 1 per block) is a row of
// the combination matrix Omega. Its null space holds the exact linear
// constraints that the event sharing imposes between the blocks (e.g. every
// block sums to the same total), along which the statistical covariance
// vanishes. A global chi2 over all blocks must be evaluated on the range of
// Omega^T Omega (in reconstructed space), with ndof = rank.
//
// Two versions are provided:
//   - structural: from the binning geometry alone (BinSchemeBase blocks),
//     enumerating the cells of the variables shared between blocks; it
//     predicts the rank before any simulation exists;
//   - MC: from the bin patterns that the selected events actually fill
//     (collected by UniverseMaker for every input file).
// Both are stored as Omega^T Omega over a set of distinct rows (patterns).
// Only its range (rank and null space) is meaningful: duplicating rows or
// adding rows that are combinations of others changes the entries but not
// the range.
//
// Reference: the range-projected global chi2 for blockwise measurements,
// JINST 21 (2026) 11, P11017.

// Standard library includes
#include <set>
#include <string>
#include <vector>

// ROOT includes
#include "TMatrixD.h"
#include "TMatrixDSym.h"

class BinSchemeBase;

namespace CombinationMatrix {

  // A row of Omega: the sorted global reco bin indices that one event fills
  using Pattern = std::vector< int >;

  // Geometry of one block of ordinary reco bins
  struct BlockGeometry {
    int block_index = -1;  // position in BinSchemeBase::vect_block
    int offset = 0;        // global index of the first reco bin of the block
    int num_bins = 0;
    std::string selection; // blocks with different selections are treated
                           // as disjoint event populations
    bool is_2d = false;
    // 1D: var = the binned variable, edges = its bin edges
    // 2D: var = the slice variable, edges = the slice edges; inner_var and
    //     inner_edges[s] = the variable binned within slice s and its edges
    std::string var;
    std::vector< double > edges;
    std::string inner_var;
    std::vector< std::vector< double > > inner_edges;
  };

  // The geometry of every ordinary reco block of a bin scheme, with global
  // offsets following the reco bin order that MakeConfig writes (blocks in
  // order; 2D blocks slice by slice). num_reco_bins is set to the total
  // number of reco bins, including the sideband bins appended after them.
  std::vector< BlockGeometry > reco_block_geometry( BinSchemeBase& scheme,
    int& num_reco_bins );

  // Structural Omega^T Omega over num_bins columns, accumulated over a set of
  // event patterns achievable under the binning geometry whose span equals
  // that of the full enumeration of the reference: for every cell of
  // the refined grid of the variables shared between blocks, one base
  // pattern (the first possible bin of every block); and, for every pair of
  // bins (a, b) of a block that are both possible in some cell, one pattern
  // that differs from that cell's base pattern only by b in place of a. Any
  // achievable pattern is a base pattern plus such (e_b - e_a) differences,
  // so the range is exact while the number of rows stays ~ the number of
  // shared cells. num_rows is set to the number of rows accumulated.
  TMatrixDSym structural_gram( const std::vector< BlockGeometry >& blocks,
    int num_bins, long long& num_rows, bool verbose = true );

  // Omega^T Omega for a set of patterns over num_bins columns
  TMatrixDSym pattern_gram( const std::set< Pattern >& patterns, int num_bins );

  // Number of eigenvalues of a symmetric positive semi-definite matrix above
  // rel_threshold times the largest one
  int psd_rank( const TMatrixDSym& m, double rel_threshold = 1e-9 );

  // Orthonormal basis (columns) of the range of a symmetric PSD matrix, with
  // the eigenvalue threshold of psd_rank
  TMatrixD range_basis( const TMatrixDSym& m, double rel_threshold = 1e-9 );

  // Range-projected global chi2 of an unfolded result (JINST 21 (2026) 11,
  // P11017): the residual (unfolded signal - A_C * prediction, an n x 1
  // column) and its covariance are taken back to reconstructed space with
  // the inverse of the unfolding matrix U (unfolded = U * background-
  // subtracted data), projected onto the columns of V (an orthonormal basis
  // of the range of the combination matrix, see range_basis), and compared
  // with an ordinary chi2 on the reduced, full-rank system. ndof = the
  // number of columns of V. In reconstructed space the value does not depend
  // on the regularisation.
  double global_chi2( const TMatrixD& residual, const TMatrixD& cov,
    const TMatrixD& U, const TMatrixD& V );

  // Default eigenvalue threshold, relative to the largest eigenvalue
  constexpr double DEFAULT_REL_THRESHOLD = 1e-9;

  // Names of the objects written to ROOT files
  const std::string STRUCTURAL_NAME = "structural_combination";
  const std::string MC_NAME = "combination_matrix";
  const std::string MC_NUM_PATTERNS_NAME = "num_combination_patterns";

}
