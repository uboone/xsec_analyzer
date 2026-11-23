// Combination matrices for a global goodness of fit across blocks that share
// events. See XSecAnalyzer/CombinationMatrix.hh.

// Standard library includes
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <unordered_set>

// ROOT includes
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/BinSchemeBase.hh"
#include "XSecAnalyzer/CombinationMatrix.hh"

namespace {

  constexpr double EDGE_TOLERANCE = 1e-9;

  // Index of the bin of x in the sorted edges, or -1 outside them. Lower
  // edges are inclusive and upper edges exclusive, as in the reco bin
  // definitions; the cell midpoints used here never sit on an edge.
  int find_bin( double x, const std::vector< double >& edges ) {
    if ( edges.size() < 2u || x < edges.front() || x >= edges.back() ) {
      return -1;
    }
    auto it = std::upper_bound( edges.cbegin(), edges.cend(), x );
    return static_cast< int >( it - edges.cbegin() ) - 1;
  }

  // Sorted union of edges, merging values closer than EDGE_TOLERANCE
  std::vector< double > merge_edges( std::vector< double > e ) {
    std::sort( e.begin(), e.end() );
    std::vector< double > out;
    for ( double x : e ) {
      if ( out.empty() || std::abs( x - out.back() ) > EDGE_TOLERANCE ) {
        out.push_back( x );
      }
    }
    return out;
  }

  // Inner edges of a Block2D slice. A slice binned { -DBL_MAX, DBL_MAX } is a
  // single bin over the whole slice (Block2D stores it as { 0, 1 }).
  std::vector< double > inner_edges_of( const std::vector< double >& v ) {
    if ( v.size() == 2u && v.front() == -DBL_MAX && v.back() == DBL_MAX ) {
      return { -DBL_MAX, DBL_MAX };
    }
    return v;
  }

  int count_bins( Block* b ) {
    if ( b->Is1D() ) return b->GetNBinsX();
    int n = 0;
    if ( b->Is2D() ) {
      for ( int j = 0; j < b->GetNBinsX(); ++j ) n += b->GetNBinsY( j );
      return n;
    }
    for ( int j = 0; j < b->GetNBinsX(); ++j ) {
      for ( int k = 0; k < b->GetNBinsY( j ); ++k ) n += b->GetNBinsZ( j, k );
    }
    return n;
  }

  void add_row( TMatrixDSym& g, const CombinationMatrix::Pattern& p ) {
    for ( int a : p ) for ( int b : p ) g( a, b ) += 1.;
  }

  // Structural Omega^T Omega for the blocks of one selection, added to g
  long long accumulate_group( const std::vector< CombinationMatrix::BlockGeometry >& blocks,
    TMatrixDSym& g, bool verbose )
  {
    // Which blocks use each variable
    std::map< std::string, std::set< int > > users;
    for ( size_t k = 0u; k < blocks.size(); ++k ) {
      users[ blocks[ k ].var ].insert( k );
      if ( blocks[ k ].is_2d ) users[ blocks[ k ].inner_var ].insert( k );
    }
    std::vector< std::string > shared;
    for ( const auto& u : users ) if ( u.second.size() >= 2u ) shared.push_back( u.first );

    // Refined edges of the shared variables: the union of every block's edges
    std::map< std::string, std::vector< double > > refined;
    for ( const auto& b : blocks ) {
      auto& r_outer = refined[ b.var ];
      r_outer.insert( r_outer.end(), b.edges.cbegin(), b.edges.cend() );
      if ( b.is_2d ) {
        auto& r_inner = refined[ b.inner_var ];
        for ( const auto& e : b.inner_edges ) {
          auto ie = inner_edges_of( e );
          if ( ie.front() == -DBL_MAX ) continue; // unbinned slice: no edges
          r_inner.insert( r_inner.end(), ie.cbegin(), ie.cend() );
        }
      }
    }
    std::vector< std::vector< double > > cells_edges;
    long double num_cells = 1.;
    for ( const auto& v : shared ) {
      cells_edges.push_back( merge_edges( refined.at( v ) ) );
      num_cells *= cells_edges.back().size() - 1u;
    }
    std::map< std::string, int > shared_index;
    for ( size_t i = 0u; i < shared.size(); ++i ) shared_index[ shared[ i ] ] = i;

    if ( verbose ) {
      std::cout << "  shared variables:";
      for ( size_t i = 0u; i < shared.size(); ++i ) {
        std::cout << ' ' << shared[ i ] << " (" << cells_edges[ i ].size() - 1u
          << " cells)";
      }
      std::cout << "\n  shared cells to enumerate: " << static_cast< double >( num_cells )
        << '\n';
    }
    if ( num_cells > 5e8 ) {
      throw std::runtime_error( "CombinationMatrix: too many shared cells to"
        " enumerate" );
    }

    // Per block: the slice bin offsets (2D) and whether its options do not
    // depend on the cell (no shared variable)
    std::vector< std::vector< int > > slice_offset( blocks.size() );
    std::vector< bool > constant( blocks.size() );
    for ( size_t k = 0u; k < blocks.size(); ++k ) {
      const auto& b = blocks[ k ];
      bool s = shared_index.count( b.var );
      bool i = b.is_2d && shared_index.count( b.inner_var );
      constant[ k ] = !s && !i;
      if ( b.is_2d ) {
        int o = b.offset;
        for ( const auto& e : b.inner_edges ) {
          slice_offset[ k ].push_back( o );
          o += e.size() - 1u;
        }
      }
    }

    // Possible global bins of block k in the cell with midpoints mid
    std::vector< double > mid( shared.size() );
    auto options = [ & ]( size_t k, std::vector< int >& opts ) {
      opts.clear();
      const auto& b = blocks[ k ];
      auto mid_of = [ & ]( const std::string& v, double& x ) {
        auto it = shared_index.find( v );
        if ( it == shared_index.end() ) return false;
        x = mid[ it->second ];
        return true;
      };
      double xs, xi;
      bool s = mid_of( b.var, xs );
      if ( !b.is_2d ) {
        if ( s ) {
          int j = find_bin( xs, b.edges );
          if ( j >= 0 ) opts.push_back( b.offset + j );
        }
        else for ( int j = 0; j < b.num_bins; ++j ) opts.push_back( b.offset + j );
        return;
      }
      bool i = mid_of( b.inner_var, xi );
      int ns = b.inner_edges.size();
      int s_lo = 0, s_hi = ns;
      if ( s ) {
        int js = find_bin( xs, b.edges );
        if ( js < 0 ) return;
        s_lo = js;
        s_hi = js + 1;
      }
      for ( int js = s_lo; js < s_hi; ++js ) {
        auto ie = inner_edges_of( b.inner_edges[ js ] );
        if ( i ) {
          int j = find_bin( xi, ie );
          if ( j >= 0 ) opts.push_back( slice_offset[ k ][ js ] + j );
        }
        else {
          for ( size_t j = 0u; j + 1u < b.inner_edges[ js ].size(); ++j ) {
            opts.push_back( slice_offset[ k ][ js ] + j );
          }
        }
      }
    };

    long long rows = 0;
    long long valid_cells = 0;
    std::vector< std::vector< int > > opts( blocks.size() );
    std::vector< std::unordered_set< long long > > seen_pairs( blocks.size() );
    std::vector< bool > constant_done( blocks.size(), false );
    const long long stride = g.GetNrows();
    CombinationMatrix::Pattern base( blocks.size() ), row;

    std::vector< size_t > cell( shared.size(), 0u );
    bool done = false;
    while ( !done ) {
      for ( size_t i = 0u; i < shared.size(); ++i ) {
        const auto& e = cells_edges[ i ];
        mid[ i ] = 0.5 * ( e[ cell[ i ] ] + e[ cell[ i ] + 1u ] );
      }
      bool valid = true;
      for ( size_t k = 0u; k < blocks.size() && valid; ++k ) {
        options( k, opts[ k ] );
        valid = !opts[ k ].empty();
      }
      if ( valid ) {
        ++valid_cells;
        for ( size_t k = 0u; k < blocks.size(); ++k ) base[ k ] = opts[ k ].front();
        add_row( g, base );
        ++rows;
        for ( size_t k = 0u; k < blocks.size(); ++k ) {
          if ( opts[ k ].size() < 2u ) continue;
          if ( constant[ k ] && constant_done[ k ] ) continue;
          for ( size_t j = 1u; j < opts[ k ].size(); ++j ) {
            long long key = opts[ k ].front() * stride + opts[ k ][ j ];
            if ( !seen_pairs[ k ].insert( key ).second ) continue;
            row = base;
            row[ k ] = opts[ k ][ j ];
            add_row( g, row );
            ++rows;
          }
          if ( constant[ k ] ) constant_done[ k ] = true;
        }
      }
      // Next cell (odometer)
      done = true;
      for ( size_t i = 0u; i < shared.size(); ++i ) {
        if ( ++cell[ i ] + 1u < cells_edges[ i ].size() ) { done = false; break; }
        cell[ i ] = 0u;
      }
    }
    if ( verbose ) {
      std::cout << "  valid cells: " << valid_cells << ", rows: " << rows << '\n';
    }
    return rows;
  }

}

std::vector< CombinationMatrix::BlockGeometry >
  CombinationMatrix::reco_block_geometry( BinSchemeBase& scheme, int& num_reco_bins )
{
  std::vector< BlockGeometry > out;
  int offset = 0;
  for ( size_t i = 0u; i < scheme.vect_block.size(); ++i ) {
    Block* b = scheme.vect_block[ i ].block_reco_;
    BlockGeometry g;
    g.block_index = i;
    g.offset = offset;
    g.num_bins = count_bins( b );
    g.selection = b->GetSelection();
    if ( b->Is1D() ) {
      g.var = b->GetXName();
      g.edges = b->GetVector();
    }
    else if ( b->Is2D() ) {
      auto* b2 = dynamic_cast< Block2D* >( b );
      if ( !b2 ) throw std::runtime_error( "CombinationMatrix: unexpected 2D"
        " block type" );
      g.is_2d = true;
      g.var = b->GetXName();
      g.inner_var = b->GetYName();
      const auto map = b2->GetMap();
      for ( auto it = map.cbegin(); it != map.cend(); ++it ) {
        g.edges.push_back( it->first );
        auto next = it;
        if ( ++next == map.cend() ) break;
        g.inner_edges.push_back( it->second );
      }
    }
    else {
      throw std::runtime_error( "CombinationMatrix: 3D blocks are not"
        " supported yet" );
    }
    if ( b->GetBinType() == kOrdinaryRecoBin ) out.push_back( g );
    offset += g.num_bins;
  }
  for ( auto& sb : scheme.vect_sideband ) offset += count_bins( sb.block_reco_ );
  num_reco_bins = offset;
  return out;
}

TMatrixDSym CombinationMatrix::structural_gram(
  const std::vector< BlockGeometry >& blocks, int num_bins,
  long long& num_rows, bool verbose )
{
  TMatrixDSym g( num_bins );
  num_rows = 0;

  // Blocks of different selections share no events: one group per selection
  std::map< std::string, std::vector< BlockGeometry > > groups;
  for ( const auto& b : blocks ) groups[ b.selection ].push_back( b );
  if ( groups.size() > 1u && verbose ) {
    std::cout << "CombinationMatrix: " << groups.size() << " selections;"
      " blocks of different selections are treated as disjoint event"
      " populations\n";
  }
  for ( const auto& gr : groups ) {
    if ( verbose ) {
      std::cout << "Structural combination matrix, selection \"" << gr.first
        << "\" (" << gr.second.size() << " blocks)\n";
    }
    num_rows += accumulate_group( gr.second, g, verbose );
  }
  return g;
}

TMatrixDSym CombinationMatrix::pattern_gram( const std::set< Pattern >& patterns,
  int num_bins )
{
  TMatrixDSym g( num_bins );
  for ( const auto& p : patterns ) add_row( g, p );
  return g;
}

int CombinationMatrix::psd_rank( const TMatrixDSym& m, double rel_threshold ) {
  TMatrixDSymEigen eig( m );
  const TVectorD& w = eig.GetEigenValues();
  double wmax = 0.;
  for ( int i = 0; i < w.GetNrows(); ++i ) wmax = std::max( wmax, w( i ) );
  int r = 0;
  for ( int i = 0; i < w.GetNrows(); ++i ) if ( w( i ) > rel_threshold * wmax ) ++r;
  return r;
}

TMatrixD CombinationMatrix::range_basis( const TMatrixDSym& m, double rel_threshold ) {
  TMatrixDSymEigen eig( m );
  const TVectorD& w = eig.GetEigenValues();
  const TMatrixD& v = eig.GetEigenVectors();
  double wmax = 0.;
  for ( int i = 0; i < w.GetNrows(); ++i ) wmax = std::max( wmax, w( i ) );
  std::vector< int > keep;
  for ( int i = 0; i < w.GetNrows(); ++i ) if ( w( i ) > rel_threshold * wmax ) keep.push_back( i );
  TMatrixD out( m.GetNrows(), keep.size() );
  for ( size_t c = 0u; c < keep.size(); ++c ) {
    for ( int r = 0; r < m.GetNrows(); ++r ) out( r, c ) = v( r, keep[ c ] );
  }
  return out;
}

double CombinationMatrix::global_chi2( const TMatrixD& residual,
  const TMatrixD& cov, const TMatrixD& U, const TMatrixD& V )
{
  const int n = U.GetNrows();
  if ( U.GetNcols() != n || cov.GetNrows() != n || cov.GetNcols() != n
    || residual.GetNrows() != n || V.GetNrows() != n )
  {
    throw std::runtime_error( "CombinationMatrix::global_chi2: inconsistent"
      " dimensions" );
  }
  // Back to reconstructed space
  TMatrixD Uinv( U );
  double det = 0.;
  Uinv.Invert( &det );
  if ( det == 0. ) {
    throw std::runtime_error( "CombinationMatrix::global_chi2: singular"
      " unfolding matrix" );
  }
  TMatrixD Vt( TMatrixD::kTransposed, V );
  TMatrixD W( Vt, TMatrixD::kMult, Uinv );             // r x n
  TMatrixD r( W, TMatrixD::kMult, residual );          // r x 1
  TMatrixD Wt( TMatrixD::kTransposed, W );
  TMatrixD C( W, TMatrixD::kMult, TMatrixD( cov, TMatrixD::kMult, Wt ) );
  // Symmetrise against rounding, then solve the reduced system
  TMatrixD Ct( TMatrixD::kTransposed, C );
  C += Ct;
  C *= 0.5;
  TMatrixD Cinv( C );
  Cinv.Invert( &det );
  if ( det == 0. ) {
    throw std::runtime_error( "CombinationMatrix::global_chi2: singular"
      " reduced covariance" );
  }
  TMatrixD rt( TMatrixD::kTransposed, r );
  TMatrixD chi2( rt, TMatrixD::kMult, TMatrixD( Cinv, TMatrixD::kMult, r ) );
  return chi2( 0, 0 );
}
