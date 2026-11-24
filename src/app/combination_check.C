// Checks of the combination matrices used for a global chi2 across blocks
// that share events (see XSecAnalyzer/CombinationMatrix.hh).
//
//   CombinationCheck UNIV_FILE STRUCTURAL_FILE
//     Compares the MC combination matrices that univmake saved for every
//     input file of UNIV_FILE with the structural one written by
//     "BinScheme -s" (configs/<prefix>structural_combination.root): rank of
//     each, of their sum (the union of the patterns of all files), and
//     whether the MC range lies inside the structural range.
//
//   CombinationCheck --config BINNING_TXT [OUT_TXT]
//     Structural matrix of a binning in the text format of the code that
//     accompanies JINST 21 (2026) 11, P11017 (block headers "block <idx>
//     <var> [<inner var>] [sel <idx>]", followed by the edges, or by the
//     slices "<slice idx> <low> <high>" each with its inner edges), for
//     cross-checks against its Python implementation; OUT_TXT receives
//     Omega^T Omega.

// Standard library includes
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// ROOT includes
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TKey.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TParameter.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/CombinationMatrix.hh"

namespace CM = CombinationMatrix;

namespace {

  // Blocks of a binning file in the format described above
  std::vector< CM::BlockGeometry > parse_config( const std::string& name,
    int& num_bins )
  {
    std::ifstream in( name );
    if ( !in ) throw std::runtime_error( "cannot open " + name );
    std::vector< std::string > lines;
    std::string l;
    while ( std::getline( in, l ) ) {
      auto p = l.find_first_not_of( " \t" );
      if ( p == std::string::npos || l[ p ] == '#' ) continue;
      lines.push_back( l.substr( p ) );
    }
    auto numbers = []( const std::string& s ) {
      std::istringstream is( s );
      std::vector< double > v;
      double x;
      while ( is >> x ) v.push_back( x );
      return v;
    };
    std::vector< CM::BlockGeometry > out;
    num_bins = 0;
    size_t i = 0u;
    while ( i < lines.size() ) {
      std::istringstream hs( lines[ i ] );
      std::vector< std::string > parts;
      std::string w;
      while ( hs >> w ) parts.push_back( w );
      if ( parts.empty() || parts[ 0 ] != "block" ) {
        throw std::runtime_error( "expected 'block': " + lines[ i ] );
      }
      CM::BlockGeometry g;
      g.block_index = std::stoi( parts[ 1 ] );
      g.selection = "0";
      for ( size_t k = 0u; k + 1u < parts.size(); ++k ) {
        if ( parts[ k ] == "sel" ) {
          g.selection = parts[ k + 1u ];
          parts.resize( k );
          break;
        }
      }
      g.offset = num_bins;
      ++i;
      if ( parts.size() == 3u ) {
        g.var = "v" + parts[ 2 ];
        g.edges = numbers( lines[ i++ ] );
        g.num_bins = g.edges.size() - 1u;
      }
      else if ( parts.size() == 4u ) {
        g.is_2d = true;
        g.var = "v" + parts[ 2 ];
        g.inner_var = "v" + parts[ 3 ];
        std::map< int, std::pair< std::pair< double, double >, std::vector< double > > > sl;
        while ( i < lines.size() && lines[ i ].rfind( "block", 0 ) != 0 ) {
          auto h = numbers( lines[ i++ ] );
          sl[ static_cast< int >( h.at( 0 ) ) ] = { { h.at( 1 ), h.at( 2 ) },
            numbers( lines[ i++ ] ) };
        }
        for ( const auto& s : sl ) {
          if ( g.edges.empty() ) g.edges.push_back( s.second.first.first );
          if ( std::abs( g.edges.back() - s.second.first.first ) > 1e-9 ) {
            throw std::runtime_error( "non-contiguous slices in block "
              + parts[ 1 ] );
          }
          g.edges.push_back( s.second.first.second );
          g.inner_edges.push_back( s.second.second );
          g.num_bins += s.second.second.size() - 1u;
        }
      }
      else throw std::runtime_error( "bad block header: " + lines[ i - 1 ] );
      num_bins += g.num_bins;
      out.push_back( g );
    }
    return out;
  }

  int run_config( const std::string& cfg, const std::string& out_txt ) {
    int n = 0;
    auto blocks = parse_config( cfg, n );
    long long rows = 0;
    auto g = CM::structural_gram( blocks, n, rows );
    int r = CM::psd_rank( g );
    std::cout << n << " bins, " << rows << " rows, rank " << r << ", "
      << n - r << " null directions\n";
    if ( !out_txt.empty() ) {
      std::ofstream o( out_txt );
      o << std::setprecision( 17 );
      for ( int a = 0; a < n; ++a ) {
        for ( int b = 0; b < n; ++b ) o << g( a, b ) << ( b + 1 < n ? ' ' : '\n' );
      }
      std::cout << "Omega^T Omega written to " << out_txt << '\n';
    }
    return 0;
  }

  // Fraction of the norm of m outside the range spanned by the columns of V
  double outside_fraction( const TMatrixDSym& m, const TMatrixD& V ) {
    TMatrixD M( m );
    TMatrixD Vt( TMatrixD::kTransposed, V );
    TMatrixD PM( V, TMatrixD::kMult, TMatrixD( Vt, TMatrixD::kMult, M ) );
    TMatrixD R( M, TMatrixD::kMinus, PM );
    double num = 0., den = 0.;
    for ( int a = 0; a < M.GetNrows(); ++a ) {
      for ( int b = 0; b < M.GetNcols(); ++b ) {
        num += R( a, b ) * R( a, b );
        den += M( a, b ) * M( a, b );
      }
    }
    return den > 0. ? std::sqrt( num / den ) : 0.;
  }

  int run_compare( const std::string& univ_name, const std::string& struct_name ) {
    TFile sf( struct_name.c_str(), "read" );
    TMatrixDSym* S = nullptr;
    sf.GetObject( CM::STRUCTURAL_NAME.c_str(), S );
    auto* n_ord = dynamic_cast< TParameter< int >* >( sf.Get( "num_ordinary_reco_bins" ) );
    if ( !S || !n_ord ) throw std::runtime_error( "no structural matrix in "
      + struct_name );
    const int n = S->GetNrows();
    const int no = n_ord->GetVal();
    const int rs = CM::psd_rank( *S );
    TMatrixD Vs = CM::range_basis( *S );
    std::cout << "Structural: " << no << " ordinary reco bins (" << n
      << " reco bins), rank " << rs << ", " << no - rs << " null directions\n\n";

    TFile uf( univ_name.c_str(), "read" );
    TMatrixDSym total( n );
    int num_files = 0;
    std::cout << std::left << std::setw( 70 ) << "file" << std::right
      << std::setw( 10 ) << "patterns" << std::setw( 7 ) << "rank"
      << std::setw( 14 ) << "outside S" << '\n';
    // Only the TDirectoryFiles are opened: the per-file folders hold
    // thousands of universe histograms that are not needed here
    auto subdirs = []( TDirectory* d ) {
      std::vector< TDirectory* > out;
      TIter next( d->GetListOfKeys() );
      while ( auto* key = static_cast< TKey* >( next() ) ) {
        if ( std::string( key->GetClassName() ) != "TDirectoryFile" ) continue;
        if ( auto* sd = d->GetDirectory( key->GetName() ) ) out.push_back( sd );
      }
      return out;
    };
    for ( auto* top : subdirs( &uf ) ) {
      for ( auto* sub : subdirs( top ) ) {
        TMatrixDSym* M = nullptr;
        sub->GetObject( CM::MC_NAME.c_str(), M );
        if ( !M ) continue;
        if ( M->GetNrows() != n ) throw std::runtime_error( "size mismatch in "
          + std::string( sub->GetName() ) );
        auto* np = dynamic_cast< TParameter< int >* >(
          sub->Get( CM::MC_NUM_PATTERNS_NAME.c_str() ) );
        std::string label = sub->GetName();
        if ( label.size() > 68u ) label = "..." + label.substr( label.size() - 65u );
        std::cout << std::left << std::setw( 70 ) << label << std::right
          << std::setw( 10 ) << ( np ? np->GetVal() : -1 ) << std::setw( 7 )
          << ( np && np->GetVal() > 0 ? CM::psd_rank( *M ) : 0 )
          << std::setw( 14 ) << std::scientific << std::setprecision( 2 )
          << outside_fraction( *M, Vs ) << std::defaultfloat << '\n';
        total += *M;
        ++num_files;
      }
    }
    if ( num_files == 0 ) throw std::runtime_error( "no " + CM::MC_NAME
      + " in " + univ_name + " (made before the combination matrices were"
      " added to univmake?)" );
    const int rt = CM::psd_rank( total );
    std::cout << "\nAll " << num_files << " files: rank " << rt << ", "
      << no - rt << " null directions; fraction outside the structural range "
      << std::scientific << std::setprecision( 2 ) << outside_fraction( total, Vs )
      << std::defaultfloat << '\n';
    if ( rt < rs ) {
      std::cout << "The selected events populate " << rs - rt << " fewer"
        " directions than the binning geometry allows.\n";
    }
    else if ( rt > rs ) {
      std::cout << "The selected events populate " << rt - rs << " more"
        " directions than the structural matrix: some events fill a pattern"
        " the geometry excludes (e.g. a value outside a block's range).\n";
    }
    else std::cout << "MC and structural ranks agree.\n";
    return 0;
  }

}

int main( int argc, char* argv[] ) {
  std::vector< std::string > args( argv + 1, argv + argc );
  try {
    if ( args.size() >= 2u && args[ 0 ] == "--config" ) {
      return run_config( args[ 1 ], args.size() > 2u ? args[ 2 ] : "" );
    }
    if ( args.size() == 2u ) return run_compare( args[ 0 ], args[ 1 ] );
  }
  catch ( const std::exception& e ) {
    std::cerr << "ERROR: " << e.what() << '\n';
    return 1;
  }
  std::cout << "Usage: CombinationCheck UNIV_FILE STRUCTURAL_FILE\n"
    "       CombinationCheck --config BINNING_TXT [OUT_TXT]\n";
  return 1;
}
