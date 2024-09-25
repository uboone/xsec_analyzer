// xsec_analyzer utility for automated creation of bin and slice configuration
// files by Liang Liu <liangliu@fnal.gov>

// Standard library includes
#include <cstdlib>
#include <iostream>

// Gnu Portability library (Gnulib) includes
#include <getopt.h>

// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/MakeConfig.hh"

// ROOT includes
#include "TRint.h"

static int verbose_flag;

int main( int argc, char** argv ) {

  int c;
  TRint app( "app", nullptr, nullptr );
  MakeConfig mm;
  mm.BinScheme();
  while ( true ) {

    static struct option long_options[] =
    {
      // These options set a flag
      {"config", no_argument, 0, 'c'},
      {"save", no_argument, 0, 's'},
      {"help", no_argument, 0, 'h'},

      {0, 0, 0, 0}
    };
    // getopt_long stores the option index here
    int option_index = 0;

    c = getopt_long( argc, argv, "csh", long_options, &option_index );

    if( c == -1 ) break;

    // Detect the end of the options
    switch ( c )
    {
      case 'c':
        // Print (plot) a response matrix that facilitates binning schemes
        mm.ResPlots();
        break;
      case 's':
        // Save binning configuration into text files
        mm.Print();
        break;
      case 'h':
      case '?':
      default:
        std::cout << "Usage: [options] \n";
        std::cout << "Options: \n";
        std::cout << "    -c, --config; Print (plot) a response matrix"
          << " that facilitates binning schemes.  \n";
        std::cout << "    -s, --save;   Save binning configuration into text"
          << " files \n";
        std::cout << "    -h, --help;   Print this help information. \n";
        abort ();
    }

  } // option parsing loop

  // Instead of reporting '--verbose' and '--brief' as they are encountered,
  // we report the final status resulting from them
  if ( verbose_flag ) puts( "verbose flag is set" );

  // Print any remaining command line arguments (not options)
  if ( optind < argc ) {
    printf( "non-option ARGV-elements: " );
    while ( optind < argc ) {
      printf( "%s ", argv[optind++] );
    }
    putchar( '\n' );
  }
  if ( argc <= 1 ) {
    std::cout << "Usage: [options] \n";
    std::cout << "Options: \n";
    std::cout << "    -c, --config; Print (plot) a response matrix that"
      << " facilitates binning schemes.  \n";
    std::cout << "    -s, --save;   Save binning configuration into text"
      << " files \n";
    std::cout << "    -h, --help;   Print this help information. \n";
    abort();
  }

  app.Run();
  return 0;
}
