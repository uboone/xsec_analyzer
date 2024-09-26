// Standard library includes
#include <iostream>

// XSecAnalyzer includes
#include "XSecAnalyzer/StandaloneUnfolding.hh"

int main( int argc, char** argv ) {

  if ( argc != 2 ) {
    std::cout << "Usage: standalone_unfold CONFIG_FILE\n";
    return 1;
  }

  StandaloneUnfolding my_su( argv[1] );
  my_su.run_unfolding();
  return 0;
}
