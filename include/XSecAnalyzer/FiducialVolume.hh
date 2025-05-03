#pragma once
#include "TVector3.h"
#include "Constants.hh"

class FiducialVolume {

public:

  FiducialVolume( double x_min, double x_max, double y_min, double y_max,
    double z_min, double z_max ) : x_min_( x_min ), x_max_( x_max ),
    y_min_( y_min ), y_max_( y_max ), z_min_( z_min ), z_max_( z_max ) {}

  FiducialVolume() : FiducialVolume( BOGUS, BOGUS, BOGUS,
    BOGUS, BOGUS, BOGUS ) {}

  // Use a template to allow arbitrary numerical types as input
  template < typename Number > bool is_inside( Number x, Number y, Number z )
    const
  {
    bool x_inside = ( x_min_ < x ) && ( x < x_max_ );
    bool y_inside = ( y_min_ < y ) && ( y < y_max_ );
    bool z_inside = ( z_min_ < z ) && ( z < z_max_ );

    return ( x_inside && y_inside && z_inside );
  }

  inline bool is_inside( const TVector3& pos ) const {
    return this->is_inside( pos.X(), pos.Y(), pos.Z() );
  }

  // Returns the number of argon atoms inside the volume
  double num_Ar_targets() const;

  // Returns the total BNB or NuMI neutrino flux (nu / cm^2) in the fiducial
  // volume as a function of a given beam exposure (measured in
  // protons-on-target)
  // NOTE: This is currently approximated using the flux in the *active volume*.
  // TODO: Revisit this approximation
  // TODO: Make beam and flavor configurable rather than relying on hard-coding
  double integrated_numu_flux( double pot ) const;

protected:

  double x_min_, x_max_, y_min_, y_max_, z_min_, z_max_;
};
