#pragma once
#include "TVector3.h"
#include "Constants.h"

// Definition of the fiducial volume for the analysis, together with a utility
// function template to check if a given point lies within it or not

//DB
//Fiducial volume definitions are actually different between Stephen's and Samantha's analyses. See the definition for Stephen's above
//(now being used as the default). Samantha's defintion is here: https://github.com/ssfehlberg/CC2p-Event-Selection/blob/9492ff121a2eb884f464e1c166d067f217a04900/PeLEE_ntuples/helper_funcs.h#L18
//Discusions with Afro stated that as long as a consistent definition of FV being used, it should be ok:
//https://microboone.slack.com/archives/D04A8CUB1EW/p1691423102739129 . Thus using Stephen's defintion as the default

// Use a template here so that this function can take float or double values as
// input
template <typename Number> bool point_inside_FV( Number x, Number y, Number z )
{
  bool x_inside_FV = ( FV_X_MIN < x ) && ( x < FV_X_MAX );
  bool y_inside_FV = ( FV_Y_MIN < y ) && ( y < FV_Y_MAX );
  bool z_inside_FV = ( FV_Z_MIN < z ) && ( z < FV_Z_MAX );
  return ( x_inside_FV && y_inside_FV && z_inside_FV );
}

inline bool point_inside_FV( const TVector3& pos ) {
  return point_inside_FV( pos.X(), pos.Y(), pos.Z() );
}

// Returns the number of Ar nuclei inside the fiducial volume
inline double num_Ar_targets_in_FV() {
  double volume = ( FV_X_MAX - FV_X_MIN ) * ( FV_Y_MAX - FV_Y_MIN )
    * ( FV_Z_MAX - FV_Z_MIN ); // cm^3
  constexpr double m_mol_Ar = 39.948; // g/mol
  constexpr double N_Avogadro = 6.02214076e23; // mol^(-1)
  constexpr double mass_density_LAr = 1.3836; // g/cm^3

  double num_Ar = volume * mass_density_LAr * N_Avogadro / m_mol_Ar;
  return num_Ar;
}

// Returns the total BNB muon neutrino flux (numu / cm^2) in the fiducial
// volume as a function of a given beam exposure (measured in
// protons-on-target)
// NOTE: This is currently approximated using the flux in the *active volume*.
// TODO: Revisit this approximation
inline double integrated_numu_flux_in_FV( double pot ) {
  // Obtained using the histogram hEnumu_cv (the central-value numu flux
  // in the MicroBooNE active volume as a function of neutrino energy)
  // stored in /pnfs/uboone/persistent/uboonebeam/bnb_gsimple
  // /bnb_gsimple_fluxes_01.09.2019_463_hist/. The ROOT commmands executed
  // were
  // root [3] hEnumu_cv->Scale( 1/(4997.*5e8)/(256.35*233.) )
  // root [4] hEnumu_cv->Integral()
  // (double) 7.3762291e-10
  // See the README file in that same folder for details.
  constexpr double numu_per_cm2_per_POT_in_AV = 7.3762291e-10;
  double flux = pot * numu_per_cm2_per_POT_in_AV; // numu / cm^2
  return flux;
}
