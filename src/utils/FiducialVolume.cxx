// XSecAnalyzer includes
#include "XSecAnalyzer/FiducialVolume.hh"

// Returns the number of Ar nuclei inside the fiducial volume
double FiducialVolume::num_Ar_targets() const {
  double volume = ( x_max_ - x_min_ ) * ( y_max_ - y_min_ )
    * ( z_max_ - z_min_ ); // cm^3
  constexpr double m_mol_Ar = 39.948; // g/mol
  constexpr double num_avogadro = 6.02214076e23; // mol^(-1)
  constexpr double mass_density_LAr = 1.3836; // g/cm^3

  double num_Ar = volume * mass_density_LAr * num_avogadro / m_mol_Ar;
  return num_Ar;
}

double FiducialVolume::integrated_numu_flux( double pot ) const {
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
