#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include "Constants.h"
#include "EventCategory.hh"
#include "AnalysisEvent.h"
#include "FiducialVolume.hh"

#include "TVector2.h"
#include "TVector3.h"

// Helper function that avoids NaNs when taking square roots of negative
// numbers
inline double real_sqrt( double x ) {
  if ( x < 0. ) return 0.;
  else return std::sqrt( x );
}

// Helper function that returns true if a given PDG code represents a meson or
// antimeson. Otherwise returns false. Based on points 10, 12, and 13 of the
// Particle Data Group's "Monte Carlo Particle Numbering Scheme"
// (2019 revision).

inline bool is_meson_or_antimeson( int pdg_code ) {
  // Ignore differences between mesons and antimesons for this test. Mesons
  // will have positive PDG codes, while antimesons will have negative ones.
  int abs_pdg = std::abs( pdg_code );

  // Meson PDG codes have no more than seven digits. Seven-digit
  // codes beginning with "99" are reserved for generator-specific
  // particles
  if ( abs_pdg >= 9900000 ) return false;

  // Mesons have a value of zero for $n_{q1}$, the thousands digit
  int thousands_digit = ( abs_pdg / 1000 ) % 10;
  if ( thousands_digit != 0 ) return false;

  // They also have a nonzero value for $n_{q2}$, the hundreds digit
  int hundreds_digit = ( abs_pdg / 100 ) % 10;
  if ( hundreds_digit == 0 ) return false;

  // Reserved codes for Standard Model parton distribution functions
  if ( abs_pdg >= 901 && abs_pdg <= 930 ) return false;

  // Reggeon and pomeron
  if ( abs_pdg == 110 || abs_pdg == 990 ) return false;

  // Reserved codes for GEANT tracking purposes
  if ( abs_pdg == 998 || abs_pdg == 999 ) return false;

  // Reserved code for generator-specific pseudoparticles
  if ( abs_pdg == 100 ) return false;

  // If we've passed all of the tests above, then the particle is a meson
  return true;
}

// Function that defines the track-length-dependent proton PID cut
inline double proton_pid_cut( double track_length ) {

  double cut = DEFAULT_PROTON_PID_CUT;

  // Piecewise cut removed 27 June 2021
  //// All track length values are in cm
  //if ( track_length >= 0. && track_length <= 10.5 ) {
  //  cut = -0.0034219305*std::pow( track_length, 2 );
  //  cut += 0.018436866*track_length + 0.062718401;
  //}
  //else if ( track_length > 10.5 && track_length <= 33.1776508 ) {
  //  cut = 0.014153245*( track_length - 10.5 ) - 0.12096235;
  //}
  
  return cut;
}

// Helper function for computing STVs (either reco or true)
inline void compute_stvs( const TVector3& p3mu, const TVector3& p3p, float& delta_pT,
  float& delta_phiT, float& delta_alphaT, float& delta_pL, float& pn,
  float& delta_pTx, float& delta_pTy )
{
  delta_pT = (p3mu + p3p).Perp();

  delta_phiT = std::acos( (-p3mu.X()*p3p.X() - p3mu.Y()*p3p.Y())
    / (p3mu.XYvector().Mod() * p3p.XYvector().Mod()) );

  TVector2 delta_pT_vec = (p3mu + p3p).XYvector();
  delta_alphaT = std::acos( (-p3mu.X()*delta_pT_vec.X()
    - p3mu.Y()*delta_pT_vec.Y())
    / (p3mu.XYvector().Mod() * delta_pT_vec.Mod()) );

  float Emu = std::sqrt(std::pow(MUON_MASS, 2) + p3mu.Mag2());
  float Ep = std::sqrt(std::pow(PROTON_MASS, 2) + p3p.Mag2());
  float R = TARGET_MASS + p3mu.Z() + p3p.Z() - Emu - Ep;

  // Estimated mass of the final remnant nucleus (CCQE assumption)
  float mf = TARGET_MASS - NEUTRON_MASS + BINDING_ENERGY;
  delta_pL = 0.5*R - (std::pow(mf, 2) + std::pow(delta_pT, 2)) / (2.*R);

  pn = std::sqrt( std::pow(delta_pL, 2) + std::pow(delta_pT, 2) );

  // Components of the 2D delta_pT vector (see arXiv:1910.08658)
  
  // We assume that the neutrino travels along the +z direction (also done
  // in the other expressions above)
  TVector3 zUnit( 0., 0., 1. );

  // Defines the x direction for the components of the delta_pT vector
  TVector2 xTUnit = zUnit.Cross( p3mu ).XYvector().Unit();

  delta_pTx = xTUnit.X()*delta_pT_vec.X() + xTUnit.Y()*delta_pT_vec.Y();

  // Defines the y direction for the components of the delta_T vector
  TVector2 yTUnit = ( -p3mu ).XYvector().Unit();

  delta_pTy = yTUnit.X()*delta_pT_vec.X() + yTUnit.Y()*delta_pT_vec.Y();
}


//The only insta-return values should be Unknown (i.e. data) or OOFV
inline EventCategory categorize_event(AnalysisEvent* Event, FiducialVolume FV) {

  // Real data has a bogus true neutrino PDG code that is not one of the
  // allowed values (±12, ±14, ±16)
  int abs_mc_nu_pdg = std::abs( Event->mc_nu_pdg_ );
  Event->is_mc_ = ( abs_mc_nu_pdg == ELECTRON_NEUTRINO || abs_mc_nu_pdg == MUON_NEUTRINO || abs_mc_nu_pdg == TAU_NEUTRINO );
  if ( !Event->is_mc_ ) {
    return kUnknown;
  }

  if (abs_mc_nu_pdg == TAU_NEUTRINO) {
    std::cerr << "Did not expect to be dealing with nutaus. Currently defined EventCategory as kOther" << std::endl;
    return kOther;
  }

  bool MCVertexInFV = point_inside_FV(FV, Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_);
  if ( !MCVertexInFV ) {
    return kOOFV;
  }

  bool isNC = (Event->mc_nu_ccnc_ == NEUTRAL_CURRENT);
  //DB Currently only one NC category is supported so test first. Will likely want to change this in the future
  if (isNC) return kNC;
  
  //DB According to P. Green (17/06/23), nu taus are not considered. There use this boolean as a switch for nue/numu
  bool isNumu = (Event->mc_nu_pdg_ == MUON_NEUTRINO);  

  if ( isNumu ) {
      if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCCQE; // QE
      else if ( Event->mc_nu_interaction_type_ == 10 ) return kNuMuCCMEC; // MEC
      else if ( Event->mc_nu_interaction_type_ == 1 ) return kNuMuCCRES; // RES
      //else if ( mc_nu_interaction_type_ == 2 ) // DIS
      //else if ( mc_nu_interaction_type_ == 3 ) // COH
  } else {
      if ( Event->mc_nu_interaction_type_ == 0 ) return kNuECCQE; // QE
      else if ( Event->mc_nu_interaction_type_ == 10 ) return kNuECCMEC; // MEC
      else if ( Event->mc_nu_interaction_type_ == 1 ) return kNuECCRES; // RES
      //else if ( mc_nu_interaction_type_ == 2 ) // DIS
      //else if ( mc_nu_interaction_type_ == 3 ) // COH
  }

  //Assumed that if nothing has been selected so far, thus kOther
  return kOther;
  
}

#endif
