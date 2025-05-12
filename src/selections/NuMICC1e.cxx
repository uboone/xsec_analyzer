// XSecAnalyzer includes
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/Selections/NuMICC1e.hh"

NuMICC1e::NuMICC1e() : SelectionBase( "NuMICC1e" ) {
  // FV definition as in PeLEE analysis
  // x_min, x_max, y_min, y_max, z_min, z_max
  this->define_fv( 10., 246., -105., 105., 10., 1026. );
}

std::string NuMICC1e::categorize_event( AnalysisEvent& ev ) {

  int nu_pdg, ccnc, npi0, nelec;
  float nu_x, nu_y, nu_z, Ee;
  const auto& in = ev.in();
  in.at( "nu_pdg" ) >> nu_pdg;
  in.at( "ccnc" ) >> ccnc;
  in.at( "npi0" ) >> npi0;
  in.at( "nelec" ) >> nelec;
  in.at( "elec_e" ) >> Ee;
  in.at( "true_nu_vtx_x" ) >> nu_x;
  in.at( "true_nu_vtx_y" ) >> nu_y;
  in.at( "true_nu_vtx_z" ) >> nu_z;

  // Require signal events to be inside the true fiducial volume
  bool sig_inFV = this->get_fv().is_inside( nu_x, nu_y, nu_z );

  // Require an incident electron (anti)neutrino
  bool sig_isNuE = ( std::abs( nu_pdg ) == ELECTRON_NEUTRINO );

  // Require a charged-current interaction
  bool sig_isCC = ( ccnc == CHARGED_CURRENT );

  // Require a final-state electron above threshold
  bool sig_has_fs_electron = ( nelec > 0 ); // 30 MeV threshold

  bool is_signal = sig_inFV && sig_isNuE && sig_isCC && sig_has_fs_electron;


  // Evaluate the true kinematic variables of interest
  float mc_electron_energy = BOGUS;
  // Check if there is a true final-state electron in this event. If there is
  // one, then store its energy
  bool has_true_electron = ( sig_isNuE && sig_isCC );
  if ( has_true_electron ) mc_electron_energy = Ee;

  // Save truth information to the output TTree
  auto& out = ev.out();
  out[ "mc_is_nue" ] = sig_isNuE;
  out[ "mc_is_CC" ] = sig_isCC;
  out[ "mc_vertex_in_FV" ] = sig_inFV;
  out[ "mc_has_fs_electron" ] = sig_has_fs_electron;
  out[ "mc_is_signal" ] = is_signal;
  out[ "mc_electron_energy" ] = mc_electron_energy;

  // All events outside of the true fiducial volume should be categorized
  // as "out of fiducial volume"
  if ( !sig_inFV ) return "OOFV";
  // NC categories
  else if ( !sig_isCC ) {
    if ( npi0 > 0 ) return "NC #pi^{0}";
    else return "NC Other";
  }
  // CC muon (anti)neutrinos
  else if ( std::abs( nu_pdg ) == MUON_NEUTRINO ) {
    if ( npi0 > 0 ) return "#nu_{#mu} CC #pi^{0}";
    else return "#nu_{#mu} CC Other";
  }
  // CC electron (anti)neutrinos
  else if ( sig_isNuE ) {
    // signal events
    if ( is_signal ) return "#nu_{e} CC Signal";
    // non-signal nues
    else return "#nu_{e} CC Other";
  }
  // We shouldn't ever get here, but return "Unknown" just in case
  std::cout << "Warning: Unknown event! Check the categorization logic.\n";
  return "Unknown";
}

bool NuMICC1e::is_selected( AnalysisEvent& ev ) {

  // Get access to the reco information needed to apply the selection
  const auto& in = ev.in();
  unsigned int shr_id;
  int nslice, n_showers, n_tracks;
  float nu_vx, nu_vy, nu_vz, contained_frac, topo_score,
    cosmic_ip, shr_energy_cali, shr_score, hits_ratio, shrmoliereavg,
    shr_tkfit_gap10_dedx_Y, shr_distance;

  in.at( "shr_id" ) >> shr_id;
  in.at( "nslice" ) >> nslice;
  in.at( "n_showers" ) >> n_showers;
  in.at( "n_tracks" ) >> n_tracks;
  in.at( "reco_nu_vtx_sce_x" ) >> nu_vx;
  in.at( "reco_nu_vtx_sce_y" ) >> nu_vy;
  in.at( "reco_nu_vtx_sce_z" ) >> nu_vz;
  in.at( "contained_fraction" ) >> contained_frac;
  in.at( "topological_score" ) >> topo_score;
  in.at( "CosmicIP" ) >> cosmic_ip;
  in.at( "shr_energy_cali" ) >> shr_energy_cali;
  in.at( "shr_score" ) >> shr_score;
  in.at( "hits_ratio" ) >> hits_ratio;
  in.at( "shrmoliereavg" ) >> shrmoliereavg;
  in.at( "shr_tkfit_gap10_dedx_Y" ) >> shr_tkfit_gap10_dedx_Y;
  in.at( "shr_distance" ) >> shr_distance;

  std::vector< unsigned int >* gen_vec;
  in.at( "pfp_generation_v" ) >> gen_vec;

  // PRE-SELECTION
  // neutrino slice
  bool has_nu_slice = ( nslice == 1 );
  // vertex inside FV
  bool in_fv = this->get_fv().is_inside( nu_vx, nu_vy, nu_vz );
  // at least one shower
  bool has_shower = ( n_showers >= 1 );
  // contained fraction
  bool contained_cut_ok = ( contained_frac >= 0.85 );

  bool sel_pass_preselection = has_nu_slice && in_fv
    && has_shower && contained_cut_ok;

  // COSMIC REJECTION
  // topological score
  bool topo_ok = ( topo_score >= 0.2 );
  // cosmic impact parameter
  bool cip_ok = ( cosmic_ip >= 10 );

  bool sel_pass_cosmic_rejection = topo_ok && cip_ok;

  // SHOWER IDENTIFICATION
  // valid shower ID
  bool valid_ID = ( shr_id != 0 ); // zero is the default value (not filled)
  // shower pfp generation
  // NOTE: The second expression after the && is only evaluated if we have
  // a valid shower ID, thus avoiding any problems with an invalid
  // argument to the std::vector::at() function.
  bool second_generation = valid_ID && gen_vec->at( shr_id - 1 ) == 2;
  // shower energy (threshold matches the signal definition)
  bool energy_ok = ( ( shr_energy_cali / 0.83 ) >= 0.03 );
  // shower score
  bool score_ok = ( shr_score <= 0.15 );
  // shower hits ratio
  bool hits_ratio_ok = ( hits_ratio >= 0.5 );

  bool sel_pass_shower_identification = valid_ID && second_generation
    && energy_ok && score_ok && hits_ratio_ok;

  // ELECTRON IDENTIFICATION
  // moliere average angle
  bool moliere_ok = ( shrmoliereavg <= 7. );
  // shower distance and dE/dx (default to passing the cut)
  bool dist_and_dEdx_ok = true;
  if ( n_tracks > 0 ) {
    // track present, 2D distance-dE/dx cut
    if ( shr_tkfit_gap10_dedx_Y >= 0. && shr_tkfit_gap10_dedx_Y < 1.75 ) {
      if ( shr_distance > 3.0 ) dist_and_dEdx_ok = false;
    }
    else if ( shr_tkfit_gap10_dedx_Y >= 1.75 && shr_tkfit_gap10_dedx_Y < 2.5 ) {
      if ( shr_distance > 12.0 ) dist_and_dEdx_ok = false;
    }
    else if ( shr_tkfit_gap10_dedx_Y >= 2.5 && shr_tkfit_gap10_dedx_Y < 3.5 ) {
      if ( shr_distance > 3.0 ) dist_and_dEdx_ok = false;
    }
    else if ( shr_tkfit_gap10_dedx_Y >= 3.5 && shr_tkfit_gap10_dedx_Y < 4.7 ) {
      dist_and_dEdx_ok = false;
    }
    else if ( shr_tkfit_gap10_dedx_Y >= 4.7 ) {
      if ( shr_distance > 3.0 ) dist_and_dEdx_ok = false;
    }
    else dist_and_dEdx_ok = false;
  }
  else {
    // no track, 1D dE/dx cut
    if ( shr_tkfit_gap10_dedx_Y < 1.7 ) dist_and_dEdx_ok = false;
    if ( shr_tkfit_gap10_dedx_Y > 2.7 && shr_tkfit_gap10_dedx_Y < 5.5 ) {
      dist_and_dEdx_ok = false;
    }
  }

  bool sel_pass_electron_identification = moliere_ok && dist_and_dEdx_ok;

  // Flag indicating whether the event passed the full selection
  bool sel_nu_e_cc = sel_pass_preselection && sel_pass_cosmic_rejection
    && sel_pass_shower_identification && sel_pass_electron_identification;

  // Set the reco energy of the electron candidate if we found one
  double reco_electron_energy = BOGUS;
  if ( sel_pass_shower_identification ) {
    // Apply the shower energy correction factor
    reco_electron_energy = shr_energy_cali / 0.83;
  }

  // We're done. Store the results in the output tree and return whether the
  // event passed the full selection
  auto& out = ev.out();
  out[ "sel_pass_preselection" ] = sel_pass_preselection;
  out[ "sel_pass_cosmic_rejection" ] = sel_pass_cosmic_rejection;
  out[ "sel_pass_shower_identification" ] = sel_pass_shower_identification;
  out[ "sel_nu_e_cc" ] = sel_nu_e_cc;
  out[ "reco_electron_energy" ] = reco_electron_energy;

  return sel_nu_e_cc;
}
