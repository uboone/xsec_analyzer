#include "AnalysisEvent.h"

#include "Functions.h"
#include <iostream>

// Sets the signal definition flags and returns an event category based on MC
// truth information
EventCategory AnalysisEvent::categorize_event() {

  // Real data has a bogus true neutrino PDG code that is not one of the
  // allowed values (±12, ±14, ±16)
  int abs_mc_nu_pdg = std::abs( mc_nu_pdg_ );
  is_mc_ = ( abs_mc_nu_pdg == ELECTRON_NEUTRINO
    || abs_mc_nu_pdg == MUON_NEUTRINO || abs_mc_nu_pdg == TAU_NEUTRINO );
  if ( !is_mc_ ) return kUnknown;

  mc_vertex_in_FV_ = mc_vertex_inside_FV();
  mc_neutrino_is_numu_ = ( mc_nu_pdg_ == MUON_NEUTRINO );

  if ( !mc_vertex_in_FV_ ) {
    mc_is_signal_ = false;
    return kOOFV;
  }
  else if ( mc_nu_ccnc_ == NEUTRAL_CURRENT ) {
    mc_is_signal_ = false;
    return kNC;
  }
  else if ( !mc_neutrino_is_numu_ ) {
    mc_is_signal_ = false;
    if ( mc_nu_pdg_ == ELECTRON_NEUTRINO
      && mc_nu_ccnc_ == CHARGED_CURRENT ) return kNuECC;
    else return kOther;
  }

  // Set flags to their default values here
  mc_muon_in_mom_range_ = false;
  mc_lead_p_in_mom_range_ = false;
  mc_no_fs_pi0_ = true;
  mc_no_charged_pi_above_threshold_ = true;
  mc_no_fs_mesons_ = true;

  double lead_p_mom = LOW_FLOAT;

  for ( size_t p = 0u; p < mc_nu_daughter_pdg_->size(); ++p ) {
    int pdg = mc_nu_daughter_pdg_->at( p );
    float energy = mc_nu_daughter_energy_->at( p );

    // Do the general check for (anti)mesons first before considering
    // any individual PDG codes
    if ( is_meson_or_antimeson(pdg) ) {
      mc_no_fs_mesons_ = false;
    }

    // Check that the muon has a momentum within the allowed range
    if ( pdg == MUON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(MUON_MASS, 2) );
      if ( mom >= MUON_P_MIN_MOM_CUT && mom <= MUON_P_MAX_MOM_CUT ) {
        mc_muon_in_mom_range_ = true;
      }
    }
    else if ( pdg == PROTON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PROTON_MASS, 2) );
      if ( mom > lead_p_mom ) lead_p_mom = mom;
    }
    else if ( pdg == PI_ZERO ) {
      mc_no_fs_pi0_ = false;
    }
    else if ( std::abs(pdg) == PI_PLUS ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PI_PLUS_MASS, 2) );
      if ( mom > CHARGED_PI_MOM_CUT ) {
        mc_no_charged_pi_above_threshold_ = false;
      }
    }
  }

  // Check that the leading proton has a momentum within the allowed range
  if ( lead_p_mom >= LEAD_P_MIN_MOM_CUT && lead_p_mom <= LEAD_P_MAX_MOM_CUT ) {
    mc_lead_p_in_mom_range_ = true;
  }

  mc_is_signal_ = mc_vertex_in_FV_ && mc_neutrino_is_numu_
    && mc_muon_in_mom_range_ && mc_lead_p_in_mom_range_
    && mc_no_fs_mesons_;

  // Sort signal by interaction mode
  if ( mc_is_signal_ ) {
    if ( mc_nu_interaction_type_ == 0 ) return kSignalCCQE; // QE
    else if ( mc_nu_interaction_type_ == 10 ) return kSignalCCMEC; // MEC
    else if ( mc_nu_interaction_type_ == 1 ) return kSignalCCRES; // RES
    //else if ( mc_nu_interaction_type_ == 2 ) // DIS
    //else if ( mc_nu_interaction_type_ == 3 ) // COH
    else return kSignalOther;
  }
  else if ( !mc_no_fs_pi0_ || !mc_no_charged_pi_above_threshold_ ) {
    return kNuMuCCNpi;
  }
  else if ( !mc_lead_p_in_mom_range_ ) {
    return kNuMuCC0pi0p;
  }
  else return kNuMuCCOther;
}

// Sets the index of the muon candidate in the track vectors, or BOGUS_INDEX if
// one could not be found. The sel_has_muon_candidate_ flag is also set by this
// function.c
void AnalysisEvent::find_muon_candidate() {
  /*
  std::vector<int> muon_candidate_indices;
  std::vector<int> muon_pid_scores;

  for ( int p = 0; p < num_pf_particles_; ++p ) {
    // Only direct neutrino daughters (generation == 2) will be considered as
    // possible muon candidates
    unsigned int generation = pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    float track_score = pfp_track_score_->at( p );
    float start_dist = track_start_distance_->at( p );
    float track_length = track_length_->at( p );
    float pid_score = track_llr_pid_score_->at( p );

    if ( track_score > MUON_TRACK_SCORE_CUT
      && start_dist < MUON_VTX_DISTANCE_CUT
      && track_length > MUON_LENGTH_CUT
      && pid_score > MUON_PID_CUT )
    {
      muon_candidate_indices.push_back( p );
      muon_pid_scores.push_back( pid_score );
    }
  }

  //DB Samantha's analysis explicitly cuts out events with num_candidates!=1
  size_t num_candidates = muon_candidate_indices.size();
  if ( num_candidates > 0u ) sel_has_muon_candidate_CCNP_ = true;

  if ( num_candidates == 1u ) {
    muon_candidate_idx_CCNP_ = muon_candidate_indices.front();
  }
  else if ( num_candidates > 1u ) {
    // In the case of multiple muon candidates, choose the one with the highest
    // PID score (most muon-like) as the one to use
    float highest_score = LOW_FLOAT;
    int chosen_index = BOGUS_INDEX;
    for ( size_t c = 0; c < num_candidates; ++c ) {
      float score = muon_pid_scores.at( c );
      if ( highest_score < score ) {
        highest_score = score;
        chosen_index = muon_candidate_indices.at( c );
      }
    }
    muon_candidate_idx_CCNP_ = chosen_index;
  }
  else {
    muon_candidate_idx_CCNP_ = BOGUS_INDEX;
  }
  */
}

void AnalysisEvent::find_lead_p_candidate() {
  /*
  float lead_p_track_length = LOW_FLOAT;
  size_t lead_p_index = 0u;
  for ( int p = 0; p < num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    // Skip the muon candidate reco track (this function assumes that it has
    // already been found)
    if ( p == muon_candidate_idx_CCNP_ ) continue;

    // Skip PFParticles that are shower-like (track scores near 0)
    float track_score = pfp_track_score_->at( p );
    if ( track_score <= TRACK_SCORE_CUT ) continue;

    // All non-muon-candidate reco tracks are considered proton candidates
    float track_length = track_length_->at( p );
    if ( track_length <= 0. ) continue;

    if ( track_length > lead_p_track_length ) {
      lead_p_track_length = track_length;
      lead_p_index = p;
    }
  }

  // If the leading proton track length changed from its initial
  // value, then we found one. Set the index appropriately.
  if ( lead_p_track_length != LOW_FLOAT ) lead_p_candidate_idx_ = lead_p_index;
  // Otherwise, set the index to BOGUS_INDEX
  else lead_p_candidate_idx_ = BOGUS_INDEX;
  */
}

//DB We have assumed that the CCNP muon candidate is the one being used here... but that's most likely incorrect if considering CC2P selection
void AnalysisEvent::compute_observables() {
  /*
  // First compute the MC truth observables (if this is a signal MC event)
  this->compute_mc_truth_observables();

  // In cases where we failed to find a muon candidate, check whether there are
  // at least two generation == 2 PFParticles. If there are, then compute the
  // usual observables using the longest track as the muon candidate and the
  // second-longest track as the leading proton candidate. This will enable
  // sideband studies of NC backgrounds in the STV phase space.
  if ( !sel_has_muon_candidate_CCNP_ ) {

    float max_trk_len = LOW_FLOAT;
    int max_trk_idx = BOGUS_INDEX;

    float next_to_max_trk_len = LOW_FLOAT;
    int next_to_max_trk_idx = BOGUS_INDEX;

    for ( int p = 0; p < num_pf_particles_; ++p ) {

      // Only include direct neutrino daughters (generation == 2)
      unsigned int generation = pfp_generation_->at( p );
      if ( generation != 2u ) continue;

      float trk_len = track_length_->at( p );

      if ( trk_len > next_to_max_trk_len ) {

        next_to_max_trk_len = trk_len;
        next_to_max_trk_idx = p;

        if ( next_to_max_trk_len > max_trk_len ) {

          next_to_max_trk_len = max_trk_len;
          next_to_max_trk_idx = max_trk_idx;

          max_trk_len = trk_len;
          max_trk_idx = p;
        }
      }
    }

    // If we found at least two usable PFParticles, then assign the indices to
    // be used below
    if ( max_trk_idx != BOGUS_INDEX && next_to_max_trk_idx != BOGUS_INDEX ) {
      muon_candidate_idx_CCNP_ = max_trk_idx;
      lead_p_candidate_idx_ = next_to_max_trk_idx;
    }
  }

  // Abbreviate some of the calculations below by using these handy
  // references to the muon and leading proton 3-momenta
  auto& p3mu = *p3_mu_;
  auto& p3p = *p3_lead_p_;

  // Set the reco 3-momentum of the muon candidate if we found one
  bool muon = muon_candidate_idx_CCNP_ != BOGUS_INDEX;
  if ( muon ) {
    float mu_dirx = track_dirx_->at( muon_candidate_idx_CCNP_ );
    float mu_diry = track_diry_->at( muon_candidate_idx_CCNP_ );
    float mu_dirz = track_dirz_->at( muon_candidate_idx_CCNP_ );

    // The selection flag indicating whether the muon candidate is contained
    // was already set when the selection was applied. Use it to choose the
    // best momentum estimator to use.
    float muon_mom = LOW_FLOAT;
    if ( sel_muon_contained_ ) {
      muon_mom = track_range_mom_mu_->at( muon_candidate_idx_CCNP_ );
    }
    else {
      muon_mom = track_mcs_mom_mu_->at( muon_candidate_idx_CCNP_ );
    }

    p3mu = TVector3( mu_dirx, mu_diry, mu_dirz );
    p3mu = p3mu.Unit() * muon_mom;
  }

  // Set the reco 3-momentum of the leading proton candidate if we found one
  bool lead_p = lead_p_candidate_idx_ != BOGUS_INDEX;
  if ( lead_p ) {

    float p_dirx = track_dirx_->at( lead_p_candidate_idx_ );
    float p_diry = track_diry_->at( lead_p_candidate_idx_ );
    float p_dirz = track_dirz_->at( lead_p_candidate_idx_ );
    float KEp = track_kinetic_energy_p_->at( lead_p_candidate_idx_ );
    float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

    p3p = TVector3( p_dirx, p_diry, p_dirz );
    p3p = p3p.Unit() * p_mom;
  }

  // Reset the vector of reconstructed proton candidate 3-momenta
  p3_p_vec_->clear();

  // Set the reco 3-momenta of all proton candidates (i.e., all generation == 2
  // tracks except the muon candidate) assuming we found both a muon candidate
  // and at least one proton candidate.
  if ( muon && lead_p ) {
    for ( int p = 0; p < num_pf_particles_; ++p ) {
      // Skip the muon candidate
      if ( p == muon_candidate_idx_CCNP_ ) continue;

      // Only include direct neutrino daughters (generation == 2)
      unsigned int generation = pfp_generation_->at( p );
      if ( generation != 2u ) continue;

      float p_dirx = track_dirx_->at( p );
      float p_diry = track_diry_->at( p );
      float p_dirz = track_dirz_->at( p );
      float KEp = track_kinetic_energy_p_->at( p );
      float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

      TVector3 p3_temp( p_dirx, p_diry, p_dirz );
      p3_temp = p3_temp.Unit() * p_mom;

      p3_p_vec_->push_back( p3_temp );
    }

    // TODO: reduce code duplication by just getting the leading proton
    // 3-momentum from this sorted vector
    // Sort the reco proton 3-momenta in order from highest to lowest magnitude
    std::sort( p3_p_vec_->begin(), p3_p_vec_->end(), [](const TVector3& a,
      const TVector3& b) -> bool { return a.Mag() > b.Mag(); } );
  }

  // Compute reco STVs if we have both a muon candidate
  // and a leading proton candidate in the event
  if ( muon && lead_p ) {
    compute_stvs( p3mu, p3p, delta_pT_, delta_phiT_,
      delta_alphaT_, delta_pL_, pn_, delta_pTx_, delta_pTy_ );

    theta_mu_p_ = std::acos( p3mu.Dot(p3p) / p3mu.Mag() / p3p.Mag() );
  }
  */
}

void AnalysisEvent::compute_mc_truth_observables() {

  // If this is not an MC event, then just return without doing anything
  if ( !is_mc_ ) return;

  size_t num_mc_daughters = mc_nu_daughter_pdg_->size();

  // Set the true 3-momentum of the final-state muon if there is one
  bool true_muon = ( mc_neutrino_is_numu_ && mc_nu_ccnc_ == CHARGED_CURRENT );
  if ( true_muon ) {
    // Loop over the MC neutrino daughters, find the muon, and get its
    // true 3-momentum. Note that we assume there is only one muon in
    // this loop.
    
    bool found_muon = false;
    for ( size_t d = 0u; d < num_mc_daughters; ++d ) {
      int pdg = mc_nu_daughter_pdg_->at( d );
      if ( pdg == MUON ) {
        found_muon = true;
        float px = mc_nu_daughter_px_->at( d );
        float py = mc_nu_daughter_py_->at( d );
        float pz = mc_nu_daughter_pz_->at( d );
        *mc_p3_mu_ = TVector3( px, py, pz );
        break;
      }
    }

    if ( !found_muon ) {
      std::cout << "WARNING: Missing muon in MC signal event!\n";
      return;
    }
  }

  // Reset the vector of true MC proton 3-momenta
  mc_p3_p_vec_->clear();

  // Set the true 3-momentum of the leading proton (if there is one)
  float max_mom = LOW_FLOAT;
  for ( int p = 0; p < num_mc_daughters; ++p ) {
    int pdg = mc_nu_daughter_pdg_->at( p );
    if ( pdg == PROTON )
    {
      float px = mc_nu_daughter_px_->at( p );
      float py = mc_nu_daughter_py_->at( p );
      float pz = mc_nu_daughter_pz_->at( p );
      TVector3 temp_p3 = TVector3( px, py, pz );

      mc_p3_p_vec_->push_back( temp_p3 );

      float mom = temp_p3.Mag();
      if ( mom > max_mom ) {
        max_mom = mom;
        *mc_p3_lead_p_ = temp_p3;
      }
    }
  }

  // TODO: reduce code duplication by just getting the leading proton
  // 3-momentum from this sorted vector
  // Sort the true proton 3-momenta in order from highest to lowest magnitude
  std::sort( mc_p3_p_vec_->begin(), mc_p3_p_vec_->end(), [](const TVector3& a,
    const TVector3& b) -> bool { return a.Mag() > b.Mag(); } );

  // If the event contains a leading proton, then set the 3-momentum
  // accordingly
  bool true_lead_p = ( max_mom != LOW_FLOAT );
  if ( !true_lead_p && mc_is_signal_ ) {
    // If it doesn't for a signal event, then something is wrong.
    std::cout << "WARNING: Missing leading proton in MC signal event!\n";
    return;
  }

  // Compute true STVs if the event contains both a muon and a leading
  // proton
  if ( true_muon && true_lead_p ) {
    compute_stvs( *mc_p3_mu_, *mc_p3_lead_p_, mc_delta_pT_, mc_delta_phiT_,
      mc_delta_alphaT_, mc_delta_pL_, mc_pn_, mc_delta_pTx_, mc_delta_pTy_ );

    mc_theta_mu_p_ = std::acos( mc_p3_mu_->Dot(*mc_p3_lead_p_)
      / mc_p3_mu_->Mag() / mc_p3_lead_p_->Mag() );
  }
}

bool AnalysisEvent::reco_vertex_inside_FV() {
  return point_inside_FV( nu_vx_, nu_vy_, nu_vz_ );
}

bool AnalysisEvent::mc_vertex_inside_FV() {
  return point_inside_FV( mc_nu_vx_, mc_nu_vy_, mc_nu_vz_ );
}

bool AnalysisEvent::in_proton_containment_vol( float x, float y, float z ) {
  bool x_inside_PCV = ( PCV_X_MIN < x ) && ( x < PCV_X_MAX );
  bool y_inside_PCV = ( PCV_Y_MIN < y ) && ( y < PCV_Y_MAX );
  bool z_inside_PCV = ( PCV_Z_MIN < z ) && ( z < PCV_Z_MAX );
  return ( x_inside_PCV && y_inside_PCV && z_inside_PCV );
}
