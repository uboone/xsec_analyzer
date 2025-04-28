// XSecAnalyzer includes
#include "XSecAnalyzer/Constants.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/KICalculator.hh"
#include "XSecAnalyzer/Selections/CC1muNp0pi.hh"

CC1muNp0pi::CC1muNp0pi() : SelectionBase( "CC1muNp0pi" ) {
  ki_calc_type_ = KICalculator::kOpt1;
  this->define_fv( 21.5, 234.85, -95.0, 95.0, 21.5, 966.8 );
}

bool CC1muNp0pi::is_selected( AnalysisEvent& ev ) {

  // Get access to the input PeLEE ntuple branches
  const auto& in = ev.in();

  // Get access to the output TTree branches
  auto& out = ev.out();

  FiducialVolume pcv( 10., 246.35, -106.5, 106.5, 10., 1026.8 );

  // Check whether the space-charge-corrected position of the reconstructed
  // neutrino vertex lies within the fiducial volume
  float nu_vx, nu_vy, nu_vz;
  in.at( "reco_nu_vtx_sce_x" ) >> nu_vx;
  in.at( "reco_nu_vtx_sce_y" ) >> nu_vy;
  in.at( "reco_nu_vtx_sce_z" ) >> nu_vz;

  bool sel_reco_vertex_in_FV = this->get_fv().is_inside( nu_vx, nu_vy, nu_vz );

  // Evaluate cuts on the topological score cut and the cosmic impact parameter
  float topological_score, cosmic_ip;
  in.at( "topological_score" ) >> topological_score;
  in.at( "CosmicIP" ) >> cosmic_ip;
  bool sel_topo_cut_passed = topological_score > TOPO_SCORE_CUT;
  bool sel_cosmic_ip_cut_passed = cosmic_ip > COSMIC_IP_CUT;

  // Retrieve information about the PFParticles in the event
  int num_pfps;
  std::vector< unsigned int >* gen_vec;
  std::vector< float > *track_start_x, *track_start_y, *track_start_z;
  in.at( "n_pfps" ) >> num_pfps;
  in.at( "pfp_generation_v" ) >> gen_vec;
  in.at( "trk_sce_start_x_v" ) >> track_start_x;
  in.at( "trk_sce_start_y_v" ) >> track_start_y;
  in.at( "trk_sce_start_z_v" ) >> track_start_z;

  // Evaluate the containment cut by looping over all PFParticles and checking
  // their starting positions. Pass the cut by default.
  bool sel_pfp_starts_in_PCV = true;
  for ( int p = 0; p < num_pfps; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = gen_vec->at( p );
    if ( generation != 2u ) continue;

    // Use the track reconstruction results to get the start point for
    // every PFParticle for the purpose of verifying containment. We could
    // in principle differentiate between tracks and showers here, but
    // (1) we cut out all showers later on in the selection anyway, and
    // (2) the blinded PeLEE data ntuples do not include shower information.
    // We therefore apply the track reconstruction here unconditionally.
    float x = track_start_x->at( p );
    float y = track_start_y->at( p );
    float z = track_start_z->at( p );

    // Verify that the start of the PFParticle lies within the containment
    // volume.
    //
    // See https://stackoverflow.com/a/2488507 for an explanation of the use
    // of &= here. Don't worry, it's type-safe since both operands are bool.
    sel_pfp_starts_in_PCV &= pcv.is_inside( x, y, z );
  }

  // Search for a muon candidate among the PFParticles in the event. A cut on
  // the muon momentum is handled later.
  std::vector< float > *track_score, *track_start_distance, *track_length,
    *track_llr_pid_score;
  in.at( "trk_score_v" ) >> track_score;
  in.at( "trk_distance_v" ) >> track_start_distance;
  in.at( "trk_len_v" ) >> track_length;
  in.at( "trk_llr_pid_score_v" ) >> track_llr_pid_score;

  std::vector< int > muon_candidate_indices;
  std::vector< int > muon_pid_scores;

  bool sel_has_muon_candidate = false;

  for ( int p = 0; p < num_pfps; ++p ) {
    // Only direct neutrino daughters (generation == 2) will be considered as
    // possible muon candidates
    unsigned int generation = gen_vec->at( p );
    if ( generation != 2u ) continue;

    float trk_score = track_score->at( p );
    float start_dist = track_start_distance->at( p );
    float trk_length = track_length->at( p );
    float pid_score = track_llr_pid_score->at( p );

    if ( trk_score > MUON_TRACK_SCORE_CUT
      && start_dist < MUON_VTX_DISTANCE_CUT
      && trk_length > MUON_LENGTH_CUT
      && pid_score > MUON_PID_CUT )
    {
      muon_candidate_indices.push_back( p );
      muon_pid_scores.push_back( pid_score );
    }
  }

  size_t num_candidates = muon_candidate_indices.size();
  if ( num_candidates > 0u ) sel_has_muon_candidate = true;

  // By default, set the muon candidate index to the bogus value
  int muon_candidate_idx = BOGUS_INDEX;
  // If exactly one muon candidate was found, then assign the index of its
  // corresponding PFParticle
  if ( num_candidates == 1u ) {
    muon_candidate_idx = muon_candidate_indices.front();
  }
  // In the case of multiple muon candidates, choose the one with the highest
  // PID score (most muon-like) as the one to use
  else if ( num_candidates > 1u ) {
    float highest_score = LOW_FLOAT;
    int chosen_index = BOGUS_INDEX;
    for ( size_t c = 0; c < num_candidates; ++c ) {
      float score = muon_pid_scores.at( c );
      if ( highest_score < score ) {
        highest_score = score;
        chosen_index = muon_candidate_indices.at( c );
      }
    }
    muon_candidate_idx = chosen_index;
  }

  // Combine some of the prior cuts into a numu CC preselection
  bool sel_nu_mu_cc = sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV
    && sel_has_muon_candidate && sel_topo_cut_passed;

  // Fail the shower cut if any showers were reconstructed
  // NOTE: We could do this quicker like this,
  //   sel_no_reco_showers_ = ( num_showers_ > 0 );
  // but it might be nice to be able to adjust the track score for this cut.
  // Therefore, we do it the hard way.
  int reco_shower_count = 0;
  for ( int p = 0; p < num_pfps; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = gen_vec->at( p );
    if ( generation != 2u ) continue;

    float tscore = track_score->at( p );
    if ( tscore <= TRACK_SCORE_CUT ) ++reco_shower_count;
  }
  // Check the shower cut
  bool sel_no_reco_showers = ( reco_shower_count == 0 );

  // Set default values of the next few cuts before checking them
  bool sel_passed_proton_pid_cut = true;
  bool sel_protons_contained = true;
  bool sel_muon_contained = false;
  bool sel_muon_passed_mom_cuts = false;
  bool sel_muon_quality_ok = false;
  bool sel_has_p_candidate = false;
  bool sel_lead_p_passed_mom_cuts = false;

  std::vector< float > *track_end_x, *track_end_y, *track_end_z,
    *muon_mom_range, *muon_mom_mcs;
  in.at( "trk_sce_end_x_v" ) >> track_end_x;
  in.at( "trk_sce_end_y_v" ) >> track_end_y;
  in.at( "trk_sce_end_z_v" ) >> track_end_z;
  in.at( "trk_range_muon_mom_v" ) >> muon_mom_range;
  in.at( "trk_mcs_muon_mom_v" ) >> muon_mom_mcs;

  for ( int p = 0; p < num_pfps; ++p ) {
    // Only worry about direct neutrino daughters (PFParticles considered
    // daughters of the reconstructed neutrino)
    unsigned int generation = gen_vec->at( p );
    if ( generation != 2u ) continue;

    // Check that we can find a muon candidate in the event. If more than
    // one is found, also fail the cut.
    if ( p == muon_candidate_idx ) {

      // Check whether the muon candidate is contained. Use the same
      // containment volume as the protons. TODO: revisit this as needed.
      float endx = track_end_x->at( p );
      float endy = track_end_y->at( p );
      float endz = track_end_z->at( p );
      bool end_contained = pcv.is_inside( endx, endy, endz );

      if ( end_contained ) sel_muon_contained = true;

      // Check that the muon candidate is above threshold. Use the best
      // momentum based on whether it was contained or not.
      float muon_mom = LOW_FLOAT;
      float range_muon_mom = muon_mom_range->at( p );
      float mcs_muon_mom = muon_mom_mcs->at( p );

      if ( sel_muon_contained ) muon_mom = range_muon_mom;
      else muon_mom = mcs_muon_mom;

      if ( muon_mom >= MUON_P_MIN_MOM_CUT && muon_mom <= MUON_P_MAX_MOM_CUT ) {
        sel_muon_passed_mom_cuts = true;
      }

      // Apply muon candidate quality cut by comparing MCS and range-based
      // momentum estimators.
      double frac_diff_range_mcs = std::abs( range_muon_mom - mcs_muon_mom );
      if ( range_muon_mom > 0. ) {
        frac_diff_range_mcs /= range_muon_mom;
        if ( frac_diff_range_mcs < MUON_MOM_QUALITY_CUT ) {
          sel_muon_quality_ok = true;
        }
      }
    }
    else {

      float trk_score = track_score->at( p );
      if ( trk_score <= TRACK_SCORE_CUT ) continue;

      // Bad tracks in the searchingfornues TTree can have
      // bogus track lengths. This skips those.
      float trk_length = track_length->at( p );
      if ( trk_length <= 0. ) continue;

      // We found a reco track that is not the muon candidate. All such
      // tracks are considered proton candidates.
      sel_has_p_candidate = true;

      float llr_pid_score = track_llr_pid_score->at( p );

      // Check whether the current proton candidate fails the proton PID cut
      if ( llr_pid_score > DEFAULT_PROTON_PID_CUT ) {
        sel_passed_proton_pid_cut = false;
      }

      // Check whether the current proton candidate fails the containment cut
      float endx = track_end_x->at( p );
      float endy = track_end_y->at( p );
      float endz = track_end_z->at( p );
      bool end_contained = pcv.is_inside( endx, endy, endz );
      if ( !end_contained ) sel_protons_contained = false;
    }

  }

  // All that remains is to apply the leading proton candidate cuts. We could
  // search for it above, but doing it here makes the code more readable (with
  // likely negligible impact on performance)
  // Only do the search if we actually have a proton candidate in the event.

  std::vector< float >* track_KE_proton;
  in.at( "trk_energy_proton_v" ) >> track_KE_proton;

  int lead_p_candidate_idx = BOGUS_INDEX;
  if ( sel_has_p_candidate ) {

    float lead_p_track_length = LOW_FLOAT;
    size_t lead_p_index = 0u;
    for ( int p = 0; p < num_pfps; ++p ) {

      // Only check direct neutrino daughters (generation == 2)
      unsigned int generation = gen_vec->at( p );
      if ( generation != 2u ) continue;

      // Skip the muon candidate reco track (if we have one, it was already
      // found above)
      if ( p == muon_candidate_idx ) continue;

      // Skip PFParticles that are shower-like (track scores near 0)
      float trk_score = track_score->at( p );
      if ( trk_score <= TRACK_SCORE_CUT ) continue;

      // All non-muon-candidate reco tracks are considered proton candidates
      float trk_length = track_length->at( p );
      if ( trk_length <= 0. ) continue;

      if ( trk_length > lead_p_track_length ) {
        lead_p_track_length = trk_length;
        lead_p_index = p;
      }

    }

    // If the leading proton track length changed from its initial
    // value, then we found one. Set the index appropriately.
    if ( lead_p_track_length != LOW_FLOAT ) {
      lead_p_candidate_idx = lead_p_index;
    }

    // Check the range-based reco momentum for the leading proton candidate
    float lead_p_KE = track_KE_proton->at( lead_p_candidate_idx );
    float range_mom_lead_p = real_sqrt( lead_p_KE*lead_p_KE
      + 2.*PROTON_MASS*lead_p_KE );
    if ( range_mom_lead_p >= LEAD_P_MIN_MOM_CUT
      && range_mom_lead_p <= LEAD_P_MAX_MOM_CUT )
    {
      sel_lead_p_passed_mom_cuts = true;
    }

  } // A proton candidate was found

  // All right, we've applied all selection cuts. Set the flag that indicates
  // whether all were passed (and thus the event is selected as a CCNp0pi
  // candidate)
  bool sel_CCNp0pi = sel_nu_mu_cc && sel_no_reco_showers
    && sel_muon_passed_mom_cuts && sel_muon_contained && sel_muon_quality_ok
    && sel_has_p_candidate && sel_passed_proton_pid_cut
    && sel_protons_contained && sel_lead_p_passed_mom_cuts;

  // Store the muon and leading proton candidate indices in the output
  // tree. We do this here with the actual values (which may be BOGUS_INDEX
  // if such a candidate was not found) because we will replace them with
  // other values below if needed to allow studies of an NC-enriched sideband.
  out[ "lead_p_candidate_idx" ] = lead_p_candidate_idx;
  out[ "muon_candidate_idx" ] = muon_candidate_idx;

  // In cases where we failed to find a muon candidate, check whether there are
  // at least two generation == 2 PFParticles. If there are, then compute the
  // usual observables using the longest track as the muon candidate and the
  // second-longest track as the leading proton candidate. This will enable
  // sideband studies of NC backgrounds in the CC1mu0piNp phase space.
  if ( !sel_has_muon_candidate ) {
    float max_trk_len = LOW_FLOAT;
    int max_trk_idx = BOGUS_INDEX;

    float next_to_max_trk_len = LOW_FLOAT;
    int next_to_max_trk_idx = BOGUS_INDEX;

    for ( int p = 0; p < num_pfps; ++p ) {

      // Only include direct neutrino daughters (generation == 2)
      unsigned int generation = gen_vec->at( p );
      if ( generation != 2u ) continue;

      float trk_len = track_length->at( p );

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
      muon_candidate_idx = max_trk_idx;
      lead_p_candidate_idx = next_to_max_trk_idx;
    }
  }

  // Set the reco 3-momentum of the muon candidate if we found one
  std::vector< float > *track_dir_x, *track_dir_y, *track_dir_z;
  in.at( "trk_dir_x_v" ) >> track_dir_x;
  in.at( "trk_dir_y_v" ) >> track_dir_y;
  in.at( "trk_dir_z_v" ) >> track_dir_z;

  TVector3 p3mu;
  bool muon = muon_candidate_idx != BOGUS_INDEX;
  if ( muon ) {
    float mu_dirx = track_dir_x->at( muon_candidate_idx );
    float mu_diry = track_dir_y->at( muon_candidate_idx );
    float mu_dirz = track_dir_z->at( muon_candidate_idx );

    // The selection flag indicating whether the muon candidate is contained
    // was already set when the selection was applied. Use it to choose the
    // best momentum estimator to use.
    float muon_mom = LOW_FLOAT;
    if ( sel_muon_contained ) {
      muon_mom = muon_mom_range->at( muon_candidate_idx );
    }
    else {
      muon_mom = muon_mom_mcs->at( muon_candidate_idx );
    }

    p3mu = TVector3( mu_dirx, mu_diry, mu_dirz );
    p3mu = p3mu.Unit() * muon_mom;
  }

  // Set the reco 3-momentum of the leading proton candidate if we found one
  TVector3 p3p;
  bool lead_p = lead_p_candidate_idx != BOGUS_INDEX;
  if ( lead_p ) {

    float p_dirx = track_dir_x->at( lead_p_candidate_idx );
    float p_diry = track_dir_y->at( lead_p_candidate_idx );
    float p_dirz = track_dir_z->at( lead_p_candidate_idx );
    float KEp = track_KE_proton->at( lead_p_candidate_idx );
    float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

    p3p = TVector3( p_dirx, p_diry, p_dirz );
    p3p = p3p.Unit() * p_mom;
  }

  // Set the reco 3-momenta of all proton candidates (i.e., all generation == 2
  // tracks except the muon candidate) assuming we found both a muon candidate
  // and at least one proton candidate.
  std::vector< TVector3 > p3_p_vec;
  if ( muon && lead_p ) {
    for ( int p = 0; p < num_pfps; ++p ) {
      // Skip the muon candidate
      if ( p == muon_candidate_idx ) continue;

      // Only include direct neutrino daughters (generation == 2)
      unsigned int generation = gen_vec->at( p );
      if ( generation != 2u ) continue;

      float p_dirx = track_dir_x->at( p );
      float p_diry = track_dir_y->at( p );
      float p_dirz = track_dir_z->at( p );
      float KEp = track_KE_proton->at( p );
      float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

      TVector3 p3_temp( p_dirx, p_diry, p_dirz );
      p3_temp = p3_temp.Unit() * p_mom;

      p3_p_vec.push_back( p3_temp );
    }

    // TODO: reduce code duplication by just getting the leading proton
    // 3-momentum from this sorted vector
    // Sort the reco proton 3-momenta in order from highest to lowest magnitude
    std::sort( p3_p_vec.begin(), p3_p_vec.end(), [](const TVector3& a,
      const TVector3& b) -> bool { return a.Mag() > b.Mag(); } );
  }

  // Compute reco STVs if we have both a muon candidate
  // and a leading proton candidate in the event
  double delta_pT = BOGUS, delta_phiT = BOGUS, delta_alphaT = BOGUS,
    delta_pL = BOGUS, pn = BOGUS, delta_pTx = BOGUS, delta_pTy = BOGUS,
    theta_mu_p = BOGUS;

  if ( muon && lead_p ) {

    KICalculator ki_calc( p3mu, p3p, ki_calc_type_ );

    delta_pT = ki_calc.pT();
    delta_phiT = ki_calc.delta_phiT();
    delta_alphaT = ki_calc.delta_alphaT();
    delta_pL = ki_calc.pL();
    pn = ki_calc.pn();
    delta_pTx = ki_calc.pTx();
    delta_pTy = ki_calc.pTy();
    theta_mu_p = ki_calc.theta_lep_had();
  }

  // We're done, store the results in the output tree and return
  // the boolean flag indicating whether the full selection was passed or not
  out[ "reco_vertex_in_FV" ] = sel_reco_vertex_in_FV;
  out[ "pfp_starts_in_PCV" ] = sel_pfp_starts_in_PCV;
  out[ "has_muon_candidate" ] = sel_has_muon_candidate;
  out[ "topo_cut_passed" ] = sel_topo_cut_passed;
  out[ "nu_mu_cc" ] = sel_nu_mu_cc;
  out[ "muon_contained" ] = sel_muon_contained;
  out[ "muon_passed_mom_cuts" ] = sel_muon_passed_mom_cuts;
  out[ "no_reco_showers" ] = sel_no_reco_showers;
  out[ "has_p_candidate" ] = sel_has_p_candidate;
  out[ "muon_quality_ok" ] = sel_muon_quality_ok;
  out[ "protons_contained" ] = sel_protons_contained;
  out[ "passed_proton_pid_cut" ] = sel_passed_proton_pid_cut;
  out[ "lead_p_passed_mom_cuts" ] = sel_lead_p_passed_mom_cuts;
  out[ "cosmic_ip_cut_passed" ] = sel_cosmic_ip_cut_passed;

  out[ "reco_delta_pT" ] = delta_pT;
  out[ "reco_delta_phiT" ] = delta_phiT;
  out[ "reco_delta_alphaT" ] = delta_alphaT;
  out[ "reco_delta_pL" ] = delta_pL;
  out[ "reco_pn" ] = pn;
  out[ "reco_delta_pTx" ] = delta_pTx;
  out[ "reco_delta_pTy" ] = delta_pTy;
  out[ "reco_theta_mu_p" ] = theta_mu_p;

  out[ "reco_p3_mu" ] = p3mu;
  out[ "reco_p3_lead_p" ] = p3p;
  out[ "reco_p3_p_vec" ] = p3_p_vec;

  return sel_CCNp0pi;
}

std::string CC1muNp0pi::categorize_event( AnalysisEvent& ev ) {

  // Get access to the input and output trees
  const auto& in = ev.in();
  auto& out = ev.out();

  // Check whether the true neutrino vertex lies within the fiducial volume
  float mc_nu_vx, mc_nu_vy, mc_nu_vz;
  int mc_nu_pdg, mc_ccnc, mc_interaction_type;
  in.at( "true_nu_vtx_x" ) >> mc_nu_vx;
  in.at( "true_nu_vtx_y" ) >> mc_nu_vy;
  in.at( "true_nu_vtx_z" ) >> mc_nu_vz;
  in.at( "nu_pdg" ) >> mc_nu_pdg;
  in.at( "ccnc" ) >> mc_ccnc;
  in.at( "interaction" ) >> mc_interaction_type;

  bool sig_in_fv = this->get_fv().is_inside( mc_nu_vx, mc_nu_vy, mc_nu_vz );
  bool sig_is_numu = ( mc_nu_pdg == MUON_NEUTRINO );
  bool is_nc = ( mc_ccnc == NEUTRAL_CURRENT );

  // Set some flags to their default values before checking them
  bool sig_no_fs_mesons = true;
  bool sig_mc_no_fs_pi0 = true;
  bool sig_mc_no_charged_pi_above_threshold = true;
  bool sig_mu_in_mom_range = false;
  bool sig_lead_p_mom_in_range = false;

  std::vector< int >* nu_daughter_pdg;
  std::vector< float > *nu_daughter_energy, *nu_daughter_px,
    *nu_daughter_py, *nu_daughter_pz;

  in.at( "mc_pdg" ) >> nu_daughter_pdg;
  in.at( "mc_E" ) >> nu_daughter_energy;
  in.at( "mc_px" ) >> nu_daughter_px;
  in.at( "mc_py" ) >> nu_daughter_py;
  in.at( "mc_pz" ) >> nu_daughter_pz;

  int num_p_in_mom_range = 0;
  double mom_lead_p = 0.;

  for ( size_t p = 0u; p < nu_daughter_pdg->size(); ++p ) {
    int pdg = nu_daughter_pdg->at( p );
    float energy = nu_daughter_energy->at( p );

    // Do the general check for (anti)mesons first before considering
    // any individual PDG codes
    if ( is_meson_or_antimeson( pdg ) ) {
      sig_no_fs_mesons = false;
    }

    if ( pdg == MUON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(MUON_MASS, 2) );
      if ( mom >= MUON_P_MIN_MOM_CUT && mom <= MUON_P_MAX_MOM_CUT ) {
	sig_mu_in_mom_range = true;
      }
    }
    else if ( pdg == PROTON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PROTON_MASS, 2) );
      if ( mom > mom_lead_p ) mom_lead_p = mom;
      if ( mom >= LEAD_P_MIN_MOM_CUT && mom <= LEAD_P_MAX_MOM_CUT ) {
        ++num_p_in_mom_range;
      }
    }
    // Not part of the signal definition, but the flag is recorded anyway
    else if ( pdg == PI_ZERO ) {
      sig_mc_no_fs_pi0 = false;
    }
    // Also not part of the signal definition
    else if ( std::abs(pdg) == PI_PLUS ) {
      double mom = real_sqrt( std::pow(energy, 2)
        - std::pow(PI_PLUS_MASS, 2) );
      if ( mom > CHARGED_PI_MOM_CUT ) {
        sig_mc_no_charged_pi_above_threshold = false;
      }
    }
  }

  // Check that the leading proton has a momentum within the allowed range
  if ( mom_lead_p >= LEAD_P_MIN_MOM_CUT
    && mom_lead_p <= LEAD_P_MAX_MOM_CUT )
  {
    sig_lead_p_mom_in_range = true;
  }

  bool is_signal = sig_in_fv && !is_nc && sig_is_numu && sig_mu_in_mom_range
    && sig_lead_p_mom_in_range && sig_no_fs_mesons;

  // Set the true 3-momentum of the final-state muon if there is one
  size_t num_mc_daughters = nu_daughter_pdg->size();

  TVector3 mc_p3mu;
  bool true_muon = ( sig_is_numu && mc_ccnc == CHARGED_CURRENT );
  if ( true_muon ) {
    // Loop over the MC neutrino daughters, find the muon, and get its
    // true 3-momentum. Note that we assume there is only one muon in
    // this loop.
    bool found_muon = false;
    for ( size_t d = 0u; d < num_mc_daughters; ++d ) {
      int pdg = nu_daughter_pdg->at( d );
      if ( pdg == MUON ) {
        found_muon = true;
        float px = nu_daughter_px->at( d );
        float py = nu_daughter_py->at( d );
        float pz = nu_daughter_pz->at( d );
	mc_p3mu = TVector3( px, py, pz );
        break;
      }
    }

    if ( !found_muon ) {
      std::cout << "WARNING: Missing muon in true numu CC event!\n";
    }
  }

  // Set the true 3-momentum of the final-state protons (if there are any) and
  // the leading proton specifically
  TVector3 mc_p3p;
  std::vector< TVector3 > mc_p3_p_vec;

  float max_mom = LOW_FLOAT;
  for ( int p = 0; p < num_mc_daughters; ++p ) {
    int pdg = nu_daughter_pdg->at( p );
    if ( pdg == PROTON )
    {
      float px = nu_daughter_px->at( p );
      float py = nu_daughter_py->at( p );
      float pz = nu_daughter_pz->at( p );
      TVector3 temp_p3 = TVector3( px, py, pz );

      mc_p3_p_vec.push_back( temp_p3 );

      float mom = temp_p3.Mag();
      if ( mom > max_mom ) {
        max_mom = mom;
	mc_p3p = temp_p3;
      }
    }
  }

  // TODO: reduce code duplication by just getting the leading proton
  // 3-momentum from this sorted vector
  // Sort the true proton 3-momenta in order from highest to lowest magnitude
  std::sort( mc_p3_p_vec.begin(), mc_p3_p_vec.end(), [](const TVector3& a,
    const TVector3& b) -> bool { return a.Mag() > b.Mag(); } );

  // Check for a true leading proton here and warn the user if one was
  // not found for a signal event
  bool true_lead_p = ( max_mom != LOW_FLOAT );
  if ( !true_lead_p && is_signal ) {
    std::cout << "WARNING: Missing leading proton in MC signal event!\n";
  }

  // Compute true STVs if the event contains both a muon and a leading proton
  double mc_delta_pT = BOGUS, mc_delta_phiT = BOGUS, mc_delta_alphaT = BOGUS,
    mc_delta_pL = BOGUS, mc_pn = BOGUS, mc_delta_pTx = BOGUS,
    mc_delta_pTy = BOGUS, mc_theta_mu_p = BOGUS;

  if ( true_muon && true_lead_p ) {

    KICalculator ki_calc( mc_p3mu, mc_p3p, ki_calc_type_ );

    mc_delta_pT = ki_calc.pT();
    mc_delta_phiT = ki_calc.delta_phiT();
    mc_delta_alphaT = ki_calc.delta_alphaT();
    mc_delta_pL = ki_calc.pL();
    mc_pn = ki_calc.pn();
    mc_delta_pTx = ki_calc.pTx();
    mc_delta_pTy = ki_calc.pTy();
    mc_theta_mu_p = ki_calc.theta_lep_had();
  }

  // We're done. Store MC truth information for this event and return the
  // string indicating the category to which it belongs.
  out[ "mc_is_numu" ] = sig_is_numu;
  out[ "mc_vertex_in_FV" ] = sig_in_fv;
  out[ "mc_lead_p_in_range" ] = sig_lead_p_mom_in_range;
  out[ "mc_muon_in_mom_range" ] = sig_mu_in_mom_range;
  out[ "mc_no_FS_mesons" ] = sig_no_fs_mesons;
  out[ "mc_no_charged_pions_above_threshold" ]
    = sig_mc_no_charged_pi_above_threshold;
  out[ "mc_no_pi0s" ] = sig_mc_no_fs_pi0;
  out[ "num_p_in_mom_range" ] = num_p_in_mom_range;

  out[ "true_delta_pT" ] = mc_delta_pT;
  out[ "true_delta_phiT" ] = mc_delta_phiT;
  out[ "true_delta_alphaT" ] = mc_delta_alphaT;
  out[ "true_delta_pL" ] = mc_delta_pL;
  out[ "true_pn" ] = mc_pn;
  out[ "true_delta_pTx" ] = mc_delta_pTx;
  out[ "true_delta_pTy" ] = mc_delta_pTy;
  out[ "true_theta_mu_p" ] = mc_theta_mu_p;

  out[ "true_p3_mu" ] = mc_p3mu;
  out[ "true_p3_lead_p" ] = mc_p3p;
  out[ "true_p3_p_vec" ] = mc_p3_p_vec;

  if ( !sig_in_fv ) return "Out FV";
  else if ( is_nc ) return "NC";
  else if ( mc_nu_pdg == ELECTRON_NEUTRINO ) return "#nu_{e} CC";
  else if ( mc_nu_pdg != MUON_NEUTRINO ) return "Other";
  else if ( is_signal ) {
    switch ( mc_interaction_type ) {
      case QE_INTERACTION: return "Signal (CCQE)"; break;
      case MEC_INTERACTION: return "Signal (CC2p2h)"; break;
      case RES_INTERACTION: return "Signal (CCRES)"; break;
      default: return "Signal (CC other)";
    }
  }
  else if ( !sig_mc_no_fs_pi0 || !sig_mc_no_charged_pi_above_threshold ) {
    return "#nu_{#mu} CCN#pi";
  }
  else return "Other #nu_{#mu} CC";
}
