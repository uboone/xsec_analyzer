// XSecAnalyzer includes
#include "XSecAnalyzer/TreeUtils.hh"
#include "XSecAnalyzer/Functions.hh"

#include "XSecAnalyzer/Selections/CC1muNp0pi.hh"
#include "XSecAnalyzer/Selections/EventCategoriesXp.hh"

CC1muNp0pi::CC1muNp0pi() : SelectionBase("CC1muNp0pi") {
  CalcType = kOpt1;
}

void CC1muNp0pi::DefineConstants() {
  DefineTrueFV(21.5,234.85,-95.0,95.0,21.5,966.8);
  DefineRecoFV(21.5,234.85,-95.0,95.0,21.5,966.8);
}

void CC1muNp0pi::ComputeTrueObservables(AnalysisEvent* Event) {
  size_t num_mc_daughters = Event->mc_nu_daughter_pdg_->size();

  // Set the true 3-momentum of the final-state muon if there is one
  bool true_muon = ( sig_isNuMu_ && Event->mc_nu_ccnc_ == CHARGED_CURRENT );
  if ( true_muon ) {
    // Loop over the MC neutrino daughters, find the muon, and get its
    // true 3-momentum. Note that we assume there is only one muon in
    // this loop.
    bool found_muon = false;
    for ( size_t d = 0u; d < num_mc_daughters; ++d ) {
      int pdg = Event->mc_nu_daughter_pdg_->at( d );
      if ( pdg == MUON ) {
        found_muon = true;
        float px = Event->mc_nu_daughter_px_->at( d );
        float py = Event->mc_nu_daughter_py_->at( d );
        float pz = Event->mc_nu_daughter_pz_->at( d );
	*mc_p3mu = TVector3( px, py, pz );
        break;
      }
    }

    if ( !found_muon ) {
      std::cout << "WARNING: Missing muon in MC signal event!\n";
      return;
    }
  }

  // Set the true 3-momentum of the leading proton (if there is one)
  float max_mom = LOW_FLOAT;
  for ( int p = 0; p < num_mc_daughters; ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    if ( pdg == PROTON )
    {
      float px = Event->mc_nu_daughter_px_->at( p );
      float py = Event->mc_nu_daughter_py_->at( p );
      float pz = Event->mc_nu_daughter_pz_->at( p );
      TVector3 temp_p3 = TVector3( px, py, pz );

      mc_p3_p_vec_->push_back( temp_p3 );

      float mom = temp_p3.Mag();
      if ( mom > max_mom ) {
        max_mom = mom;
	*mc_p3p = temp_p3;
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
  if ( !true_lead_p && IsEventMCSignal() ) {
    // If it doesn't for a signal event, then something is wrong.
    std::cout << "WARNING: Missing leading proton in MC signal event!\n";
    return;
  }

  // Compute true STVs if the event contains both a muon and a leading
  // proton
  if ( true_muon && true_lead_p ) {
    double MuonEnergy = real_sqrt(mc_p3mu->Mag()*mc_p3mu->Mag() + MUON_MASS*MUON_MASS);
    double ProtonEnergy = real_sqrt(mc_p3p->Mag()*mc_p3p->Mag() + PROTON_MASS*PROTON_MASS);
    STVTools.CalculateSTVs(*(mc_p3mu.get()), *(mc_p3p.get()), MuonEnergy, ProtonEnergy, CalcType);

    mc_delta_pT_ = STVTools.ReturnPt();
    mc_delta_phiT_ = STVTools.ReturnDeltaPhiT() * TMath::Pi()/180.;
    mc_delta_alphaT_ = STVTools.ReturnDeltaAlphaT() * TMath::Pi()/180.;
    mc_delta_pL_ = STVTools.ReturnPL();
    mc_pn_ = STVTools.ReturnPn();
    mc_delta_pTx_ = STVTools.ReturnPtx();
    mc_delta_pTy_ = STVTools.ReturnPty();

    mc_theta_mu_p_ = std::acos( mc_p3mu->Dot(*mc_p3p)
      / mc_p3mu->Mag() / mc_p3p->Mag() );
  }
}

void CC1muNp0pi::ComputeRecoObservables(AnalysisEvent* Event) {

  // In cases where we failed to find a muon candidate, check whether there are
  // at least two generation == 2 PFParticles. If there are, then compute the
  // usual observables using the longest track as the muon candidate and the
  // second-longest track as the leading proton candidate. This will enable
  // sideband studies of NC backgrounds in the STV phase space.
  if ( !sel_has_muon_candidate_ ) {
    float max_trk_len = LOW_FLOAT;
    int max_trk_idx = BOGUS_INDEX;

    float next_to_max_trk_len = LOW_FLOAT;
    int next_to_max_trk_idx = BOGUS_INDEX;

    for ( int p = 0; p < Event->num_pf_particles_; ++p ) {

      // Only include direct neutrino daughters (generation == 2)
      unsigned int generation = Event->pfp_generation_->at( p );
      if ( generation != 2u ) continue;

      float trk_len = Event->track_length_->at( p );

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
      muon_candidate_idx_ = max_trk_idx;
      lead_p_candidate_idx_ = next_to_max_trk_idx;
    }
  }

  // Set the reco 3-momentum of the muon candidate if we found one
  bool muon = muon_candidate_idx_ != BOGUS_INDEX;
  if ( muon ) {
    float mu_dirx = Event->track_dirx_->at( muon_candidate_idx_ );
    float mu_diry = Event->track_diry_->at( muon_candidate_idx_ );
    float mu_dirz = Event->track_dirz_->at( muon_candidate_idx_ );

    // The selection flag indicating whether the muon candidate is contained
    // was already set when the selection was applied. Use it to choose the
    // best momentum estimator to use.
    float muon_mom = LOW_FLOAT;
    if ( sel_muon_contained_ ) {
      muon_mom = Event->track_range_mom_mu_->at( muon_candidate_idx_ );
    }
    else {
      muon_mom = Event->track_mcs_mom_mu_->at( muon_candidate_idx_ );
    }

    *p3mu = TVector3( mu_dirx, mu_diry, mu_dirz );
    *p3mu = p3mu->Unit() * muon_mom;
  }

    // Set the reco 3-momentum of the leading proton candidate if we found one
  bool lead_p = lead_p_candidate_idx_ != BOGUS_INDEX;
  if ( lead_p ) {

    float p_dirx = Event->track_dirx_->at( lead_p_candidate_idx_ );
    float p_diry = Event->track_diry_->at( lead_p_candidate_idx_ );
    float p_dirz = Event->track_dirz_->at( lead_p_candidate_idx_ );
    float KEp = Event->track_kinetic_energy_p_->at( lead_p_candidate_idx_ );
    float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

    *p3p = TVector3( p_dirx, p_diry, p_dirz );
    *p3p = p3p->Unit() * p_mom;
  }

  // Set the reco 3-momenta of all proton candidates (i.e., all generation == 2
  // tracks except the muon candidate) assuming we found both a muon candidate
  // and at least one proton candidate.
  if ( muon && lead_p ) {
    for ( int p = 0; p < Event->num_pf_particles_; ++p ) {
      // Skip the muon candidate
      if ( p == muon_candidate_idx_ ) continue;

      // Only include direct neutrino daughters (generation == 2)
      unsigned int generation = Event->pfp_generation_->at( p );
      if ( generation != 2u ) continue;

      float p_dirx = Event->track_dirx_->at( p );
      float p_diry = Event->track_diry_->at( p );
      float p_dirz = Event->track_dirz_->at( p );
      float KEp = Event->track_kinetic_energy_p_->at( p );
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
    double MuonEnergy = real_sqrt(p3mu->Mag()*p3mu->Mag() + MUON_MASS*MUON_MASS);
    double ProtonEnergy	= real_sqrt(p3p->Mag()*p3p->Mag() + PROTON_MASS*PROTON_MASS);
    STVTools.CalculateSTVs(*(p3mu), *(p3p), MuonEnergy, ProtonEnergy, CalcType);

    delta_pT_ = STVTools.ReturnPt();
    delta_phiT_ = STVTools.ReturnDeltaPhiT() * TMath::Pi()/180.;
    delta_alphaT_ = STVTools.ReturnDeltaAlphaT() * TMath::Pi()/180.;
    delta_pL_ = STVTools.ReturnPL();
    pn_ = STVTools.ReturnPn();
    delta_pTx_ = STVTools.ReturnPtx();
    delta_pTy_ = STVTools.ReturnPty();

    theta_mu_p_ = std::acos( p3mu->Dot(*p3p) / p3mu->Mag() / p3p->Mag() );
  }

}

bool CC1muNp0pi::DefineSignal(AnalysisEvent* Event) {
  sig_inFV_ = point_inside_FV(ReturnTrueFV(), Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_);
  sig_isNuMu_ = (Event->mc_nu_pdg_ == MUON_NEUTRINO);
  bool IsNC = (Event->mc_nu_ccnc_ == NEUTRAL_CURRENT);

  sig_noFSMesons_= true;
  sig_mc_no_fs_pi0_ = true;
  sig_mc_no_charged_pi_above_threshold_ = true;

  sig_muonInMomRange_ = false;

  sig_nProtons_in_Momentum_range = 0;

  double LeadProtonMomentum = 0.;

  for ( size_t p = 0u; p < Event->mc_nu_daughter_pdg_->size(); ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    float energy = Event->mc_nu_daughter_energy_->at( p );

    // Do the general check for (anti)mesons first before considering
    // any individual PDG codes
    if ( is_meson_or_antimeson(pdg) ) {
      sig_noFSMesons_ = false;
    }

    if ( pdg == MUON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(MUON_MASS, 2) );
      if ( mom >= MUON_P_MIN_MOM_CUT && mom <= MUON_P_MAX_MOM_CUT ) {
	sig_muonInMomRange_= true;
      }
    }
    else if ( pdg == PROTON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PROTON_MASS, 2) );
      if ( mom > LeadProtonMomentum ) LeadProtonMomentum = mom;
      if ( mom >= LEAD_P_MIN_MOM_CUT && mom <= LEAD_P_MAX_MOM_CUT ) sig_nProtons_in_Momentum_range++;
    }
    //Not used for selection purposes
    else if ( pdg == PI_ZERO ) {
      sig_mc_no_fs_pi0_ = false;
    }
    //Not used for selection purposes
    else if ( std::abs(pdg) == PI_PLUS ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PI_PLUS_MASS, 2) );
      if ( mom > CHARGED_PI_MOM_CUT ) {
        sig_mc_no_charged_pi_above_threshold_ = false;
      }
    }
  }

  sig_leadProtonMomInRange_ = false;
  // Check that the leading proton has a momentum within the allowed range
  if ( LeadProtonMomentum >= LEAD_P_MIN_MOM_CUT && LeadProtonMomentum <= LEAD_P_MAX_MOM_CUT ) {
    sig_leadProtonMomInRange_ = true;
  }

  bool ReturnVal = sig_inFV_ && !IsNC && sig_isNuMu_ && sig_muonInMomRange_ && sig_leadProtonMomInRange_ && sig_noFSMesons_;
  return ReturnVal;
}

bool CC1muNp0pi::Selection(AnalysisEvent* Event) {
  FiducialVolume PCV;
  PCV.X_Min = 10.;
  PCV.X_Max = 246.35;
  PCV.Y_Min = -106.5;
  PCV.Y_Max = 106.5;
  PCV.Z_Min = 10.;
  PCV.Z_Max = 1026.8;

  sel_reco_vertex_in_FV_ = point_inside_FV(ReturnRecoFV(), Event->nu_vx_, Event->nu_vy_, Event->nu_vz_);
  //std::cout << ReturnRecoFV().X_Min << " " << ReturnRecoFV().X_Max << " " << Event->nu_vx_ << " " << sel_reco_vertex_in_FV_ << std::endl;

  sel_topo_cut_passed_ = Event->topological_score_ > TOPO_SCORE_CUT;
  sel_cosmic_ip_cut_passed_ = Event->cosmic_impact_parameter_ > COSMIC_IP_CUT;

  // Apply the containment cut to the starting positions of all
  // reconstructed tracks and showers. Pass this cut by default.
  sel_pfp_starts_in_PCV_ = true;

  // Loop over each PFParticle in the event
  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    // Use the track reconstruction results to get the start point for
    // every PFParticle for the purpose of verifying containment. We could
    // in principle differentiate between tracks and showers here, but
    // (1) we cut out all showers later on in the selection anyway, and
    // (2) the blinded PeLEE data ntuples do not include shower information.
    // We therefore apply the track reconstruction here unconditionally.
    float x = Event->track_startx_->at( p );
    float y = Event->track_starty_->at( p );
    float z = Event->track_startz_->at( p );

    // Verify that the start of the PFParticle lies within the containment
    // volume.
    // TODO: revisit which containment volume to use for PFParticle start
    // positions. See https://stackoverflow.com/a/2488507 for an explanation
    // of the use of &= here. Don't worry, it's type-safe since both operands
    // are bool.
    sel_pfp_starts_in_PCV_ &= point_inside_FV(PCV, x, y, z );
  }

  // Sets the sel_has_muon_candidate_ flag as appropriate. The threshold check
  // is handled later.
  std::vector<int> muon_candidate_indices;
  std::vector<int> muon_pid_scores;

  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {
    // Only direct neutrino daughters (generation == 2) will be considered as
    // possible muon candidates
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    float track_score = Event->pfp_track_score_->at( p );
    float start_dist = Event->track_start_distance_->at( p );
    float track_length = Event->track_length_->at( p );
    float pid_score = Event->track_llr_pid_score_->at( p );

    if ( track_score > MUON_TRACK_SCORE_CUT
      && start_dist < MUON_VTX_DISTANCE_CUT
      && track_length > MUON_LENGTH_CUT
      && pid_score > MUON_PID_CUT )
    {
      muon_candidate_indices.push_back( p );
      muon_pid_scores.push_back( pid_score );
    }
  }

  size_t num_candidates = muon_candidate_indices.size();
  if ( num_candidates > 0u ) sel_has_muon_candidate_ = true;

  if ( num_candidates == 1u ) {
    muon_candidate_idx_ = muon_candidate_indices.front();
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
    muon_candidate_idx_ = chosen_index;
  }
  else {
    muon_candidate_idx_ = BOGUS_INDEX;
  }

  sel_nu_mu_cc_ = sel_reco_vertex_in_FV_ && sel_pfp_starts_in_PCV_
    && sel_has_muon_candidate_ && sel_topo_cut_passed_;

  // Fail the shower cut if any showers were reconstructed
  // NOTE: We could do this quicker like this,
  //   sel_no_reco_showers_ = ( num_showers_ > 0 );
  // but it might be nice to be able to adjust the track score for this cut.
  // Thus, we do it the hard way.
  int reco_shower_count = 0;
  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    float tscore = Event->pfp_track_score_->at( p );
    if ( tscore <= TRACK_SCORE_CUT ) ++reco_shower_count;
  }
  // Check the shower cut
  sel_no_reco_showers_ = ( reco_shower_count == 0 );

  // Set flags that default to true here
  sel_passed_proton_pid_cut_ = true;
  sel_protons_contained_ = true;

  // Set flags that default to false here
  sel_muon_contained_ = false;

  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {
    // Only worry about direct neutrino daughters (PFParticles considered
    // daughters of the reconstructed neutrino)
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    // Check that we can find a muon candidate in the event. If more than
    // one is found, also fail the cut.
    if ( p == muon_candidate_idx_ ) {

      // Check whether the muon candidate is contained. Use the same
      // containment volume as the protons. TODO: revisit this as needed.
      float endx = Event->track_endx_->at( p );
      float endy = Event->track_endy_->at( p );
      float endz = Event->track_endz_->at( p );
      bool end_contained = point_inside_FV(PCV, endx, endy, endz );

      if ( end_contained ) sel_muon_contained_ = true;

      // Check that the muon candidate is above threshold. Use the best
      // momentum based on whether it was contained or not.

      float muon_mom = LOW_FLOAT;
      float range_muon_mom = Event->track_range_mom_mu_->at( p );
      float mcs_muon_mom = Event->track_mcs_mom_mu_->at( p );

      if ( sel_muon_contained_ ) muon_mom = range_muon_mom;
      else muon_mom = mcs_muon_mom;

      if ( muon_mom >= MUON_P_MIN_MOM_CUT && muon_mom <= MUON_P_MAX_MOM_CUT ) {
        sel_muon_passed_mom_cuts_ = true;
      }

      // Apply muon candidate quality cut by comparing MCS and range-based
      // momentum estimators. Default to failing the cut.
      sel_muon_quality_ok_ = false;

      double frac_diff_range_mcs = std::abs( range_muon_mom - mcs_muon_mom );
      if ( range_muon_mom > 0. ) {
        frac_diff_range_mcs /= range_muon_mom;
        if ( frac_diff_range_mcs < MUON_MOM_QUALITY_CUT ) {
          sel_muon_quality_ok_ = true;
        }
      }
    }
    else {

      float track_score = Event->pfp_track_score_->at( p );
      if ( track_score <= TRACK_SCORE_CUT ) continue;

      // Bad tracks in the searchingfornues TTree can have
      // bogus track lengths. This skips those.
      float track_length = Event->track_length_->at( p );
      if ( track_length <= 0. ) continue;

      // We found a reco track that is not the muon candidate. All such
      // tracks are considered proton candidates.
      sel_has_p_candidate_ = true;

      float llr_pid_score = Event->track_llr_pid_score_->at( p );

      // Check whether the current proton candidate fails the proton PID cut
      if ( llr_pid_score > proton_pid_cut(track_length) ) {
        sel_passed_proton_pid_cut_ = false;
      }

      // Check whether the current proton candidate fails the containment cut
      float endx = Event->track_endx_->at( p );
      float endy = Event->track_endy_->at( p );
      float endz = Event->track_endz_->at( p );
      bool end_contained = point_inside_FV(PCV, endx, endy, endz );
      if ( !end_contained ) sel_protons_contained_ = false;
    }

  }

  // Don't bother to apply the cuts that involve the leading
  // proton candidate if we don't have one
  if ( !sel_has_p_candidate_ ) {
    return false;
  }

  // All that remains is to apply the leading proton candidate cuts. We could
  // search for it above, but doing it here makes the code more readable (with
  // likely negligible impact on performance)
  float lead_p_track_length = LOW_FLOAT;
  size_t lead_p_index = 0u;
  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    // Skip the muon candidate reco track (this function assumes that it has
    // already been found)
    if ( p == muon_candidate_idx_ ) continue;

    // Skip PFParticles that are shower-like (track scores near 0)
    float track_score = Event->pfp_track_score_->at( p );
    if ( track_score <= TRACK_SCORE_CUT ) continue;

    // All non-muon-candidate reco tracks are considered proton candidates
    float track_length = Event->track_length_->at( p );
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

  // Check the range-based reco momentum for the leading proton candidate
  float lead_p_KE = Event->track_kinetic_energy_p_->at( lead_p_candidate_idx_ );
  float range_mom_lead_p = real_sqrt( lead_p_KE*lead_p_KE
    + 2.*PROTON_MASS*lead_p_KE );
  if ( range_mom_lead_p >= LEAD_P_MIN_MOM_CUT
    && range_mom_lead_p <= LEAD_P_MAX_MOM_CUT )
  {
    sel_lead_p_passed_mom_cuts_ = true;
  }

  // All right, we've applied all selection cuts. Set the flag that indicates
  // whether all were passed (and thus the event is selected as a CCNp0pi
  // candidate)
  bool sel_CCNp0pi_ = sel_nu_mu_cc_ && sel_no_reco_showers_
    && sel_muon_passed_mom_cuts_ && sel_muon_contained_ && sel_muon_quality_ok_
    && sel_has_p_candidate_ && sel_passed_proton_pid_cut_
    && sel_protons_contained_ && sel_lead_p_passed_mom_cuts_;

  return sel_CCNp0pi_;
}

int CC1muNp0pi::CategorizeEvent(AnalysisEvent* Event) {
  // Real data has a bogus true neutrino PDG code that is not one of the
  // allowed values (±12, ±14, ±16)
  int abs_mc_nu_pdg = std::abs( Event->mc_nu_pdg_ );
  Event->is_mc_ = ( abs_mc_nu_pdg == ELECTRON_NEUTRINO || abs_mc_nu_pdg == MUON_NEUTRINO || abs_mc_nu_pdg == TAU_NEUTRINO );
  if ( !Event->is_mc_ ) {
    return kUnknown;
  }

  bool MCVertexInFV = point_inside_FV(ReturnTrueFV(), Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_);
  if ( !MCVertexInFV ) {
    return kOOFV;
  }

  bool isNC = (Event->mc_nu_ccnc_ == NEUTRAL_CURRENT);
  //DB Currently only one NC category is supported so test first. Will likely want to change this in the future
  if (isNC) return kNC;

  if (Event->mc_nu_pdg_ == ELECTRON_NEUTRINO) {
    return kNuECC;
  }
  if (!(Event->mc_nu_pdg_ == MUON_NEUTRINO)) {
    return kOther;
  }

  if ( IsEventMCSignal() ) {
    if (sig_nProtons_in_Momentum_range == 1) {
      if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCC1p0pi_CCQE; // QE
      else if ( Event->mc_nu_interaction_type_ == 10 ) return kNuMuCC1p0pi_CCMEC; // MEC
      else if ( Event->mc_nu_interaction_type_ == 1 ) return kNuMuCC1p0pi_CCRES; // RES
      else return kNuMuCCMp0pi_Other;
    } else if (sig_nProtons_in_Momentum_range == 2) {
      if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCC2p0pi_CCQE; // QE
      else if ( Event->mc_nu_interaction_type_ == 10 ) return kNuMuCC2p0pi_CCMEC; // MEC
      else if ( Event->mc_nu_interaction_type_ == 1 ) return kNuMuCC2p0pi_CCRES; // RES
      else return kNuMuCCMp0pi_Other;
    } else { // i.e. >=3
      if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCCMp0pi_CCQE; // QE
      else if ( Event->mc_nu_interaction_type_ == 10 ) return kNuMuCCMp0pi_CCMEC; // MEC
      else if ( Event->mc_nu_interaction_type_ == 1 ) return kNuMuCCMp0pi_CCRES; // RES
      else return kNuMuCCMp0pi_Other;
    }
  }
  else if (!sig_mc_no_fs_pi0_ || !sig_mc_no_charged_pi_above_threshold_) {
    return kNuMuCCNpi;
  } else if (!sig_leadProtonMomInRange_) {
    if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCC0p0pi_CCQE; // QE
    else if ( Event->mc_nu_interaction_type_ == 10 ) return kNuMuCC0p0pi_CCMEC; // MEC
    else if ( Event->mc_nu_interaction_type_ == 1 ) return kNuMuCC0p0pi_CCRES; // RES
    else return kNuMuCC0p0pi_Other;
  }
  return kNuMuCCOther;
}

void CC1muNp0pi::DefineOutputBranches() {
  SetBranch(&sig_isNuMu_,"mc_is_numu",kBool);
  SetBranch(&sig_inFV_,"mc_vertex_in_FV",kBool);
  SetBranch(&sig_leadProtonMomInRange_,"mc_lead_p_in_range",kBool);
  SetBranch(&sig_muonInMomRange_,"mc_muon_in_mom_range",kBool);
  SetBranch(&sig_noFSMesons_,"mc_no_FS_mesons",kBool);
  SetBranch(&sig_mc_no_charged_pi_above_threshold_,"mc_no_charged_pions_above_thres",kBool);
  SetBranch(&sig_mc_no_fs_pi0_,"mc_no_pi0s",kBool);
  SetBranch(&sig_nProtons_in_Momentum_range,"nProtons_in_Momentum_range",kInteger);

  SetBranch(&sel_reco_vertex_in_FV_,"reco_vertex_in_FV",kBool);
  SetBranch(&sel_pfp_starts_in_PCV_,"pfp_starts_in_PCV",kBool);
  SetBranch(&sel_has_muon_candidate_,"has_muon_candidate",kBool);
  SetBranch(&sel_topo_cut_passed_,"topo_cut_passed",kBool);
  SetBranch(&sel_nu_mu_cc_,"nu_mu_cc",kBool);
  SetBranch(&sel_muon_contained_,"muon_contained",kBool);
  SetBranch(&sel_muon_passed_mom_cuts_,"muon_passed_mom_cuts",kBool);
  SetBranch(&sel_no_reco_showers_,"no_reco_showers",kBool);
  SetBranch(&sel_has_p_candidate_,"has_p_candidate",kBool);
  SetBranch(&sel_muon_quality_ok_,"muon_quality_ok",kBool);
  SetBranch(&sel_protons_contained_,"protons_contained",kBool);
  SetBranch(&sel_passed_proton_pid_cut_,"passed_proton_pid_cut",kBool);
  SetBranch(&sel_lead_p_passed_mom_cuts_,"lead_p_passed_mom_cuts",kBool);
  SetBranch(&sel_cosmic_ip_cut_passed_,"cosmic_ip_cut_passed",kBool);

  SetBranch(&lead_p_candidate_idx_,"lead_p_candidate_idx",kInteger);
  SetBranch(&muon_candidate_idx_,"muon_candidate_idx",kInteger);

  SetBranch(&delta_pT_,"reco_delta_pT",kDouble);
  SetBranch(&delta_phiT_,"reco_delta_phiT",kDouble);
  SetBranch(&delta_alphaT_,"reco_delta_alphaT",kDouble);
  SetBranch(&delta_pL_,"reco_delta_pL",kDouble);
  SetBranch(&pn_,"reco_pn",kDouble);
  SetBranch(&delta_pTx_,"reco_delta_pTx",kDouble);
  SetBranch(&delta_pTy_,"reco_delta_pTy",kDouble);
  SetBranch(&theta_mu_p_,"reco_theta_mu_p",kDouble);
  SetBranch(p3mu,"reco_p3_mu",kTVector);
  SetBranch(p3p,"reco_p3_lead_p",kTVector);
  SetBranch(p3_p_vec_,"reco_p3_p_vec",kSTDVector);

  SetBranch(&mc_delta_pT_,"true_delta_pT",kDouble);
  SetBranch(&mc_delta_phiT_,"true_delta_phiT",kDouble);
  SetBranch(&mc_delta_alphaT_,"true_delta_alphaT",kDouble);
  SetBranch(&mc_delta_pL_,"true_delta_pL",kDouble);
  SetBranch(&mc_pn_,"true_pn",kDouble);
  SetBranch(&mc_delta_pTx_,"true_delta_pTx",kDouble);
  SetBranch(&mc_delta_pTy_,"true_delta_pTy",kDouble);
  SetBranch(&mc_theta_mu_p_,"true_theta_mu_p",kDouble);
  SetBranch(mc_p3mu,"true_p3_mu",kTVector);
  SetBranch(mc_p3p,"true_p3_lead_p",kTVector);
  SetBranch(mc_p3_p_vec_,"true_p3_p_vec",kSTDVector);
}

void CC1muNp0pi::DefineCategoryMap() {
  // Use the shared category map for 1p/2p/Np/Xp
  categ_map_ = CC1muXp_MAP;
}
