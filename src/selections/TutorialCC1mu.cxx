// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/EventCategoriesXp.hh"
#include "XSecAnalyzer/Selections/TutorialCC1mu.hh"

TutorialCC1mu::TutorialCC1mu() : SelectionBase( "TutorialCC1mu" ) {
}

void TutorialCC1mu::DefineConstants() {
  // Define reco & true fiducial volumes, alongside any other constants used
  // within selection cuts

  // FV definitions as in Gardiner's CC0piNp analysis
  // (see https://arxiv.org/abs/2403.19574 and the CC1muNp0pi class)
  // x_min, x_max, y_min, y_max, z_min, z_max
  this->DefineTrueFV( 21.5, 234.85, -95.0, 95.0, 21.5, 966.8 );
  this->DefineRecoFV( 21.5, 234.85, -95.0, 95.0, 21.5, 966.8 );
}

void TutorialCC1mu::ComputeRecoObservables( AnalysisEvent* Event ) {
  // Evaluate the reconstructed kinematic variables of interest for the xsec
  // measurement

  // Set the reco 3-momentum of the muon candidate if we found one
  // during the call to Selection()
  bool has_muon_candidate = muon_candidate_idx_ != BOGUS_INDEX;
  if ( !has_muon_candidate ) return;

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

  *reco_p3mu_ = TVector3( mu_dirx, mu_diry, mu_dirz );
  *reco_p3mu_ = reco_p3mu_->Unit() * muon_mom;
}

void TutorialCC1mu::ComputeTrueObservables( AnalysisEvent* Event ) {
  // Evaluate the true kinematic variables of interest

  // Check if there is a true final-state muon in this event
  bool has_true_muon = ( sig_isNuMu_ && Event->mc_nu_ccnc_ == CHARGED_CURRENT );

  // If there isn't one, we don't need to do anything else
  if ( !has_true_muon ) return;

  // Loop over the true final-state particles of the input event
  size_t num_fs_particles = Event->mc_nu_daughter_pdg_->size();

  // Find the final-state muon with the highest momentum, and store its true
  // 3-momentum in the owned TVector3 object.
  mc_p3mu_->SetXYZ( 0., 0., 0. );

  bool found_muon = false;
  for ( size_t f = 0u; f < num_fs_particles; ++f ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( f );
    if ( pdg == MUON ) {

      found_muon = true;
      float px = Event->mc_nu_daughter_px_->at( f );
      float py = Event->mc_nu_daughter_py_->at( f );
      float pz = Event->mc_nu_daughter_pz_->at( f );

      // Replace the stored 3-momentum with the one for the current final-state
      // muon if the latter is larger
      TVector3 temp_p3mu( px, py, pz );
      if ( temp_p3mu.Mag() > mc_p3mu_->Mag() ) {
        *mc_p3mu_ = temp_p3mu;
      }

    } // final-state muon
  } // loop over final-state particles

  if ( !found_muon ) {
    std::cout << "WARNING: Missing muon in MC signal event!\n";
  }

}

int TutorialCC1mu::CategorizeEvent( AnalysisEvent* Event ) {
  // Identify the event category of the selected event

  // Real data has a bogus true neutrino PDG code that is not one of the
  // allowed values (±12, ±14, ±16)
  int abs_mc_nu_pdg = std::abs( Event->mc_nu_pdg_ );
  Event->is_mc_ = ( abs_mc_nu_pdg == ELECTRON_NEUTRINO ||
    abs_mc_nu_pdg == MUON_NEUTRINO || abs_mc_nu_pdg == TAU_NEUTRINO );
  if ( !Event->is_mc_ ) {
    return kUnknown;
  }

  // All events outside of the true fiducial volume should be categorized
  // as "out of fiducial volume"
  bool mcVertexInFV = point_inside_FV( this->ReturnTrueFV(),
    Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_ );
  if ( !mcVertexInFV ) {
    return kOOFV;
  }

  bool isNC = ( Event->mc_nu_ccnc_ == NEUTRAL_CURRENT );
  if ( isNC ) return kNC;

  if ( Event->mc_nu_pdg_ == ELECTRON_NEUTRINO ) {
    return kNuECC;
  }
  if ( !(Event->mc_nu_pdg_ == MUON_NEUTRINO) ) {
    return kOther;
  }

  if ( this->IsEventMCSignal() ) {
    // Categorize all numuCC events as "CC other" since we don't look at the
    // hadronic content in this selection
    // TODO: revisit this!
    return kNuMuCCOther;
  }

  // We shouldn't ever get here, but return "unknown" just in case
  return kUnknown;
}

bool TutorialCC1mu::DefineSignal( AnalysisEvent* Event ) {
  // Determine whether or not the input event matches the signal definition
  // for this selection. This determination should be done only with MC truth
  // information.

  // Require signal events to be inside the true fiducial volume
  sig_inFV_ = point_inside_FV( this->ReturnTrueFV(), Event->mc_nu_vx_,
    Event->mc_nu_vy_, Event->mc_nu_vz_ );

  // Require an incident muon neutrino
  sig_isNuMu_ = ( Event->mc_nu_pdg_ == MUON_NEUTRINO );

  // Require a final-state muon
  for ( size_t p = 0u; p < Event->mc_nu_daughter_pdg_->size(); ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    if ( pdg == MUON ) {
      sig_has_fs_muon_ = true;
    }
  }

  sig_is_signal_ = sig_inFV_ && sig_isNuMu_ && sig_has_fs_muon_;
  return sig_is_signal_;
}

bool TutorialCC1mu::Selection(AnalysisEvent* Event) {
  // Apply the selection cuts on reco variables only to determine whether this
  // event is a candidate signal event or not.

  // "Proton containment volume" from https://arxiv.org/abs/2403.19574. We keep
  // this definition so that this toy selection matches the CC inclusive
  // portion of the selection used in Gardiner's CC0piNp analysis.
  FiducialVolume PCV;
  PCV.X_Min = 10.;
  PCV.X_Max = 246.35;
  PCV.Y_Min = -106.5;
  PCV.Y_Max = 106.5;
  PCV.Z_Min = 10.;
  PCV.Z_Max = 1026.8;

  // Require selected events to have a reco vertex within the reco fiducial
  // volume
  sel_reco_vertex_in_FV_ = point_inside_FV( this->ReturnRecoFV(),
    Event->nu_vx_, Event->nu_vy_, Event->nu_vz_ );

  // Require the event to have a topological score that passes the cut
  sel_topo_cut_passed_ = Event->topological_score_ > TOPO_SCORE_CUT;

  // Apply the containment cut to the starting positions of all
  // reconstructed tracks and showers. Pass this cut by default.
  sel_pfp_starts_in_PCV_ = true;

  // Loop over each PFParticle in the event
  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    // Use the track reconstruction results to get the start point for every
    // PFParticle for the purpose of verifying containment. We could in
    // principle differentiate between tracks and showers here, but track
    // information was used in the CC0piNp analysis (which rejected all showers
    // downstream of this cut anyway). For consistency, we therefore apply the
    // track reconstruction here unconditionally.
    float x = Event->track_startx_->at( p );
    float y = Event->track_starty_->at( p );
    float z = Event->track_startz_->at( p );

    // Verify that the start of the PFParticle lies within the containment
    // volume. See https://stackoverflow.com/a/2488507 for an explanation of
    // the use of &= here. Don't worry, it's type-safe since both operands are
    // bool.
    sel_pfp_starts_in_PCV_ &= point_inside_FV( PCV, x, y, z );
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

  // Check whether the muon candidate is contained. Use the PCV as the
  // containment volume.
  sel_muon_contained_ = false;
  if ( muon_candidate_idx_ != BOGUS_INDEX ) {
    float endx = Event->track_endx_->at( muon_candidate_idx_ );
    float endy = Event->track_endy_->at( muon_candidate_idx_ );
    float endz = Event->track_endz_->at( muon_candidate_idx_ );
    bool end_contained = point_inside_FV( PCV, endx, endy, endz );
    if ( end_contained ) sel_muon_contained_ = true;
  }

  sel_nu_mu_cc_ = sel_reco_vertex_in_FV_ && sel_pfp_starts_in_PCV_
    && sel_has_muon_candidate_ && sel_topo_cut_passed_;

  return sel_nu_mu_cc_;
}

void TutorialCC1mu::DefineOutputBranches() {
  // Save any additional variables to the output TTree
  this->SetBranch( &sig_isNuMu_, "mc_is_numu", kBool );
  this->SetBranch( &sig_inFV_, "mc_vertex_in_FV", kBool );
  this->SetBranch( &sig_has_fs_muon_, "mc_has_fs_muon", kBool );
  this->SetBranch( &sig_is_signal_, "mc_is_signal", kBool );
  this->SetBranch( &sel_reco_vertex_in_FV_, "sel_reco_vertex_in_FV", kBool);
  this->SetBranch( &sel_pfp_starts_in_PCV_, "sel_pfp_starts_in_PCV", kBool);
  this->SetBranch( &sel_has_muon_candidate_, "sel_has_muon_candidate", kBool);
  this->SetBranch( &sel_topo_cut_passed_, "sel_topo_cut_passed", kBool);
  this->SetBranch( &sel_nu_mu_cc_, "sel_nu_mu_cc", kBool);
  this->SetBranch( &sel_muon_contained_, "sel_muon_contained", kBool);
  this->SetBranch( &muon_candidate_idx_, "muon_candidate_idx", kInteger);
  this->SetBranch( reco_p3mu_, "reco_p3_mu", kTVector );
  this->SetBranch( mc_p3mu_, "true_p3_mu", kTVector );
}

void TutorialCC1mu::DefineCategoryMap() {
  // Use the shared category map for 1p/2p/Np/Xp
  categ_map_ = CC1muXp_MAP;
}
