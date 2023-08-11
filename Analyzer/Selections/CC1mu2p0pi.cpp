#include "CC1mu2p0pi.h"

#include "TreeUtils.hh"

CC1mu2p0pi::CC1mu2p0pi() : SelectionBase("CC1mu2p0pi") {
}

void CC1mu2p0pi::ComputeObservables() {

}

bool CC1mu2p0pi::Selection(AnalysisEvent* Event) {
  float_t xmin;
  float_t xmax;
  float_t ymin;
  float_t ymax;
  float_t zmin;
  float_t zmax;
  
  //==============================================================================================================================
  //DB This has a difference between Stephen's and Samanatha's original analyses
  //sel_reco_vertex_in_FV_ = this->reco_vertex_inside_FV();
  
  //DB Let's fix this and perform our own FV search for comparisons to Samantha's original selection code
  float_t x = Event->nu_vx_;
  float_t y = Event->nu_vy_;
  float_t z = Event->nu_vz_;

  xmin = 10.;
  xmax = 246.35;
  ymin = -106.5;
  ymax = 106.5;
  zmin = 10.;
  zmax = 1026.8;

  if((x <= xmin || x >= xmax) || (y <= ymin || y >= ymax) || (z <= zmin || z >= zmax)){
    sel_reco_vertex_in_FV_ = false;
  } else{
    sel_reco_vertex_in_FV_ = true;
  }

  //==============================================================================================================================
  //DB Samantha's analysis explicitly cuts out events with num_candidates!=1 (n_muons) 
  int n_muons = 0;
  int chosen_index = 0;

  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {
    // Only direct neutrino daughters (generation == 2) will be considered as
    // possible muon candidates
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    float pid_score = Event->track_llr_pid_score_->at( p );
    if ( pid_score >= MUON_PID_CUT && pid_score > -1 && pid_score < 1) {
      n_muons += 1;

      //Gets overwritten if multuple muon candidates, but that's fine because we require exactly one muon candidate
      chosen_index = p;
    }
  }

  if ( n_muons == 1u ) {
    sel_has_muon_candidate_ = true;
    muon_candidate_idx_ = chosen_index;
  } else {
    muon_candidate_idx_ = BOGUS_INDEX;
  }

  //==============================================================================================================================
  //DB Does the event pass the numuCC0pi selection?
  sel_nu_mu_cc_ = sel_reco_vertex_in_FV_ && sel_has_muon_candidate_;

  //==============================================================================================================================
  //DB Require exactly 3 PFP's
  if (Event->num_pf_particles_ == 3) {
    sel_npfps_eq_3 = true;
  }

  //==============================================================================================================================
  //DB Require 3 tracks (track_score > 0.8 [MUON_TRACK_SCORE_CUT]) whose "vertex distance attachment is less than 4 cm"
  int nTracks = 0;

  TVector3 nu_vtx_reco(Event->nu_vx_,Event->nu_vy_,Event->nu_vz_);
  for(int i = 0; i < Event->num_pf_particles_; i ++){
    float track_score = Event->pfp_track_score_->at(i);
    float track_distance = Event->track_start_distance_->at(i);
    float track_pid = Event->track_llr_pid_score_->at(i);
    TVector3 track_end(Event->track_endx_->at(i),Event->track_endy_->at(i),Event->track_endz_->at(i)); //leading track end reco
    track_end -= nu_vtx_reco;
    double track_end_distance = track_end.Mag();

    //DB I think this cut has a different logic than Stephen's
    if (track_score >= MUON_TRACK_SCORE_CUT && (track_distance <= MUON_VTX_DISTANCE_CUT || track_end_distance <= MUON_VTX_DISTANCE_CUT)){
      nTracks++;
    }
  }

  if (nTracks == 3) {
    sel_ntracks_eq_3 = true;
  }

  //==============================================================================================================================
  //DB Make sure there's two Protons and one Muon.
  //Muon selection already performed in apply_numu_CC_selection()
  //      -> sel_has_muon_candidate_ == true if only 1 muon candidate
  int nProtons = 0;
  for(int i = 0; i < Event->num_pf_particles_; i ++){
    float track_pid = Event->track_llr_pid_score_->at(i);
    if(track_pid < DEFAULT_PROTON_PID_CUT && track_pid < 1 && track_pid > -1){
      nProtons += 1;
    }
  }

  if (sel_has_muon_candidate_ && nProtons == 2) {
    sel_correctparticles = true;
  }

  //==============================================================================================================================
  //DB Now check that all 3PFPs are contained (start and end) within the FV
  
  //DB Initially samantha precalculates which index corresponds to the muon, leading p proton and recoil proton
  //But then applies the same containment cut to all three. The function used is:
  //https://github.com/ssfehlberg/CC2p-Event-Selection/blob/9492ff121a2eb884f464e1c166d067f217a04900/PeLEE_ntuples/helper_funcs.h#L18-L24
  //Where the start of the vertex requires an addition 10cm tighter FV cut (i.e. the argument variables), and the end of the vertex
  //uses the FV cut provided (i.e. the argument variables are 0cm). I don't think there's any need to precalculate the indexs...
  //because we've already selected events which only have 3 PFPs within them (thus 1muon and 2protons)...
  //but this could be wrong
  bool StartContained = true;
  bool EndContained = true;

  for(int i = 0; i < Event->num_pf_particles_; i ++){
    TVector3 track_start(Event->track_startx_->at(i),Event->track_starty_->at(i),Event->track_startz_->at(i));

    xmin = 10.;
    xmax = 246.35;
    ymin = -106.5;
    ymax = 106.5;
    zmin = 10.;
    zmax = 1026.8;

    if((track_start.x() <= xmin || track_start.x() >= xmax) || (track_start.y() <= ymin || track_start.y() >= ymax) || (track_start.z() <= zmin || track_start.z() >= zmax)){
      StartContained = false;
      break;
    }

    TVector3 track_end(Event->track_endx_->at(i),Event->track_endy_->at(i),Event->track_endz_->at(i));

    xmin = 0.0;
    xmax = 256.35;
    ymin = -116.5;
    ymax = 116.5;
    zmin = 0.0;
    zmax = 1036.8;

    if((track_end.x() <= xmin || track_end.x() >= xmax) || (track_end.y() <= ymin || track_end.y() >= ymax) || (track_end.z() <= zmin || track_end.z() >= zmax)){
      EndContained = false;
      break;
    }

  }
  if (StartContained && EndContained) {
    sel_containedparticles = true;
  }

  //==============================================================================================================================
  //DB Now ensure that the muon and proton candidates pass the momentum threshold requirements of
  // 0.1 <= MuonMomentum <= 1.2
  // 0.3 <= ProtonMomentum <= 1.0
  
  bool MomentumThresholdPassed = true;
  for(int i = 0; i < Event->num_pf_particles_; i ++){
    if (i == muon_candidate_idx_) {
      if ( Event->track_range_mom_mu_->at(i) < MUON_P_MIN_MOM_CUT || Event->track_range_mom_mu_->at(i) > MUON_P_MAX_MOM_CUT ) {
        MomentumThresholdPassed = false;
      }
    } else {
      float ProtonMomentum = std::sqrt(std::pow(Event->track_kinetic_energy_p_->at(i) + PROTON_MASS,2) - std::pow(PROTON_MASS,2));
      if ( ProtonMomentum < PROTON_MIN_MOM_CUT || ProtonMomentum > PROTON_MAX_MOM_CUT ) {
        MomentumThresholdPassed = false;
      }
    }
  }

  if (MomentumThresholdPassed == true) sel_momentum_threshold_passed_ = true;

  //==============================================================================================================================
  //Does everything pass selection?

  bool Passed = sel_nu_mu_cc_ && sel_npfps_eq_3 && sel_ntracks_eq_3 && sel_correctparticles
    && sel_containedparticles && sel_momentum_threshold_passed_;

  return Passed;
}

void CC1mu2p0pi::DefineBranches() {
  SetBranch(&sel_reco_vertex_in_FV_,"reco_vertex_in_FV",kBool);
  SetBranch(&sel_has_muon_candidate_,"has_muon_candidate",kBool);
  SetBranch(&sel_nu_mu_cc_,"nu_mu_cc",kBool);
  SetBranch(&sel_npfps_eq_3,"npfps_eq_3",kBool);
  SetBranch(&sel_ntracks_eq_3,"ntracks_eq_3",kBool);
  SetBranch(&sel_containedparticles,"containedparticles",kBool);
  SetBranch(&sel_correctparticles,"correctparticles",kBool);
  SetBranch(&sel_momentum_threshold_passed_,"momentum_threshold_passed",kBool);
}
