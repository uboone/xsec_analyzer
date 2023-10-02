#include "CC1mu2p0pi.h"

#include "TreeUtils.hh"
#include "Functions.h"
#include "FiducialVolume.hh"

CC1mu2p0pi::CC1mu2p0pi() : SelectionBase("CC1mu2p0pi") {
}

void CC1mu2p0pi::DefineConstants() {
  TrueFV.X_Min = 10.;
  TrueFV.X_Max = 246.35;
  TrueFV.Y_Min = -106.5;
  TrueFV.Y_Max = 106.5;
  TrueFV.Z_Min = 10.;
  TrueFV.Z_Max = 1026.8;
}

void CC1mu2p0pi::ComputeRecoObservables(AnalysisEvent* Event) {
}

void CC1mu2p0pi::ComputeTrueObservables(AnalysisEvent* Event) {
}

EventCategory CC1mu2p0pi::CategorizeEvent(AnalysisEvent* Event)	{
  return kUnknown;
}

//Taken from https://github.com/ssfehlberg/CC2p-Event-Selection/blob/9492ff121a2eb884f464e1c166d067f217a04900/PeLEE_ntuples/mc_efficiency.C#L109-L112
bool CC1mu2p0pi::DefineSignal(AnalysisEvent* Event) {

  //==============================================================================================================================
  //DB Calculate the values which we need
  int mc_n_threshold_muon = 0;
  int mc_n_threshold_proton = 0;
  int mc_n_threshold_pion0 = 0;
  int mc_n_threshold_pionpm = 0;
  
  for ( size_t p = 0u; p < Event->mc_nu_daughter_pdg_->size(); ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    float energy = Event->mc_nu_daughter_energy_->at( p );
    if ( std::abs(pdg) == MUON) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(MUON_MASS, 2) );
      if ( mom > 0.1 && mom < 1.2){
	mc_n_threshold_muon++;
      }
    } else if (std::abs(pdg) == PROTON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PROTON_MASS, 2) );
      if ( mom > 0.3 && mom < 1.0){
	mc_n_threshold_proton++;
      }
    } else if ( pdg == PI_ZERO ) {
      mc_n_threshold_pion0++;
      
    } else if (std::abs(pdg) == PI_PLUS ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PI_PLUS_MASS, 2) );
      if ( mom > 0.065 ) {
	mc_n_threshold_pionpm++;
      }
    }
  }  

  //==============================================================================================================================
  //Calculate the booleans related to the different signal cuts
  
  //DB Discussions (https://microboone.slack.com/archives/C05TCS17EHL/p1695988699125549) - Afro says we should not be using Space Charge Effects (SCE) in the true FV definition
  //Currently included for validation purposes
  //sig_truevertex_in_fv_ = point_inside_FV(TrueFV, Event->mc_nu_sce_vx_, Event->mc_nu_sce_vy_, Event->mc_nu_sce_vz_);
  sig_truevertex_in_fv_ = point_inside_FV(TrueFV, Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_);
  
  sig_ccnc_ = (Event->mc_nu_ccnc_ == CHARGED_CURRENT);
  sig_is_numu_ = (Event->mc_nu_pdg_ == MUON_NEUTRINO);
  sig_two_protons_above_thresh_ = (mc_n_threshold_proton == 2);
  sig_one_muon_above_thres_ = (mc_n_threshold_muon == 1);
  sig_no_pions_ = ((mc_n_threshold_pion0 == 0) && (mc_n_threshold_pionpm == 0));

  //==============================================================================================================================
  //Is the event signal?
  
  bool IsSignal = sig_ccnc_ && sig_is_numu_ && sig_two_protons_above_thresh_ && sig_one_muon_above_thres_ && sig_no_pions_ && sig_truevertex_in_fv_;
  return IsSignal;
}

bool CC1mu2p0pi::Selection(AnalysisEvent* Event) {
  FiducialVolume FV;
  FV.X_Min = 10.;
  FV.X_Max = 246.35;
  FV.Y_Min = -106.5;
  FV.Y_Max = 106.5;
  FV.Z_Min = 10.;
  FV.Z_Max = 1026.8;

  FiducialVolume FV_noBorder;
  FV_noBorder.X_Min = 0.;
  FV_noBorder.X_Max = 256.35;
  FV_noBorder.Y_Min = -116.5;
  FV_noBorder.Y_Max = 116.5;
  FV_noBorder.Z_Min = 0.;
  FV_noBorder.Z_Max = 1036.8;
  
  //==============================================================================================================================
  //Vertex in FV?
  float_t x = Event->nu_vx_;
  float_t y = Event->nu_vy_;
  float_t z = Event->nu_vz_;

  sel_reco_vertex_in_FV_ = point_inside_FV(FV,x,y,z);

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

  //DB Discussion (https://microboone.slack.com/archives/D04A8CUB1EW/p1691519899965369)
  //   - Logic fault in Samantha's original code: https://github.com/ssfehlberg/CC2p-Event-Selection/blob/9492ff121a2eb884f464e1c166d067f217a04900/PeLEE_ntuples/twoproton_pelee_bnb.C#L138
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
  bool Contained = true;

  for(int i = 0; i < Event->num_pf_particles_; i ++){
    bool StartContained_i = point_inside_FV(FV,Event->track_startx_->at(i),Event->track_starty_->at(i),Event->track_startz_->at(i));
    bool EndContained_i = point_inside_FV(FV_noBorder,Event->track_endx_->at(i),Event->track_endy_->at(i),Event->track_endz_->at(i));

    if (!StartContained_i || !EndContained_i) {
      Contained = false;
    }
  }
  if (Contained) sel_containedparticles = true;

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

void CC1mu2p0pi::DefineOutputBranches() {
  SetBranch(&sel_reco_vertex_in_FV_,"sel_reco_vertex_in_FV",kBool);
  SetBranch(&sel_has_muon_candidate_,"sel_has_muon_candidate",kBool);
  SetBranch(&sel_nu_mu_cc_,"sel_nu_mu_cc",kBool);
  SetBranch(&sel_npfps_eq_3,"sel_npfps_eq_3",kBool);
  SetBranch(&sel_ntracks_eq_3,"sel_ntracks_eq_3",kBool);
  SetBranch(&sel_containedparticles,"sel_containedparticles",kBool);
  SetBranch(&sel_correctparticles,"sel_correctparticles",kBool);
  SetBranch(&sel_momentum_threshold_passed_,"sel_momentum_threshold_passed",kBool);

  SetBranch(&sig_truevertex_in_fv_,"sig_truevertex_in_fv",kBool);
  SetBranch(&sig_ccnc_,"sig_ccnc_",kBool);
  SetBranch(&sig_is_numu_,"sig_is_numu_",kBool);
  SetBranch(&sig_two_protons_above_thresh_,"sig_two_protons_above_thresh_",kBool);
  SetBranch(&sig_one_muon_above_thres_,"sig_one_muon_above_thres_",kBool);
  SetBranch(&sig_no_pions_,"sig_no_pions_",kBool);
}
