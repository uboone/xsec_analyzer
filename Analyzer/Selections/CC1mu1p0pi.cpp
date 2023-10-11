#include "CC1mu1p0pi.h"

#include "TreeUtils.hh"
#include "FiducialVolume.hh"
#include "EventCategory.hh"
#include "Functions.h"

CC1mu1p0pi::CC1mu1p0pi() : SelectionBase("CC1mu1p0pi") {
}

void CC1mu1p0pi::DefineConstants() {
  DefineTrueFV(10.,246.,-105.,105.,10.,1026.);
  DefineRecoFV(10.,246.,-105.,105.,10.,1026.);
}

void CC1mu1p0pi::ComputeRecoObservables(AnalysisEvent* Event) {
}

void CC1mu1p0pi::ComputeTrueObservables(AnalysisEvent* Event) {
}

EventCategory CC1mu1p0pi::CategorizeEvent(AnalysisEvent* Event) {
  return kUnknown;
}

bool CC1mu1p0pi::DefineSignal(AnalysisEvent* Event) {
  //==============================================================================================================================
  //DB Calculate the values which we need
  int mc_n_threshold_muon = 0;
  int mc_n_threshold_proton = 0;
  int mc_n_threshold_pion0 = 0;
  int mc_n_threshold_pionpm = 0;
  int mc_n_heaviermeson = 0;
  
  for ( size_t p = 0u; p < Event->mc_nu_daughter_pdg_->size(); ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    TVector3 MCParticle(Event->mc_nu_daughter_px_->at(p),Event->mc_nu_daughter_py_->at(p),Event->mc_nu_daughter_pz_->at(p));
    double ParticleMomentum = MCParticle.Mag();
    
    if ( pdg == MUON && ParticleMomentum >= 0.1 ) {mc_n_threshold_muon++;}
    else if ( pdg == PROTON && ParticleMomentum >= 0.3 ) {mc_n_threshold_proton++;}
    else if ( pdg == PI_ZERO ) {mc_n_threshold_pion0++;}
    else if ( std::abs(pdg) == PI_PLUS && ParticleMomentum >= 0.07) {mc_n_threshold_pionpm++;}
    else if ( pdg != PI_ZERO && std::abs(pdg) != PI_PLUS && is_meson_or_antimeson(pdg) ) {mc_n_heaviermeson++;}
  }

  //==============================================================================================================================
  //Calculate the booleans related to the different signal cuts
  
  sig_truevertex_in_fv_ = point_inside_FV(ReturnTrueFV(), Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_);
  
  sig_ccnc_= (Event->mc_nu_ccnc_ == CHARGED_CURRENT);
  sig_is_numu_ = (Event->mc_nu_pdg_ == MUON_NEUTRINO);
  sig_one_muon_above_thresh_ = (mc_n_threshold_muon == 1);
  sig_one_proton_above_thresh_ = (mc_n_threshold_proton == 1);
  sig_no_pions_ = ((mc_n_threshold_pion0 == 0) && (mc_n_threshold_pionpm == 0));
  sig_no_heavy_mesons_ = (mc_n_heaviermeson == 0);
  
  //==============================================================================================================================
  //Is the event signal?

  bool isSignal = sig_truevertex_in_fv_ && sig_ccnc_ && sig_is_numu_ && sig_one_muon_above_thresh_ && sig_one_proton_above_thresh_ && sig_no_pions_ && sig_no_heavy_mesons_;
  return isSignal;
}

bool CC1mu1p0pi::Selection(AnalysisEvent* Event) {
  //==============================================================================================================================
  // Requirement for exactly one neutrino slice
  if (Event->nslice_ == 1) sel_nslice_eq_1_ = true;

  //==============================================================================================================================
  // Requirement for exactly 2 tracks and 0 showers for a CC1p0pi selection

  int reco_shower_count = 0;
  int reco_track_count = 0;
  std::vector<int> CandidateIndex;
  
  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {    
    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2u ) continue;
    
    float tscore = Event->pfp_track_score_->at( p );
    if ( tscore <= TRACK_SCORE_CUT ) {
      ++reco_shower_count;
    } else {
      ++reco_track_count;
      CandidateIndex.push_back(p);
    }
    
  }

  if (reco_shower_count == 0) sel_nshower_eq_0_ = true;
  if (reco_track_count == 2) sel_ntrack_eq_2_ = true;

  //==============================================================================================================================
  // Identify candidate muon & proton
  // Muon = the one with the highest LLR PID Score
  // Proton = the one with the lowest LLR PID Score

  int CandidateMuonIndex = -1.;
  int CandidateProtonIndex = -1.;
  
  if (sel_ntrack_eq_2_) {
    float first_pid_score = (Event->track_llr_pid_score_)->at( CandidateIndex.at(0) );
    float second_pid_score = (Event->track_llr_pid_score_)->at( CandidateIndex.at(1) );
    
    if (first_pid_score > second_pid_score) {
      CandidateMuonIndex = CandidateIndex.at(0);
      CandidateProtonIndex = CandidateIndex.at(1);
    } else {
      CandidateMuonIndex = CandidateIndex.at(1);
      CandidateProtonIndex = CandidateIndex.at(0);
    }
    
    if ((Event->pfp_reco_pdg_)->at(CandidateMuonIndex) == MUON) {
      sel_muoncandidate_tracklike_ = true;
    }
    if ((Event->pfp_reco_pdg_)->at(CandidateProtonIndex) == MUON) {
      sel_protoncandidate_tracklike_ = true;
    }
  }

  //==============================================================================================================================
  //Nuetrino vertex in FV?

  sel_nuvertex_contained_ = point_inside_FV(ReturnRecoFV(),Event->nu_vx_,Event->nu_vy_,Event->nu_vz_);

  //==============================================================================================================================
  //Containment check on the muon and proton

  if (sel_ntrack_eq_2_) {

    //Muon variables
    double MuonTrackStartX = Event->track_startx_->at(CandidateMuonIndex);
    double MuonTrackStartY = Event->track_starty_->at(CandidateMuonIndex);
    double MuonTrackStartZ = Event->track_startz_->at(CandidateMuonIndex);
    
    double MuonTrackEndX = Event->track_endx_->at(CandidateMuonIndex);
    double MuonTrackEndY = Event->track_endy_->at(CandidateMuonIndex);
    double MuonTrackEndZ = Event->track_endz_->at(CandidateMuonIndex);
    
    bool CandidateMuonTrackStartContainment = point_inside_FV(ReturnRecoFV(),MuonTrackStartX,MuonTrackStartY,MuonTrackStartZ);
    bool CandidateMuonTrackEndContainment = point_inside_FV(ReturnRecoFV(),MuonTrackEndX,MuonTrackEndY,MuonTrackEndZ);
    
    double CandidateMuMom = Event->track_range_mom_mu_->at(CandidateMuonIndex); // GeV/c
    double CandidateMuE_GeV = TMath::Sqrt(TMath::Power(CandidateMuMom,2.) + TMath::Power(MUON_MASS,2.)); // GeV
    
    // If exiting muon, switch to MCS and recalculate the energy
    if (!CandidateMuonTrackEndContainment) {
      CandidateMuMom = Event->track_mcs_mom_mu_->at(CandidateMuonIndex); // GeV/c
      CandidateMuE_GeV = TMath::Sqrt( TMath::Power(CandidateMuMom,2.) + TMath::Power(MUON_MASS,2.)); // GeV    
    }
    
    //Proton variables
    double ProtonTrackStartX = Event->track_startx_->at(CandidateProtonIndex);
    double ProtonTrackStartY = Event->track_starty_->at(CandidateProtonIndex);
    double ProtonTrackStartZ = Event->track_startz_->at(CandidateProtonIndex);
    
    double ProtonTrackEndX = Event->track_endx_->at(CandidateProtonIndex);
    double ProtonTrackEndY = Event->track_endy_->at(CandidateProtonIndex);
    double ProtonTrackEndZ = Event->track_endz_->at(CandidateProtonIndex);
    
    bool CandidateProtonTrackStartContainment = point_inside_FV(ReturnRecoFV(),ProtonTrackStartX,ProtonTrackStartY,ProtonTrackStartZ);
    bool CandidateProtonTrackEndContainment = point_inside_FV(ReturnRecoFV(),ProtonTrackEndX,ProtonTrackEndY,ProtonTrackEndZ);
    
    double CandidatePKE_GeV = Event->track_kinetic_energy_p_->at(CandidateProtonIndex); // GeV // Watch out, kinetic energy not energy
    double CandidatePE_GeV = CandidatePKE_GeV + PROTON_MASS; // GeV
    double CandidatePMom = TMath::Sqrt(TMath::Power(CandidatePE_GeV,2.) - TMath::Power(PROTON_MASS,2.)); // GeV/c

    // Comes from /uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Constants.h:L1762
    if (CandidateMuMom >= MUON_P_MIN_MOM_CUT) sel_muoncandidate_above_p_thresh = true;

    // Comes from /uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Constants.h:L1775
    if (CandidatePMom >= 0.3) sel_protoncandidate_above_p_thresh = true;

    if (CandidateMuonTrackStartContainment && CandidateMuonTrackEndContainment) sel_muoncandidate_contained = true;
    if (CandidateProtonTrackStartContainment && CandidateProtonTrackEndContainment) sel_protoncandidate_contained = true;
  }

  //==============================================================================================================================
  //Check that the muon momemnt estimators agree within 25%

  if (CandidateMuonIndex != -1) {
    double MuonMom_MCS = Event->track_mcs_mom_mu_->at(CandidateMuonIndex);
    double MuonMom_Range = Event->track_range_mom_mu_->at(CandidateMuonIndex);
  
    double Reso =  TMath::Abs(MuonMom_MCS - MuonMom_Range) / MuonMom_Range;
    sel_muon_momentum_quality = (Reso <= 0.25);
  }

  //==============================================================================================================================
  //Check for flipped tracks

  if (CandidateMuonIndex != -1 && CandidateProtonIndex != -1) {
    TVector3 VertexLocation(Event->nu_vx_,Event->nu_vy_,Event->nu_vz_);
    TVector3 Candidate_MuonTrack_Start(Event->track_startx_->at(CandidateMuonIndex),Event->track_starty_->at(CandidateMuonIndex),Event->track_startz_->at(CandidateMuonIndex));
    TVector3 Candidate_MuonTrack_End(Event->track_endx_->at(CandidateMuonIndex),Event->track_endy_->at(CandidateMuonIndex),Event->track_endz_->at(CandidateMuonIndex));
    TVector3 Candidate_ProtonTrack_Start(Event->track_startx_->at(CandidateProtonIndex),Event->track_starty_->at(CandidateProtonIndex),Event->track_startz_->at(CandidateProtonIndex));
    TVector3 Candidate_ProtonTrack_End(Event->track_endx_->at(CandidateProtonIndex),Event->track_endy_->at(CandidateProtonIndex),Event->track_endz_->at(CandidateProtonIndex));
    
    double Vertex_MuonTrackStart_Mag = (VertexLocation - Candidate_MuonTrack_Start).Mag();
    double Vertex_MuonTrackEnd_Mag = (VertexLocation - Candidate_MuonTrack_End).Mag();
    double Vertex_ProtonTrackStart_Mag = (VertexLocation - Candidate_ProtonTrack_Start).Mag();
    double Vertex_ProtonTrackEnd_Mag = (VertexLocation - Candidate_ProtonTrack_End).Mag();
    
    double MuonTrackStart_to_ProtonTrackStart_Mag = (Candidate_MuonTrack_Start - Candidate_ProtonTrack_Start).Mag();
    double MuonTrackEnd_to_ProtonTrackEnd_Mag = (Candidate_MuonTrack_End - Candidate_ProtonTrack_End).Mag();
    
    if ( !( (Vertex_MuonTrackStart_Mag > Vertex_MuonTrackEnd_Mag) || (Vertex_ProtonTrackStart_Mag > Vertex_ProtonTrackEnd_Mag) || (MuonTrackStart_to_ProtonTrackStart_Mag > MuonTrackEnd_to_ProtonTrackEnd_Mag) ) ) {
      sel_no_flipped_tracks_ = true;
    }
  }

  //==============================================================================================================================
  //Check Proton Candidate's LLH to be a proton

  if (CandidateMuonIndex != -1 && CandidateProtonIndex != -1) {
    double Candidate_Proton_LLR = Event->track_llr_pid_score_->at(CandidateProtonIndex);
    sel_proton_cand_passed_LLRCut = (Candidate_Proton_LLR < 0.05); 
  }

  //==============================================================================================================================
  //Apply kinematic cuts

  if (CandidateMuonIndex != -1) {
    //This is essentially duplicate of sel_muoncandidate_above_p_thresh cut
    double MuonMom = Event->track_range_mom_mu_->at(CandidateMuonIndex);
    if (! (MuonMom < 0.1 || MuonMom > 1.2) ) {sel_muon_momentum_in_range = true;}

    double MuonCosTheta = cos(Event->track_theta_->at(CandidateMuonIndex));
    if (! (MuonCosTheta < -1. || MuonCosTheta > 1.) ) {sel_muon_costheta_in_range = true;}

    double MuonPhi = Event->track_phi_->at(CandidateMuonIndex) * 180./ TMath::Pi();
    if (! (MuonPhi < -180. || MuonPhi > 180.) ) {sel_muon_phi_in_range = true;}
  }
  if (CandidateProtonIndex != -1) {
    //This is essentially duplicate of sel_protoncandidate_above_p_thresh cut
    double ProtonMomentum = TMath::Sqrt( TMath::Power(Event->track_kinetic_energy_p_->at(CandidateProtonIndex) + PROTON_MASS,2.) - TMath::Power(PROTON_MASS,2.));
    if (! (ProtonMomentum < 0.3 || ProtonMomentum > 1.) ) {sel_proton_momentum_in_range = true;}

    double ProtonCosTheta = cos(Event->track_theta_->at(CandidateProtonIndex));
    if (! (ProtonCosTheta < -1. || ProtonCosTheta > 1.) ) {sel_proton_costheta_in_range = true;}

    double ProtonPhi = Event->track_phi_->at(CandidateProtonIndex) * 180./ TMath::Pi();
    if (! (ProtonPhi < -180. || ProtonPhi > 180.) ) {sel_proton_phi_in_range = true;}
  }
  
  
  //==============================================================================================================================
  //Does everything pass selection?
  bool Passed = sel_nslice_eq_1_ && sel_nshower_eq_0_ && sel_ntrack_eq_2_
    && sel_muoncandidate_tracklike_ && sel_protoncandidate_tracklike_ && sel_nuvertex_contained_
    && sel_muoncandidate_above_p_thresh && sel_protoncandidate_above_p_thresh
    && sel_muoncandidate_contained && sel_protoncandidate_contained && sel_muon_momentum_quality
    && sel_no_flipped_tracks_ && sel_proton_cand_passed_LLRCut
    && sel_muon_momentum_in_range && sel_muon_costheta_in_range && sel_muon_phi_in_range
    && sel_proton_momentum_in_range && sel_proton_costheta_in_range && sel_proton_phi_in_range;

  return Passed;
}

void CC1mu1p0pi::DefineOutputBranches() {
  SetBranch(&sel_nslice_eq_1_,"nslice_eq_1",kBool);
  SetBranch(&sel_nshower_eq_0_,"nshower_eq_0_",kBool);
  SetBranch(&sel_ntrack_eq_2_,"ntrack_eq_2_",kBool);
  SetBranch(&sel_muoncandidate_tracklike_,"muoncandidate_tracklike",kBool);
  SetBranch(&sel_protoncandidate_tracklike_,"protoncandidate_tracklike",kBool);
  SetBranch(&sel_nuvertex_contained_,"nuvertex_contained_",kBool);
  SetBranch(&sel_muoncandidate_above_p_thresh,"muoncandidate_above_p_thresh",kBool);
  SetBranch(&sel_protoncandidate_above_p_thresh,"protoncandidate_above_p_thresh",kBool);
  SetBranch(&sel_muoncandidate_contained,"muoncandidate_contained",kBool);
  SetBranch(&sel_protoncandidate_contained,"protoncandidate_contained",kBool);
  SetBranch(&sel_muon_momentum_quality,"sel_muon_momentum_quality",kBool);
  SetBranch(&sel_no_flipped_tracks_,"sel_no_flipped_tracks",kBool);
  SetBranch(&sel_proton_cand_passed_LLRCut,"sel_proton_cand_passed_LLRCut",kBool);
  SetBranch(&sel_muon_momentum_in_range,"sel_muon_momentum_in_range",kBool);
  SetBranch(&sel_muon_costheta_in_range,"sel_muon_costheta_in_range",kBool);
  SetBranch(&sel_muon_phi_in_range,"sel_muon_phi_in_range",kBool);
  SetBranch(&sel_proton_momentum_in_range,"sel_proton_momentum_in_range",kBool);
  SetBranch(&sel_proton_costheta_in_range,"sel_proton_costheta_in_range",kBool);
  SetBranch(&sel_proton_phi_in_range,"sel_proton_phi_in_range",kBool);

  SetBranch(&sig_truevertex_in_fv_,"sig_truevertex_in_fv",kBool);
  SetBranch(&sig_ccnc_,"sig_ccnc",kBool);
  SetBranch(&sig_is_numu_,"sig_is_numu",kBool);
  SetBranch(&sig_one_muon_above_thresh_,"sig_one_muon_above_thresh",kBool);
  SetBranch(&sig_one_proton_above_thresh_,"sig_one_proton_above_thresh",kBool);
  SetBranch(&sig_no_pions_,"sig_no_pions",kBool);
  SetBranch(&sig_no_heavy_mesons_,"sig_no_heavy_mesons",kBool);

  SetBranch(&CandidateMuonIndex,"CandidateMuonIndex",kInteger);
  SetBranch(&CandidateProtonIndex,"CandidateProtonIndex",kInteger);
}
