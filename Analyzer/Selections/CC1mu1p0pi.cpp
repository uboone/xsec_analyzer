#include "CC1mu1p0pi.h"

#include "TreeUtils.hh"

CC1mu1p0pi::CC1mu1p0pi() : SelectionBase("CC1mu1p0pi") {
}

void CC1mu1p0pi::ComputeObservables() {

}

bool CC1mu1p0pi::InFV(double x, double y, double z) {
  double FVx = 256.;
  double FVy = 230;
  double FVz = 1036.;
  double borderx = 10.;
  double bordery = 10.;
  double borderz = 10.;

  if (!(x < (FVx - borderx) && (x > borderx))) {
    return false;
  }
  if (!(y < (FVy/2. - bordery) && (y > (-FVy/2. + bordery)))) {
    return false;
  }
  if (!(z < (FVz - borderz) && (z > borderz))) {
    return false;
  }
  return true;
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

  sel_nuvertex_contained_ = InFV(Event->nu_vx_,Event->nu_vy_,Event->nu_vz_);

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
    
    bool CandidateMuonTrackStartContainment = InFV(MuonTrackStartX,MuonTrackStartY,MuonTrackStartZ);
    bool CandidateMuonTrackEndContainment = InFV(MuonTrackEndX,MuonTrackEndY,MuonTrackEndZ);
    
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
    
    bool CandidateProtonTrackStartContainment = InFV(ProtonTrackStartX,ProtonTrackStartY,ProtonTrackStartZ);
    bool CandidateProtonTrackEndContainment = InFV(ProtonTrackEndX,ProtonTrackEndY,ProtonTrackEndZ);
    
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
  //Does everything pass selection?
  bool Passed = sel_nslice_eq_1_ && sel_nshower_eq_0_ && sel_ntrack_eq_2_
    && sel_muoncandidate_tracklike_ && sel_protoncandidate_tracklike_ && sel_nuvertex_contained_
    && sel_muoncandidate_above_p_thresh && sel_protoncandidate_above_p_thresh
    && sel_muoncandidate_contained && sel_protoncandidate_contained;

  return Passed;
}

void CC1mu1p0pi::DefineBranches() {
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
}
