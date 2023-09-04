#include "CC1mu1p0pi.h"

#include "TreeUtils.hh"
#include "FiducialVolume.hh"
#include "EventCategory.hh"
#include "Functions.h"

CC1mu1p0pi::CC1mu1p0pi() : SelectionBase("CC1mu1p0pi") {
}

void CC1mu1p0pi::DefineConstants() {
  TrueFV.X_Min = 10.;
  TrueFV.X_Max = 246.35;
  TrueFV.Y_Min = -106.5;
  TrueFV.Y_Max = 106.5;
  TrueFV.Z_Min = 10.;
  TrueFV.Z_Max = 1026.8;
}

void CC1mu1p0pi::ComputeObservables(AnalysisEvent* Event) {
}

bool CC1mu1p0pi::DefineSignal(AnalysisEvent* Event) {
  bool inFV = point_inside_FV(TrueFV, Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_);
  bool IsNuMu = (Event->mc_nu_pdg_ == MUON_NEUTRINO );

  bool NoFSMesons = true;
  bool NoChargedPiAboveThres = true;
  bool NoFSPi0s = true;
  bool MuonInRange = false;
  double ProtonMomentum = 0.;

  for ( size_t p = 0u; p < Event->mc_nu_daughter_pdg_->size(); ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    float energy = Event->mc_nu_daughter_energy_->at( p );

    // Do the general check for (anti)mesons first before considering
    // any individual PDG codes
    if ( is_meson_or_antimeson(pdg) ) {
      NoFSMesons = false;
    }

    if ( pdg == MUON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(MUON_MASS, 2) );
      if ( mom >= MUON_P_MIN_MOM_CUT && mom <= MUON_P_MAX_MOM_CUT ) {
        MuonInRange = true;
      }
    }
    else if ( pdg == PROTON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PROTON_MASS, 2) );
      if ( mom > ProtonMomentum ) ProtonMomentum = mom;
    }
    else if ( pdg == PI_ZERO ) {
      NoFSPi0s = false;
    }
    else if ( std::abs(pdg) == PI_PLUS ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PI_PLUS_MASS, 2) );
      if ( mom > CHARGED_PI_MOM_CUT ) {
        NoChargedPiAboveThres = false;
      }
    }
  }

  bool LeadProtonMomInRange = false;
  // Check that the leading proton has a momentum within the allowed range
  if ( ProtonMomentum >= LEAD_P_MIN_MOM_CUT && ProtonMomentum <= LEAD_P_MAX_MOM_CUT ) {
    LeadProtonMomInRange = true;
  }

  bool ReturnVal = inFV && IsNuMu && MuonInRange && LeadProtonMomInRange && NoFSMesons;
  return ReturnVal;

}

bool CC1mu1p0pi::Selection(AnalysisEvent* Event) {
  FiducialVolume FV;
  FV.X_Min = 10.;
  FV.X_Max = 246.;
  FV.Y_Min = -105;
  FV.Y_Max = 105;
  FV.Z_Min = 10.;
  FV.Z_Max = 1026.;
  
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

  sel_nuvertex_contained_ = point_inside_FV(FV,Event->nu_vx_,Event->nu_vy_,Event->nu_vz_);

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
    
    bool CandidateMuonTrackStartContainment = point_inside_FV(FV,MuonTrackStartX,MuonTrackStartY,MuonTrackStartZ);
    bool CandidateMuonTrackEndContainment = point_inside_FV(FV,MuonTrackEndX,MuonTrackEndY,MuonTrackEndZ);
    
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
    
    bool CandidateProtonTrackStartContainment = point_inside_FV(FV,ProtonTrackStartX,ProtonTrackStartY,ProtonTrackStartZ);
    bool CandidateProtonTrackEndContainment = point_inside_FV(FV,ProtonTrackEndX,ProtonTrackEndY,ProtonTrackEndZ);
    
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
