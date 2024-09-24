// XSecAnalyzer includes
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/TreeUtils.hh"

#include "XSecAnalyzer/Selections/CC1mu1p0pi.hh"
#include "XSecAnalyzer/Selections/EventCategoriesXp.hh"

CC1mu1p0pi::CC1mu1p0pi() : SelectionBase("CC1mu1p0pi") {
  CalcType = kOpt1;
  fPP = new TF1("fPP","29.354172*x-14.674918",0.3,0.5);
}

void CC1mu1p0pi::DefineConstants() {
  DefineTrueFV(10.,246.,-105.,105.,10.,1026.);
  DefineRecoFV(10.,246.,-105.,105.,10.,1026.);
}

void CC1mu1p0pi::ComputeRecoObservables(AnalysisEvent* Event) {
  if (CandidateMuonIndex != BOGUS_INDEX && CandidateProtonIndex != BOGUS_INDEX) {
    TVector3 TVector3CandidateMuon(-1,-1,-1);
    TVector3 TVector3CandidateProton(-1,-1,-1);

    double CandidateMuMom = Event->track_range_mom_mu_->at(CandidateMuonIndex);
    double CandidateMuE_GeV = TMath::Sqrt( TMath::Power(CandidateMuMom,2.) + TMath::Power(MUON_MASS,2.)); // GeV

    //Adjust for non-contained muons

    double MuonTrackEndX = Event->track_endx_->at(CandidateMuonIndex);
    double MuonTrackEndY = Event->track_endy_->at(CandidateMuonIndex);
    double MuonTrackEndZ = Event->track_endz_->at(CandidateMuonIndex);
    bool CandidateMuonTrackEndContainment = point_inside_FV(ReturnRecoFV(),MuonTrackEndX,MuonTrackEndY,MuonTrackEndZ);
    // If exiting muon, switch to MCS and recalculate the energy
    if (!CandidateMuonTrackEndContainment) {
      CandidateMuMom = Event->track_mcs_mom_mu_->at(CandidateMuonIndex); // GeV/c
      CandidateMuE_GeV = TMath::Sqrt( TMath::Power(CandidateMuMom,2.) + TMath::Power(MUON_MASS,2.)); // GeV
    }

    double reco_Pmu_cos_theta = cos(Event->track_theta_->at(CandidateMuonIndex));
    double reco_Pmu_phi = Event->track_phi_->at(CandidateMuonIndex);
    double reco_Emu = TMath::Sqrt(CandidateMuMom*CandidateMuMom + MUON_MASS*MUON_MASS);

    TVector3CandidateMuon.SetMag(CandidateMuMom);
    TVector3CandidateMuon.SetTheta(TMath::ACos(reco_Pmu_cos_theta));
    TVector3CandidateMuon.SetPhi(reco_Pmu_phi);

    double CandidatePKE_GeV = Event->track_kinetic_energy_p_->at(CandidateProtonIndex); // GeV // Watch out, kinetic energy not energy
    double CandidatePE_GeV = CandidatePKE_GeV + PROTON_MASS; // GeV
    double CandidatePMom = TMath::Sqrt(TMath::Power(CandidatePE_GeV,2.) - TMath::Power(PROTON_MASS,2.)); // GeV/c

    // STV redefinition if P_p < 0.5 where the biases have been observed
    if (CandidatePMom < 0.5) {
      CandidatePMom = ( 1.-0.01*fPP->Eval(CandidatePMom) ) * CandidatePMom;
      CandidatePE_GeV = TMath::Sqrt(TMath::Power(CandidatePMom,2) + TMath::Power(PROTON_MASS,2.)); // GeV/c
    }

    double CandidateProtonTrackTheta = Event->track_theta_->at(CandidateProtonIndex);
    double CandidateProtonTrackPhi = Event->track_phi_->at(CandidateProtonIndex);

    TVector3CandidateProton.SetMag(CandidatePMom);
    TVector3CandidateProton.SetTheta(CandidateProtonTrackTheta); // rad
    TVector3CandidateProton.SetPhi(CandidateProtonTrackPhi); // rad

    STVTools.CalculateSTVs(TVector3CandidateMuon,TVector3CandidateProton,CandidateMuE_GeV,CandidatePE_GeV,CalcType);

    Reco_Pt = STVTools.ReturnPt();
    Reco_Ptx = STVTools.ReturnPtx();
    Reco_Pty = STVTools.ReturnPty();
    Reco_PL = STVTools.ReturnPL();
    Reco_Pn = STVTools.ReturnPn();
    Reco_PnPerp = STVTools.ReturnPnPerp();
    Reco_PnPerpx = STVTools.ReturnPnPerpx();
    Reco_PnPerpy = STVTools.ReturnPnPerpy();
    Reco_PnPar = STVTools.ReturnPnPar();
    Reco_DeltaAlphaT = STVTools.ReturnDeltaAlphaT();
    Reco_DeltaAlpha3Dq = STVTools.ReturnDeltaAlpha3Dq();
    Reco_DeltaAlpha3DMu = STVTools.ReturnDeltaAlpha3DMu();
    Reco_DeltaPhiT = STVTools.ReturnDeltaPhiT();
    Reco_DeltaPhi3D = STVTools.ReturnDeltaPhi3D();
    Reco_ECal = STVTools.ReturnECal();
    Reco_EQE = STVTools.ReturnEQE();
    Reco_Q2 = STVTools.ReturnQ2();
    Reco_A = STVTools.ReturnA();
    Reco_EMiss = STVTools.ReturnEMiss();
    Reco_kMiss = STVTools.ReturnkMiss();
    Reco_PMiss = STVTools.ReturnPMiss();
    Reco_PMissMinus = STVTools.ReturnPMissMinus();

    //Now deal with the backtracking

    double CandidateMuonPx = Event->pfp_true_px_->at(CandidateMuonIndex);
    double CandidateMuonPy = Event->pfp_true_py_->at(CandidateMuonIndex);
    double CandidateMuonPz = Event->pfp_true_pz_->at(CandidateMuonIndex);
    TVector3 BackTrackCandidateMuonP(CandidateMuonPx,CandidateMuonPy,CandidateMuonPz);
    double BackTrackCandidateMuonTrackMomentum_GeV = BackTrackCandidateMuonP.Mag(); // GeV
    double BackTrackCandidateMuonTrack_E_GeV = TMath::Sqrt( TMath::Power(BackTrackCandidateMuonTrackMomentum_GeV,2.) + TMath::Power(MUON_MASS,2.) ); // GeV

    double CandidateProtonPx = Event->pfp_true_px_->at(CandidateProtonIndex);
    double CandidateProtonPy = Event->pfp_true_py_->at(CandidateProtonIndex);
    double CandidateProtonPz = Event->pfp_true_pz_->at(CandidateProtonIndex);
    TVector3 BackTrackCandidateProtonP(CandidateProtonPx,CandidateProtonPy,CandidateProtonPz);
    double BackTrackCandidateProtonTrackMomentum_GeV = BackTrackCandidateProtonP.Mag(); // GeV
    double BackTrackCandidateProtonTrack_E_GeV = TMath::Sqrt( TMath::Power(BackTrackCandidateProtonTrackMomentum_GeV,2.) + TMath::Power(PROTON_MASS,2.) ); // GeV

    STVTools.CalculateSTVs(BackTrackCandidateMuonP,BackTrackCandidateProtonP,BackTrackCandidateMuonTrack_E_GeV,BackTrackCandidateProtonTrack_E_GeV,CalcType);

    BackTrack_Pt = STVTools.ReturnPt();
    BackTrack_Ptx = STVTools.ReturnPtx();
    BackTrack_Pty = STVTools.ReturnPty();
    BackTrack_PL = STVTools.ReturnPL();
    BackTrack_Pn = STVTools.ReturnPn();
    BackTrack_PnPerp = STVTools.ReturnPnPerp();
    BackTrack_PnPerpx = STVTools.ReturnPnPerpx();
    BackTrack_PnPerpy = STVTools.ReturnPnPerpy();
    BackTrack_PnPar = STVTools.ReturnPnPar();
    BackTrack_DeltaAlphaT = STVTools.ReturnDeltaAlphaT();
    BackTrack_DeltaAlpha3Dq = STVTools.ReturnDeltaAlpha3Dq();
    BackTrack_DeltaAlpha3DMu = STVTools.ReturnDeltaAlpha3DMu();
    BackTrack_DeltaPhiT = STVTools.ReturnDeltaPhiT();
    BackTrack_DeltaPhi3D = STVTools.ReturnDeltaPhi3D();
    BackTrack_ECal = STVTools.ReturnECal();
    BackTrack_EQE = STVTools.ReturnEQE();
    BackTrack_Q2 = STVTools.ReturnQ2();
    BackTrack_A = STVTools.ReturnA();
    BackTrack_EMiss = STVTools.ReturnEMiss();
    BackTrack_kMiss = STVTools.ReturnkMiss();
    BackTrack_PMiss = STVTools.ReturnPMiss();
    BackTrack_PMissMinus = STVTools.ReturnPMissMinus();
  }

}

void CC1mu1p0pi::ComputeTrueObservables(AnalysisEvent* Event) {
  bool isInteresting = sig_one_muon_above_thresh_ && sig_one_proton_above_thresh_ && sig_no_pions_ && sig_no_heavy_mesons_;
  if (isInteresting) {
    double Muon_MCParticlePx = Event->mc_nu_daughter_px_->at(truemuonindex);
    double Muon_MCParticlePy = Event->mc_nu_daughter_px_->at(truemuonindex);
    double Muon_MCParticlePz = Event->mc_nu_daughter_px_->at(truemuonindex);
    TVector3 Muon_TVector3True(Muon_MCParticlePx,Muon_MCParticlePy,Muon_MCParticlePz);
    double Muon_TrueMomentum_GeV = Muon_TVector3True.Mag(); // GeV
    double Muon_TrueE_GeV = TMath::Sqrt( TMath::Power(Muon_TrueMomentum_GeV,2.) + TMath::Power(MUON_MASS,2.) ); // GeV

    double Proton_MCParticlePx = Event->mc_nu_daughter_px_->at(trueprotonindex);
    double Proton_MCParticlePy = Event->mc_nu_daughter_px_->at(trueprotonindex);
    double Proton_MCParticlePz = Event->mc_nu_daughter_px_->at(trueprotonindex);
    TVector3 Proton_TVector3True(Proton_MCParticlePx,Proton_MCParticlePy,Proton_MCParticlePz);
    double Proton_TrueMomentum_GeV = Proton_TVector3True.Mag(); // GeV
    double Proton_TrueE_GeV = TMath::Sqrt( TMath::Power(Proton_TrueMomentum_GeV,2.) + TMath::Power(PROTON_MASS,2.) ); // GeV

    STVTools.CalculateSTVs(Muon_TVector3True,Proton_TVector3True,Muon_TrueE_GeV,Proton_TrueE_GeV,CalcType);

    True_Pt = STVTools.ReturnPt();
    True_Ptx = STVTools.ReturnPtx();
    True_Pty = STVTools.ReturnPty();
    True_PL = STVTools.ReturnPL();
    True_Pn = STVTools.ReturnPn();
    True_PnPerp = STVTools.ReturnPnPerp();
    True_PnPerpx = STVTools.ReturnPnPerpx();
    True_PnPerpy = STVTools.ReturnPnPerpy();
    True_PnPar = STVTools.ReturnPnPar();
    True_DeltaAlphaT = STVTools.ReturnDeltaAlphaT();
    True_DeltaAlpha3Dq = STVTools.ReturnDeltaAlpha3Dq();
    True_DeltaAlpha3DMu = STVTools.ReturnDeltaAlpha3DMu();
    True_DeltaPhiT = STVTools.ReturnDeltaPhiT();
    True_DeltaPhi3D = STVTools.ReturnDeltaPhi3D();
    True_ECal = STVTools.ReturnECal();
    True_EQE = STVTools.ReturnEQE();
    True_Q2 = STVTools.ReturnQ2();
    True_A = STVTools.ReturnA();
    True_EMiss = STVTools.ReturnEMiss();
    True_kMiss = STVTools.ReturnkMiss();
    True_PMiss = STVTools.ReturnPMiss();
    True_PMissMinus = STVTools.ReturnPMissMinus();
  }
}

int CC1mu1p0pi::CategorizeEvent(AnalysisEvent* Event) {
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

  //Boolean which basically MC Signal selection without requesting a particular number of protons (N >= 1)
  bool Is_CC1muNp0pi_Event = (sig_mc_n_threshold_proton >= 1) && sig_no_pions_ && sig_one_muon_above_thresh_ && sig_no_heavy_mesons_;

  if ( Is_CC1muNp0pi_Event ) {
    if (sig_mc_n_threshold_proton == 1) {
      if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCC1p0pi_CCQE; // QE
      else if ( Event->mc_nu_interaction_type_ == 10 ) return kNuMuCC1p0pi_CCMEC; // MEC
      else if ( Event->mc_nu_interaction_type_ == 1 ) return kNuMuCC1p0pi_CCRES; // RES
      else return kNuMuCCMp0pi_Other;
    } else if (sig_mc_n_threshold_proton == 2) {
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
  else if (!sig_no_pions_) {
    return kNuMuCCNpi;
  } else if (sig_mc_n_threshold_proton == 0) {
    if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCC0p0pi_CCQE; // QE
    else if ( Event->mc_nu_interaction_type_ == 10 ) return kNuMuCC0p0pi_CCMEC; // MEC
    else if ( Event->mc_nu_interaction_type_ == 1 ) return kNuMuCC0p0pi_CCRES; // RES
    else return kNuMuCC0p0pi_Other;
  }
  return kNuMuCCOther;
}

bool CC1mu1p0pi::DefineSignal(AnalysisEvent* Event) {
  //==============================================================================================================================
  //DB Calculate the values which we need
  sig_mc_n_threshold_muon = 0;
  sig_mc_n_threshold_proton = 0;
  sig_mc_n_threshold_pion0 = 0;
  sig_mc_n_threshold_pionpm = 0;
  sig_mc_n_heaviermeson = 0;

  for ( size_t p = 0u; p < Event->mc_nu_daughter_pdg_->size(); ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    TVector3 MCParticle(Event->mc_nu_daughter_px_->at(p),Event->mc_nu_daughter_py_->at(p),Event->mc_nu_daughter_pz_->at(p));
    double ParticleMomentum = MCParticle.Mag();

    if ( pdg == MUON && ParticleMomentum >= 0.1 ) {sig_mc_n_threshold_muon++; truemuonindex = p;}
    else if ( pdg == PROTON && ParticleMomentum >= 0.3 ) {sig_mc_n_threshold_proton++; trueprotonindex = p;}
    else if ( pdg == PI_ZERO ) {sig_mc_n_threshold_pion0++;}
    else if ( std::abs(pdg) == PI_PLUS && ParticleMomentum >= 0.07) {sig_mc_n_threshold_pionpm++;}
    else if ( pdg != PI_ZERO && std::abs(pdg) != PI_PLUS && is_meson_or_antimeson(pdg) ) {sig_mc_n_heaviermeson++;}
  }

  //==============================================================================================================================
  //Calculate the booleans related to the different signal cuts

  sig_truevertex_in_fv_ = point_inside_FV(ReturnTrueFV(), Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_);

  sig_ccnc_= (Event->mc_nu_ccnc_ == CHARGED_CURRENT);
  sig_is_numu_ = (Event->mc_nu_pdg_ == MUON_NEUTRINO);
  sig_one_muon_above_thresh_ = (sig_mc_n_threshold_muon == 1);
  sig_one_proton_above_thresh_ = (sig_mc_n_threshold_proton == 1);
  sig_no_pions_ = ((sig_mc_n_threshold_pion0 == 0) && (sig_mc_n_threshold_pionpm == 0));
  sig_no_heavy_mesons_ = (sig_mc_n_heaviermeson == 0);

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

    // STV redefinition if P_p < 0.5 where the biases have been observed
    if (CandidatePMom < 0.5) {
      CandidatePMom = ( 1.-0.01*fPP->Eval(CandidatePMom) ) * CandidatePMom;
      CandidatePE_GeV = TMath::Sqrt(TMath::Power(CandidatePMom,2) + TMath::Power(PROTON_MASS,2.)); // GeV/c
    }

    // Comes from /uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Constants.hh:L1762
    if (CandidateMuMom >= MUON_P_MIN_MOM_CUT) sel_muoncandidate_above_p_thresh = true;

    // Comes from /uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Constants.hh:L1775
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
  SetBranch(&CandidateMuonIndex,"CandidateMuonIndex",kInteger);
  SetBranch(&CandidateProtonIndex,"CandidateProtonIndex",kInteger);

  SetBranch(&sig_truevertex_in_fv_,"sig_truevertex_in_fv",kBool);
  SetBranch(&sig_ccnc_,"sig_ccnc",kBool);
  SetBranch(&sig_is_numu_,"sig_is_numu",kBool);
  SetBranch(&sig_one_muon_above_thresh_,"sig_one_muon_above_thresh",kBool);
  SetBranch(&sig_one_proton_above_thresh_,"sig_one_proton_above_thresh",kBool);
  SetBranch(&sig_no_pions_,"sig_no_pions",kBool);
  SetBranch(&sig_no_heavy_mesons_,"sig_no_heavy_mesons",kBool);
  SetBranch(&sig_mc_n_threshold_muon,"mc_n_threshold_muon",kInteger);
  SetBranch(&sig_mc_n_threshold_proton,"mc_n_threshold_proton",kInteger);
  SetBranch(&sig_mc_n_threshold_pion0,"mc_n_threshold_pion0",kInteger);
  SetBranch(&sig_mc_n_threshold_pionpm,"mc_n_threshold_pionpm",kInteger);
  SetBranch(&sig_mc_n_heaviermeson,"mc_n_heaviermeson",kInteger);
  SetBranch(&truemuonindex,"truemuonindex",kInteger);
  SetBranch(&trueprotonindex,"trueprotonindex",kInteger);

  SetBranch(&Reco_Pt,"Reco_Pt",kDouble);
  SetBranch(&Reco_Ptx,"Reco_Ptx",kDouble);
  SetBranch(&Reco_Pty,"Reco_Pty",kDouble);
  SetBranch(&Reco_PL,"Reco_PL",kDouble);
  SetBranch(&Reco_Pn,"Reco_Pn",kDouble);
  SetBranch(&Reco_PnPerp,"Reco_PnPerp",kDouble);
  SetBranch(&Reco_PnPerpx,"Reco_PnPerpx",kDouble);
  SetBranch(&Reco_PnPerpy,"Reco_PnPerpy",kDouble);
  SetBranch(&Reco_PnPar,"Reco_PnPar",kDouble);
  SetBranch(&Reco_DeltaAlphaT,"Reco_DeltaAlphaT",kDouble);
  SetBranch(&Reco_DeltaAlpha3Dq,"Reco_DeltaAlpha3Dq",kDouble);
  SetBranch(&Reco_DeltaAlpha3DMu,"Reco_DeltaAlpha3DMu",kDouble);
  SetBranch(&Reco_DeltaPhiT,"Reco_DeltaPhiT",kDouble);
  SetBranch(&Reco_DeltaPhi3D,"Reco_DeltaPhi3D",kDouble);
  SetBranch(&Reco_ECal,"Reco_ECal",kDouble);
  SetBranch(&Reco_EQE,"Reco_EQE",kDouble);
  SetBranch(&Reco_Q2,"Reco_Q2",kDouble);
  SetBranch(&Reco_A,"Reco_A",kDouble);
  SetBranch(&Reco_EMiss,"Reco_EMiss",kDouble);
  SetBranch(&Reco_kMiss,"Reco_kMiss",kDouble);
  SetBranch(&Reco_PMiss,"Reco_PMiss",kDouble);
  SetBranch(&Reco_PMissMinus,"Reco_PMissMinus",kDouble);
  SetBranch(&BackTrack_Pt,"BackTrack_Pt",kDouble);
  SetBranch(&BackTrack_Ptx,"BackTrack_Ptx",kDouble);
  SetBranch(&BackTrack_Pty,"BackTrack_Pty",kDouble);
  SetBranch(&BackTrack_PL,"BackTrack_PL",kDouble);
  SetBranch(&BackTrack_Pn,"BackTrack_Pn",kDouble);
  SetBranch(&BackTrack_PnPerp,"BackTrack_PnPerp",kDouble);
  SetBranch(&BackTrack_PnPerpx,"BackTrack_PnPerpx",kDouble);
  SetBranch(&BackTrack_PnPerpy,"BackTrack_PnPerpy",kDouble);
  SetBranch(&BackTrack_PnPar,"BackTrack_PnPar",kDouble);
  SetBranch(&BackTrack_DeltaAlphaT,"BackTrack_DeltaAlphaT",kDouble);
  SetBranch(&BackTrack_DeltaAlpha3Dq,"BackTrack_DeltaAlpha3Dq",kDouble);
  SetBranch(&BackTrack_DeltaAlpha3DMu,"BackTrack_DeltaAlpha3DMu",kDouble);
  SetBranch(&BackTrack_DeltaPhiT,"BackTrack_DeltaPhiT",kDouble);
  SetBranch(&BackTrack_DeltaPhi3D,"BackTrack_DeltaPhi3D",kDouble);
  SetBranch(&BackTrack_ECal,"BackTrack_ECal",kDouble);
  SetBranch(&BackTrack_EQE,"BackTrack_EQE",kDouble);
  SetBranch(&BackTrack_Q2,"BackTrack_Q2",kDouble);
  SetBranch(&BackTrack_A,"BackTrack_A",kDouble);
  SetBranch(&BackTrack_EMiss,"BackTrack_EMiss",kDouble);
  SetBranch(&BackTrack_kMiss,"BackTrack_kMiss",kDouble);
  SetBranch(&BackTrack_PMiss,"BackTrack_PMiss",kDouble);
  SetBranch(&BackTrack_PMissMinus,"BackTrack_PMissMinus",kDouble);
  SetBranch(&True_Pt,"True_Pt",kDouble);
  SetBranch(&True_Ptx,"True_Ptx",kDouble);
  SetBranch(&True_Pty,"True_Pty",kDouble);
  SetBranch(&True_PL,"True_PL",kDouble);
  SetBranch(&True_Pn,"True_Pn",kDouble);
  SetBranch(&True_PnPerp,"True_PnPerp",kDouble);
  SetBranch(&True_PnPerpx,"True_PnPerpx",kDouble);
  SetBranch(&True_PnPerpy,"True_PnPerpy",kDouble);
  SetBranch(&True_PnPar,"True_PnPar",kDouble);
  SetBranch(&True_DeltaAlphaT,"True_DeltaAlphaT",kDouble);
  SetBranch(&True_DeltaAlpha3Dq,"True_DeltaAlpha3Dq",kDouble);
  SetBranch(&True_DeltaAlpha3DMu,"True_DeltaAlpha3DMu",kDouble);
  SetBranch(&True_DeltaPhiT,"True_DeltaPhiT",kDouble);
  SetBranch(&True_DeltaPhi3D,"True_DeltaPhi3D",kDouble);
  SetBranch(&True_ECal,"True_ECal",kDouble);
  SetBranch(&True_EQE,"True_EQE",kDouble);
  SetBranch(&True_Q2,"True_Q2",kDouble);
  SetBranch(&True_A,"True_A",kDouble);
  SetBranch(&True_EMiss,"True_EMiss",kDouble);
  SetBranch(&True_kMiss,"True_kMiss",kDouble);
  SetBranch(&True_PMiss,"True_PMiss",kDouble);
  SetBranch(&True_PMissMinus,"True_PMissMinus",kDouble);
}


void CC1mu1p0pi::DefineCategoryMap() {
  // Use the shared category map for 1p/2p/Np/Xp
  categ_map_ = CC1muXp_MAP;
}
