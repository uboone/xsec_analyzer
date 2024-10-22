// XSecAnalyzer includes
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/TreeUtils.hh"

#include "XSecAnalyzer/Selections/CC1mu1p0pi.hh"
#include "XSecAnalyzer/Selections/EventCategoriesXp.hh"

CC1mu1p0pi::CC1mu1p0pi() : SelectionBase( "CC1mu1p0pi" ) {
  CalcType = kOpt1;
  fPP = new TF1( "fPP", "29.354172*x-14.674918", 0.3, 0.5 );
}

void CC1mu1p0pi::define_constants() {
  this->define_true_FV( 10., 246., -105., 105., 10., 1026. );
  this->define_reco_FV( 10., 246., -105., 105., 10., 1026. );
}

void CC1mu1p0pi::compute_reco_observables( AnalysisEvent* Event ) {

  if ( CandidateMuonIndex != BOGUS_INDEX
    && CandidateProtonIndex != BOGUS_INDEX)
  {
    TVector3 TVector3CandidateMuon( -1, -1, -1 );
    TVector3 TVector3CandidateProton( -1, -1, -1 );

    double CandidateMuMom = Event->track_range_mom_mu_->at(CandidateMuonIndex);
    double CandidateMuE_GeV = TMath::Sqrt( TMath::Power(CandidateMuMom,2.)
      + TMath::Power(MUON_MASS,2.)); // GeV

    //Adjust for non-contained muons

    double MuonTrackEndX = Event->track_endx_->at(CandidateMuonIndex);
    double MuonTrackEndY = Event->track_endy_->at(CandidateMuonIndex);
    double MuonTrackEndZ = Event->track_endz_->at(CandidateMuonIndex);
    bool CandidateMuonTrackEndContainment = point_inside_FV( this->reco_FV(),
      MuonTrackEndX, MuonTrackEndY, MuonTrackEndZ );
    // If exiting muon, switch to MCS and recalculate the energy
    if ( !CandidateMuonTrackEndContainment ) {
      // GeV/c
      CandidateMuMom = Event->track_mcs_mom_mu_->at(CandidateMuonIndex);
      CandidateMuE_GeV = TMath::Sqrt( TMath::Power(CandidateMuMom,2.)
        + TMath::Power(MUON_MASS,2.)); // GeV
    }

    double reco_Pmu_cos_theta = std::cos(
      Event->track_theta_->at(CandidateMuonIndex) );
    double reco_Pmu_phi = Event->track_phi_->at(CandidateMuonIndex);
    double reco_Emu = TMath::Sqrt( CandidateMuMom*CandidateMuMom
      + MUON_MASS*MUON_MASS );

    TVector3CandidateMuon.SetMag(CandidateMuMom);
    TVector3CandidateMuon.SetTheta(TMath::ACos(reco_Pmu_cos_theta));
    TVector3CandidateMuon.SetPhi(reco_Pmu_phi);

    // GeV, watch out, kinetic energy not energy
    double CandidatePKE_GeV = Event->track_kinetic_energy_p_
      ->at( CandidateProtonIndex );
    double CandidatePE_GeV = CandidatePKE_GeV + PROTON_MASS; // GeV
    double CandidatePMom = TMath::Sqrt(TMath::Power(CandidatePE_GeV,2.)
      - TMath::Power(PROTON_MASS,2.)); // GeV/c

    // STV redefinition if P_p < 0.5 where the biases have been observed
    if (CandidatePMom < 0.5) {
      CandidatePMom = ( 1.-0.01*fPP->Eval(CandidatePMom) ) * CandidatePMom;
      CandidatePE_GeV = TMath::Sqrt(TMath::Power(CandidatePMom,2)
        + TMath::Power(PROTON_MASS,2.)); // GeV/c
    }

    double CandidateProtonTrackTheta = Event->track_theta_
      ->at( CandidateProtonIndex );
    double CandidateProtonTrackPhi = Event->track_phi_
      ->at( CandidateProtonIndex );

    TVector3CandidateProton.SetMag(CandidatePMom);
    TVector3CandidateProton.SetTheta(CandidateProtonTrackTheta); // rad
    TVector3CandidateProton.SetPhi(CandidateProtonTrackPhi); // rad

    STVTools stv_tools;
    stv_tools.CalculateSTVs( TVector3CandidateMuon, TVector3CandidateProton,
      CandidateMuE_GeV, CandidatePE_GeV, CalcType);

    Reco_Pt = stv_tools.ReturnPt();
    Reco_Ptx = stv_tools.ReturnPtx();
    Reco_Pty = stv_tools.ReturnPty();
    Reco_PL = stv_tools.ReturnPL();
    Reco_Pn = stv_tools.ReturnPn();
    Reco_PnPerp = stv_tools.ReturnPnPerp();
    Reco_PnPerpx = stv_tools.ReturnPnPerpx();
    Reco_PnPerpy = stv_tools.ReturnPnPerpy();
    Reco_PnPar = stv_tools.ReturnPnPar();
    Reco_DeltaAlphaT = stv_tools.ReturnDeltaAlphaT();
    Reco_DeltaAlpha3Dq = stv_tools.ReturnDeltaAlpha3Dq();
    Reco_DeltaAlpha3DMu = stv_tools.ReturnDeltaAlpha3DMu();
    Reco_DeltaPhiT = stv_tools.ReturnDeltaPhiT();
    Reco_DeltaPhi3D = stv_tools.ReturnDeltaPhi3D();
    Reco_ECal = stv_tools.ReturnECal();
    Reco_EQE = stv_tools.ReturnEQE();
    Reco_Q2 = stv_tools.ReturnQ2();
    Reco_A = stv_tools.ReturnA();
    Reco_EMiss = stv_tools.ReturnEMiss();
    Reco_kMiss = stv_tools.ReturnkMiss();
    Reco_PMiss = stv_tools.ReturnPMiss();
    Reco_PMissMinus = stv_tools.ReturnPMissMinus();

    //Now deal with the backtracking

    double CandidateMuonPx = Event->pfp_true_px_->at(CandidateMuonIndex);
    double CandidateMuonPy = Event->pfp_true_py_->at(CandidateMuonIndex);
    double CandidateMuonPz = Event->pfp_true_pz_->at(CandidateMuonIndex);
    TVector3 BackTrackCandidateMuonP(CandidateMuonPx, CandidateMuonPy,
      CandidateMuonPz );
    double BackTrackCandidateMuonTrackMomentum_GeV
      = BackTrackCandidateMuonP.Mag(); // GeV
    double BackTrackCandidateMuonTrack_E_GeV = TMath::Sqrt(
      TMath::Power(BackTrackCandidateMuonTrackMomentum_GeV, 2.)
      + TMath::Power(MUON_MASS,2.) ); // GeV

    double CandidateProtonPx = Event->pfp_true_px_->at(CandidateProtonIndex);
    double CandidateProtonPy = Event->pfp_true_py_->at(CandidateProtonIndex);
    double CandidateProtonPz = Event->pfp_true_pz_->at(CandidateProtonIndex);
    TVector3 BackTrackCandidateProtonP(CandidateProtonPx, CandidateProtonPy,
      CandidateProtonPz );
    double BackTrackCandidateProtonTrackMomentum_GeV
      = BackTrackCandidateProtonP.Mag(); // GeV
    double BackTrackCandidateProtonTrack_E_GeV = TMath::Sqrt(
      TMath::Power(BackTrackCandidateProtonTrackMomentum_GeV,2.)
      + TMath::Power(PROTON_MASS,2.) ); // GeV

    stv_tools.CalculateSTVs( BackTrackCandidateMuonP, BackTrackCandidateProtonP,
      BackTrackCandidateMuonTrack_E_GeV, BackTrackCandidateProtonTrack_E_GeV,
      CalcType );

    BackTrack_Pt = stv_tools.ReturnPt();
    BackTrack_Ptx = stv_tools.ReturnPtx();
    BackTrack_Pty = stv_tools.ReturnPty();
    BackTrack_PL = stv_tools.ReturnPL();
    BackTrack_Pn = stv_tools.ReturnPn();
    BackTrack_PnPerp = stv_tools.ReturnPnPerp();
    BackTrack_PnPerpx = stv_tools.ReturnPnPerpx();
    BackTrack_PnPerpy = stv_tools.ReturnPnPerpy();
    BackTrack_PnPar = stv_tools.ReturnPnPar();
    BackTrack_DeltaAlphaT = stv_tools.ReturnDeltaAlphaT();
    BackTrack_DeltaAlpha3Dq = stv_tools.ReturnDeltaAlpha3Dq();
    BackTrack_DeltaAlpha3DMu = stv_tools.ReturnDeltaAlpha3DMu();
    BackTrack_DeltaPhiT = stv_tools.ReturnDeltaPhiT();
    BackTrack_DeltaPhi3D = stv_tools.ReturnDeltaPhi3D();
    BackTrack_ECal = stv_tools.ReturnECal();
    BackTrack_EQE = stv_tools.ReturnEQE();
    BackTrack_Q2 = stv_tools.ReturnQ2();
    BackTrack_A = stv_tools.ReturnA();
    BackTrack_EMiss = stv_tools.ReturnEMiss();
    BackTrack_kMiss = stv_tools.ReturnkMiss();
    BackTrack_PMiss = stv_tools.ReturnPMiss();
    BackTrack_PMissMinus = stv_tools.ReturnPMissMinus();
  }

}

void CC1mu1p0pi::compute_true_observables( AnalysisEvent* Event ) {

  bool isInteresting = sig_one_muon_above_thresh_
    && sig_one_proton_above_thresh_ && sig_no_pions_ && sig_no_heavy_mesons_;
  if ( isInteresting ) {
    double Muon_MCParticlePx = Event->mc_nu_daughter_px_->at(truemuonindex);
    double Muon_MCParticlePy = Event->mc_nu_daughter_px_->at(truemuonindex);
    double Muon_MCParticlePz = Event->mc_nu_daughter_px_->at(truemuonindex);
    TVector3 Muon_TVector3True( Muon_MCParticlePx, Muon_MCParticlePy,
      Muon_MCParticlePz );
    double Muon_TrueMomentum_GeV = Muon_TVector3True.Mag(); // GeV
    double Muon_TrueE_GeV = TMath::Sqrt( TMath::Power(Muon_TrueMomentum_GeV,2.)
      + TMath::Power(MUON_MASS,2.) ); // GeV

    double Proton_MCParticlePx = Event
      ->mc_nu_daughter_px_->at( trueprotonindex );
    double Proton_MCParticlePy = Event
      ->mc_nu_daughter_px_->at( trueprotonindex );
    double Proton_MCParticlePz = Event
      ->mc_nu_daughter_px_->at( trueprotonindex );
    TVector3 Proton_TVector3True( Proton_MCParticlePx, Proton_MCParticlePy,
      Proton_MCParticlePz );
    double Proton_TrueMomentum_GeV = Proton_TVector3True.Mag(); // GeV
    double Proton_TrueE_GeV = TMath::Sqrt(
      TMath::Power(Proton_TrueMomentum_GeV, 2.)
      + TMath::Power(PROTON_MASS,2.) ); // GeV

    STVTools stv_tools;
    stv_tools.CalculateSTVs( Muon_TVector3True, Proton_TVector3True,
      Muon_TrueE_GeV, Proton_TrueE_GeV, CalcType );

    True_Pt = stv_tools.ReturnPt();
    True_Ptx = stv_tools.ReturnPtx();
    True_Pty = stv_tools.ReturnPty();
    True_PL = stv_tools.ReturnPL();
    True_Pn = stv_tools.ReturnPn();
    True_PnPerp = stv_tools.ReturnPnPerp();
    True_PnPerpx = stv_tools.ReturnPnPerpx();
    True_PnPerpy = stv_tools.ReturnPnPerpy();
    True_PnPar = stv_tools.ReturnPnPar();
    True_DeltaAlphaT = stv_tools.ReturnDeltaAlphaT();
    True_DeltaAlpha3Dq = stv_tools.ReturnDeltaAlpha3Dq();
    True_DeltaAlpha3DMu = stv_tools.ReturnDeltaAlpha3DMu();
    True_DeltaPhiT = stv_tools.ReturnDeltaPhiT();
    True_DeltaPhi3D = stv_tools.ReturnDeltaPhi3D();
    True_ECal = stv_tools.ReturnECal();
    True_EQE = stv_tools.ReturnEQE();
    True_Q2 = stv_tools.ReturnQ2();
    True_A = stv_tools.ReturnA();
    True_EMiss = stv_tools.ReturnEMiss();
    True_kMiss = stv_tools.ReturnkMiss();
    True_PMiss = stv_tools.ReturnPMiss();
    True_PMissMinus = stv_tools.ReturnPMissMinus();
  }

}

int CC1mu1p0pi::categorize_event( AnalysisEvent* Event ) {

  // Real data has a bogus true neutrino PDG code that is not one of the
  // allowed values (±12, ±14, ±16)
  int abs_mc_nu_pdg = std::abs( Event->mc_nu_pdg_ );
  Event->is_mc_ = ( abs_mc_nu_pdg == ELECTRON_NEUTRINO
    || abs_mc_nu_pdg == MUON_NEUTRINO || abs_mc_nu_pdg == TAU_NEUTRINO );
  if ( !Event->is_mc_ ) {
    return kUnknown;
  }

  bool MCVertexInFV = point_inside_FV( this->true_FV(),
    Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_ );
  if ( !MCVertexInFV ) {
    return kOOFV;
  }

  bool isNC = ( Event->mc_nu_ccnc_ == NEUTRAL_CURRENT );
  // DB Currently only one NC category is supported so test first.
  // Will likely want to change this in the future
  if (isNC) return kNC;

  if ( Event->mc_nu_pdg_ == ELECTRON_NEUTRINO ) {
    return kNuECC;
  }
  if ( !(Event->mc_nu_pdg_ == MUON_NEUTRINO) ) {
    return kOther;
  }

  // Boolean which basically MC Signal selection without requesting a
  // particular number of protons (N >= 1)
  bool Is_CC1muNp0pi_Event = (sig_mc_n_threshold_proton >= 1)
    && sig_no_pions_ && sig_one_muon_above_thresh_ && sig_no_heavy_mesons_;

  if ( Is_CC1muNp0pi_Event ) {
    if (sig_mc_n_threshold_proton == 1) {
      if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCC1p0pi_CCQE; //QE
      else if ( Event->mc_nu_interaction_type_ == 10 ) {
        return kNuMuCC1p0pi_CCMEC; // MEC
      }
      else if ( Event->mc_nu_interaction_type_ == 1 ) {
        return kNuMuCC1p0pi_CCRES; // RES
      }
      else return kNuMuCCMp0pi_Other;
    } else if (sig_mc_n_threshold_proton == 2) {
      if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCC2p0pi_CCQE; //QE
      else if ( Event->mc_nu_interaction_type_ == 10 ) {
        return kNuMuCC2p0pi_CCMEC; // MEC
      }
      else if ( Event->mc_nu_interaction_type_ == 1 ) {
        return kNuMuCC2p0pi_CCRES; // RES
      }
      else return kNuMuCCMp0pi_Other;
    } else { // i.e. >=3
      if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCCMp0pi_CCQE; //QE
      else if ( Event->mc_nu_interaction_type_ == 10 ) {
        return kNuMuCCMp0pi_CCMEC; // MEC
      }
      else if ( Event->mc_nu_interaction_type_ == 1 ) {
        return kNuMuCCMp0pi_CCRES; // RES
      }
      else return kNuMuCCMp0pi_Other;
    }
  }
  else if (!sig_no_pions_) {
    return kNuMuCCNpi;
  } else if (sig_mc_n_threshold_proton == 0) {
    if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCC0p0pi_CCQE; // QE
    else if ( Event->mc_nu_interaction_type_ == 10 ) {
      return kNuMuCC0p0pi_CCMEC; // MEC
    }
    else if ( Event->mc_nu_interaction_type_ == 1 ) {
      return kNuMuCC0p0pi_CCRES; // RES
    }
    else return kNuMuCC0p0pi_Other;
  }
  return kNuMuCCOther;

}

bool CC1mu1p0pi::define_signal( AnalysisEvent* Event ) {

  // =====================================
  // DB Calculate the values which we need
  sig_mc_n_threshold_muon = 0;
  sig_mc_n_threshold_proton = 0;
  sig_mc_n_threshold_pion0 = 0;
  sig_mc_n_threshold_pionpm = 0;
  sig_mc_n_heaviermeson = 0;

  for ( size_t p = 0u; p < Event->mc_nu_daughter_pdg_->size(); ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    TVector3 MCParticle( Event->mc_nu_daughter_px_->at(p),
      Event->mc_nu_daughter_py_->at(p), Event->mc_nu_daughter_pz_->at(p) );
    double ParticleMomentum = MCParticle.Mag();

    if ( pdg == MUON && ParticleMomentum >= 0.1 ) {
      sig_mc_n_threshold_muon++;
      truemuonindex = p;
    }
    else if ( pdg == PROTON && ParticleMomentum >= 0.3 ) {
      sig_mc_n_threshold_proton++;
      trueprotonindex = p;
    }
    else if ( pdg == PI_ZERO ) {
      sig_mc_n_threshold_pion0++;
    }
    else if ( std::abs(pdg) == PI_PLUS && ParticleMomentum >= 0.07 ) {
      sig_mc_n_threshold_pionpm++;
    }
    else if ( pdg != PI_ZERO && std::abs(pdg) != PI_PLUS
      && is_meson_or_antimeson(pdg) )
    {
      sig_mc_n_heaviermeson++;
    }
  }

  // ===========================================================
  // Calculate the booleans related to the different signal cuts

  sig_truevertex_in_fv_ = point_inside_FV( this->true_FV(), Event->mc_nu_vx_,
    Event->mc_nu_vy_, Event->mc_nu_vz_ );

  sig_ccnc_= ( Event->mc_nu_ccnc_ == CHARGED_CURRENT );
  sig_is_numu_ = ( Event->mc_nu_pdg_ == MUON_NEUTRINO );
  sig_one_muon_above_thresh_ = ( sig_mc_n_threshold_muon == 1 );
  sig_one_proton_above_thresh_ = ( sig_mc_n_threshold_proton == 1 );
  sig_no_pions_ = ( (sig_mc_n_threshold_pion0 == 0)
    && (sig_mc_n_threshold_pionpm == 0) );
  sig_no_heavy_mesons_ = ( sig_mc_n_heaviermeson == 0 );

  // ====================
  // Is the event signal?

  bool isSignal = sig_truevertex_in_fv_ && sig_ccnc_ && sig_is_numu_
    && sig_one_muon_above_thresh_ && sig_one_proton_above_thresh_
    && sig_no_pions_ && sig_no_heavy_mesons_;

  return isSignal;

}

bool CC1mu1p0pi::selection( AnalysisEvent* Event) {

  // ==========================================
  // Requirement for exactly one neutrino slice
  if (Event->nslice_ == 1) sel_nslice_eq_1_ = true;

  // ======================================================================
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

  // ==============================================
  // Identify candidate muon & proton
  // Muon = the one with the highest LLR PID Score
  // Proton = the one with the lowest LLR PID Score

  int CandidateMuonIndex = -1.;
  int CandidateProtonIndex = -1.;

  if (sel_ntrack_eq_2_) {
    float first_pid_score = Event->track_llr_pid_score_
      ->at( CandidateIndex.at(0) );
    float second_pid_score = Event->track_llr_pid_score_
      ->at( CandidateIndex.at(1) );

    if (first_pid_score > second_pid_score) {
      CandidateMuonIndex = CandidateIndex.at(0);
      CandidateProtonIndex = CandidateIndex.at(1);
    } else {
      CandidateMuonIndex = CandidateIndex.at(1);
      CandidateProtonIndex = CandidateIndex.at(0);
    }

    if ( Event->pfp_reco_pdg_->at(CandidateMuonIndex) == MUON ) {
      sel_muoncandidate_tracklike_ = true;
    }
    if ( Event->pfp_reco_pdg_->at(CandidateProtonIndex) == MUON ) {
      sel_protoncandidate_tracklike_ = true;
    }
  }

  // ======================
  // Neutrino vertex in FV?

  sel_nuvertex_contained_ = point_inside_FV( this->reco_FV(),
    Event->nu_vx_,Event->nu_vy_,Event->nu_vz_ );

  // ========================================
  // Containment check on the muon and proton

  if ( sel_ntrack_eq_2_ ) {

    // Muon variables
    double MuonTrackStartX = Event->track_startx_->at( CandidateMuonIndex );
    double MuonTrackStartY = Event->track_starty_->at( CandidateMuonIndex );
    double MuonTrackStartZ = Event->track_startz_->at( CandidateMuonIndex );

    double MuonTrackEndX = Event->track_endx_->at( CandidateMuonIndex );
    double MuonTrackEndY = Event->track_endy_->at( CandidateMuonIndex );
    double MuonTrackEndZ = Event->track_endz_->at( CandidateMuonIndex );

    bool CandidateMuonTrackStartContainment = point_inside_FV(
      this->reco_FV(),MuonTrackStartX,MuonTrackStartY,MuonTrackStartZ );
    bool CandidateMuonTrackEndContainment = point_inside_FV( this->reco_FV(),
      MuonTrackEndX,MuonTrackEndY,MuonTrackEndZ );

    double CandidateMuMom = Event->track_range_mom_mu_
      ->at( CandidateMuonIndex ); // GeV/c
    double CandidateMuE_GeV = TMath::Sqrt( TMath::Power(CandidateMuMom,2.)
      + TMath::Power(MUON_MASS,2.) ); // GeV

    // If exiting muon, switch to MCS and recalculate the energy
    if ( !CandidateMuonTrackEndContainment ) {
      CandidateMuMom = Event->track_mcs_mom_mu_
        ->at(CandidateMuonIndex); // GeV/c
      CandidateMuE_GeV = TMath::Sqrt( TMath::Power(CandidateMuMom,2.)
        + TMath::Power(MUON_MASS,2.)); // GeV
    }

    // Proton variables
    double ProtonTrackStartX = Event->track_startx_->at( CandidateProtonIndex );
    double ProtonTrackStartY = Event->track_starty_->at( CandidateProtonIndex );
    double ProtonTrackStartZ = Event->track_startz_->at( CandidateProtonIndex );

    double ProtonTrackEndX = Event->track_endx_->at( CandidateProtonIndex );
    double ProtonTrackEndY = Event->track_endy_->at( CandidateProtonIndex );
    double ProtonTrackEndZ = Event->track_endz_->at( CandidateProtonIndex );

    bool CandidateProtonTrackStartContainment = point_inside_FV(
      this->reco_FV(), ProtonTrackStartX, ProtonTrackStartY,
      ProtonTrackStartZ );
    bool CandidateProtonTrackEndContainment = point_inside_FV( this->reco_FV(),
      ProtonTrackEndX, ProtonTrackEndY, ProtonTrackEndZ );

    double CandidatePKE_GeV = Event->track_kinetic_energy_p_
      ->at(CandidateProtonIndex); // GeV // Watch out, kinetic energy not energy
    double CandidatePE_GeV = CandidatePKE_GeV + PROTON_MASS; // GeV
    double CandidatePMom = TMath::Sqrt( TMath::Power(CandidatePE_GeV,2.)
      - TMath::Power(PROTON_MASS,2.) ); // GeV/c

    // STV redefinition if P_p < 0.5 where the biases have been observed
    if ( CandidatePMom < 0.5 ) {
      CandidatePMom = ( 1.-0.01*fPP->Eval(CandidatePMom) ) * CandidatePMom;
      CandidatePE_GeV = TMath::Sqrt( TMath::Power(CandidatePMom,2)
        + TMath::Power(PROTON_MASS,2.) ); // GeV/c
    }

    // Comes from /uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs
    // /ubana/ubana/myClasses/Constants.hh:L1762
    if (CandidateMuMom >= MUON_P_MIN_MOM_CUT) {
      sel_muoncandidate_above_p_thresh = true;
    }

    // Comes from /uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/
    // ubana/ubana/myClasses/Constants.hh:L1775
    if (CandidatePMom >= 0.3) sel_protoncandidate_above_p_thresh = true;

    if (CandidateMuonTrackStartContainment
      && CandidateMuonTrackEndContainment) sel_muoncandidate_contained = true;
    if (CandidateProtonTrackStartContainment
      && CandidateProtonTrackEndContainment)
    {
      sel_protoncandidate_contained = true;
    }
  }

  // ========================================================
  // Check that the muon momentum estimators agree within 25%

  if (CandidateMuonIndex != -1) {
    double MuonMom_MCS = Event->track_mcs_mom_mu_->at(CandidateMuonIndex);
    double MuonMom_Range = Event->track_range_mom_mu_->at(CandidateMuonIndex);

    double Reso =  TMath::Abs(MuonMom_MCS - MuonMom_Range) / MuonMom_Range;
    sel_muon_momentum_quality = (Reso <= 0.25);
  }

  // ========================
  // Check for flipped tracks

  if (CandidateMuonIndex != -1 && CandidateProtonIndex != -1) {
    TVector3 VertexLocation(Event->nu_vx_,Event->nu_vy_,Event->nu_vz_);

    TVector3 Candidate_MuonTrack_Start(
      Event->track_startx_->at(CandidateMuonIndex),
      Event->track_starty_->at(CandidateMuonIndex),
      Event->track_startz_->at(CandidateMuonIndex)
    );

    TVector3 Candidate_MuonTrack_End(
      Event->track_endx_->at(CandidateMuonIndex),
      Event->track_endy_->at(CandidateMuonIndex),
      Event->track_endz_->at(CandidateMuonIndex)
    );

    TVector3 Candidate_ProtonTrack_Start(
      Event->track_startx_->at(CandidateProtonIndex),
      Event->track_starty_->at(CandidateProtonIndex),
      Event->track_startz_->at(CandidateProtonIndex)
    );

    TVector3 Candidate_ProtonTrack_End(
      Event->track_endx_->at(CandidateProtonIndex),
      Event->track_endy_->at(CandidateProtonIndex),
      Event->track_endz_->at(CandidateProtonIndex)
    );

    double Vertex_MuonTrackStart_Mag
      = (VertexLocation - Candidate_MuonTrack_Start).Mag();
    double Vertex_MuonTrackEnd_Mag
      = (VertexLocation - Candidate_MuonTrack_End).Mag();
    double Vertex_ProtonTrackStart_Mag
      = (VertexLocation - Candidate_ProtonTrack_Start).Mag();
    double Vertex_ProtonTrackEnd_Mag
      = (VertexLocation - Candidate_ProtonTrack_End).Mag();

    double MuonTrackStart_to_ProtonTrackStart_Mag
      = (Candidate_MuonTrack_Start - Candidate_ProtonTrack_Start).Mag();
    double MuonTrackEnd_to_ProtonTrackEnd_Mag
      = (Candidate_MuonTrack_End - Candidate_ProtonTrack_End).Mag();

    if ( !( (Vertex_MuonTrackStart_Mag > Vertex_MuonTrackEnd_Mag)
      || (Vertex_ProtonTrackStart_Mag > Vertex_ProtonTrackEnd_Mag)
      || (MuonTrackStart_to_ProtonTrackStart_Mag
        > MuonTrackEnd_to_ProtonTrackEnd_Mag) ) )
    {
      sel_no_flipped_tracks_ = true;
    }
  }

  // ===========================================
  // Check Proton Candidate's LLH to be a proton

  if ( CandidateMuonIndex != -1 && CandidateProtonIndex != -1 ) {
    double Candidate_Proton_LLR
      = Event->track_llr_pid_score_->at(CandidateProtonIndex);
    sel_proton_cand_passed_LLRCut = (Candidate_Proton_LLR < 0.05);
  }

  // ====================
  // Apply kinematic cuts

  if ( CandidateMuonIndex != -1 ) {
    // This is essentially duplicate of sel_muoncandidate_above_p_thresh cut
    double MuonMom = Event->track_range_mom_mu_->at(CandidateMuonIndex);
    if ( ! (MuonMom < 0.1 || MuonMom > 1.2) ) {
      sel_muon_momentum_in_range = true;
    }

    double MuonCosTheta = std::cos(
      Event->track_theta_->at(CandidateMuonIndex) );
    if ( ! (MuonCosTheta < -1. || MuonCosTheta > 1.) ) {
      sel_muon_costheta_in_range = true;
    }

    double MuonPhi = Event->track_phi_->at(CandidateMuonIndex)
      * 180./ TMath::Pi();
    if ( ! (MuonPhi < -180. || MuonPhi > 180.) ) {
      sel_muon_phi_in_range = true;
    }
  }

  if ( CandidateProtonIndex != -1 ) {
    // This is essentially duplicate of sel_protoncandidate_above_p_thresh cut
    double ProtonMomentum = TMath::Sqrt(
      TMath::Power(Event->track_kinetic_energy_p_->at(CandidateProtonIndex)
      + PROTON_MASS,2.) - TMath::Power(PROTON_MASS,2.)
    );

    if ( ! (ProtonMomentum < 0.3 || ProtonMomentum > 1.) ) {
      sel_proton_momentum_in_range = true;
    }

    double ProtonCosTheta = std::cos(
      Event->track_theta_->at(CandidateProtonIndex)
    );

    if ( ! (ProtonCosTheta < -1. || ProtonCosTheta > 1.) ) {
      sel_proton_costheta_in_range = true;
    }

    double ProtonPhi = Event->track_phi_->at(CandidateProtonIndex)
      * 180./ TMath::Pi();
    if ( ! (ProtonPhi < -180. || ProtonPhi > 180.) ) {
      sel_proton_phi_in_range = true;
    }
  }


  // ================================
  // Does everything pass selection?
  bool Passed = sel_nslice_eq_1_ && sel_nshower_eq_0_ && sel_ntrack_eq_2_
    && sel_muoncandidate_tracklike_ && sel_protoncandidate_tracklike_
    && sel_nuvertex_contained_ && sel_muoncandidate_above_p_thresh
    && sel_protoncandidate_above_p_thresh && sel_muoncandidate_contained
    && sel_protoncandidate_contained && sel_muon_momentum_quality
    && sel_no_flipped_tracks_ && sel_proton_cand_passed_LLRCut
    && sel_muon_momentum_in_range && sel_muon_costheta_in_range
    && sel_muon_phi_in_range && sel_proton_momentum_in_range
    && sel_proton_costheta_in_range && sel_proton_phi_in_range;

  return Passed;
}

void CC1mu1p0pi::define_output_branches() {

  set_branch( &sel_nslice_eq_1_, "nslice_eq_1" );
  set_branch( &sel_nshower_eq_0_, "nshower_eq_0_" );
  set_branch( &sel_ntrack_eq_2_, "ntrack_eq_2_" );
  set_branch( &sel_muoncandidate_tracklike_, "muoncandidate_tracklike" );
  set_branch( &sel_protoncandidate_tracklike_, "protoncandidate_tracklike" );
  set_branch( &sel_nuvertex_contained_, "nuvertex_contained_" );

  set_branch( &sel_muoncandidate_above_p_thresh,
    "muoncandidate_above_p_thresh" );

  set_branch( &sel_protoncandidate_above_p_thresh,
    "protoncandidate_above_p_thresh" );

  set_branch( &sel_muoncandidate_contained, "muoncandidate_contained" );
  set_branch( &sel_protoncandidate_contained, "protoncandidate_contained" );
  set_branch( &sel_muon_momentum_quality, "sel_muon_momentum_quality" );
  set_branch( &sel_no_flipped_tracks_, "sel_no_flipped_tracks" );
  set_branch( &sel_proton_cand_passed_LLRCut, "sel_proton_cand_passed_LLRCut" );
  set_branch( &sel_muon_momentum_in_range, "sel_muon_momentum_in_range" );
  set_branch( &sel_muon_costheta_in_range, "sel_muon_costheta_in_range" );
  set_branch( &sel_muon_phi_in_range, "sel_muon_phi_in_range" );
  set_branch( &sel_proton_momentum_in_range, "sel_proton_momentum_in_range" );
  set_branch( &sel_proton_costheta_in_range, "sel_proton_costheta_in_range" );
  set_branch( &sel_proton_phi_in_range, "sel_proton_phi_in_range" );
  set_branch( &CandidateMuonIndex, "CandidateMuonIndex" );
  set_branch( &CandidateProtonIndex, "CandidateProtonIndex" );

  set_branch( &sig_truevertex_in_fv_, "sig_truevertex_in_fv" );
  set_branch( &sig_ccnc_, "sig_ccnc" );
  set_branch( &sig_is_numu_, "sig_is_numu" );
  set_branch( &sig_one_muon_above_thresh_, "sig_one_muon_above_thresh" );
  set_branch( &sig_one_proton_above_thresh_, "sig_one_proton_above_thresh" );
  set_branch( &sig_no_pions_, "sig_no_pions" );
  set_branch( &sig_no_heavy_mesons_, "sig_no_heavy_mesons" );
  set_branch( &sig_mc_n_threshold_muon, "mc_n_threshold_muon" );
  set_branch( &sig_mc_n_threshold_proton, "mc_n_threshold_proton" );
  set_branch( &sig_mc_n_threshold_pion0, "mc_n_threshold_pion0" );
  set_branch( &sig_mc_n_threshold_pionpm, "mc_n_threshold_pionpm" );
  set_branch( &sig_mc_n_heaviermeson, "mc_n_heaviermeson" );
  set_branch( &truemuonindex, "truemuonindex" );
  set_branch( &trueprotonindex, "trueprotonindex" );

  set_branch( &Reco_Pt, "Reco_Pt" );
  set_branch( &Reco_Ptx, "Reco_Ptx" );
  set_branch( &Reco_Pty, "Reco_Pty" );
  set_branch( &Reco_PL, "Reco_PL" );
  set_branch( &Reco_Pn, "Reco_Pn" );
  set_branch( &Reco_PnPerp, "Reco_PnPerp" );
  set_branch( &Reco_PnPerpx, "Reco_PnPerpx" );
  set_branch( &Reco_PnPerpy, "Reco_PnPerpy" );
  set_branch( &Reco_PnPar, "Reco_PnPar" );
  set_branch( &Reco_DeltaAlphaT, "Reco_DeltaAlphaT" );
  set_branch( &Reco_DeltaAlpha3Dq, "Reco_DeltaAlpha3Dq" );
  set_branch( &Reco_DeltaAlpha3DMu, "Reco_DeltaAlpha3DMu" );
  set_branch( &Reco_DeltaPhiT, "Reco_DeltaPhiT" );
  set_branch( &Reco_DeltaPhi3D, "Reco_DeltaPhi3D" );
  set_branch( &Reco_ECal, "Reco_ECal" );
  set_branch( &Reco_EQE, "Reco_EQE" );
  set_branch( &Reco_Q2, "Reco_Q2" );
  set_branch( &Reco_A, "Reco_A" );
  set_branch( &Reco_EMiss, "Reco_EMiss" );
  set_branch( &Reco_kMiss, "Reco_kMiss" );
  set_branch( &Reco_PMiss, "Reco_PMiss" );
  set_branch( &Reco_PMissMinus, "Reco_PMissMinus" );
  set_branch( &BackTrack_Pt, "BackTrack_Pt" );
  set_branch( &BackTrack_Ptx, "BackTrack_Ptx" );
  set_branch( &BackTrack_Pty, "BackTrack_Pty" );
  set_branch( &BackTrack_PL, "BackTrack_PL" );
  set_branch( &BackTrack_Pn, "BackTrack_Pn" );
  set_branch( &BackTrack_PnPerp, "BackTrack_PnPerp" );
  set_branch( &BackTrack_PnPerpx, "BackTrack_PnPerpx" );
  set_branch( &BackTrack_PnPerpy, "BackTrack_PnPerpy" );
  set_branch( &BackTrack_PnPar, "BackTrack_PnPar" );
  set_branch( &BackTrack_DeltaAlphaT, "BackTrack_DeltaAlphaT" );
  set_branch( &BackTrack_DeltaAlpha3Dq, "BackTrack_DeltaAlpha3Dq" );
  set_branch( &BackTrack_DeltaAlpha3DMu, "BackTrack_DeltaAlpha3DMu" );
  set_branch( &BackTrack_DeltaPhiT, "BackTrack_DeltaPhiT" );
  set_branch( &BackTrack_DeltaPhi3D, "BackTrack_DeltaPhi3D" );
  set_branch( &BackTrack_ECal, "BackTrack_ECal" );
  set_branch( &BackTrack_EQE, "BackTrack_EQE" );
  set_branch( &BackTrack_Q2, "BackTrack_Q2" );
  set_branch( &BackTrack_A, "BackTrack_A" );
  set_branch( &BackTrack_EMiss, "BackTrack_EMiss" );
  set_branch( &BackTrack_kMiss, "BackTrack_kMiss" );
  set_branch( &BackTrack_PMiss, "BackTrack_PMiss" );
  set_branch( &BackTrack_PMissMinus, "BackTrack_PMissMinus" );
  set_branch( &True_Pt, "True_Pt" );
  set_branch( &True_Ptx, "True_Ptx" );
  set_branch( &True_Pty, "True_Pty" );
  set_branch( &True_PL, "True_PL" );
  set_branch( &True_Pn, "True_Pn" );
  set_branch( &True_PnPerp, "True_PnPerp" );
  set_branch( &True_PnPerpx, "True_PnPerpx" );
  set_branch( &True_PnPerpy, "True_PnPerpy" );
  set_branch( &True_PnPar, "True_PnPar" );
  set_branch( &True_DeltaAlphaT, "True_DeltaAlphaT" );
  set_branch( &True_DeltaAlpha3Dq, "True_DeltaAlpha3Dq" );
  set_branch( &True_DeltaAlpha3DMu, "True_DeltaAlpha3DMu" );
  set_branch( &True_DeltaPhiT, "True_DeltaPhiT" );
  set_branch( &True_DeltaPhi3D, "True_DeltaPhi3D" );
  set_branch( &True_ECal, "True_ECal" );
  set_branch( &True_EQE, "True_EQE" );
  set_branch( &True_Q2, "True_Q2" );
  set_branch( &True_A, "True_A" );
  set_branch( &True_EMiss, "True_EMiss" );
  set_branch( &True_kMiss, "True_kMiss" );
  set_branch( &True_PMiss, "True_PMiss" );
  set_branch( &True_PMissMinus, "True_PMissMinus" );
}

void CC1mu1p0pi::reset() {

  sel_nslice_eq_1_ = false;
  sel_nshower_eq_0_ = false;
  sel_ntrack_eq_2_ = false;
  sel_muoncandidate_tracklike_ = false;
  sel_protoncandidate_tracklike_= false;
  sel_nuvertex_contained_ = false;
  sel_muoncandidate_above_p_thresh = false;
  sel_protoncandidate_above_p_thresh = false;
  sel_muoncandidate_contained = false;
  sel_protoncandidate_contained = false;
  sel_muon_momentum_quality = false;
  sel_no_flipped_tracks_ = false;
  sel_proton_cand_passed_LLRCut = false;
  sel_muon_momentum_in_range = false;
  sel_muon_costheta_in_range = false;
  sel_muon_phi_in_range = false;
  sel_proton_momentum_in_range = false;
  sel_proton_costheta_in_range = false;
  sel_proton_phi_in_range = false;
  CandidateMuonIndex = BOGUS_INDEX;
  CandidateProtonIndex = BOGUS_INDEX;

  sig_truevertex_in_fv_ = false;
  sig_ccnc_ = false;
  sig_is_numu_ = false;
  sig_one_muon_above_thresh_ = false;
  sig_one_proton_above_thresh_ = false;
  sig_no_pions_ = false;
  sig_no_heavy_mesons_ = false;
  sig_mc_n_threshold_muon = BOGUS_INDEX;
  sig_mc_n_threshold_proton = BOGUS_INDEX;
  sig_mc_n_threshold_pion0 = BOGUS_INDEX;
  sig_mc_n_threshold_pionpm = BOGUS_INDEX;
  sig_mc_n_heaviermeson = BOGUS_INDEX;
  truemuonindex = BOGUS_INDEX;
  trueprotonindex = BOGUS_INDEX;

  Reco_Pt = BOGUS;
  Reco_Ptx = BOGUS;
  Reco_Pty = BOGUS;
  Reco_PL = BOGUS;
  Reco_Pn = BOGUS;
  Reco_PnPerp = BOGUS;
  Reco_PnPerpx = BOGUS;
  Reco_PnPerpy = BOGUS;
  Reco_PnPar = BOGUS;
  Reco_DeltaAlphaT = BOGUS;
  Reco_DeltaAlpha3Dq = BOGUS;
  Reco_DeltaAlpha3DMu = BOGUS;
  Reco_DeltaPhiT = BOGUS;
  Reco_DeltaPhi3D = BOGUS;
  Reco_ECal = BOGUS;
  Reco_EQE = BOGUS;
  Reco_Q2 = BOGUS;
  Reco_A = BOGUS;
  Reco_EMiss = BOGUS;
  Reco_kMiss = BOGUS;
  Reco_PMiss = BOGUS;
  Reco_PMissMinus = BOGUS;
  BackTrack_Pt = BOGUS;
  BackTrack_Ptx = BOGUS;
  BackTrack_Pty = BOGUS;
  BackTrack_PL = BOGUS;
  BackTrack_Pn = BOGUS;
  BackTrack_PnPerp = BOGUS;
  BackTrack_PnPerpx = BOGUS;
  BackTrack_PnPerpy = BOGUS;
  BackTrack_PnPar = BOGUS;
  BackTrack_DeltaAlphaT = BOGUS;
  BackTrack_DeltaAlpha3Dq = BOGUS;
  BackTrack_DeltaAlpha3DMu = BOGUS;
  BackTrack_DeltaPhiT = BOGUS;
  BackTrack_DeltaPhi3D = BOGUS;
  BackTrack_ECal = BOGUS;
  BackTrack_EQE = BOGUS;
  BackTrack_Q2 = BOGUS;
  BackTrack_A = BOGUS;
  BackTrack_EMiss = BOGUS;
  BackTrack_kMiss = BOGUS;
  BackTrack_PMiss = BOGUS;
  BackTrack_PMissMinus = BOGUS;
  True_Pt = BOGUS;
  True_Ptx = BOGUS;
  True_Pty = BOGUS;
  True_PL = BOGUS;
  True_Pn = BOGUS;
  True_PnPerp = BOGUS;
  True_PnPerpx = BOGUS;
  True_PnPerpy = BOGUS;
  True_PnPar = BOGUS;
  True_DeltaAlphaT = BOGUS;
  True_DeltaAlpha3Dq = BOGUS;
  True_DeltaAlpha3DMu = BOGUS;
  True_DeltaPhiT = BOGUS;
  True_DeltaPhi3D = BOGUS;
  True_ECal = BOGUS;
  True_EQE = BOGUS;
  True_Q2 = BOGUS;
  True_A = BOGUS;
  True_EMiss = BOGUS;
  True_kMiss = BOGUS;
  True_PMiss = BOGUS;
  True_PMissMinus = BOGUS;
}

void CC1mu1p0pi::define_category_map() {
  // Use the shared category map for 1p/2p/Np/Xp
  categ_map_ = CC1muXp_MAP;
}
