// XSecAnalyzer includes
#include "XSecAnalyzer/TreeUtils.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/FiducialVolume.hh"

#include "XSecAnalyzer/Selections/CC1mu2p0pi.hh"
#include "XSecAnalyzer/Selections/EventCategoriesXp.hh"

CC1mu2p0pi::CC1mu2p0pi() : SelectionBase( "CC1mu2p0pi" ) {
  CalcType = kOpt1;
}

void CC1mu2p0pi::define_constants() {
  this->define_true_FV( 10., 246.35, -106.5, 106.5, 10., 1026.8 );
  this->define_reco_FV( 10., 246.35, -106.5, 106.5, 10., 1026.8 );
}

void CC1mu2p0pi::compute_reco_observables( AnalysisEvent* Event ) {

  if ( LeadingProtonIndex != BOGUS_INDEX && RecoilProtonIndex != BOGUS_INDEX
    && muon_candidate_idx_ != BOGUS_INDEX)
  {
    TVector3 NuVertex(Event->nu_vx_,Event->nu_vy_,Event->nu_vz_);

    //=====================================================================

    float MuonMomentum = Event->track_range_mom_mu_->at( muon_candidate_idx_ );
    float MuonEnergy = std::sqrt( MuonMomentum*MuonMomentum
      + MUON_MASS*MUON_MASS);

    float MuonTrackStartDistance = Event->track_start_distance_
      ->at( muon_candidate_idx_ );

    TVector3 MuonTrackEndDistanceVector(
      Event->track_endx_->at(muon_candidate_idx_),
      Event->track_endy_->at(muon_candidate_idx_),
      Event->track_endz_->at(muon_candidate_idx_)
    );

    MuonTrackEndDistanceVector -= NuVertex;
    float MuonTrackEndDistance = MuonTrackEndDistanceVector.Mag();

    float Muon_X = Event->track_dirx_->at(muon_candidate_idx_);
    float Muon_Y = Event->track_diry_->at(muon_candidate_idx_);
    float Muon_Z = Event->track_dirz_->at(muon_candidate_idx_);
    TVector3 MuonMomentumVector = TVector3(Muon_X, Muon_Y, Muon_Z);
    MuonMomentumVector = MuonMomentumVector.Unit() * MuonMomentum;

    if (MuonTrackStartDistance > MuonTrackEndDistance) {
      MuonMomentumVector *= -1.0;
    }

    //================================================================

    float LeadingProtonMomentum = std::sqrt(
      std::pow( Event->track_kinetic_energy_p_->at(LeadingProtonIndex)
      + PROTON_MASS, 2 ) - std::pow( PROTON_MASS, 2 )
    );

    float LeadingProtonEnergy = Event->track_kinetic_energy_p_
      ->at( LeadingProtonIndex ) + PROTON_MASS;

    float LeadingProtonTrackStartDistance
      = Event->track_start_distance_->at(LeadingProtonIndex);
    TVector3 LeadingProtonTrackEndDistanceVector(
      Event->track_endx_->at(LeadingProtonIndex),
      Event->track_endy_->at(LeadingProtonIndex),
      Event->track_endz_->at(LeadingProtonIndex)
    );
    LeadingProtonTrackEndDistanceVector -= NuVertex;
    float LeadingProtonTrackEndDistance
      = LeadingProtonTrackEndDistanceVector.Mag();

    float LeadingProton_X = Event->track_dirx_->at( LeadingProtonIndex );
    float LeadingProton_Y = Event->track_diry_->at( LeadingProtonIndex );
    float LeadingProton_Z = Event->track_dirz_->at( LeadingProtonIndex );
    TVector3 LeadingProtonMomentumVector = TVector3( LeadingProton_X,
      LeadingProton_Y, LeadingProton_Z );
    LeadingProtonMomentumVector
      = LeadingProtonMomentumVector.Unit() * LeadingProtonMomentum;

    if ( LeadingProtonTrackStartDistance > LeadingProtonTrackEndDistance ) {
      LeadingProtonMomentumVector *= -1.0;
    }

    // ============================================================

    float RecoilProtonMomentum = std::sqrt(
      std::pow( Event->track_kinetic_energy_p_->at(RecoilProtonIndex)
        + PROTON_MASS, 2 ) - std::pow( PROTON_MASS, 2 )
    );
    float RecoilProtonEnergy
      = Event->track_kinetic_energy_p_->at(RecoilProtonIndex) + PROTON_MASS;

    float RecoilProtonTrackStartDistance
      = Event->track_start_distance_->at(RecoilProtonIndex);
    TVector3 RecoilProtonTrackEndDistanceVector(
      Event->track_endx_->at(RecoilProtonIndex),
      Event->track_endy_->at(RecoilProtonIndex),
      Event->track_endz_->at(RecoilProtonIndex)
    );

    RecoilProtonTrackEndDistanceVector -= NuVertex;
    float RecoilProtonTrackEndDistance
      = RecoilProtonTrackEndDistanceVector.Mag();

    float RecoilProton_X = Event->track_dirx_->at( RecoilProtonIndex );
    float RecoilProton_Y = Event->track_diry_->at( RecoilProtonIndex );
    float RecoilProton_Z = Event->track_dirz_->at( RecoilProtonIndex );
    TVector3 RecoilProtonMomentumVector = TVector3(
      RecoilProton_X, RecoilProton_Y, RecoilProton_Z );
    RecoilProtonMomentumVector = RecoilProtonMomentumVector.Unit()
      * RecoilProtonMomentum;

    if ( RecoilProtonTrackStartDistance > RecoilProtonTrackEndDistance ) {
      RecoilProtonMomentumVector *= -1.0;
    }

    // ======================================================

    TVector3 ProtonSummedMomentumVector = LeadingProtonMomentumVector
      + RecoilProtonMomentumVector;
    float ProtonSummedEnergy = LeadingProtonEnergy + RecoilProtonEnergy;

    Reco_CosPlPr = LeadingProtonMomentumVector
      .Angle( RecoilProtonMomentumVector );
    Reco_CosMuPsum = MuonMomentumVector
      .Angle( ProtonSummedMomentumVector );

    STVTools stv_tools;
    stv_tools.CalculateSTVs( MuonMomentumVector, ProtonSummedMomentumVector,
      MuonEnergy, ProtonSummedEnergy );

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
  }

}

void CC1mu2p0pi::compute_true_observables( AnalysisEvent* Event ) {
  if ( sig_two_protons_above_thresh_ && sig_one_muon_above_thres_
    && sig_no_pions_)
  {
    double Muon_MCParticlePx = Event->mc_nu_daughter_px_->at(TrueMuonIndex);
    double Muon_MCParticlePy = Event->mc_nu_daughter_px_->at(TrueMuonIndex);
    double Muon_MCParticlePz = Event->mc_nu_daughter_px_->at(TrueMuonIndex);
    TVector3 Muon_TVector3True( Muon_MCParticlePx, Muon_MCParticlePy,
      Muon_MCParticlePz );
    double Muon_TrueMomentum_GeV = Muon_TVector3True.Mag(); // GeV
    double Muon_TrueE_GeV = TMath::Sqrt(
     TMath::Power(Muon_TrueMomentum_GeV,2.) // GeV
     + TMath::Power(MUON_MASS,2.)
    );

    double LeadingProton_MCParticlePx
      = Event->mc_nu_daughter_px_->at(TrueLeadingProtonIndex);
    double LeadingProton_MCParticlePy
      = Event->mc_nu_daughter_px_->at(TrueLeadingProtonIndex);
    double LeadingProton_MCParticlePz
      = Event->mc_nu_daughter_px_->at(TrueLeadingProtonIndex);

    TVector3 LeadingProton_TVector3True(
      LeadingProton_MCParticlePx,
      LeadingProton_MCParticlePy,
      LeadingProton_MCParticlePz
    );

    double LeadingProton_TrueMomentum_GeV
      = LeadingProton_TVector3True.Mag(); // GeV
    double LeadingProton_TrueE_GeV = // GeV
      TMath::Sqrt( TMath::Power(LeadingProton_TrueMomentum_GeV,2.)
      + TMath::Power(PROTON_MASS,2.)
    );

    double RecoilProton_MCParticlePx
      = Event->mc_nu_daughter_px_->at(TrueRecoilProtonIndex);
    double RecoilProton_MCParticlePy
      = Event->mc_nu_daughter_px_->at(TrueRecoilProtonIndex);
    double RecoilProton_MCParticlePz
      = Event->mc_nu_daughter_px_->at(TrueRecoilProtonIndex);
    TVector3 RecoilProton_TVector3True( RecoilProton_MCParticlePx,
      RecoilProton_MCParticlePy,
      RecoilProton_MCParticlePz
    );
    double RecoilProton_TrueMomentum_GeV
      = RecoilProton_TVector3True.Mag(); // GeV
    double RecoilProton_TrueE_GeV = TMath::Sqrt(
      TMath::Power(RecoilProton_TrueMomentum_GeV, 2.)
      + TMath::Power(PROTON_MASS,2.)
    ); // GeV

    TVector3 ProtonSum_TVector3True
      = LeadingProton_TVector3True+RecoilProton_TVector3True;

    double ProtonSum_TrueE_GeV
      = LeadingProton_TrueE_GeV+RecoilProton_TrueE_GeV;

    STVTools stv_tools;
    stv_tools.CalculateSTVs( Muon_TVector3True, ProtonSum_TVector3True,
      Muon_TrueE_GeV, ProtonSum_TrueE_GeV,CalcType );

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

int CC1mu2p0pi::categorize_event( AnalysisEvent* Event ) {
  // Real data has a bogus true neutrino PDG code that is not one of the
  // allowed values (±12, ±14, ±16)
  int abs_mc_nu_pdg = std::abs( Event->mc_nu_pdg_ );
  Event->is_mc_ = ( abs_mc_nu_pdg == ELECTRON_NEUTRINO
    || abs_mc_nu_pdg == MUON_NEUTRINO || abs_mc_nu_pdg == TAU_NEUTRINO );
  if ( !Event->is_mc_ ) {
    return kUnknown;
  }

  bool MCVertexInFV = point_inside_FV( this->true_FV(), Event->mc_nu_vx_,
    Event->mc_nu_vy_, Event->mc_nu_vz_ );
  if ( !MCVertexInFV ) {
    return kOOFV;
  }

  bool isNC = ( Event->mc_nu_ccnc_ == NEUTRAL_CURRENT );
  // DB Currently only one NC category is supported so test first.
  // Will likely want to change this in the future
  if ( isNC ) return kNC;

  if (Event->mc_nu_pdg_ == ELECTRON_NEUTRINO) {
    return kNuECC;
  }
  if (!(Event->mc_nu_pdg_ == MUON_NEUTRINO)) {
    return kOther;
  }

  // Boolean which basically MC Signal selection without requesting a
  // particular number of protons (N >= 1)
  bool Is_CC1muNp0pi_Event = (sig_mc_n_threshold_proton >= 1) && sig_no_pions_
    && sig_one_muon_above_thres_;

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
    if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCC0p0pi_CCQE; //QE
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

// Taken from https://github.com/ssfehlberg/CC2p-Event-Selection/blob
// /9492ff121a2eb884f464e1c166d067f217a04900/PeLEE_ntuples
// /mc_efficiency.C#L109-L112
bool CC1mu2p0pi::define_signal( AnalysisEvent* Event ) {

  // ============================================================
  // DB Calculate the values which we need
  sig_mc_n_threshold_muon = 0;
  sig_mc_n_threshold_proton = 0;
  sig_mc_n_threshold_pion0 = 0;
  sig_mc_n_threshold_pionpm = 0;

  std::vector<double> TrueProtonMomenta = std::vector<double>();
  std::vector<int> TrueProtonIndices = std::vector<int>();

  int counter = 0;
  for ( size_t p = 0u; p < Event->mc_nu_daughter_pdg_->size(); ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    float energy = Event->mc_nu_daughter_energy_->at( p );
    if ( std::abs(pdg) == MUON) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(MUON_MASS, 2) );
      if ( mom > 0.1 && mom < 1.2){
	sig_mc_n_threshold_muon++;

	TrueMuonIndex = p;
      }
    } else if (std::abs(pdg) == PROTON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PROTON_MASS, 2) );
      if ( mom > 0.3 && mom < 1.0){
	sig_mc_n_threshold_proton++;

	TrueProtonMomenta.push_back(mom);
	TrueProtonIndices.push_back(p);
	counter++;
      }
    } else if ( pdg == PI_ZERO ) {
      sig_mc_n_threshold_pion0++;

    } else if (std::abs(pdg) == PI_PLUS ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PI_PLUS_MASS, 2) );
      if ( mom > 0.065 ) {
	sig_mc_n_threshold_pionpm++;
      }
    }
  }

  if (TrueProtonMomenta.size() == 2) {
    if (TrueProtonMomenta[0] > TrueProtonMomenta[1]) {
      TrueLeadingProtonIndex = TrueProtonIndices[0];
      TrueRecoilProtonIndex = TrueProtonIndices[1];
    } else {
      TrueLeadingProtonIndex = TrueProtonIndices[1];
      TrueRecoilProtonIndex = TrueProtonIndices[0];
    }
  }

  // ===========================================================
  // Calculate the booleans related to the different signal cuts

  // DB Discussions (https://microboone.slack.com/archives/C05TCS17EHL
  // /p1695988699125549) - Afro says we should not be using Space Charge
  // Effects (SCE) in the true FV definition
  // Currently included for validation purposes
  // sig_truevertex_in_fv_ = point_inside_FV( this->true_FV(),
  //   Event->mc_nu_sce_vx_, Event->mc_nu_sce_vy_, Event->mc_nu_sce_vz_ );
  sig_truevertex_in_fv_ = point_inside_FV( this->true_FV(), Event->mc_nu_vx_,
    Event->mc_nu_vy_, Event->mc_nu_vz_ );

  sig_ccnc_ = (Event->mc_nu_ccnc_ == CHARGED_CURRENT);
  sig_is_numu_ = (Event->mc_nu_pdg_ == MUON_NEUTRINO);
  sig_two_protons_above_thresh_ = (sig_mc_n_threshold_proton == 2);
  sig_one_muon_above_thres_ = (sig_mc_n_threshold_muon == 1);
  sig_no_pions_ = ((sig_mc_n_threshold_pion0 == 0)
    && (sig_mc_n_threshold_pionpm == 0));

  // ====================
  // Is the event signal?

  bool IsSignal = sig_ccnc_ && sig_is_numu_ && sig_two_protons_above_thresh_ &&
    sig_one_muon_above_thres_ && sig_no_pions_ && sig_truevertex_in_fv_;
  return IsSignal;
}

bool CC1mu2p0pi::selection( AnalysisEvent* Event ) {

  FiducialVolume FV_noBorder;
  FV_noBorder.X_Min = 0.;
  FV_noBorder.X_Max = 256.35;
  FV_noBorder.Y_Min = -116.5;
  FV_noBorder.Y_Max = 116.5;
  FV_noBorder.Z_Min = 0.;
  FV_noBorder.Z_Max = 1036.8;

  // =============
  // Vertex in FV?
  float_t x = Event->nu_vx_;
  float_t y = Event->nu_vy_;
  float_t z = Event->nu_vz_;

  sel_reco_vertex_in_FV_ = point_inside_FV( this->reco_FV(), x, y, z );

  // =======================================================================
  // DB Samantha's analysis explicitly cuts out events with num_candidates!=1
  // (n_muons)
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

      // Gets overwritten if multiple muon candidates, but that's fine because
      // we require exactly one muon candidate
      chosen_index = p;
    }
  }

  if ( n_muons == 1u ) {
    sel_has_muon_candidate_ = true;
    muon_candidate_idx_ = chosen_index;
  } else {
    muon_candidate_idx_ = BOGUS_INDEX;
  }

  // ===============================================
  // DB Does the event pass the numuCC0pi selection?
  sel_nu_mu_cc_ = sel_reco_vertex_in_FV_ && sel_has_muon_candidate_;

  // ==========================
  // DB Require exactly 3 PFP's
  if (Event->num_pf_particles_ == 3) {
    sel_npfps_eq_3 = true;
  }

  // ============================================
  // DB Require 3 tracks (track_score > 0.8 [MUON_TRACK_SCORE_CUT]) whose
  // "vertex distance attachment is less than 4 cm"
  int nTracks = 0;

  TVector3 nu_vtx_reco(Event->nu_vx_,Event->nu_vy_,Event->nu_vz_);
  for(int i = 0; i < Event->num_pf_particles_; i ++){
    float track_score = Event->pfp_track_score_->at(i);
    float track_distance = Event->track_start_distance_->at(i);
    float track_pid = Event->track_llr_pid_score_->at(i);

    // leading track end reco
    TVector3 track_end( Event->track_endx_->at(i),
      Event->track_endy_->at(i), Event->track_endz_->at(i) );
    track_end -= nu_vtx_reco;

    double track_end_distance = track_end.Mag();

    // DB I think this cut has a different logic than Stephen's
    if ( track_score >= MUON_TRACK_SCORE_CUT
      && (track_distance <= MUON_VTX_DISTANCE_CUT
        || track_end_distance <= MUON_VTX_DISTANCE_CUT) )
    {
      nTracks++;
    }
  }

  if (nTracks == 3) {
    sel_ntracks_eq_3 = true;
  }

  // ==============================================
  // DB Make sure there's two Protons and one Muon.
  // Muon selection already performed in apply_numu_CC_selection()
  // -> sel_has_muon_candidate_ == true if only 1 muon candidate
  int nProtons = 0;
  for(int i = 0; i < Event->num_pf_particles_; i ++){
    float track_pid = Event->track_llr_pid_score_->at(i);
    if(track_pid < DEFAULT_PROTON_PID_CUT && track_pid < 1 && track_pid > -1) {
      nProtons += 1;
    }
  }

  //DB Discussion (https://microboone.slack.com/archives
  // /D04A8CUB1EW/p1691519899965369)
  // - Logic fault in Samantha's original code: https://github.com/ssfehlberg
  // /CC2p-Event-Selection/blob/9492ff121a2eb884f464e1c166d067f217a04900
  // /PeLEE_ntuples/twoproton_pelee_bnb.C#L138
  if (sel_has_muon_candidate_ && nProtons == 2) {
    sel_correctparticles = true;
  }

  // =======================================================================
  // DB Now check that all 3PFPs are contained (start and end) within the FV

  // DB Initially samantha precalculates which index corresponds to the muon,
  // leading p proton and recoil proton But then applies the same containment
  // cut to all three. The function used is:
  // https://github.com/ssfehlberg/CC2p-Event-Selection/blob
  // /9492ff121a2eb884f464e1c166d067f217a04900/PeLEE_ntuples
  // /helper_funcs.h#L18-L24
  // Where the start of the vertex requires an addition 10cm tighter FV cut
  // (i.e. the argument variables), and the end of the vertex uses the FV cut
  // provided (i.e. the argument variables are 0cm). I don't think there's any
  // need to precalculate the indexs... because we've already selected events
  // which only have 3 PFPs within them (thus 1muon and 2protons)... but this
  // could be wrong
  bool Contained = true;

  for(int i = 0; i < Event->num_pf_particles_; i ++){
    bool StartContained_i = point_inside_FV( this->reco_FV(),
      Event->track_startx_->at(i),
      Event->track_starty_->at(i),
      Event->track_startz_->at(i)
    );

    bool EndContained_i = point_inside_FV( FV_noBorder,
      Event->track_endx_->at(i), Event->track_endy_->at(i),
      Event->track_endz_->at(i)
    );

    if (!StartContained_i || !EndContained_i) {
      Contained = false;
    }
  }
  if (Contained) sel_containedparticles = true;

  // ==============================================================
  // DB Now ensure that the muon and proton candidates pass the momentum
  // threshold requirements of
  // 0.1 <= MuonMomentum <= 1.2
  // 0.3 <= ProtonMomentum <= 1.0

  std::vector<double> ProtonMomenta = std::vector<double>();
  std::vector<int> ProtonIndices = std::vector<int>();

  int counter = 0;
  bool MomentumThresholdPassed = true;
  for(int i = 0; i < Event->num_pf_particles_; i ++){
    if (i == muon_candidate_idx_) {
      if ( Event->track_range_mom_mu_->at(i) < MUON_P_MIN_MOM_CUT
        || Event->track_range_mom_mu_->at(i) > MUON_P_MAX_MOM_CUT )
      {
        MomentumThresholdPassed = false;
      }
    } else {
      float ProtonMomentum = std::sqrt(
        std::pow(Event->track_kinetic_energy_p_->at(i) + PROTON_MASS,2)
          - std::pow(PROTON_MASS,2));
      if ( ProtonMomentum < PROTON_MIN_MOM_CUT
        || ProtonMomentum > PROTON_MAX_MOM_CUT )
      {
        MomentumThresholdPassed = false;
      }
      ProtonMomenta.push_back(ProtonMomentum);
      ProtonIndices.push_back(i);
      counter++;
    }
  }

  if (ProtonMomenta.size() == 2) {
    if (ProtonMomenta[0] > ProtonMomenta[1]) {
      LeadingProtonIndex = ProtonIndices[0];
      RecoilProtonIndex = ProtonIndices[1];
    } else {
      LeadingProtonIndex = ProtonIndices[1];
      RecoilProtonIndex = ProtonIndices[0];
    }
  }

  if (MomentumThresholdPassed == true) sel_momentum_threshold_passed_ = true;

  // ==============================
  //Does everything pass selection?

  bool passed = sel_nu_mu_cc_ && sel_npfps_eq_3 && sel_ntracks_eq_3
    && sel_correctparticles && sel_containedparticles
    && sel_momentum_threshold_passed_;

  return passed;
}

void CC1mu2p0pi::define_output_branches() {

  set_branch( &sel_reco_vertex_in_FV_, "sel_reco_vertex_in_FV" );
  set_branch( &sel_has_muon_candidate_, "sel_has_muon_candidate" );
  set_branch( &sel_nu_mu_cc_, "sel_nu_mu_cc" );
  set_branch( &sel_npfps_eq_3, "sel_npfps_eq_3" );
  set_branch( &sel_ntracks_eq_3, "sel_ntracks_eq_3" );
  set_branch( &sel_containedparticles, "sel_containedparticles" );
  set_branch( &sel_correctparticles, "sel_correctparticles" );
  set_branch( &sel_momentum_threshold_passed_,
    "sel_momentum_threshold_passed" );

  set_branch( &sig_truevertex_in_fv_, "sig_truevertex_in_fv" );
  set_branch( &sig_ccnc_, "sig_ccnc_" );
  set_branch( &sig_is_numu_, "sig_is_numu_" );
  set_branch( &sig_two_protons_above_thresh_,
    "sig_two_protons_above_thresh_" );
  set_branch( &sig_one_muon_above_thres_, "sig_one_muon_above_thres_" );
  set_branch( &sig_no_pions_, "sig_no_pions_" );
  set_branch( &sig_mc_n_threshold_muon, "mc_n_threshold_muon" );
  set_branch( &sig_mc_n_threshold_proton, "mc_n_threshold_proton" );
  set_branch( &sig_mc_n_threshold_pion0, "mc_n_threshold_pion0" );
  set_branch( &sig_mc_n_threshold_pionpm, "mc_n_threshold_pionpm" );

  set_branch( &LeadingProtonIndex, "LeadingProtonIndex" );
  set_branch( &RecoilProtonIndex, "RecoilProtonIndex" );

  set_branch( &Reco_CosPlPr, "Reco_CosPlPr" );
  set_branch( &Reco_CosMuPsum, "Reco_CosMuPsum" );
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

void CC1mu2p0pi::reset() {

  sel_reco_vertex_in_FV_ = false;
  sel_has_muon_candidate_ = false;
  sel_nu_mu_cc_ = false;
  sel_npfps_eq_3 = false;
  sel_ntracks_eq_3 = false;
  sel_containedparticles = false;
  sel_correctparticles = false;
  sel_momentum_threshold_passed_ = false;

  sig_truevertex_in_fv_ = false;
  sig_ccnc_ = false;
  sig_is_numu_ = false;
  sig_two_protons_above_thresh_ = false;
  sig_one_muon_above_thres_ = false;
  sig_no_pions_ = false;
  sig_mc_n_threshold_muon = BOGUS_INDEX;
  sig_mc_n_threshold_proton = BOGUS_INDEX;
  sig_mc_n_threshold_pion0 = BOGUS_INDEX;
  sig_mc_n_threshold_pionpm = BOGUS_INDEX;

  LeadingProtonIndex = BOGUS_INDEX;
  RecoilProtonIndex = BOGUS_INDEX;

  Reco_CosPlPr = BOGUS;
  Reco_CosMuPsum = BOGUS;
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

void CC1mu2p0pi::define_category_map() {
  // Use the shared category map for 1p/2p/Np/Xp
  categ_map_ = CC1muXp_MAP;
}
