#include "XSecAnalyzer/Selections/CC1muXp0pi.hh"
#include "XSecAnalyzer/Selections/EventCategoriesXp.hh"
#include "XSecAnalyzer/TreeUtils.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/Constants.hh"

CC1muXp0pi::CC1muXp0pi() : SelectionBase("CC1muXp0pi") {
  CalcType = kOpt4;
  DefineCategoryMap();
}

void CC1muXp0pi::DefineConstants() {
  DefineTrueFV(21.5,234.85,-95.0,95.0,21.5,966.8);
  DefineRecoFV(21.5,234.85,-95.0,95.0,21.5,966.8);
}

void CC1muXp0pi::ComputeTrueObservables(AnalysisEvent* Event) {

  if(!IsEventMCSignal()) return;

  size_t num_mc_daughters = Event->mc_nu_daughter_pdg_->size();

  TVector3 v3muon, v3hadron(0, 0, 0), v3leadproton;
  double   muon_energy = 0, hadron_energy = 0, lead_proton_energy = 0; 

  TLorentzVector p4Muon;
  std::vector<TLorentzVector> p4Hadron;
  std::vector<TLorentzVector> p4leadproton;

  double lead_proton_mom = -999;
  int num_p_candidates = 0;

  for ( size_t d = 0u; d < num_mc_daughters; ++d ) {
      int pdg = Event->mc_nu_daughter_pdg_->at( d );
      float energy = Event->mc_nu_daughter_energy_->at( d );
      float px = Event->mc_nu_daughter_px_->at( d );
      float py = Event->mc_nu_daughter_py_->at( d );
      float pz = Event->mc_nu_daughter_pz_->at( d );

      if(pdg == MUON){
        double mom = real_sqrt( std::pow(energy, 2) - std::pow(MUON_MASS, 2) );
        v3muon.SetXYZ(px, py, pz);
        p4Muon.SetPxPyPzE(px, py, pz, energy);
        muon_energy = energy;
        if(fabs(mom - v3muon.Mag()) < 1e-12)  std::cout << "Weird case, the mom is inconsistent with two method! : " << 
          mom << " != " << v3muon.Mag()  << std::endl;
      }
      else if(pdg == PROTON){
        double mom = real_sqrt( std::pow(energy, 2) - std::pow(PROTON_MASS, 2) );
        if ( mom >= LEAD_P_MIN_MOM_CUT && mom <= LEAD_P_MAX_MOM_CUT ){
          num_p_candidates++;
          TVector3 v3p(px, py, pz);
          v3hadron += v3p;
          hadron_energy+=energy;
          TLorentzVector p4proton(px, py, pz, energy);
          p4Hadron.push_back(p4proton);

          if(mom > lead_proton_mom){
            lead_proton_mom = mom;
            lead_proton_energy = energy;
            v3leadproton.SetXYZ(px, py, pz);
            mc_lead_proton_energy_ = energy;
            mc_lead_proton_costh_ = v3p.CosTheta();
            mc_lead_proton_p_ = v3p.Mag();
          }
        }
      }
  }

  mc_muon_energy_ = muon_energy;
  mc_muon_costh_ = v3muon.CosTheta();
  mc_muon_p_ = v3muon.Mag();
  mc_hadron_energy_ = hadron_energy;
  mc_hadron_costh_ = v3hadron.CosTheta();
  mc_hadron_p_ = v3hadron.Mag();

  // Compute true STVs if the event contains both a muon and a leading
  // proton
  //
  if(num_p_candidates != 0){
    STVTools.CalculateSTVs(p4Muon, p4Hadron, CalcType);

    mc_hadron_delta_pT_ = STVTools.ReturnPt();
    mc_hadron_delta_phiT_ = STVTools.ReturnDeltaPhiT() * TMath::Pi()/180.;
    mc_hadron_delta_alphaT_ = STVTools.ReturnDeltaAlphaT() * TMath::Pi()/180.;
    mc_hadron_delta_alpha3Dq_ = STVTools.ReturnDeltaAlpha3Dq() * TMath::Pi()/180.;
    mc_hadron_delta_alpha3DMu_ = STVTools.ReturnDeltaAlpha3DMu() * TMath::Pi()/180.;
    mc_hadron_delta_phi3D_ = STVTools.ReturnDeltaPhi3D() * TMath::Pi()/180.;
    mc_hadron_delta_pL_ = STVTools.ReturnPL();
    mc_hadron_pn_ = STVTools.ReturnPn();
    mc_hadron_delta_pTx_ = STVTools.ReturnPtx();
    mc_hadron_delta_pTy_ = STVTools.ReturnPty();
    mc_hadron_theta_mu_p_ = STVTools.ReturnThetaMuHadron();

    STVTools.CalculateSTVs(v3muon, v3leadproton, muon_energy, lead_proton_energy, CalcType);

    mc_lead_proton_delta_pT_ = STVTools.ReturnPt();
    mc_lead_proton_delta_phiT_ = STVTools.ReturnDeltaPhiT() * TMath::Pi()/180.;
    mc_lead_proton_delta_alphaT_ = STVTools.ReturnDeltaAlphaT() * TMath::Pi()/180.;
    mc_lead_proton_delta_alpha3Dq_ = STVTools.ReturnDeltaAlpha3Dq() * TMath::Pi()/180.;
    mc_lead_proton_delta_alpha3DMu_ = STVTools.ReturnDeltaAlpha3DMu() * TMath::Pi()/180.;
    mc_lead_proton_delta_phi3D_ = STVTools.ReturnDeltaPhi3D() * TMath::Pi()/180.;
    mc_lead_proton_delta_pL_ = STVTools.ReturnPL();
    mc_lead_proton_pn_ = STVTools.ReturnPn();
    mc_lead_proton_delta_pTx_ = STVTools.ReturnPtx();
    mc_lead_proton_delta_pTy_ = STVTools.ReturnPty();
    mc_lead_proton_theta_mu_p_ = STVTools.ReturnThetaMuHadron();

  }

}

void CC1muXp0pi::ComputeRecoObservables(AnalysisEvent* Event) {
  
 
  TLorentzVector p4Muon;
  std::vector<TLorentzVector> p4Hadron;
  std::vector<TLorentzVector> p4leadproton;


  if(muon_candidate_idx_ != BOGUS_INDEX){

    float p_dirx = Event->track_dirx_->at( muon_candidate_idx_ );
    float p_diry = Event->track_diry_->at( muon_candidate_idx_ );
    float p_dirz = Event->track_dirz_->at( muon_candidate_idx_ );

    float range_muon_mom = Event->track_range_mom_mu_->at( muon_candidate_idx_ );
    float mcs_muon_mom = Event->track_mcs_mom_mu_->at( muon_candidate_idx_ );

    TVector3 direction(p_dirx, p_diry, p_dirz);
    double MuonEnergy = real_sqrt(range_muon_mom*range_muon_mom + MUON_MASS*MUON_MASS);
    p4Muon.SetPxPyPzE(direction.Unit().X() * range_muon_mom,
        direction.Unit().Y() * range_muon_mom,
        direction.Unit().Z() * range_muon_mom,
        MuonEnergy);
  }

  muon_energy_ = p4Muon.E();
  muon_costh_ = p4Muon.CosTheta();
  muon_p_ = p4Muon.Rho();
  

  if(proton_index.size() > 0){
    for(auto & itp: proton_index){
      TLorentzVector p4p(proton_px.at(itp.first), proton_py.at(itp.first), 
          proton_pz.at(itp.first), proton_e.at(itp.first));
      p4Hadron.push_back(p4p);
    }
    p4leadproton.push_back(p4Hadron[0]);
  }



  // Compute reco STVs if we have both a muon candidate
  // and a leading proton candidate in the event
  if ( muon_candidate_idx_ != BOGUS_INDEX && p4Hadron.size() > 0) {

    STVTools.CalculateSTVs(p4Muon, p4Hadron, CalcType);
    hadron_delta_pT_ = STVTools.ReturnPt();
    hadron_delta_phiT_ = STVTools.ReturnDeltaPhiT() * TMath::Pi()/180.;
    hadron_delta_alphaT_ = STVTools.ReturnDeltaAlphaT() * TMath::Pi()/180.;
    hadron_delta_alpha3Dq_ = STVTools.ReturnDeltaAlpha3Dq() * TMath::Pi()/180.;
    hadron_delta_alpha3DMu_ = STVTools.ReturnDeltaAlpha3DMu() * TMath::Pi()/180;
    hadron_delta_phi3D_ = STVTools.ReturnDeltaPhi3D() * TMath::Pi()/180.;
    hadron_delta_pL_ = STVTools.ReturnPL();
    hadron_pn_ = STVTools.ReturnPn();
    hadron_delta_pTx_ = STVTools.ReturnPtx();
    hadron_delta_pTy_ = STVTools.ReturnPty();
    hadron_theta_mu_p_ = STVTools.ReturnThetaMuHadron();

    STVTools.CalculateSTVs(p4Muon, p4leadproton, CalcType);
    lead_proton_delta_pT_ = STVTools.ReturnPt();
    lead_proton_delta_phiT_ = STVTools.ReturnDeltaPhiT() * TMath::Pi()/180.;
    lead_proton_delta_alphaT_ = STVTools.ReturnDeltaAlphaT() * TMath::Pi()/180.;
    lead_proton_delta_alpha3Dq_ = STVTools.ReturnDeltaAlpha3Dq() * TMath::Pi()/180.;
    lead_proton_delta_alpha3DMu_ = STVTools.ReturnDeltaAlpha3DMu() * TMath::Pi()/180;
    lead_proton_delta_phi3D_ = STVTools.ReturnDeltaPhi3D() * TMath::Pi()/180.;
    lead_proton_delta_pL_ = STVTools.ReturnPL();
    lead_proton_pn_ = STVTools.ReturnPn();
    lead_proton_delta_pTx_ = STVTools.ReturnPtx();
    lead_proton_delta_pTy_ = STVTools.ReturnPty();
    lead_proton_theta_mu_p_ = STVTools.ReturnThetaMuHadron();

  }

}

// Select signal from the inclusive MC

bool CC1muXp0pi::DefineSignal(AnalysisEvent* Event) {

  // there are five criteria:
  // 1. A muon neutrino undergoes a CC interaction with an argon nucleus. Interaction 
  //    vertex should be within FV.
  // 2. There are no mesons in final states.
  // 3. The momentum of the outgoing muon lies within the interval [0.1, 1.2] GeV/c.

  // vertex inside FV
  sig_inFV_ = point_inside_FV(ReturnTrueFV(), Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_);

  // initial state is a muon neutrion
  sig_isNuMu_ = (Event->mc_nu_pdg_ == MUON_NEUTRINO);

  // CC or NC process
  bool IsNC = (Event->mc_nu_ccnc_ == NEUTRAL_CURRENT);

  sig_noFSMesons_= true;
  sig_mc_no_fs_pi0_ = true;
  sig_mc_no_fs_eta_ = true;
  sig_mc_no_charged_pi_above_threshold_ = true;

  sig_muonInMomRange_ = false;

  sig_nProtons_in_Momentum_range = 0;

  double LeadProtonMomentum = 0.;

  TLorentzVector p4g, p4mu, p4p, p4n;
  std::vector<TLorentzVector> p4gv, p4muv, p4pv, p4nv;


  p4gv.clear();
  p4muv.clear();
  p4pv.clear();
  p4nv.clear();


  sig_num_gamma = 0;
  sig_num_muplus = 0;
  sig_num_muminus = 0;
  sig_num_eplus = 0;
  sig_num_eminus = 0;
  sig_num_nu = 0;
  sig_num_antinu = 0;
  sig_num_proton = 0;
  sig_num_neutron = 0;

  // loop all the final state in mc truth
  for ( size_t p = 0u; p < Event->mc_nu_daughter_pdg_->size(); ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    float energy = Event->mc_nu_daughter_energy_->at( p );
    float px = Event->mc_nu_daughter_px_->at( p );
    float py = Event->mc_nu_daughter_py_->at( p );
    float pz = Event->mc_nu_daughter_pz_->at( p );

    // Do the general check for (anti)mesons first before considering
    // any individual PDG codes
    if ( is_meson_or_antimeson(pdg) ) {
      sig_noFSMesons_ = false;
    }


    if ( pdg == MUON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(MUON_MASS, 2) );
      if ( mom >= MUON_P_MIN_MOM_CUT && mom <= MUON_P_MAX_MOM_CUT ) {
        sig_muonInMomRange_= true;
        p4mu.SetPxPyPzE(px, py, pz, energy);
        p4muv.push_back(p4mu);
        sig_num_muminus++;
      }
    }
    else if(pdg == ANTI_MUON){
      sig_num_muplus++;
    }
    else if(pdg == ELECTRON){
      sig_num_eminus++;
    }
    else if(pdg == ANTI_ELECTRON){
      sig_num_eplus++;
    }
    else if(pdg == MUON_NEUTRINO){
      sig_num_nu++;
    }
    else if(pdg == ANTI_MUON_NEUTRINO){
      sig_num_antinu++;
    }
    else if ( pdg == PROTON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PROTON_MASS, 2) );
      if ( mom > LeadProtonMomentum ) LeadProtonMomentum = mom;
      if ( mom >= LEAD_P_MIN_MOM_CUT && mom <= LEAD_P_MAX_MOM_CUT ) sig_nProtons_in_Momentum_range++;
      p4p.SetPxPyPzE(px, py, pz, energy);
      p4pv.push_back(p4p);
      sig_num_proton++;
    }
    else if( pdg == NEUTRON){
      p4n.SetPxPyPzE(px, py, pz, energy);
      p4nv.push_back(p4n);
      sig_num_neutron++;
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
    else if(std::abs(pdg) == GAMMA) {
      p4g.SetPxPyPzE(px, py, pz, energy);
      p4gv.push_back(p4g);
      sig_num_gamma++;
    }
  }

  if(p4gv.size() == 2){
    TLorentzVector p4eta = p4gv[0] + p4gv[1];
    if(fabs(p4eta.M() - ETA_MASS) < 0.00001){
      sig_noFSMesons_ = false;
      sig_mc_no_fs_eta_ = false;
    }
  }



  sig_leadProtonMomInRange_ = false;
  // Check that the leading proton has a momentum within the allowed range
  if ( LeadProtonMomentum >= LEAD_P_MIN_MOM_CUT && LeadProtonMomentum <= LEAD_P_MAX_MOM_CUT ) {
    sig_leadProtonMomInRange_ = true;
  }

  //  bool ReturnVal = sig_inFV_ && !IsNC && sig_isNuMu_ && sig_muonInMomRange_ && sig_leadProtonMomInRange_ && sig_noFSMesons_;
  //
  //  1. 
  //  2. 
  //  3. 
  bool ReturnVal = sig_inFV_ && !IsNC && sig_isNuMu_ && sig_muonInMomRange_ && sig_noFSMesons_ && !sig_num_gamma && sig_num_muminus == 1
    && !sig_num_muplus && !sig_num_eplus && !sig_num_eminus && !sig_num_nu && !sig_num_antinu;
  return ReturnVal;

}

bool CC1muXp0pi::Selection(AnalysisEvent* Event) {
  FiducialVolume PCV;
  PCV.X_Min = 10.;
  PCV.X_Max = 246.35;
  PCV.Y_Min = -106.5;
  PCV.Y_Max = 106.5;
  PCV.Z_Min = 10.;
  PCV.Z_Max = 1026.8;

  //  The basic selection criteria 
  //  1. fiducial volume
  //  2. topological score
  //  3. 
  sel_reco_vertex_in_FV_ = point_inside_FV(ReturnRecoFV(), Event->nu_vx_, Event->nu_vy_, Event->nu_vz_);

  sel_cosmic_ip_cut_passed_ = Event->cosmic_impact_parameter_ > COSMIC_IP_CUT;



  std::vector<int> muon_index, muon_bt_pdg, rest_bt_pdg;
  std::vector<double> muon_range_mom, muon_mcs_mom, muon_px, muon_py, muon_pz;
  muon_bt_pdg.clear();
  rest_bt_pdg.clear();
  muon_index.clear();
  muon_range_mom.clear();
  muon_mcs_mom.clear();
  muon_px.clear();
  muon_py.clear();
  muon_pz.clear();

  proton_index.clear();
  proton_pid_score.clear();
  proton_e.clear();
  proton_px.clear();
  proton_py.clear();
  proton_pz.clear();

  std::vector<double> rest_llr_pid, rest_trk_score, rest_trk_length, rest_proton_mom;
  rest_llr_pid.clear();
  rest_trk_score.clear();
  rest_trk_length.clear();
  rest_proton_mom.clear();


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

    float p_dirx = Event->track_dirx_->at( p );
    float p_diry = Event->track_diry_->at( p );
    float p_dirz = Event->track_dirz_->at( p );

    float range_muon_mom = Event->track_range_mom_mu_->at( p );
    float mcs_muon_mom = Event->track_mcs_mom_mu_->at( p );

    TVector3 direction(p_dirx, p_diry, p_dirz);

    if ( track_score > MUON_TRACK_SCORE_CUT
        && start_dist < MUON_VTX_DISTANCE_CUT
        && track_length > MUON_LENGTH_CUT
        && pid_score > MUON_PID_CUT )
    {
      muon_candidate_indices.push_back( p );

      if(Event->pfp_true_pdg_->size() > 0)
        muon_bt_pdg.push_back(Event->pfp_true_pdg_->at(p));
      else
        muon_bt_pdg.push_back(0);

      muon_pid_scores.push_back( pid_score );
      muon_index.push_back(p);
      muon_range_mom.push_back(range_muon_mom);
      muon_mcs_mom.push_back(mcs_muon_mom);
      muon_px.push_back(direction.x());
      muon_py.push_back(direction.y());
      muon_pz.push_back(direction.z());
    }
  }

  size_t num_candidates = muon_candidate_indices.size();
  sel_num_muon_candidates = num_candidates;

  for(int num_mu = 0; num_mu < num_candidates; num_mu++){
    sel_muon_candidate_index[num_mu] = muon_index[num_mu];
    sel_muon_bt_pdg[num_mu] = muon_bt_pdg[num_mu];
    sel_muon_candidate_range_mom[num_mu] = muon_range_mom[num_mu];
    sel_muon_candidate_mcs_mom[num_mu] = muon_mcs_mom[num_mu];
    sel_muon_quality[num_mu] = (muon_range_mom[num_mu] - muon_mcs_mom[num_mu])/ muon_range_mom[num_mu];
    sel_muon_candidate_px[num_mu] = muon_px[num_mu];
    sel_muon_candidate_py[num_mu] = muon_py[num_mu];
    sel_muon_candidate_pz[num_mu] = muon_pz[num_mu];
  }

  if ( num_candidates > 0u ) sel_has_muon_candidate_ = true;
   
  muon_candidate_filter_idx_ = -1;

  if ( num_candidates == 1u ) {
    muon_candidate_idx_ = muon_candidate_indices.front();
    muon_candidate_filter_idx_ = 0;
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
        muon_candidate_filter_idx_  = c;
      }
    }
    muon_candidate_idx_ = chosen_index;
  }
  else {
    muon_candidate_idx_ = BOGUS_INDEX;
  }

  sel_nu_mu_cc_ = sel_reco_vertex_in_FV_ && sel_pfp_starts_in_PCV_
    && sel_has_muon_candidate_;


  // Set flags that default to true here
  sel_has_non_proton_particles_ = false;
  sel_no_reco_showers_ = true;
  sel_passed_proton_pid_cut_ = true;
  sel_protons_contained_ = true;

  // Set flags that default to false here
  sel_muon_contained_ = false;
  sel_muon_passed_mom_cuts_ = false;
  sel_muon_quality_ok_ = false;

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

      double frac_diff_range_mcs = std::abs( range_muon_mom - mcs_muon_mom );
      if ( range_muon_mom > 0. ) {
        frac_diff_range_mcs /= range_muon_mom;
        if ( frac_diff_range_mcs < MUON_MOM_QUALITY_CUT ) {
          sel_muon_quality_ok_ = true;
        }
      }
    }
    else {

      // instead of select a lead proton candidate, 
      // we remove the shower and non-proton tracks
      // if we find any shower or non-proton tracks 
      // in an event, than we kill this event.

      // if the events have a shower, than we set sel_has_non_proton_particles_ to be true
      float track_score = Event->pfp_track_score_->at( p );
      if ( track_score <= TRACK_SCORE_CUT ){
        sel_no_reco_showers_ = false;
        sel_has_non_proton_particles_ = true;
      }

      float llr_pid_score = Event->track_llr_pid_score_->at( p );
      float track_length = Event->track_length_->at( p );


      if(Event->pfp_true_pdg_->size() > 0)
        rest_bt_pdg.push_back(Event->pfp_true_pdg_->at(p));
      else
        rest_bt_pdg.push_back(0);

      rest_llr_pid.push_back(llr_pid_score);
      rest_trk_score.push_back(track_score);
      rest_trk_length.push_back(track_length);

      // Check whether the current proton candidate fails the proton PID cut
      if ( llr_pid_score > proton_pid_cut(track_length) ) {
        sel_passed_proton_pid_cut_ = true;
        sel_has_non_proton_particles_ = true;
      }     

      // Check whether the current proton candidate fails the containment cut
      float endx = Event->track_endx_->at( p );
      float endy = Event->track_endy_->at( p );
      float endz = Event->track_endz_->at( p );
      bool end_contained = point_inside_FV(PCV, endx, endy, endz );
      // Bad tracks in the searchingfornues TTree can have
      // bogus track lengths. This skips those.
      //
      //
        float p_dirx = Event->track_dirx_->at( p );
        float p_diry = Event->track_diry_->at( p );
        float p_dirz = Event->track_dirz_->at( p );
        float KEp = Event->track_kinetic_energy_p_->at( p );
        float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );
        rest_proton_mom.push_back(p_mom);
      if ( track_length > 0. && end_contained){


        TVector3 v3p(p_dirx, p_diry, p_dirz);
        v3p = v3p.Unit() * p_mom;
        if(p_mom > 0.25 && p_mom < 1.0){
          proton_index.insert(std::pair<double, int>(track_length, p));
          proton_e.insert(std::pair<double, double>(track_length, KEp + PROTON_MASS));
          proton_px.insert(std::pair<double, double>(track_length, v3p.X()));
          proton_py.insert(std::pair<double, double>(track_length, v3p.Y()));
          proton_pz.insert(std::pair<double, double>(track_length, v3p.Z()));
        }
      }
    }
  }

  size_t num_proton_candidates = proton_index.size();
  sel_num_proton_candidates = num_proton_candidates;

  int i_num_p = 0;
  for(auto & itp: proton_index){
    sel_proton_candidate_index[i_num_p] = itp.second;
    sel_proton_candidate_e[i_num_p] = proton_e.at(itp.first);
    sel_proton_candidate_px[i_num_p] = proton_px.at(itp.first);
    sel_proton_candidate_py[i_num_p] = proton_py.at(itp.first);
    sel_proton_candidate_pz[i_num_p] = proton_pz.at(itp.first);
    i_num_p++;
  }

  sel_num_rest_particles = rest_bt_pdg.size();
  for(int i = 0; i < sel_num_rest_particles; i++){
    sel_rest_llr_pid_score[i] = rest_llr_pid[i];
    sel_rest_track_score[i] = rest_trk_score[i];
    sel_rest_track_length[i] = rest_trk_length[i];
    sel_rest_bt_pdg[i] = rest_bt_pdg[i];
    sel_rest_proton_mom[i] = rest_proton_mom[i];
  }

  // If there are protons in the final state, the one with longest 
  // track length is set to be the leading proton track.

  lead_p_candidate_idx_ = BOGUS_INDEX;
  if(num_proton_candidates > 0u){
    lead_p_candidate_idx_ = proton_index.cbegin()->second;
  }

  // We set different requirement of topological cut according to the 
  // proton multiplicity. Cosmic rays are likely to survive the previous 
  // event selection as 0p (or 1p) events.

  if(num_proton_candidates == 0u){
    sel_topo_cut_passed_ = (Event->topological_score_ > TOPO_SCORE_CUT * 2) ? true : false; // FIXME: need to optimize this selection
  }
  else if(num_proton_candidates == 1u){
    sel_topo_cut_passed_ = (Event->topological_score_ > TOPO_SCORE_CUT * 2) ? true : false; // FIXME: need to optimize this selection
  }
  else{
    sel_topo_cut_passed_ = (Event->topological_score_ > TOPO_SCORE_CUT * 0.0) ? true : false;
  }

  // All right, we've applied all selection cuts. Set the flag that indicates
  // whether all were passed (and thus the event is selected as a CCXp0pi
  // candidate)

  bool sel_CCXp0pi_ = sel_nu_mu_cc_ && sel_muon_passed_mom_cuts_ && sel_muon_contained_ && sel_muon_quality_ok_ 
    && sel_topo_cut_passed_ &&  !sel_has_non_proton_particles_;

  return sel_CCXp0pi_;
}

int CC1muXp0pi::CategorizeEvent(AnalysisEvent* Event) {
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
    if(sig_nProtons_in_Momentum_range == 0){
      if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCC0p0pi_CCQE; // QE
      else if ( Event->mc_nu_interaction_type_ == 10 ) return kNuMuCC0p0pi_CCMEC; // MEC
      else if ( Event->mc_nu_interaction_type_ == 1 ) return kNuMuCC0p0pi_CCRES; // RES
      else return kNuMuCC0p0pi_Other;
    }
    else if (sig_nProtons_in_Momentum_range == 1) {
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
  } 
  //  else if (!sig_leadProtonMomInRange_) {
  //    if ( Event->mc_nu_interaction_type_ == 0 ) return kNuMuCC0p0pi_CCQE; // QE
  //    else if ( Event->mc_nu_interaction_type_ == 10 ) return kNuMuCC0p0pi_CCMEC; // MEC
  //    else if ( Event->mc_nu_interaction_type_ == 1 ) return kNuMuCC0p0pi_CCRES; // RES
  //    else return kNuMuCC0p0pi_Other;
  //  }
  return kNuMuCCOther;
}

void CC1muXp0pi::DefineOutputBranches() {
  SetBranch(&sig_isNuMu_,"mc_is_numu",kBool);
  SetBranch(&sig_inFV_,"mc_vertex_in_FV",kBool);
  SetBranch(&sig_leadProtonMomInRange_,"mc_lead_p_in_range",kBool);
  SetBranch(&sig_muonInMomRange_,"mc_muon_in_mom_range",kBool);
  SetBranch(&sig_noFSMesons_,"mc_no_FS_mesons",kBool);
  SetBranch(&sig_mc_no_charged_pi_above_threshold_,"mc_no_charged_pions_above_thres",kBool);
  SetBranch(&sig_mc_no_fs_pi0_,"mc_no_pi0s",kBool);
  SetBranch(&sig_mc_no_fs_eta_,"mc_no_etas",kBool);
  SetBranch(&sig_nProtons_in_Momentum_range,"nProtons_in_Momentum_range",kInteger);
  SetBranch(&sig_num_gamma,"num_gamma",kInteger);
  SetBranch(&sig_num_muplus,"num_muplus",kInteger);
  SetBranch(&sig_num_muminus,"num_muminus",kInteger);
  SetBranch(&sig_num_eplus,"num_eplus",kInteger);
  SetBranch(&sig_num_eminus,"num_eminus",kInteger);
  SetBranch(&sig_num_nu,"num_nu",kInteger);
  SetBranch(&sig_num_antinu,"num_antinu",kInteger);
  SetBranch(&sig_num_proton,"num_proton",kInteger);
  SetBranch(&sig_num_neutron,"num_neutron",kInteger);

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
  SetBranch(&sel_has_non_proton_particles_,"has_non_proton_particles",kBool);

  Tree->Branch("CC1muXp0pi_num_muon_candidates", &sel_num_muon_candidates, "CC1muXp0pi_num_muon_candidates/I");
  Tree->Branch("CC1muXp0pi_muon_candidate_index", sel_muon_candidate_index, "CC1muXp0pi_muon_candidate_index[CC1muXp0pi_num_muon_candidates]/I");
  Tree->Branch("CC1muXp0pi_muon_candidate_range_mom", sel_muon_candidate_range_mom, "CC1muXp0pi_muon_candidate_range_mom[CC1muXp0pi_num_muon_candidates]/D");
  Tree->Branch("CC1muXp0pi_muon_candidate_mcs_mom", sel_muon_candidate_mcs_mom, "CC1muXp0pi_muon_candidate_mcs_mom[CC1muXp0pi_num_muon_candidates]/D");
  Tree->Branch("CC1muXp0pi_muon_candidate_px", sel_muon_candidate_px, "CC1muXp0pi_muon_candidate_px[CC1muXp0pi_num_muon_candidates]/D");
  Tree->Branch("CC1muXp0pi_muon_candidate_py", sel_muon_candidate_py, "CC1muXp0pi_muon_candidate_py[CC1muXp0pi_num_muon_candidates]/D");
  Tree->Branch("CC1muXp0pi_muon_candidate_pz", sel_muon_candidate_pz, "CC1muXp0pi_muon_candidate_pz[CC1muXp0pi_num_muon_candidates]/D");
  Tree->Branch("CC1muXp0pi_muon_quality", sel_muon_quality, "CC1muXp0pi_muon_quality[CC1muXp0pi_num_muon_candidates]/D");
  Tree->Branch("CC1muXp0pi_muon_bt_pdg", sel_muon_bt_pdg, "CC1muXp0pi_muon_bt_pdg[CC1muXp0pi_num_muon_candidates]/I");
  // SetBranch(&sel_num_muon_candidates,"num_muon_candidates",kInteger);
  Tree->Branch("CC1muXp0pi_num_proton_candidates", &sel_num_proton_candidates, "CC1muXp0pi_num_proton_candidates/I");
  Tree->Branch("CC1muXp0pi_proton_candidate_index", sel_proton_candidate_index, "CC1muXp0pi_proton_candidate_index[CC1muXp0pi_num_proton_candidates]/I");
  Tree->Branch("CC1muXp0pi_proton_candidate_e", sel_proton_candidate_e, "CC1muXp0pi_proton_candidate_e[CC1muXp0pi_num_proton_candidates]/D");
  Tree->Branch("CC1muXp0pi_proton_candidate_px", sel_proton_candidate_px, "CC1muXp0pi_proton_candidate_px[CC1muXp0pi_num_proton_candidates]/D");
  Tree->Branch("CC1muXp0pi_proton_candidate_py", sel_proton_candidate_py, "CC1muXp0pi_proton_candidate_py[CC1muXp0pi_num_proton_candidates]/D");
  Tree->Branch("CC1muXp0pi_proton_candidate_pz", sel_proton_candidate_pz, "CC1muXp0pi_proton_candidate_pz[CC1muXp0pi_num_proton_candidates]/D");


  Tree->Branch("CC1muXp0pi_num_rest_particles", &sel_num_rest_particles, "CC1muXp0pi_num_rest_particles/I");
  Tree->Branch("CC1muXp0pi_rest_bt_pdg", sel_rest_bt_pdg, "CC1muXp0pi_rest_bt_pdg[CC1muXp0pi_num_rest_particles]/I");
  Tree->Branch("CC1muXp0pi_rest_proton_mom", sel_rest_proton_mom, "CC1muXp0pi_rest_proton_mom[CC1muXp0pi_num_rest_particles]/D");
  Tree->Branch("CC1muXp0pi_rest_llr_pid_score", sel_rest_llr_pid_score, "CC1muXp0pi_rest_llr_pid_score[CC1muXp0pi_num_rest_particles]/D");
  Tree->Branch("CC1muXp0pi_rest_track_score", sel_rest_track_score, "CC1muXp0pi_rest_track_score[CC1muXp0pi_num_rest_particles]/D");
  Tree->Branch("CC1muXp0pi_rest_track_length", sel_rest_track_length, "CC1muXp0pi_rest_track_length[CC1muXp0pi_num_rest_particles]/D");

  SetBranch(&lead_p_candidate_idx_,"lead_p_candidate_idx",kInteger);
  SetBranch(&muon_candidate_idx_,"muon_candidate_idx",kInteger);
  SetBranch(&muon_candidate_filter_idx_,"muon_candidate_filter_idx",kInteger);


  SetBranch(&muon_energy_,"reco_muon_energy",kDouble);
  SetBranch(&muon_costh_,"reco_muon_costh",kDouble);
  SetBranch(&muon_p_,"reco_muon_p",kDouble);

  SetBranch(&hadron_delta_pT_,"reco_hadron_delta_pT",kDouble);
  SetBranch(&hadron_delta_phiT_,"reco_hadron_delta_phiT",kDouble);
  SetBranch(&hadron_delta_alphaT_,"reco_hadron_delta_alphaT",kDouble);
  SetBranch(&hadron_delta_alpha3Dq_,"reco_hadron_delta_alpha3Dq",kDouble);
  SetBranch(&hadron_delta_alpha3DMu_,"reco_hadron_delta_alpha3DMu",kDouble);
  SetBranch(&hadron_delta_phi3D_,"reco_hadron_delta_phi3D",kDouble);
  SetBranch(&hadron_delta_pL_,"reco_hadron_delta_pL",kDouble);
  SetBranch(&hadron_pn_,"reco_hadron_pn",kDouble);
  SetBranch(&hadron_delta_pTx_,"reco_hadron_delta_pTx",kDouble);
  SetBranch(&hadron_delta_pTy_,"reco_hadron_delta_pTy",kDouble);
  SetBranch(&hadron_theta_mu_p_,"reco_hadron_theta_mu_p",kDouble);


  SetBranch(&lead_proton_delta_pT_,"reco_lead_proton_delta_pT",kDouble);
  SetBranch(&lead_proton_delta_phiT_,"reco_lead_proton_delta_phiT",kDouble);
  SetBranch(&lead_proton_delta_alphaT_,"reco_lead_proton_delta_alphaT",kDouble);
  SetBranch(&lead_proton_delta_alpha3Dq_,"reco_lead_proton_delta_alpha3Dq",kDouble);
  SetBranch(&lead_proton_delta_alpha3DMu_,"reco_lead_proton_delta_alpha3DMu",kDouble);
  SetBranch(&lead_proton_delta_phi3D_,"reco_lead_proton_delta_phi3D",kDouble);
  SetBranch(&lead_proton_delta_pL_,"reco_lead_proton_delta_pL",kDouble);
  SetBranch(&lead_proton_pn_,"reco_lead_proton_pn",kDouble);
  SetBranch(&lead_proton_delta_pTx_,"reco_lead_proton_delta_pTx",kDouble);
  SetBranch(&lead_proton_delta_pTy_,"reco_lead_proton_delta_pTy",kDouble);
  SetBranch(&lead_proton_theta_mu_p_,"reco_lead_proton_theta_mu_p",kDouble);




  SetBranch(p3mu,"reco_p3_mu",kTVector);
  SetBranch(p3p,"reco_p3_lead_p",kTVector);
  SetBranch(p3_p_vec_,"reco_p3_p_vec",kSTDVector);

  SetBranch(&mc_muon_energy_,"true_muon_energy",kDouble);
  SetBranch(&mc_muon_costh_,"true_muon_costh",kDouble);
  SetBranch(&mc_muon_p_,"true_muon_p",kDouble);

  SetBranch(&mc_lead_proton_energy_,"true_lead_proton_energy",kDouble);
  SetBranch(&mc_lead_proton_costh_,"true_lead_proton_costh",kDouble);
  SetBranch(&mc_lead_proton_p_,"true_lead_proton_p",kDouble);

  SetBranch(&mc_hadron_energy_,"true_hadron_energy",kDouble);
  SetBranch(&mc_hadron_costh_,"true_hadron_costh",kDouble);
  SetBranch(&mc_hadron_p_,"true_hadron_p",kDouble);

  SetBranch(&mc_hadron_delta_pT_,"true_hadron_delta_pT",kDouble);
  SetBranch(&mc_hadron_delta_phiT_,"true_hadron_delta_phiT",kDouble);
  SetBranch(&mc_hadron_delta_alphaT_,"true_hadron_delta_alphaT",kDouble);
  SetBranch(&mc_hadron_delta_alpha3Dq_,"true_hadron_delta_alpha3Dq",kDouble);
  SetBranch(&mc_hadron_delta_alpha3DMu_,"true_hadron_delta_alpha3DMu",kDouble);
  SetBranch(&mc_hadron_delta_phi3D_,"true_hadron_delta_phi3D",kDouble);
  SetBranch(&mc_hadron_delta_pL_,"true_hadron_delta_pL",kDouble);
  SetBranch(&mc_hadron_pn_,"true_hadron_pn",kDouble);
  SetBranch(&mc_hadron_delta_pTx_,"true_hadron_delta_pTx",kDouble);
  SetBranch(&mc_hadron_delta_pTy_,"true_hadron_delta_pTy",kDouble);
  SetBranch(&mc_hadron_theta_mu_p_,"true_hadron_theta_mu_p",kDouble);

  SetBranch(&mc_lead_proton_delta_pT_,"true_lead_proton_delta_pT",kDouble);
  SetBranch(&mc_lead_proton_delta_phiT_,"true_lead_proton_delta_phiT",kDouble);
  SetBranch(&mc_lead_proton_delta_alphaT_,"true_lead_proton_delta_alphaT",kDouble);
  SetBranch(&mc_lead_proton_delta_alpha3Dq_,"true_lead_proton_delta_alpha3Dq",kDouble);
  SetBranch(&mc_lead_proton_delta_alpha3DMu_,"true_lead_proton_delta_alpha3DMu",kDouble);
  SetBranch(&mc_lead_proton_delta_phi3D_,"true_lead_proton_delta_phi3D",kDouble);
  SetBranch(&mc_lead_proton_delta_pL_,"true_lead_proton_delta_pL",kDouble);
  SetBranch(&mc_lead_proton_pn_,"true_lead_proton_pn",kDouble);
  SetBranch(&mc_lead_proton_delta_pTx_,"true_lead_proton_delta_pTx",kDouble);
  SetBranch(&mc_lead_proton_delta_pTy_,"true_lead_proton_delta_pTy",kDouble);
  SetBranch(&mc_lead_proton_theta_mu_p_,"true_lead_proton_theta_mu_p",kDouble);

  SetBranch(mc_p3mu,"true_p3_mu",kTVector);
  SetBranch(mc_p3p,"true_p3_lead_p",kTVector);
  SetBranch(mc_p3_p_vec_,"true_p3_p_vec",kSTDVector);
}

void CC1muXp0pi::DefineCategoryMap() {
  // Use the shared category map for 1p/2p/Np/Xp
  categ_map_ = CC1muXp_MAP;
}
