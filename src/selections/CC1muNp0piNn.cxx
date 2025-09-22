// XSecAnalyzer includes
#include "XSecAnalyzer/TreeUtils.hh"
#include "XSecAnalyzer/Functions.hh"
#include "TGraphAsymmErrors.h"
#include "XSecAnalyzer/neutron_stepwise_reweight.hh"
#include "XSecAnalyzer/Selections/CC1muNp0piNn.hh"
#include "XSecAnalyzer/Selections/EventCategoriesXpXn.hh"
#include <numeric>
#include <functional>

CC1muNp0piNn::CC1muNp0piNn() : SelectionBase( "CC1muNp0piNn" ) {
  this->DefineCategoryMap();
  this->DefineConstants();
}

void CC1muNp0piNn::DefineConstants() {
  DefineTrueFV(21.5,234.85,-95.0,95.0,21.5,966.8);
  DefineRecoFV(21.5,234.85,-95.0,95.0,21.5,966.8);
}

void CC1muNp0piNn::ComputeTrueObservables(AnalysisEvent* Event) {

  //Compute reinteration weights for neutron argon scatters

  //double neutron_event_weight = 1.0;
  //initialize neutron multisim vector with 1000 universes and the random number generator used for the multisim
  std::vector<double> neutron_event_weights(1000, 1.0);
  std::vector<double> neutron_argon_xsec_weight(1, 1.0);

  std::vector<TVector3> starts, ends;

  // refactoring variables
  std::vector<double> all_neutron_track_weights;
  std::vector<double> all_neutron_track_lengths;
  all_neutron_track_weights.clear();
  all_neutron_track_lengths.clear();
  all_neutron_track_weights.reserve(Event->all_mc_pdg_->size());
  all_neutron_track_lengths.reserve(Event->all_mc_pdg_->size());
  std::vector<std::vector<double>> ms_neutron_track_weights(1000);
  for (auto& v : ms_neutron_track_weights) v.reserve(Event->all_mc_pdg_->size());

  // Struct that is used to compile segment lengths so that we do not need to rerun TrackNeutron 1000*number_of_neutrons times for the multisim
  struct NeutronEntry {
    bool   use_weight;   // false -> 1.0 factor
    Segment seg;         // valid only if use_weight==true
  };
  std::vector<NeutronEntry> neutron_entries;
  neutron_entries.reserve(Event->all_mc_pdg_->size());

  CrossSectionData xsec_data = LoadXSecHistogramOnce("/exp/uboone/app/users/birwin/searchingfornues/neutron_event_selection/Burke_Scripts/neutron_inelastic_cross_section.root");
  const TH1D* scaled_xsec_hist = xsec_data.scaled_hist;
  const TH1D* nominal_xsec_hist = xsec_data.nominal_hist;
  const TGraphAsymmErrors* miniCaptainGraph = xsec_data.miniCaptainGraph;

  for (size_t i = 0; i < Event->all_mc_pdg_->size(); ++i) {
    starts.emplace_back(Event->all_mc_vx_->at(i), Event->all_mc_vy_->at(i), Event->all_mc_vz_->at(i));
    ends.emplace_back(Event->all_mc_endx_->at(i), Event->all_mc_endy_->at(i), Event->all_mc_endz_->at(i));
  }

  for (size_t j = 0; j < Event->all_mc_pdg_->size(); ++j) {
    if (Event->all_mc_pdg_->at(j) != 2112) continue;

    // Build single neutron segment (will be empty if neutron starts and ends outside of TPC)
    NeutronTrackResult res = TrackNeutron(
      j,
      *Event->all_mc_pdg_,
      starts,
      ends,
      *Event->all_mc_E_,
      *Event->all_mc_end_process_
    );

    // Save the in-detector length for this neutron (cm); 0.0 if no segment
    double seg_len_cm = 0.0;
    if (!res.segments.empty()) seg_len_cm = res.segments.front().length;
    all_neutron_track_lengths.push_back(seg_len_cm);

    const double KE = Event->all_mc_E_->at(j) - 0.939565; // GeV
    if (KE < 0.1) {
      // Below 100 MeV: identity weight
      all_neutron_track_weights.push_back(1.0);
      neutron_entries.push_back({false, Segment{}});
      continue;
    }

    if (res.segments.empty()) {
      // no in-detector segment -> identity
      all_neutron_track_weights.push_back(1.0);
      neutron_entries.push_back({false, Segment{}});
      continue;
    }

    // One segment per neutron by construction
    const Segment& seg = res.segments.front();

    // scaled central value event weights for each neutron track
    double cv_w = neutron_cv_scale_weight(res.segments, nominal_xsec_hist, scaled_xsec_hist);
    if (std::isnan(cv_w)) {
      std::cout << "we have a NAN weight value. The KE of the neutron is: " << KE << std::endl;
    }
    all_neutron_track_weights.push_back(cv_w);

    neutron_entries.push_back({true, seg});
  }

  double event_scaled_cv_weight = std::accumulate(all_neutron_track_weights.begin(), all_neutron_track_weights.end(), 1.0, std::multiplies<double>());
  neutron_argon_xsec_weight[0] = event_scaled_cv_weight;

  // Multisim for neutron reinteraction uncertainty
  // First we need to initialize the array of 1000 random numbers that will serve as the scale factors for which we scale the neutron-argon xsec
  // in each unique universe. Then we need to cache that array so that we can use the exact same array for subsequent events
  constexpr int kNUniv = 1000;
  // one-time init of the per-universe z values (same for every event)
  static bool z_inited = false;
  static std::array<double, kNUniv> z_for_uni;
  if (!z_inited) {
    std::normal_distribution<> gauss(0.0, 1.0);
    for (int u = 0; u < kNUniv; ++u) {
      std::mt19937 rng(1337 + u);      // fixed seed per universe
      z_for_uni[u] = gauss(rng);       // same z[u] every event
    }
    z_inited = true;
  }

  for (int u = 0; u < kNUniv; ++u) {
    auto& per_neutron = ms_neutron_track_weights[u];

    for (const auto& entry : neutron_entries) {
      if (!entry.use_weight) {
        per_neutron.push_back(1.0);
        continue;
      }
      // single-segment iterative reweight
      const double w = stepwise_weight(entry.seg, scaled_xsec_hist, miniCaptainGraph, z_for_uni[u]);
      per_neutron.push_back(w);
    }
  }

  // Need to take the product of all track weights to get an event weight for each universe in the multisim
  for (int u = 0; u < 1000; ++u) {
  neutron_event_weights[u] =
      std::accumulate(ms_neutron_track_weights[u].begin(),
                      ms_neutron_track_weights[u].end(),
                      1.0, std::multiplies<double>());
  }

  *weight_neutron_reint = neutron_event_weights;
  *weight_neutron_argon_xsec = neutron_argon_xsec_weight;
  *neutron_track_weights = all_neutron_track_weights;


  size_t num_mc_daughters = Event->mc_nu_daughter_pdg_->size();

  // Set the true 3-momentum of the final-state muon if there is one
  bool true_muon = ( sig_isNuMu_ && Event->mc_nu_ccnc_ == CHARGED_CURRENT );
  if ( true_muon ) {
    // Loop over the MC neutrino daughters, find the muon, and get its
    // true 3-momentum. Note that we assume there is only one muon in
    // this loop.
    bool found_muon = false;
    for ( size_t d = 0u; d < num_mc_daughters; ++d ) {
      int pdg = Event->mc_nu_daughter_pdg_->at( d );
      if ( pdg == MUON ) {
        found_muon = true;
        float px = Event->mc_nu_daughter_px_->at( d );
        float py = Event->mc_nu_daughter_py_->at( d );
        float pz = Event->mc_nu_daughter_pz_->at( d );
	*mc_p3mu = TVector3( px, py, pz );
        break;
      }
    }

    if ( !found_muon ) {
      std::cout << "WARNING: Missing muon in MC signal event!\n";
      return;
    }
  }

  // Set the true 3-momentum of the leading proton (if there is one)
  float max_mom = LOW_FLOAT;
  for ( int p = 0; p < num_mc_daughters; ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    if ( pdg == PROTON )
    {
      float px = Event->mc_nu_daughter_px_->at( p );
      float py = Event->mc_nu_daughter_py_->at( p );
      float pz = Event->mc_nu_daughter_pz_->at( p );
      TVector3 temp_p3 = TVector3( px, py, pz );

      mc_p3_p_vec_->push_back( temp_p3 );

      float mom = temp_p3.Mag();
      if ( mom > max_mom ) {
        max_mom = mom;
	*mc_p3p = temp_p3;
      }
    }
  }

  // TODO: reduce code duplication by just getting the leading proton
  // 3-momentum from this sorted vector
  // Sort the true proton 3-momenta in order from highest to lowest magnitude
  std::sort( mc_p3_p_vec_->begin(), mc_p3_p_vec_->end(), [](const TVector3& a,
    const TVector3& b) -> bool { return a.Mag() > b.Mag(); } );

  // If the event contains a leading proton, then set the 3-momentum
  // accordingly
  bool true_lead_p = ( max_mom != LOW_FLOAT );
  if ( !true_lead_p && IsEventMCSignal() ) {
    // If it doesn't for a signal event, then something is wrong.
    std::cout << "WARNING: Missing leading proton in MC signal event!\n";
    return;
  }

  // Compute true STVs if the event contains both a muon and a leading
  // proton
  if ( true_muon && true_lead_p ) {
    double MuonEnergy = real_sqrt(mc_p3mu->Mag()*mc_p3mu->Mag() + MUON_MASS*MUON_MASS);
    double ProtonEnergy = real_sqrt(mc_p3p->Mag()*mc_p3p->Mag() + PROTON_MASS*PROTON_MASS);
    STVTools.CalculateSTVs(*(mc_p3mu.get()), *(mc_p3p.get()), MuonEnergy, ProtonEnergy, CalcType);

    mc_delta_pT_ = STVTools.ReturnPt();
    mc_delta_phiT_ = STVTools.ReturnDeltaPhiT() * TMath::Pi()/180.;
    mc_delta_alphaT_ = STVTools.ReturnDeltaAlphaT() * TMath::Pi()/180.;
    mc_delta_pL_ = STVTools.ReturnPL();
    mc_pn_ = STVTools.ReturnPn();
    mc_delta_pTx_ = STVTools.ReturnPtx();
    mc_delta_pTy_ = STVTools.ReturnPty();

    mc_theta_mu_p_ = std::acos( mc_p3mu->Dot(*mc_p3p)
      / mc_p3mu->Mag() / mc_p3p->Mag() );
  }
}

void CC1muNp0piNn::ComputeRecoObservables(AnalysisEvent* Event) {

  // In cases where we failed to find a muon candidate, check whether there are
  // at least two generation == 2 PFParticles. If there are, then compute the
  // usual observables using the longest track as the muon candidate and the
  // second-longest track as the leading proton candidate. This will enable
  // sideband studies of NC backgrounds in the STV phase space.
  if ( !sel_has_muon_candidate_ ) {
    float max_trk_len = LOW_FLOAT;
    int max_trk_idx = BOGUS_INDEX;

    float next_to_max_trk_len = LOW_FLOAT;
    int next_to_max_trk_idx = BOGUS_INDEX;

    for ( int p = 0; p < Event->num_pf_particles_; ++p ) {

      // Only include direct neutrino daughters (generation == 2)
      unsigned int generation = Event->pfp_generation_->at( p );
      if ( generation != 2u ) continue;

      float trk_len = Event->track_length_->at( p );

      if ( trk_len > next_to_max_trk_len ) {
        next_to_max_trk_len = trk_len;
        next_to_max_trk_idx = p;

        if ( next_to_max_trk_len > max_trk_len ) {
          next_to_max_trk_len = max_trk_len;
          next_to_max_trk_idx = max_trk_idx;

          max_trk_len = trk_len;
          max_trk_idx = p;
        }
      }
    }

    // If we found at least two usable PFParticles, then assign the indices to
    // be used below
    if ( max_trk_idx != BOGUS_INDEX && next_to_max_trk_idx != BOGUS_INDEX ) {
      muon_candidate_idx_ = max_trk_idx;
      lead_p_candidate_idx_ = next_to_max_trk_idx;
    }
  }

  // Set the reco 3-momentum of the muon candidate if we found one
  bool muon = muon_candidate_idx_ != BOGUS_INDEX;
  if ( muon ) {
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

    *p3mu = TVector3( mu_dirx, mu_diry, mu_dirz );
    *p3mu = p3mu->Unit() * muon_mom;
  }

    // Set the reco 3-momentum of the leading proton candidate if we found one
  bool lead_p = lead_p_candidate_idx_ != BOGUS_INDEX;
  if ( lead_p ) {

    float p_dirx = Event->track_dirx_->at( lead_p_candidate_idx_ );
    float p_diry = Event->track_diry_->at( lead_p_candidate_idx_ );
    float p_dirz = Event->track_dirz_->at( lead_p_candidate_idx_ );
    float KEp = Event->track_kinetic_energy_p_->at( lead_p_candidate_idx_ );
    float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

    *p3p = TVector3( p_dirx, p_diry, p_dirz );
    *p3p = p3p->Unit() * p_mom;
  }

  // Set the reco 3-momenta of all proton candidates (i.e., all generation == 2
  // tracks except the muon candidate) assuming we found both a muon candidate
  // and at least one proton candidate.
  if ( muon && lead_p ) {
    for ( int p = 0; p < Event->num_pf_particles_; ++p ) {
      // Skip the muon candidate
      if ( p == muon_candidate_idx_ ) continue;

      // Only include direct neutrino daughters (generation == 2)
      unsigned int generation = Event->pfp_generation_->at( p );
      if ( generation != 2u ) continue;

      float p_dirx = Event->track_dirx_->at( p );
      float p_diry = Event->track_diry_->at( p );
      float p_dirz = Event->track_dirz_->at( p );
      float KEp = Event->track_kinetic_energy_p_->at( p );
      float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

      TVector3 p3_temp( p_dirx, p_diry, p_dirz );
      p3_temp = p3_temp.Unit() * p_mom;

      p3_p_vec_->push_back( p3_temp );
    }

    // TODO: reduce code duplication by just getting the leading proton
    // 3-momentum from this sorted vector
    // Sort the reco proton 3-momenta in order from highest to lowest magnitude
    std::sort( p3_p_vec_->begin(), p3_p_vec_->end(), [](const TVector3& a,
      const TVector3& b) -> bool { return a.Mag() > b.Mag(); } );
  }

  // Compute reco STVs if we have both a muon candidate
  // and a leading proton candidate in the event
  if ( muon && lead_p ) {
    double MuonEnergy = real_sqrt(p3mu->Mag()*p3mu->Mag() + MUON_MASS*MUON_MASS);
    double ProtonEnergy	= real_sqrt(p3p->Mag()*p3p->Mag() + PROTON_MASS*PROTON_MASS);
    STVTools.CalculateSTVs(*(p3mu), *(p3p), MuonEnergy, ProtonEnergy, CalcType);

    delta_pT_ = STVTools.ReturnPt();
    delta_phiT_ = STVTools.ReturnDeltaPhiT() * TMath::Pi()/180.;
    delta_alphaT_ = STVTools.ReturnDeltaAlphaT() * TMath::Pi()/180.;
    delta_pL_ = STVTools.ReturnPL();
    pn_ = STVTools.ReturnPn();
    delta_pTx_ = STVTools.ReturnPtx();
    delta_pTy_ = STVTools.ReturnPty();

    theta_mu_p_ = std::acos( p3mu->Dot(*p3p) / p3mu->Mag() / p3p->Mag() );
  }

}

bool CC1muNp0piNn::DefineSignal(AnalysisEvent* Event) {
  sig_inFV_ = point_inside_FV(ReturnTrueFV(), Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_);
  sig_isNuMu_ = (Event->mc_nu_pdg_ == MUON_NEUTRINO);
  bool IsNC = (Event->mc_nu_ccnc_ == NEUTRAL_CURRENT);

  sig_noFSMesons_= true;
  sig_mc_no_fs_pi0_ = true;
  sig_mc_no_charged_pi_above_threshold_ = true;

  sig_muonInMomRange_ = false;

  sig_nProtons_in_Momentum_range = 0;

  sig_num_neutron_above_100MeV = 0;

  double LeadProtonMomentum = 0.;

  for ( size_t p = 0u; p < Event->mc_nu_daughter_pdg_->size(); ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    float energy = Event->mc_nu_daughter_energy_->at( p );

    // Do the general check for (anti)mesons first before considering
    // any individual PDG codes
    if ( is_meson_or_antimeson(pdg) ) {
      sig_noFSMesons_ = false;
    }

    if ( pdg == MUON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(MUON_MASS, 2) );
      if ( mom >= MUON_P_MIN_MOM_CUT && mom <= MUON_P_MAX_MOM_CUT ) {
	sig_muonInMomRange_= true;
      }
    }
    else if ( pdg == PROTON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PROTON_MASS, 2) );
      if ( mom > LeadProtonMomentum ) LeadProtonMomentum = mom;
      if ( mom >= LEAD_P_MIN_MOM_CUT && mom <= LEAD_P_MAX_MOM_CUT ) sig_nProtons_in_Momentum_range++;
    }
    else if ( pdg == NEUTRON && energy - NEUTRON_MASS > 0.1 ) { 
	sig_num_neutron_above_100MeV++; 
	sig_mc_FS_neutron_ = true;
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
  }

  sig_leadProtonMomInRange_ = false;
  // Check that the leading proton has a momentum within the allowed range
  if ( LeadProtonMomentum >= LEAD_P_MIN_MOM_CUT && LeadProtonMomentum <= LEAD_P_MAX_MOM_CUT ) {
    sig_leadProtonMomInRange_ = true;
  }

  bool ReturnVal = sig_inFV_ && !IsNC && sig_isNuMu_ && sig_muonInMomRange_ && sig_leadProtonMomInRange_ && sig_noFSMesons_ && sig_mc_FS_neutron_ ;
  return ReturnVal;
}

bool CC1muNp0piNn::Selection(AnalysisEvent* Event) {
  FiducialVolume PCV;
  PCV.X_Min = 10.;
  PCV.X_Max = 246.35;
  PCV.Y_Min = -106.5;
  PCV.Y_Max = 106.5;
  PCV.Z_Min = 10.;
  PCV.Z_Max = 1026.8;

  sel_reco_vertex_in_FV_ = point_inside_FV(ReturnRecoFV(), Event->nu_vx_, Event->nu_vy_, Event->nu_vz_);
  //std::cout << ReturnRecoFV().X_Min << " " << ReturnRecoFV().X_Max << " " << Event->nu_vx_ << " " << sel_reco_vertex_in_FV_ << std::endl;

  sel_topo_cut_passed_ = Event->topological_score_ > TOPO_SCORE_CUT;
  sel_cosmic_ip_cut_passed_ = Event->cosmic_impact_parameter_ > COSMIC_IP_CUT;

  // Apply the containment cut to the starting positions of all
  // reconstructed tracks and showers. Pass this cut by default.
  sel_pfp_starts_in_PCV_ = true;

  //Begin Burke Addition
  int num_secondary_proton_cands = 0;
  TVector3 event_momentum(0.0, 0.0, 0.0);
  TVector3 event_momentum_leading_p(0.0, 0.0, 0.0);
  float event_calo_energy = 0;
  float event_calo_energy_leading_p = 0;
  std::vector<int> secondary_proton_index;
  pfp_proximity_map.clear();
  pfp_direction_map.clear();
  pfp_len_map.clear();
  pfp_vtx_disp_map.clear();
  pfp_cos_theta_map.clear();
  pfp_phi_map.clear();
  pfp_start_x_map.clear();
  pfp_start_y_map.clear();
  pfp_start_z_map.clear();
  pfp_end_x_map.clear();
  pfp_end_y_map.clear();
  pfp_end_z_map.clear();
  secondary_proton_index.clear();
  sel_has_secondary_proton_candidate_ = false;
  //End Burke Addition

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

  sel_nu_mu_cc_ = sel_reco_vertex_in_FV_ && sel_pfp_starts_in_PCV_
    && sel_has_muon_candidate_ && sel_topo_cut_passed_;

  // Fail the shower cut if any showers were reconstructed
  // NOTE: We could do this quicker like this,
  //   sel_no_reco_showers_ = ( num_showers_ > 0 );
  // but it might be nice to be able to adjust the track score for this cut.
  // Thus, we do it the hard way.
  int reco_shower_count = 0;
  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    float tscore = Event->pfp_track_score_->at( p );
    if ( tscore <= TRACK_SCORE_CUT ) ++reco_shower_count;
  }
  // Check the shower cut
  sel_no_reco_showers_ = ( reco_shower_count == 0 );

  // Set flags that default to true here
  sel_passed_proton_pid_cut_ = true;
  sel_protons_contained_ = true;

  // Set flags that default to false here
  sel_muon_contained_ = false;

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


      float p_dirx = Event->track_dirx_->at( p );
      float p_diry = Event->track_diry_->at( p );
      float p_dirz = Event->track_dirz_->at( p );
      TVector3 direction(p_dirx, p_diry, p_dirz);

      // Check that the muon candidate is above threshold. Use the best
      // momentum based on whether it was contained or not.

      float muon_mom = LOW_FLOAT;
      float range_muon_mom = Event->track_range_mom_mu_->at( p );
      float mcs_muon_mom = Event->track_mcs_mom_mu_->at( p );

      if ( sel_muon_contained_ ) muon_mom = range_muon_mom;
      else muon_mom = mcs_muon_mom;

      float muon_cand_px = direction.x()*muon_mom;
      float muon_cand_py = direction.y()*muon_mom;
      float muon_cand_pz = direction.z()*muon_mom;

      event_momentum.SetX(event_momentum.X()+muon_cand_px);
      event_momentum.SetY(event_momentum.Y()+muon_cand_py);
      event_momentum.SetZ(event_momentum.Z()+muon_cand_pz);
      event_calo_energy += std::sqrt(muon_mom*muon_mom+MUON_MASS*MUON_MASS);
      event_momentum_leading_p.SetX(event_momentum_leading_p.X()+muon_cand_px);
      event_momentum_leading_p.SetY(event_momentum_leading_p.Y()+muon_cand_py);
      event_momentum_leading_p.SetZ(event_momentum_leading_p.Z()+muon_cand_pz);
      event_calo_energy_leading_p += std::sqrt(muon_mom*muon_mom+MUON_MASS*MUON_MASS);

      if ( muon_mom >= MUON_P_MIN_MOM_CUT && muon_mom <= MUON_P_MAX_MOM_CUT ) {
        sel_muon_passed_mom_cuts_ = true;
      }

      // Apply muon candidate quality cut by comparing MCS and range-based
      // momentum estimators. Default to failing the cut.
      sel_muon_quality_ok_ = false;

      double frac_diff_range_mcs = std::abs( range_muon_mom - mcs_muon_mom );
      if ( range_muon_mom > 0. ) {
        frac_diff_range_mcs /= range_muon_mom;
        if ( frac_diff_range_mcs < MUON_MOM_QUALITY_CUT ) {
          sel_muon_quality_ok_ = true;
        }
      }
    }
    else {

      float track_score = Event->pfp_track_score_->at( p );
      if ( track_score <= TRACK_SCORE_CUT ) continue;

      // Bad tracks in the searchingfornues TTree can have
      // bogus track lengths. This skips those.
      float track_length = Event->track_length_->at( p );
      if ( track_length <= 0. ) continue;

      float llr_pid_score = Event->track_llr_pid_score_->at( p );

      //Begin Burke Edits : add proton candidate transverse momentum to delta pt vector
      float p_dirx = Event->track_dirx_->at( p );
      float p_diry = Event->track_diry_->at( p );
      float p_dirz = Event->track_dirz_->at( p );
      float KEp = Event->track_kinetic_energy_p_->at( p );
      float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );
      //End Burke Edits

      //Neutron Tagging
      float sec_p_cand_startx = Event->track_startx_->at( p );
      float sec_p_cand_starty = Event->track_starty_->at( p );
      float sec_p_cand_startz = Event->track_startz_->at( p );
      float sec_p_cand_endx = Event->track_endx_->at( p );
      float sec_p_cand_endy = Event->track_endy_->at( p );
      float sec_p_cand_endz = Event->track_endz_->at( p );
      float calc_track_length = std::sqrt(pow(sec_p_cand_endx-sec_p_cand_startx,2)+pow(sec_p_cand_endy-sec_p_cand_starty,2)+pow(sec_p_cand_endz-sec_p_cand_startz,2));
      float sec_p_cand_vtx_disp = std::sqrt(pow(Event->nu_vx_-sec_p_cand_startx,2)+pow(Event->nu_vy_-sec_p_cand_starty,2)+pow(Event->nu_vz_-sec_p_cand_startz,2));
      float sec_p_cand_direction = ((sec_p_cand_startx-Event->nu_vx_)*(sec_p_cand_endx-sec_p_cand_startx)+(sec_p_cand_starty-Event->nu_vy_)*(sec_p_cand_endy-sec_p_cand_starty)+(sec_p_cand_startz-Event->nu_vz_)*(sec_p_cand_endz-sec_p_cand_startz))/std::sqrt((pow(sec_p_cand_startx-Event->nu_vx_,2)+pow(sec_p_cand_starty-Event->nu_vy_,2)+pow(sec_p_cand_startz-Event->nu_vz_,2))*(pow(sec_p_cand_endx-sec_p_cand_startx,2)+pow(sec_p_cand_endy-sec_p_cand_starty,2)+pow(sec_p_cand_endz-sec_p_cand_startz,2)));
      float sec_p_cand_costheta = (sec_p_cand_startz-Event->nu_vz_)/std::sqrt(pow(sec_p_cand_startx-Event->nu_vx_,2)+pow(sec_p_cand_starty-Event->nu_vy_,2)+pow(sec_p_cand_startz-Event->nu_vz_,2));
      float sec_p_cand_theta = atan2(std::sqrt(pow(sec_p_cand_starty-Event->nu_vy_,2)+pow(sec_p_cand_startx-Event->nu_vx_,2)),sec_p_cand_startz-Event->nu_vz_);
      float sec_p_cand_phi = atan2(sec_p_cand_starty-Event->nu_vy_,sec_p_cand_startx-Event->nu_vx_);

      bool secondary_proton_cand_contained = point_inside_FV(PCV, sec_p_cand_startx, sec_p_cand_starty, sec_p_cand_startz ) && point_inside_FV(PCV, sec_p_cand_endx, sec_p_cand_endy, sec_p_cand_endz );

      float pfp_proximity = 9000.0;
      // Loop through each track
      for ( int z = 0; z < Event->nonprim_trk_score_v_->size(); ++z ) {
	// Access the track data once per iteration
	float pfp_start_x = Event->nonprim_trk_sce_start_x_v_->at(z);
	float pfp_start_y = Event->nonprim_trk_sce_start_y_v_->at(z);
	float pfp_start_z = Event->nonprim_trk_sce_start_z_v_->at(z);
	float pfp_end_x = Event->nonprim_trk_sce_end_x_v_->at(z);
	float pfp_end_y = Event->nonprim_trk_sce_end_y_v_->at(z);
	float pfp_end_z = Event->nonprim_trk_sce_end_z_v_->at(z);
	// Only proceed if the start point is valid and not the same as the candidate start position
	if ( pfp_start_x > -500 && pfp_start_y > -500 && pfp_start_z > -500 &&
         sec_p_cand_startx != pfp_start_x && sec_p_cand_starty != pfp_start_y && sec_p_cand_startz != pfp_start_z) {
	  float dx_start_start = pfp_start_x - sec_p_cand_startx;
          float dy_start_start = pfp_start_y - sec_p_cand_starty;
          float dz_start_start = pfp_start_z - sec_p_cand_startz;
          float pfp_proximity_start_start = std::sqrt(dx_start_start * dx_start_start + dy_start_start * dy_start_start + dz_start_start * dz_start_start);

          float dx_end_start = pfp_end_x - sec_p_cand_startx;
          float dy_end_start = pfp_end_y - sec_p_cand_starty;
          float dz_end_start = pfp_end_z - sec_p_cand_startz;
          float pfp_proximity_end_start = std::sqrt(dx_end_start * dx_end_start + dy_end_start * dy_end_start + dz_end_start * dz_end_start);

          float dx_start_end = pfp_start_x - sec_p_cand_endx;
          float dy_start_end = pfp_start_y - sec_p_cand_endy;
          float dz_start_end = pfp_start_z - sec_p_cand_endz;
          float pfp_proximity_start_end = std::sqrt(dx_start_end * dx_start_end + dy_start_end * dy_start_end + dz_start_end * dz_start_end);

          float dx_end_end = pfp_end_x - sec_p_cand_endx;
          float dy_end_end = pfp_end_y - sec_p_cand_endy;
          float dz_end_end = pfp_end_z - sec_p_cand_endz;
          float pfp_proximity_end_end = std::sqrt(dx_end_end * dx_end_end + dy_end_end * dy_end_end + dz_end_end * dz_end_end);

	  // Find the minimum proximity value
	  float pfp_proximity_temp = std::min({pfp_proximity_start_start, pfp_proximity_end_start, pfp_proximity_start_end, pfp_proximity_end_end});
	  // Update the minimum proximity if needed
	  if ( pfp_proximity_temp < pfp_proximity ) {
            pfp_proximity = pfp_proximity_temp;
          }	  
	}
      }//All PFP loop
      pfp_proximity_map.insert(std::pair<int, double>(p,pfp_proximity));
      pfp_direction_map.insert(std::pair<int, double>(p,sec_p_cand_direction));
      pfp_len_map.insert(std::pair<int, double>(p,track_length));
      pfp_vtx_disp_map.insert(std::pair<int, double>(p,sec_p_cand_vtx_disp));
      pfp_cos_theta_map.insert(std::pair<int, double>(p,sec_p_cand_costheta));
      pfp_theta_map.insert(std::pair<int, double>(p,sec_p_cand_theta));
      pfp_phi_map.insert(std::pair<int, double>(p,sec_p_cand_phi));
      pfp_start_x_map.insert(std::pair<int, double>(p,sec_p_cand_startx));
      pfp_start_y_map.insert(std::pair<int, double>(p,sec_p_cand_starty));
      pfp_start_z_map.insert(std::pair<int, double>(p,sec_p_cand_startz));
      pfp_end_x_map.insert(std::pair<int, double>(p,sec_p_cand_endx));
      pfp_end_y_map.insert(std::pair<int, double>(p,sec_p_cand_endy));
      pfp_end_z_map.insert(std::pair<int, double>(p,sec_p_cand_endz));
    
      //Begin Neutron Event Selection Edits
      if ( secondary_proton_cand_contained){
        if ( track_score > NEUTRON_TRACK_SCORE_CUT ){
          if ( llr_pid_score < NEUTRON_PROTON_LLRPID ) {
            if ( NEUTRON_UPPER_DISP > sec_p_cand_vtx_disp && sec_p_cand_vtx_disp > NEUTRON_LOWER_DISP ) {
              if ( sec_p_cand_direction > NEUTRON_DIRECTION ) {
                if ( pfp_proximity > NEUTRON_PROXIMITY ) {
                  sel_has_secondary_proton_candidate_ = true;
                  secondary_proton_index.push_back( p );
                  num_secondary_proton_cands++;
		  // If this pfp is a neutron candidate, skip to the next pfp in the list to not allow for the same pfp to be a primary proton candidate
		  continue;
                }
              } //Secondary proton direction
	    } //Secondary proton vertex displacement
          } //Secondary proton llr pid
        } //Secondary proton track score
      } //Secondary proton containment
      if (std::count(secondary_proton_index.begin(),secondary_proton_index.end(),p) == 0)
        {
          event_momentum.SetX(event_momentum.X()+(p_dirx*p_mom));
          event_momentum.SetY(event_momentum.Y()+(p_diry*p_mom));
          event_momentum.SetZ(event_momentum.Z()+(p_dirz*p_mom));
          event_calo_energy += std::sqrt(p_mom*p_mom+PROTON_MASS*PROTON_MASS) - PROTON_MASS + .0309; //Proton KE + Proton binding energy to argon
        }
  
      // All tracks that start above the neutron induced proton candidate vertex separation point are not primary proton candidates
      if ( sec_p_cand_vtx_disp > 10.0 ) continue;

      // We found a reco track that is not the muon candidate and is also not a neutron candidate. All such
      // tracks are considered proton candidates
      sel_has_p_candidate_ = true;

      // Check whether the current proton candidate fails the proton PID cut
      if ( llr_pid_score > proton_pid_cut(track_length) ) {
        sel_passed_proton_pid_cut_ = false;
      }

      // Check whether the current proton candidate fails the containment cut
      float endx = Event->track_endx_->at( p );
      float endy = Event->track_endy_->at( p );
      float endz = Event->track_endz_->at( p );
      bool end_contained = point_inside_FV(PCV, endx, endy, endz );
      if ( !end_contained ) sel_protons_contained_ = false;

    }//End of non-muon candidate else statement
  }//End of Gen 2 neutrino slice pfp loop
  
  // Don't bother to apply the cuts that involve the leading
  // proton candidate if we don't have one
  if ( !sel_has_p_candidate_ ) {
    return false;
  }

  // All that remains is to apply the leading proton candidate cuts. We could
  // search for it above, but doing it here makes the code more readable (with
  // likely negligible impact on performance)
  float lead_p_track_length = LOW_FLOAT;
  size_t lead_p_index = 0u;
  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    // Skip the muon candidate reco track (this function assumes that it has
    // already been found)
    if ( p == muon_candidate_idx_ ) continue;

    // Skip PFParticles that are shower-like (track scores near 0)
    float track_score = Event->pfp_track_score_->at( p );
    if ( track_score <= TRACK_SCORE_CUT ) continue;

    // All non-muon-candidate reco tracks are considered proton candidates
    float track_length = Event->track_length_->at( p );
    if ( track_length <= 0. ) continue;

    // All tracks that start above the neutron induced proton candidate vertex separation point are not primary proton candidates
    float startx = Event->track_startx_->at( p );
    float starty = Event->track_starty_->at( p );
    float startz = Event->track_startz_->at( p );
    float nu_vx = Event->nu_vx_;
    float vtx_disp = std::sqrt(pow(Event->nu_vx_-startx,2)+pow(Event->nu_vy_-starty,2)+pow(Event->nu_vz_-startz,2));
    if ( vtx_disp > 10.0 ) continue;

    if ( track_length > lead_p_track_length ) {
      lead_p_track_length = track_length;
      lead_p_index = p;
    }
  }

  // If the leading proton track length changed from its initial
  // value, then we found one. Set the index appropriately.
  if ( lead_p_track_length != LOW_FLOAT ) lead_p_candidate_idx_ = lead_p_index;
  // Otherwise, set the index to BOGUS_INDEX
  else lead_p_candidate_idx_ = BOGUS_INDEX;

  // Check the range-based reco momentum for the leading proton candidate
  float lead_p_KE = Event->track_kinetic_energy_p_->at( lead_p_candidate_idx_ );
  float range_mom_lead_p = real_sqrt( lead_p_KE*lead_p_KE
    + 2.*PROTON_MASS*lead_p_KE );
  // Readout direction variables for lead proton track for the "lead proton only" kinematic imbalance picture
  float lead_p_dirx = Event->track_dirx_->at( lead_p_candidate_idx_ );
  float lead_p_diry = Event->track_diry_->at( lead_p_candidate_idx_ );
  float lead_p_dirz = Event->track_dirz_->at( lead_p_candidate_idx_ );
  if (std::count(secondary_proton_index.begin(),secondary_proton_index.end(),lead_p_candidate_idx_) == 0)
  {
    event_momentum_leading_p.SetX(event_momentum_leading_p.X()+(lead_p_dirx*range_mom_lead_p));
    event_momentum_leading_p.SetY(event_momentum_leading_p.Y()+(lead_p_diry*range_mom_lead_p));
    event_momentum_leading_p.SetZ(event_momentum_leading_p.Z()+(lead_p_dirz*range_mom_lead_p));
    event_calo_energy_leading_p += lead_p_KE - PROTON_MASS + .0309; //Proton KE + Proton binding energy to argon
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  if ( range_mom_lead_p >= LEAD_P_MIN_MOM_CUT
    && range_mom_lead_p <= LEAD_P_MAX_MOM_CUT )
  {
    sel_lead_p_passed_mom_cuts_ = true;
  }

  //Save missing momentum for the event
  event_missing_trans_px = -1*event_momentum.X();
  event_missing_trans_py = -1*event_momentum.Y();
  event_missing_p_L = event_calo_energy - event_momentum.Z();
  event_missing_trans_p_mag = sqrt(event_missing_trans_px*event_missing_trans_px+event_missing_trans_py*event_missing_trans_py);
  event_missing_p_mag = std::sqrt(event_missing_p_L*event_missing_p_L+event_missing_trans_p_mag*event_missing_trans_p_mag);
  //Save missing momentum for the event in the leading proton picture
  event_missing_trans_px_leading_p = -1*event_momentum_leading_p.X();
  event_missing_trans_py_leading_p = -1*event_momentum_leading_p.Y();
  event_missing_p_L_leading_p = event_calo_energy_leading_p - event_momentum_leading_p.Z();
  event_missing_trans_p_mag_leading_p = sqrt(event_missing_trans_px_leading_p*event_missing_trans_px_leading_p+event_missing_trans_py_leading_p*event_missing_trans_py_leading_p);
  event_missing_p_mag_leading_p = std::sqrt(event_missing_p_L_leading_p*event_missing_p_L_leading_p+event_missing_trans_p_mag_leading_p*event_missing_trans_p_mag_leading_p);

  sel_num_secondary_proton_cands = num_secondary_proton_cands;
  if(num_secondary_proton_cands != 0){
    for(int i = 0; i < num_secondary_proton_cands; i++){
      int proton_index = secondary_proton_index.at(i);
      secondary_proton_trk_candidate_indices[i] = proton_index;
      // Use precomputed maps to look up values for each proton
      if (pfp_proximity_map.count(proton_index))
                    secondary_proton_proximity[i] = pfp_proximity_map[proton_index];

      if (pfp_direction_map.count(proton_index))
                    secondary_proton_direction[i] = pfp_direction_map[proton_index];

      if (pfp_vtx_disp_map.count(proton_index))
                    secondary_proton_vtx_disp[i] = pfp_vtx_disp_map[proton_index];

      if (pfp_len_map.count(proton_index)){
                    secondary_proton_len[i] = pfp_len_map[proton_index];}

      if (pfp_cos_theta_map.count(proton_index)){
                    secondary_proton_costheta[i] = pfp_cos_theta_map[proton_index];}

      if (pfp_theta_map.count(proton_index))
                    secondary_proton_theta[i] = pfp_theta_map[proton_index];

      if (pfp_phi_map.count(proton_index))
                    secondary_proton_phi[i] = pfp_phi_map[proton_index];

      if (pfp_start_x_map.count(proton_index))
                    secondary_proton_startx[i] = pfp_start_x_map[proton_index];

      if (pfp_start_y_map.count(proton_index))
                    secondary_proton_starty[i] = pfp_start_y_map[proton_index];

      if (pfp_start_z_map.count(proton_index))
                    secondary_proton_startz[i] = pfp_start_z_map[proton_index];

      if (pfp_end_x_map.count(proton_index))
                    secondary_proton_endx[i] = pfp_end_x_map[proton_index];

      if (pfp_end_y_map.count(proton_index))
                    secondary_proton_endy[i] = pfp_end_y_map[proton_index];

      if (pfp_end_z_map.count(proton_index))
                    secondary_proton_endz[i] = pfp_end_z_map[proton_index];
    }
  }
  //Calculate cosine of alpha and save
  TVector3 neutron_cand_position(-999.0,-999.0,-999.0);
  if (num_secondary_proton_cands == 1)
  {
    neutron_cand_position.SetX(secondary_proton_startx[0]-Event->nu_vx_);
    neutron_cand_position.SetY(secondary_proton_starty[0]-Event->nu_vy_);
    neutron_cand_position.SetZ(secondary_proton_startz[0]-Event->nu_vz_);
    
    float neutron_cand_trans_position_mag = sqrt(pow(neutron_cand_position.X(),2)+pow(neutron_cand_position.Y(),2));
    float neutron_cand_position_mag = std::sqrt(pow(neutron_cand_position.X(),2)+pow(neutron_cand_position.Y(),2)+pow(neutron_cand_position.Z(),2));
    sel_cos_alpha_trans = (neutron_cand_position.X()*event_missing_trans_px+neutron_cand_position.Y()*event_missing_trans_py)/(neutron_cand_trans_position_mag*event_missing_trans_p_mag);
    sel_cos_alpha_3D = (neutron_cand_position.X()*event_missing_trans_px+neutron_cand_position.Y()*event_missing_trans_py+neutron_cand_position.Z()*event_missing_p_L)/(neutron_cand_position_mag*event_missing_p_mag);
    sel_alpha_trans = acos(sel_cos_alpha_trans)*(180/M_PI);
    sel_alpha_3D = acos(sel_cos_alpha_3D)*(180/M_PI);
    // Save kinematic imbalance variables for the leading proton picture
    sel_cos_alpha_trans_leading_p = (neutron_cand_position.X()*event_missing_trans_px_leading_p+neutron_cand_position.Y()*event_missing_trans_py_leading_p)/(neutron_cand_trans_position_mag*event_missing_trans_p_mag_leading_p);
    sel_cos_alpha_3D_leading_p = (neutron_cand_position.X()*event_missing_trans_px_leading_p+neutron_cand_position.Y()*event_missing_trans_py_leading_p+neutron_cand_position.Z()*event_missing_p_L_leading_p)/(neutron_cand_position_mag*event_missing_p_mag_leading_p);
    sel_alpha_trans_leading_p = acos(sel_cos_alpha_trans_leading_p)*(180/M_PI);
    sel_alpha_3D_leading_p = acos(sel_cos_alpha_3D_leading_p)*(180/M_PI);
  }


  // All right, we've applied all selection cuts. Set the flag that indicates
  // whether all were passed (and thus the event is selected as a CCNp0pi
  // candidate)
  bool sel_CCNp0piXn_ = sel_nu_mu_cc_ && sel_no_reco_showers_
    && sel_muon_passed_mom_cuts_ && sel_muon_contained_ && sel_muon_quality_ok_
    && sel_has_p_candidate_ && sel_passed_proton_pid_cut_
    && sel_protons_contained_ && sel_lead_p_passed_mom_cuts_
    && sel_has_secondary_proton_candidate_;

  return sel_CCNp0piXn_;
}

int CC1muNp0piNn::CategorizeEvent(AnalysisEvent* Event) {
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

  if ( IsEventMCSignal() ) { //CCnumu >=1 proton, 0 pions, >= 1 neutron
    if ( sig_num_neutron_above_100MeV == 1 ) { // 1 neutron
      if (sig_nProtons_in_Momentum_range == 1) {return kNuMuCC1p0pi1n;} // 1 proton
      else {return kNuMuCCNp0pi1n;} //N protons
    }
    else{ // N neutrons
      if (sig_nProtons_in_Momentum_range == 1) {return kNuMuCC1p0piNn;} // 1 proton
      else {return kNuMuCCNp0piNn;} // N protons
    }
  }
  else if (!sig_mc_no_fs_pi0_ || !sig_mc_no_charged_pi_above_threshold_) { // N pions
    return kNuMuCCNpi;
  }
  else if (sig_num_neutron_above_100MeV == 0){ // 0 neutrons
    if (sig_nProtons_in_Momentum_range == 0){return kNuMuCC0p0pi0n;} // 0 protons
    else if (sig_nProtons_in_Momentum_range == 1){return kNuMuCC1p0pi0n;} // 1 proton
    else {return kNuMuCCNp0pi0n;} // N protons
  }
  else if (sig_num_neutron_above_100MeV == 1 && sig_nProtons_in_Momentum_range == 0) {return kNuMuCC0p0pi1n;} //1 neutron 0 protons
  else if (sig_num_neutron_above_100MeV > 1 && sig_nProtons_in_Momentum_range == 0) {return kNuMuCC0p0piNn;} //N neutrons 0 protons
  else {return kNuMuCCOther;}
}

void CC1muNp0piNn::DefineOutputBranches() {
  SetBranch(&sig_isNuMu_,"mc_is_numu",kBool);
  SetBranch(&sig_inFV_,"mc_vertex_in_FV",kBool);
  SetBranch(&sig_leadProtonMomInRange_,"mc_lead_p_in_range",kBool);
  SetBranch(&sig_muonInMomRange_,"mc_muon_in_mom_range",kBool);
  SetBranch(&sig_noFSMesons_,"mc_no_FS_mesons",kBool);
  SetBranch(&sig_mc_no_charged_pi_above_threshold_,"mc_no_charged_pions_above_thres",kBool);
  SetBranch(&sig_mc_no_fs_pi0_,"mc_no_pi0s",kBool);
  SetBranch(&sig_nProtons_in_Momentum_range,"nProtons_in_Momentum_range",kInteger);

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

  SetBranch(&lead_p_candidate_idx_,"lead_p_candidate_idx",kInteger);
  SetBranch(&muon_candidate_idx_,"muon_candidate_idx",kInteger);

  SetBranch(&delta_pT_,"reco_delta_pT",kDouble);
  SetBranch(&delta_phiT_,"reco_delta_phiT",kDouble);
  SetBranch(&delta_alphaT_,"reco_delta_alphaT",kDouble);
  SetBranch(&delta_pL_,"reco_delta_pL",kDouble);
  SetBranch(&pn_,"reco_pn",kDouble);
  SetBranch(&delta_pTx_,"reco_delta_pTx",kDouble);
  SetBranch(&delta_pTy_,"reco_delta_pTy",kDouble);
  SetBranch(&theta_mu_p_,"reco_theta_mu_p",kDouble);
  SetBranch(p3mu,"reco_p3_mu",kTVector);
  SetBranch(p3p,"reco_p3_lead_p",kTVector);
  SetBranch(p3_p_vec_,"reco_p3_p_vec",kSTDVector);

  SetBranch(&mc_delta_pT_,"true_delta_pT",kDouble);
  SetBranch(&mc_delta_phiT_,"true_delta_phiT",kDouble);
  SetBranch(&mc_delta_alphaT_,"true_delta_alphaT",kDouble);
  SetBranch(&mc_delta_pL_,"true_delta_pL",kDouble);
  SetBranch(&mc_pn_,"true_pn",kDouble);
  SetBranch(&mc_delta_pTx_,"true_delta_pTx",kDouble);
  SetBranch(&mc_delta_pTy_,"true_delta_pTy",kDouble);
  SetBranch(&mc_theta_mu_p_,"true_theta_mu_p",kDouble);
  SetBranch(mc_p3mu,"true_p3_mu",kTVector);
  SetBranch(mc_p3p,"true_p3_lead_p",kTVector);
  SetBranch(mc_p3_p_vec_,"true_p3_p_vec",kSTDVector);

  //Begin Burke Additions
  SetBranch(weight_neutron_argon_xsec, "weight_neutron_argon_xsec", kSTDVector);
  SetBranch(weight_neutron_reint, "weight_neutron_reint", kSTDVector); 
  SetBranch(neutron_track_weights, "neutron_track_weights", kSTDVector);
  SetBranch(&sig_mc_FS_neutron_,"mc_FS_neutron",kBool);
  SetBranch(&sig_num_neutron_above_100MeV,"nNeutrons_above_100MeV",kInteger);
  SetBranch(&sel_has_secondary_proton_candidate_,"has_secondary_proton_cand",kBool);
  Tree->Branch("CC1muNp0piNn_num_secondary_proton_cands", &sel_num_secondary_proton_cands, "CC1muNp0piNn_num_secondary_proton_cands/I");
  Tree->Branch("CC1muNp0piNn_secondary_proton_trk_candidate_indices", secondary_proton_trk_candidate_indices, "CC1muNp0piNn_secondary_proton_trk_candidate_indices[CC1muNp0piNn_num_secondary_proton_cands]/I");
  Tree->Branch("CC1muNp0piNn_secondary_proton_proximity", secondary_proton_proximity, "CC1muNp0piNn_secondary_proton_proximity[CC1muNp0piNn_num_secondary_proton_cands]/F");
  Tree->Branch("CC1muNp0piNn_secondary_proton_direction", secondary_proton_direction, "CC1muNp0piNn_secondary_proton_direction[CC1muNp0piNn_num_secondary_proton_cands]/F");
  Tree->Branch("CC1muNp0piNn_secondary_proton_costheta", secondary_proton_costheta,"CC1muNp0piNn_secondary_proton_costheta[CC1muNp0piNn_num_secondary_proton_cands]/F");
  Tree->Branch("CC1muNp0piNn_secondary_proton_len", secondary_proton_len,"CC1muNp0piNn_secondary_proton_len[CC1muNp0piNn_num_secondary_proton_cands]/F");
  Tree->Branch("CC1muNp0piNn_secondary_proton_vtx_disp", secondary_proton_vtx_disp,"CC1muNp0piNn_secondary_proton_vtx_disp[CC1muNp0piNn_num_secondary_proton_cands]/F");
  Tree->Branch("CC1muNp0piNn_secondary_proton_theta", secondary_proton_theta,"CC1muNp0piNn_secondary_proton_theta[CC1muNp0piNn_num_secondary_proton_cands]/F");
  Tree->Branch("CC1muNp0piNn_secondary_proton_phi", secondary_proton_phi,"CC1muNp0piNn_secondary_proton_phi[CC1muNp0piNn_num_secondary_proton_cands]/F");
  Tree->Branch("CC1muNp0piNn_secondary_proton_startx", secondary_proton_startx,"CC1muNp0piNn_secondary_proton_startx[CC1muNp0piNn_num_secondary_proton_cands]/F");
  Tree->Branch("CC1muNp0piNn_secondary_proton_starty", secondary_proton_starty,"CC1muNp0piNn_secondary_proton_starty[CC1muNp0piNn_num_secondary_proton_cands]/F");
  Tree->Branch("CC1muNp0piNn_secondary_proton_startz", secondary_proton_startz,"CC1muNp0piNn_secondary_proton_startz[CC1muNp0piNn_num_secondary_proton_cands]/F");
  Tree->Branch("CC1muNp0piNn_secondary_proton_endx", secondary_proton_endx,"CC1muNp0piNn_secondary_proton_endx[CC1muNp0piNn_num_secondary_proton_cands]/F");
  Tree->Branch("CC1muNp0piNn_secondary_proton_endy", secondary_proton_endy,"CC1muNp0piNn_secondary_proton_endy[CC1muNp0piNn_num_secondary_proton_cands]/F");
  Tree->Branch("CC1muNp0piNn_secondary_proton_endz", secondary_proton_endz,"CC1muNp0piNn_secondary_proton_endz[CC1muNp0piNn_num_secondary_proton_cands]/F");

  SetBranch(&event_missing_trans_px,"missing_trans_px",kDouble);
  SetBranch(&event_missing_trans_py,"missing_trans_py",kDouble);
  SetBranch(&event_missing_trans_p_mag,"missing_trans_p_mag",kDouble);
  SetBranch(&event_missing_p_L,"event_missing_p_L",kDouble);
  SetBranch(&event_missing_p_mag,"event_missing_p_mag",kDouble);
  SetBranch(&sel_cos_alpha_trans,"cos_alpha",kDouble);
  SetBranch(&sel_alpha_trans,"alpha",kDouble);
  SetBranch(&sel_cos_alpha_3D,"cos_alpha_3D",kDouble);
  SetBranch(&sel_alpha_3D,"alpha_3D",kDouble);

  SetBranch(&event_missing_trans_px_leading_p,"missing_trans_px_leading_p",kDouble);
  SetBranch(&event_missing_trans_py_leading_p,"missing_trans_py_leading_p",kDouble);
  SetBranch(&event_missing_trans_p_mag_leading_p,"missing_trans_p_mag_leading_p",kDouble);
  SetBranch(&event_missing_p_L_leading_p,"event_missing_p_L_leading_p",kDouble);
  SetBranch(&event_missing_p_mag_leading_p,"event_missing_p_mag_leading_p",kDouble);
  SetBranch(&sel_cos_alpha_trans_leading_p,"cos_alpha_leading_p",kDouble);
  SetBranch(&sel_alpha_trans_leading_p,"alpha_leading_p",kDouble);
  SetBranch(&sel_cos_alpha_3D_leading_p,"cos_alpha_3D_leading_p",kDouble);
  SetBranch(&sel_alpha_3D_leading_p,"alpha_3D_leading_p",kDouble);
}

void CC1muNp0piNn::DefineCategoryMap() {
  // Use the shared category map for 1p/2p/Np/Xp
  categ_map_ = CC1muXpXn_MAP;
}
