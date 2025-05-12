// XSecAnalyzer includes
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/Selections/CC1mu2p0pi.hh"

namespace {
  constexpr double CC2P_MIN_P_MUON = 0.1; // GeV
  constexpr double CC2P_MAX_P_MUON = 1.2; // GeV

  constexpr double CC2P_MIN_P_PROTON = 0.3; // GeV
  constexpr double CC2P_MAX_P_PROTON = 1.0; // GeV

  constexpr double CC2P_MIN_P_PI_PLUS = 0.065; // GeV

}

CC1mu2p0pi::CC1mu2p0pi() : SelectionBase( "CC1mu2p0pi" ) {
  ki_calc_type_ = KICalculator::kOpt1;
  this->define_fv( 10., 246.35, -106.5, 106.5, 10., 1026.8 );
}

std::string CC1mu2p0pi::categorize_event( AnalysisEvent& ev ) {

  const auto& in = ev.in();
  auto& out = ev.out();

  float mc_nu_vx, mc_nu_vy, mc_nu_vz;
  int ccnc, nu_pdg, interaction_type;
  std::vector< int >* nu_daughter_pdg;
  std::vector< float > *nu_daughter_energy, *nu_daughter_px,
    *nu_daughter_py, *nu_daughter_pz;
  in.at( "true_nu_vtx_x" ) >> mc_nu_vx;
  in.at( "true_nu_vtx_y" ) >> mc_nu_vy;
  in.at( "true_nu_vtx_z" ) >> mc_nu_vz;
  in.at( "ccnc" ) >> ccnc;
  in.at( "nu_pdg" ) >> nu_pdg;
  in.at( "interaction" ) >> interaction_type;
  in.at( "mc_pdg" ) >> nu_daughter_pdg;
  in.at( "mc_E" ) >> nu_daughter_energy;
  in.at( "mc_px" ) >> nu_daughter_px;
  in.at( "mc_py" ) >> nu_daughter_py;
  in.at( "mc_pz" ) >> nu_daughter_pz;

  // Signal definition criteria taken from
  // https://tinyurl.com/uboone-samantha-cc2p-sig-def
  int num_thresh_muon = 0;
  int num_thresh_proton = 0;
  int num_thresh_pion0 = 0;
  int num_thresh_charged_pion = 0;

  std::vector< double > true_p_mom_vec;
  std::vector< int > true_p_indices;

  int true_mu_idx = BOGUS_INDEX;

  for ( size_t p = 0u; p < nu_daughter_pdg->size(); ++p ) {
    int pdg = nu_daughter_pdg->at( p );
    float energy = nu_daughter_energy->at( p );
    if ( std::abs( pdg ) == MUON ) {
      double mom = real_sqrt( energy*energy - MUON_MASS*MUON_MASS );
      if ( mom > CC2P_MIN_P_MUON && mom < CC2P_MAX_P_MUON ) {
        ++num_thresh_muon;
        true_mu_idx = p;
      }
    }
    else if ( pdg == PROTON ) {
      double mom = real_sqrt( energy*energy - PROTON_MASS*PROTON_MASS );
      if ( mom > CC2P_MIN_P_PROTON && mom < CC2P_MAX_P_PROTON ) {
        ++num_thresh_proton;
        true_p_mom_vec.push_back( mom );
        true_p_indices.push_back( p );
      }
    }
    else if ( pdg == PI_ZERO ) {
      ++num_thresh_pion0;
    }
    else if ( std::abs( pdg ) == PI_PLUS ) {
      double mom = real_sqrt( energy*energy - PI_PLUS_MASS*PI_PLUS_MASS );
      if ( mom > CC2P_MIN_P_PI_PLUS ) {
        ++num_thresh_charged_pion;
      }
    }
  }

  // Get the indices of the true leading and recoil protons if both are present
  int true_lead_p_idx = BOGUS_INDEX, true_recoil_p_idx = BOGUS_INDEX;

  if ( true_p_mom_vec.size() == 2 ) {
    if ( true_p_mom_vec[0] > true_p_mom_vec[1] ) {
      true_lead_p_idx = true_p_indices[0];
      true_recoil_p_idx = true_p_indices[1];
    }
    else {
      true_lead_p_idx = true_p_indices[1];
      true_recoil_p_idx = true_p_indices[0];
    }
  }

  // ===========================================================
  // Calculate the booleans related to the different signal cuts

  // NOTE: Original CC2p code from Samantha incorrectly used the true vertices
  // with a correction for space charge effects (SCE) in the signal definition.
  // See https://tinyurl.com/uboone-cc2p-samantha-sce-bug for details.
  // Samantha's version used the "true_nu_vtx_sce_x" (and similar for y and z)
  // branches from the PeLEE ntuples. We use the true vertex without this
  // adjustment here.
  bool sig_true_vertex_in_fv = this->get_fv().is_inside( mc_nu_vx,
    mc_nu_vy, mc_nu_vz );
  bool sig_is_cc = ( ccnc == CHARGED_CURRENT);
  bool sig_is_numu = ( nu_pdg == MUON_NEUTRINO );
  bool sig_two_protons_above_thresh = ( num_thresh_proton == 2 );
  bool sig_one_muon_above_thresh = ( num_thresh_muon == 1 );
  bool sig_no_pions = ( num_thresh_pion0 == 0 )
    && ( num_thresh_charged_pion == 0 );

  bool is_signal = sig_is_cc && sig_is_numu && sig_two_protons_above_thresh
    && sig_one_muon_above_thresh && sig_no_pions && sig_true_vertex_in_fv;

  ///////// categorize_event
  bool mc_vtx_in_fv = this->get_fv().is_inside( mc_nu_vx, mc_nu_vy, mc_nu_vz );
  bool is_nc = ( ccnc == NEUTRAL_CURRENT );

  // Boolean which is the same as the signal criteria except that any number
  // of final-state protons is allowed
  bool is_ccXp0pi = sig_is_cc && sig_is_numu && sig_true_vertex_in_fv
    && sig_no_pions && sig_one_muon_above_thresh;

  // Compute true observables to store in the output tree
  TVector3 p3_mu, p3_lead_p, p3_recoil_p;

  double true_Pt = BOGUS, true_Ptx = BOGUS, true_Pty = BOGUS, true_PL = BOGUS,
    true_Pn = BOGUS, true_PnPerp = BOGUS, true_PnPerpx = BOGUS,
    true_PnPerpy = BOGUS, true_PnPar = BOGUS, true_DeltaAlphaT = BOGUS,
    true_DeltaAlpha3Dq = BOGUS, true_DeltaAlpha3DMu = BOGUS,
    true_DeltaPhiT = BOGUS, true_DeltaPhi3D = BOGUS, true_ECal = BOGUS,
    true_EQE = BOGUS, true_Q2 = BOGUS, true_A = BOGUS, true_EMiss = BOGUS,
    true_kMiss = BOGUS, true_PMiss = BOGUS, true_PMissMinus = BOGUS;

  if ( sig_two_protons_above_thresh && sig_one_muon_above_thresh
    && sig_no_pions )
  {
    double px_mu = nu_daughter_px->at( true_mu_idx );
    double py_mu = nu_daughter_py->at( true_mu_idx );
    double pz_mu = nu_daughter_pz->at( true_mu_idx );
    p3_mu = TVector3( px_mu, py_mu, pz_mu );

    double px_lead_p = nu_daughter_px->at( true_lead_p_idx );
    double py_lead_p = nu_daughter_py->at( true_lead_p_idx );
    double pz_lead_p = nu_daughter_pz->at( true_lead_p_idx );
    p3_lead_p = TVector3( px_lead_p, py_lead_p, pz_lead_p );

    double E_lead_p = real_sqrt( p3_lead_p.Mag2() + PROTON_MASS*PROTON_MASS );
    TLorentzVector p4_lead_p( p3_lead_p, E_lead_p );

    double px_recoil_p = nu_daughter_px->at( true_recoil_p_idx );
    double py_recoil_p = nu_daughter_py->at( true_recoil_p_idx );
    double pz_recoil_p = nu_daughter_pz->at( true_recoil_p_idx );
    p3_recoil_p = TVector3( px_recoil_p, py_recoil_p, pz_recoil_p );

    double E_recoil_p = real_sqrt( p3_recoil_p.Mag2()
      + PROTON_MASS*PROTON_MASS );
    TLorentzVector p4_recoil_p( p3_recoil_p, E_recoil_p );

    TVector3 p3_p_sum = p3_lead_p + p3_recoil_p;
    TLorentzVector p4_p_sum = p4_lead_p + p4_recoil_p;

    // Hadronic invariant mass of the proton pair
    double W = p4_p_sum.M();

    KICalculator ki_calc( p3_mu, p3_p_sum, ki_calc_type_, MUON_MASS, W );

    true_Pt = ki_calc.pT();
    true_Ptx = ki_calc.pTx();
    true_Pty = ki_calc.pTy();
    true_PL = ki_calc.pL();
    true_Pn = ki_calc.pn();
    true_PnPerp = ki_calc.pn_perp();
    true_PnPerpx = ki_calc.pn_perpx();
    true_PnPerpy = ki_calc.pn_perpy();
    true_PnPar = ki_calc.pn_par();
    true_DeltaAlphaT = ki_calc.delta_alphaT();
    true_DeltaAlpha3Dq = ki_calc.delta_alpha3D_q();
    true_DeltaAlpha3DMu = ki_calc.delta_alpha3D_mu();
    true_DeltaPhiT = ki_calc.delta_phiT();
    true_DeltaPhi3D = ki_calc.delta_phi3D();
    true_ECal = ki_calc.Ecal();
    true_EQE = ki_calc.E_QE();
    true_Q2 = ki_calc.Q2();
    true_A = ki_calc.A();
    true_EMiss = ki_calc.E_miss();
    true_kMiss = ki_calc.k_miss();
    true_PMiss = ki_calc.p_miss();
    true_PMissMinus = ki_calc.p_miss_minus();
  }

  // We're done. Store the results in the output tree and return the
  // label for the current event category
  out[ "sig_truevertex_in_fv" ] = sig_true_vertex_in_fv;
  out[ "sig_ccnc" ] = sig_is_cc;
  out[ "sig_is_numu" ] = sig_is_numu;
  out[ "sig_two_protons_above_thresh" ] = sig_two_protons_above_thresh;
  out[ "sig_one_muon_above_thresh" ] = sig_one_muon_above_thresh;
  out[ "sig_no_pions" ] = sig_no_pions;

  out[ "mc_n_threshold_muon" ] = num_thresh_muon;
  out[ "mc_n_threshold_proton" ] = num_thresh_proton;
  out[ "mc_n_threshold_pion0" ] = num_thresh_pion0;
  out[ "mc_n_threshold_pionpm" ] = num_thresh_charged_pion;

  out[ "True_Pt" ] = true_Pt;
  out[ "True_Ptx" ] = true_Ptx;
  out[ "True_Pty" ] = true_Pty;
  out[ "True_PL" ] = true_PL;
  out[ "True_Pn" ] = true_Pn;
  out[ "True_PnPerp" ] = true_PnPerp;
  out[ "True_PnPerpx" ] = true_PnPerpx;
  out[ "True_PnPerpy" ] = true_PnPerpy;
  out[ "True_PnPar" ] = true_PnPar;
  out[ "True_DeltaAlphaT" ] = true_DeltaAlphaT;
  out[ "True_DeltaAlpha3Dq" ] = true_DeltaAlpha3Dq;
  out[ "True_DeltaAlpha3DMu" ] = true_DeltaAlpha3DMu;
  out[ "True_DeltaPhiT" ] = true_DeltaPhiT;
  out[ "True_DeltaPhi3D" ] = true_DeltaPhi3D;
  out[ "True_ECal" ] = true_ECal;
  out[ "True_EQE" ] = true_EQE;
  out[ "True_Q2" ] = true_Q2;
  out[ "True_A" ] = true_A;
  out[ "True_EMiss" ] = true_EMiss;
  out[ "True_kMiss" ] = true_kMiss;
  out[ "True_PMiss" ] = true_PMiss;
  out[ "True_PMissMinus" ] = true_PMissMinus;

  // Categories used here are based on https://arxiv.org/abs/2211.03734
  if ( !sig_true_vertex_in_fv ) return "OOFV";
  else if ( is_signal ) {
    if ( interaction_type == QE_INTERACTION ) return "CC1#mu2p QE";
    else if ( interaction_type == MEC_INTERACTION ) return "CC1#mu2p MEC";
    else if ( interaction_type == RES_INTERACTION ) return "CC1#mu2p RES";
    else return "CC1#mu2p Other";
  }
  else if ( !is_signal && is_ccXp0pi ) {
    if ( num_thresh_proton == 0 ) return "CC1#mu0p";
    else if ( num_thresh_proton == 1 ) return "CC1#mu1p";
    else return "CC1#mu(M>2)p";
  }
  else if ( sig_is_numu && sig_is_cc && sig_one_muon_above_thresh ) {
    if ( !sig_no_pions ) return "CC1#muN#pi";
    else return "CC1#mu Other";
  }
  else return "Other";
}

bool CC1mu2p0pi::is_selected( AnalysisEvent& ev ) {

  // Get access to the input and output ntuples
  const auto& in = ev.in();
  auto& out = ev.out();

  // Reconstructed neutrino vertex in the fiducial volume?
  float x, y, z;
  in.at( "reco_nu_vtx_sce_x" ) >> x;
  in.at( "reco_nu_vtx_sce_y" ) >> y;
  in.at( "reco_nu_vtx_sce_z" ) >> z;

  bool sel_reco_vertex_in_FV = this->get_fv().is_inside( x, y, z );

  int num_pf_particles;
  std::vector< unsigned int >* pfp_generation;
  std::vector< float > *track_llr_pid_score, *pfp_track_score,
    *track_start_distance, *track_range_mom_mu, *track_kinetic_energy_p,
    *track_startx, *track_starty, *track_startz, *track_endx, *track_endy,
    *track_endz, *track_dirx, *track_diry, *track_dirz;

  in.at( "n_pfps" ) >> num_pf_particles;
  in.at( "pfp_generation_v" ) >> pfp_generation;
  in.at( "trk_llr_pid_score_v" ) >> track_llr_pid_score;
  in.at( "trk_score_v" ) >> pfp_track_score;
  in.at( "trk_distance_v" ) >> track_start_distance;
  in.at( "trk_range_muon_mom_v" ) >> track_range_mom_mu;
  in.at( "trk_energy_proton_v" ) >> track_kinetic_energy_p;

  in.at( "trk_sce_start_x_v" ) >> track_startx;
  in.at( "trk_sce_start_y_v" ) >> track_starty;
  in.at( "trk_sce_start_z_v" ) >> track_startz;

  in.at( "trk_sce_end_x_v" ) >> track_endx;
  in.at( "trk_sce_end_y_v" ) >> track_endy;
  in.at( "trk_sce_end_z_v" ) >> track_endz;

  in.at( "trk_dir_x_v" ) >> track_dirx;
  in.at( "trk_dir_y_v" ) >> track_diry;
  in.at( "trk_dir_z_v" ) >> track_dirz;

  // Events with anything other than one muon candidate are excluded
  int n_muons = 0;
  int chosen_index = 0;

  for ( int p = 0; p < num_pf_particles; ++p ) {

    // Only direct neutrino daughters (generation == 2) will be considered as
    // possible muon candidates
    unsigned int generation = pfp_generation->at( p );
    if ( generation != 2u ) continue;

    float pid_score = track_llr_pid_score->at( p );
    if ( pid_score >= MUON_PID_CUT && pid_score > -1. && pid_score < 1. ) {
      ++n_muons;

      // Gets overwritten if there are multiple muon candidates, but that's
      // fine because we require exactly one muon candidate
      chosen_index = p;
    }
  }

  int muon_candidate_idx = BOGUS_INDEX;
  bool sel_has_muon_candidate = false;
  if ( n_muons == 1 ) {
    sel_has_muon_candidate = true;
    muon_candidate_idx = chosen_index;
  }

  // Does the event pass the numu CC pre-selection?
  bool sel_nu_mu_cc = sel_reco_vertex_in_FV && sel_has_muon_candidate;

  // Require exactly 3 PFParticles
  bool sel_npfps_eq_3 = ( num_pf_particles == 3 );

  // Require 3 tracks (track_score > 0.8 [MUON_TRACK_SCORE_CUT]) whose
  // "vertex distance attachment is less than 4 cm"
  int n_tracks = 0;

  TVector3 nu_vtx_reco( x, y, z );
  for ( int i = 0; i < num_pf_particles; ++i ) {
    float track_score = pfp_track_score->at( i );
    float track_distance = track_start_distance->at( i );
    float track_pid = track_llr_pid_score->at( i );

    // leading track end reco
    TVector3 track_end( track_endx->at(i), track_endy->at(i),
      track_endz->at(i) );
    track_end -= nu_vtx_reco;

    double track_end_distance = track_end.Mag();

    // DB: I think this cut has a different logic than Steven's
    if ( track_score >= MUON_TRACK_SCORE_CUT
      && ( track_distance <= MUON_VTX_DISTANCE_CUT
        || track_end_distance <= MUON_VTX_DISTANCE_CUT ) )
    {
      ++n_tracks;
    }
  }

  bool sel_ntracks_eq_3 = ( n_tracks == 3 );

  // Make sure that there are two reconstructed protons
  int n_protons = 0;
  for ( int i = 0; i < num_pf_particles; ++i ) {
    float track_pid = track_llr_pid_score->at( i );
    if ( track_pid < DEFAULT_PROTON_PID_CUT && track_pid < 1.
      && track_pid > -1.)
    {
      ++n_protons;
    }
  }

  // Note that there was a bug for this cut in Samantha's original code:
  // https://tinyurl.com/uboone-samantha-cc2p-cut-bug
  bool sel_correctparticles = ( sel_has_muon_candidate && n_protons == 2 );

  // Now we check that all 3 PFParticles are contained (both start and end)

  // DB: Initially Samantha precalculates which index corresponds to the muon,
  // leading proton, and recoil proton, but then she applies the same
  // containment cut to all three. The function used is given here:
  // https://tinyurl.com/uboone-cc2p-samantha-helper.  The start of the track
  // requires a tighter FV cut by an additional 10 cm (i.e. the argument
  // variables), and the end of the vertex uses the FV cut provided (i.e. the
  // argument variables are 0 cm). I don't think there's any need to
  // precalculate the indices because we've already selected events which only
  // have 3 PFParticles within them (thus 1 muon and 2 protons), but this could
  // be wrong
  bool sel_contained_particles = true;

  FiducialVolume fv_no_border( 0., 256.35, -116.5, 116.5, 0., 1036.8 );

  for ( int i = 0; i < num_pf_particles; ++i ) {
    bool start_contained_i = this->get_fv().is_inside( track_startx->at( i ),
      track_starty->at( i ), track_startz->at( i ) );

    bool end_contained_i = fv_no_border.is_inside( track_endx->at( i ),
      track_endy->at( i ), track_endz->at( i ) );

    if ( !start_contained_i || !end_contained_i ) {
      sel_contained_particles = false;
      break;
    }
  }

  // Now ensure that the muon and proton candidates pass the momentum
  // threshold requirements of
  // 0.1 GeV <= p_mu <= 1.2 GeV
  // 0.3 GeV <= p_p <= 1.0 GeV

  std::vector< double > p_mom_vec;
  std::vector< int > p_indices;

  bool sel_momenta_ok = true;
  for ( int i = 0; i < num_pf_particles; ++i ) {
    if ( i == muon_candidate_idx ) {
      if ( track_range_mom_mu->at( i ) < MUON_P_MIN_MOM_CUT
        || track_range_mom_mu->at( i ) > MUON_P_MAX_MOM_CUT )
      {
        sel_momenta_ok = false;
      }
    }
    else {
      float E_p = track_kinetic_energy_p->at( i ) + PROTON_MASS;
      float mom_p = real_sqrt( E_p*E_p - PROTON_MASS*PROTON_MASS );
      if ( mom_p < PROTON_MIN_MOM_CUT || mom_p > PROTON_MAX_MOM_CUT ) {
        sel_momenta_ok = false;
      }
      p_mom_vec.push_back( mom_p );
      p_indices.push_back( i );
    }
  }

  int lead_p_idx = BOGUS_INDEX;
  int recoil_p_idx = BOGUS_INDEX;

  if ( p_mom_vec.size() == 2 ) {
    if ( p_mom_vec[ 0 ] > p_mom_vec[ 1 ] ) {
      lead_p_idx = p_indices[ 0 ];
      recoil_p_idx = p_indices[ 1 ];
    } else {
      lead_p_idx = p_indices[ 1 ];
      recoil_p_idx = p_indices[ 0 ];
    }
  }

  // Flag that indicates whether the event passed the full CC1mu2p0pi selection
  bool passed = sel_nu_mu_cc && sel_npfps_eq_3 && sel_ntracks_eq_3
    && sel_correctparticles && sel_contained_particles && sel_momenta_ok;

  // Compute some reco observables before finishing up
  TVector3 p3_mu, p3_lead_p, p3_recoil_p;

  double reco_Pt = BOGUS, reco_Ptx = BOGUS, reco_Pty = BOGUS, reco_PL = BOGUS,
    reco_Pn = BOGUS, reco_PnPerp = BOGUS, reco_PnPerpx = BOGUS,
    reco_PnPerpy = BOGUS, reco_PnPar = BOGUS, reco_DeltaAlphaT = BOGUS,
    reco_DeltaAlpha3Dq = BOGUS, reco_DeltaAlpha3DMu = BOGUS,
    reco_DeltaPhiT = BOGUS, reco_DeltaPhi3D = BOGUS, reco_ECal = BOGUS,
    reco_EQE = BOGUS, reco_Q2 = BOGUS, reco_A = BOGUS, reco_EMiss = BOGUS,
    reco_kMiss = BOGUS, reco_PMiss = BOGUS, reco_PMissMinus = BOGUS,
    reco_CosPlPr = BOGUS, reco_CosMuPsum = BOGUS;

  if ( lead_p_idx != BOGUS_INDEX && recoil_p_idx != BOGUS_INDEX
    && muon_candidate_idx != BOGUS_INDEX )
  {
    TVector3 x3_nu_vertex( x, y, z );

    float p_mu = track_range_mom_mu->at( muon_candidate_idx );

    float mu_start_dist = track_start_distance->at( muon_candidate_idx );

    TVector3 x3_mu_end( track_endx->at( muon_candidate_idx ),
      track_endy->at( muon_candidate_idx ),
      track_endz->at( muon_candidate_idx )
    );

    x3_mu_end -= x3_nu_vertex;
    float mu_end_dist = x3_mu_end.Mag();

    float mu_dir_x = track_dirx->at( muon_candidate_idx );
    float mu_dir_y = track_diry->at( muon_candidate_idx );
    float mu_dir_z = track_dirz->at( muon_candidate_idx );
    p3_mu = TVector3( mu_dir_x, mu_dir_y, mu_dir_z );
    p3_mu = p3_mu.Unit() * p_mu;

    if ( mu_start_dist > mu_end_dist ) {
      p3_mu *= -1.0;
    }

    float E_lead_p = track_kinetic_energy_p->at( lead_p_idx ) + PROTON_MASS;
    float p_lead_p = std::sqrt( E_lead_p*E_lead_p - PROTON_MASS*PROTON_MASS );

    float lead_p_start_dist = track_start_distance->at( lead_p_idx );
    TVector3 x3_lead_p_end_dist( track_endx->at( lead_p_idx ),
      track_endy->at( lead_p_idx ), track_endz->at( lead_p_idx ) );

    x3_lead_p_end_dist -= x3_nu_vertex;
    float lead_p_end_dist = x3_lead_p_end_dist.Mag();

    float lead_p_dir_x = track_dirx->at( lead_p_idx );
    float lead_p_dir_y = track_diry->at( lead_p_idx );
    float lead_p_dir_z = track_dirz->at( lead_p_idx );
    p3_lead_p = TVector3( lead_p_dir_x,
      lead_p_dir_y, lead_p_dir_z );
    p3_lead_p = p3_lead_p.Unit() * p_lead_p;

    if ( lead_p_start_dist > lead_p_end_dist ) {
      p3_lead_p *= -1.0;
    }

    float E_recoil_p = track_kinetic_energy_p->at( recoil_p_idx ) + PROTON_MASS;
    float p_recoil_p = std::sqrt( E_recoil_p*E_recoil_p
      - PROTON_MASS*PROTON_MASS );

    float recoil_p_start_dist = track_start_distance->at( recoil_p_idx );
    TVector3 x3_recoil_p_end_dist( track_endx->at( recoil_p_idx ),
      track_endy->at( recoil_p_idx ), track_endz->at( recoil_p_idx ) );

    x3_recoil_p_end_dist -= x3_nu_vertex;
    float recoil_p_end_dist = x3_recoil_p_end_dist.Mag();

    float recoil_p_dir_x = track_dirx->at( recoil_p_idx );
    float recoil_p_dir_y = track_diry->at( recoil_p_idx );
    float recoil_p_dir_z = track_dirz->at( recoil_p_idx );
    p3_recoil_p = TVector3( recoil_p_dir_x, recoil_p_dir_y,
      recoil_p_dir_z );
    p3_recoil_p = p3_recoil_p.Unit() * p_recoil_p;

    if ( recoil_p_start_dist > recoil_p_end_dist ) {
      p3_recoil_p *= -1.0;
    }

    TVector3 p3_p_sum = p3_lead_p + p3_recoil_p;
    float E_p_sum = E_lead_p + E_recoil_p;
    TLorentzVector p4_p_sum( p3_p_sum, E_p_sum );

    // Hadronic invariant mass of the proton pair
    double W = p4_p_sum.M();

    KICalculator ki_calc( p3_mu, p3_p_sum, ki_calc_type_, MUON_MASS, W );

    reco_CosPlPr = p3_lead_p.Angle( p3_recoil_p );
    reco_CosMuPsum = p3_mu.Angle( p3_p_sum );

    reco_Pt = ki_calc.pT();
    reco_Ptx = ki_calc.pTx();
    reco_Pty = ki_calc.pTy();
    reco_PL = ki_calc.pL();
    reco_Pn = ki_calc.pn();
    reco_PnPerp = ki_calc.pn_perp();
    reco_PnPerpx = ki_calc.pn_perpx();
    reco_PnPerpy = ki_calc.pn_perpy();
    reco_PnPar = ki_calc.pn_par();
    reco_DeltaAlphaT = ki_calc.delta_alphaT();
    reco_DeltaAlpha3Dq = ki_calc.delta_alpha3D_q();
    reco_DeltaAlpha3DMu = ki_calc.delta_alpha3D_mu();
    reco_DeltaPhiT = ki_calc.delta_phiT();
    reco_DeltaPhi3D = ki_calc.delta_phi3D();
    reco_ECal = ki_calc.Ecal();
    reco_EQE = ki_calc.E_QE();
    reco_Q2 = ki_calc.Q2();
    reco_A = ki_calc.A();
    reco_EMiss = ki_calc.E_miss();
    reco_kMiss = ki_calc.k_miss();
    reco_PMiss = ki_calc.p_miss();
    reco_PMissMinus = ki_calc.p_miss_minus();
  }

  out[ "sel_reco_vertex_in_FV" ] = sel_reco_vertex_in_FV;
  out[ "sel_has_muon_candidate" ] = sel_has_muon_candidate;
  out[ "sel_nu_mu_cc" ] = sel_nu_mu_cc;
  out[ "sel_npfps_eq_3" ] = sel_npfps_eq_3;
  out[ "sel_ntracks_eq_3" ] = sel_ntracks_eq_3;
  out[ "sel_containedparticles" ] = sel_contained_particles;
  out[ "sel_correctparticles" ] = sel_correctparticles;
  out[ "sel_momentum_threshold_passed" ] = sel_momenta_ok;

  out[ "lead_p_idx" ] = lead_p_idx;
  out[ "recoil_p_idx" ] = recoil_p_idx;

  out[ "Reco_CosPlPr" ] = reco_CosPlPr;
  out[ "Reco_CosMuPsum" ] = reco_CosMuPsum;
  out[ "Reco_Pt" ] = reco_Pt;
  out[ "Reco_Ptx" ] = reco_Ptx;
  out[ "Reco_Pty" ] = reco_Pty;
  out[ "Reco_PL" ] = reco_PL;
  out[ "Reco_Pn" ] = reco_Pn;
  out[ "Reco_PnPerp" ] = reco_PnPerp;
  out[ "Reco_PnPerpx" ] = reco_PnPerpx;
  out[ "Reco_PnPerpy" ] = reco_PnPerpy;
  out[ "Reco_PnPar" ] = reco_PnPar;
  out[ "Reco_DeltaAlphaT" ] = reco_DeltaAlphaT;
  out[ "Reco_DeltaAlpha3Dq" ] = reco_DeltaAlpha3Dq;
  out[ "Reco_DeltaAlpha3DMu" ] = reco_DeltaAlpha3DMu;
  out[ "Reco_DeltaPhiT" ] = reco_DeltaPhiT;
  out[ "Reco_DeltaPhi3D" ] = reco_DeltaPhi3D;
  out[ "Reco_ECal" ] = reco_ECal;
  out[ "Reco_EQE" ] = reco_EQE;
  out[ "Reco_Q2" ] = reco_Q2;
  out[ "Reco_A" ] = reco_A;
  out[ "Reco_EMiss" ] = reco_EMiss;
  out[ "Reco_kMiss" ] = reco_kMiss;
  out[ "Reco_PMiss" ] = reco_PMiss;
  out[ "Reco_PMissMinus" ] = reco_PMissMinus;

  return passed;
}
