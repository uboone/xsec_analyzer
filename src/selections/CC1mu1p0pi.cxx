// XSecAnalyzer includes
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/KICalculator.hh"
#include "XSecAnalyzer/Selections/CC1mu1p0pi.hh"
#include "XSecAnalyzer/TreeUtils.hh"

CC1mu1p0pi::CC1mu1p0pi() : SelectionBase( "CC1mu1p0pi" ) {
  ki_calc_type_ = KICalculator::kOpt1;
  this->define_fv( 10., 246., -105., 105., 10., 1026. );
}

std::string CC1mu1p0pi::categorize_event( AnalysisEvent& ev ) {

  // Get access to the input and output trees
  const auto& in = ev.in();
  auto& out = ev.out();

  // Retrieve some truth information from the PeLEE ntuple
  float mc_nu_vx, mc_nu_vy, mc_nu_vz;
  int mc_nu_pdg, mc_ccnc, mc_interaction_type;
  in.at( "true_nu_vtx_x" ) >> mc_nu_vx;
  in.at( "true_nu_vtx_y" ) >> mc_nu_vy;
  in.at( "true_nu_vtx_z" ) >> mc_nu_vz;
  in.at( "nu_pdg" ) >> mc_nu_pdg;
  in.at( "ccnc" ) >> mc_ccnc;
  in.at( "interaction" ) >> mc_interaction_type;

  std::vector< int >* nu_daughter_pdg;
  std::vector< float > *nu_daughter_energy, *nu_daughter_px,
    *nu_daughter_py, *nu_daughter_pz;

  in.at( "mc_pdg" ) >> nu_daughter_pdg;
  in.at( "mc_E" ) >> nu_daughter_energy;
  in.at( "mc_px" ) >> nu_daughter_px;
  in.at( "mc_py" ) >> nu_daughter_py;
  in.at( "mc_pz" ) >> nu_daughter_pz;

  // Count various kinds of final-state particles taking thresholds into
  // account
  int sig_mc_n_threshold_muon = 0;
  int sig_mc_n_threshold_proton = 0;
  int sig_mc_n_threshold_pion0 = 0;
  int sig_mc_n_threshold_pionpm = 0;
  int sig_mc_n_heaviermeson = 0;

  int true_proton_index = BOGUS_INDEX;
  int true_muon_index = BOGUS_INDEX;

  for ( size_t p = 0u; p < nu_daughter_pdg->size(); ++p ) {
    int pdg = nu_daughter_pdg->at( p );
    TVector3 p3( nu_daughter_px->at( p ),
      nu_daughter_py->at( p ), nu_daughter_pz->at( p ) );
    double mom = p3.Mag();

    if ( pdg == MUON && mom >= 0.1 ) {
      ++sig_mc_n_threshold_muon;
      true_muon_index = p;
    }
    else if ( pdg == PROTON && mom >= 0.3 ) {
      ++sig_mc_n_threshold_proton;
      true_proton_index = p;
    }
    else if ( pdg == PI_ZERO ) {
      ++sig_mc_n_threshold_pion0;
    }
    else if ( std::abs( pdg ) == PI_PLUS && mom >= 0.07 ) {
      ++sig_mc_n_threshold_pionpm;
    }
    else if ( pdg != PI_ZERO && std::abs( pdg ) != PI_PLUS
      && is_meson_or_antimeson( pdg ) )
    {
      ++sig_mc_n_heaviermeson;
    }
  }

  // Calculate the booleans related to the different signal cuts
  bool sig_true_vertex_in_fv = this->get_fv().is_inside( mc_nu_vx,
    mc_nu_vy, mc_nu_vz );
  bool sig_is_cc = ( mc_ccnc == CHARGED_CURRENT );
  bool sig_is_numu = ( mc_nu_pdg == MUON_NEUTRINO );
  bool sig_one_muon_above_thresh = ( sig_mc_n_threshold_muon == 1 );
  bool sig_one_proton_above_thresh = ( sig_mc_n_threshold_proton == 1 );
  bool sig_no_pions = ( sig_mc_n_threshold_pion0 == 0
    && sig_mc_n_threshold_pionpm == 0 );
  bool sig_no_heavy_mesons = ( sig_mc_n_heaviermeson == 0 );

  // Is the event signal?
  bool is_signal = sig_true_vertex_in_fv && sig_is_cc && sig_is_numu
    && sig_one_muon_above_thresh && sig_one_proton_above_thresh
    && sig_no_pions && sig_no_heavy_mesons;

  // Store some true values of observables for interesting events
  bool is_interesting = sig_one_muon_above_thresh
    && sig_one_proton_above_thresh && sig_no_pions && sig_no_heavy_mesons;

  double true_pT = BOGUS, true_pTx = BOGUS, true_pTy = BOGUS, true_pL = BOGUS,
    true_pn = BOGUS, true_pn_perp = BOGUS, true_pn_perpx = BOGUS,
    true_pn_perpy = BOGUS, true_pn_par = BOGUS, true_delta_alphaT = BOGUS,
    true_delta_alpha3D_q = BOGUS, true_delta_alpha3D_mu = BOGUS,
    true_delta_phiT = BOGUS, true_delta_phi3D = BOGUS, true_Ecal = BOGUS,
    true_E_QE = BOGUS, true_Q2 = BOGUS, true_A = BOGUS, true_E_miss = BOGUS,
    true_k_miss = BOGUS, true_p_miss = BOGUS, true_p_miss_minus = BOGUS;

  if ( is_interesting ) {

    double mu_px = nu_daughter_px->at( true_muon_index );
    double mu_py = nu_daughter_py->at( true_muon_index );
    double mu_pz = nu_daughter_pz->at( true_muon_index );
    TVector3 p3_mu( mu_px, mu_py, mu_pz );

    double p_px = nu_daughter_px->at( true_proton_index );
    double p_py = nu_daughter_py->at( true_proton_index );
    double p_pz = nu_daughter_pz->at( true_proton_index );
    TVector3 p3_p( p_px, p_py, p_pz );

    KICalculator ki_calc( p3_mu, p3_p, ki_calc_type_ );

    true_pT = ki_calc.pT();
    true_pTx = ki_calc.pTx();
    true_pTy = ki_calc.pTy();
    true_pL = ki_calc.pL();
    true_pn = ki_calc.pn();
    true_pn_perp = ki_calc.pn_perp();
    true_pn_perpx = ki_calc.pn_perpx();
    true_pn_perpy = ki_calc.pn_perpy();
    true_pn_par = ki_calc.pn_par();
    true_delta_alphaT = ki_calc.delta_alphaT();
    true_delta_alpha3D_q = ki_calc.delta_alpha3D_q();
    true_delta_alpha3D_mu = ki_calc.delta_alpha3D_mu();
    true_delta_phiT = ki_calc.delta_phiT();
    true_delta_phi3D = ki_calc.delta_phi3D();
    true_Ecal = ki_calc.Ecal();
    true_E_QE = ki_calc.E_QE();
    true_Q2 = ki_calc.Q2();
    true_A = ki_calc.A();
    true_E_miss = ki_calc.E_miss();
    true_k_miss = ki_calc.k_miss();
    true_p_miss = ki_calc.p_miss();
    true_p_miss_minus = ki_calc.p_miss_minus();
  }

  out[ "sig_true_vertex_in_fv" ] = sig_true_vertex_in_fv;
  out[ "sig_is_cc" ] = sig_is_cc;
  out[ "sig_is_numu" ] = sig_is_numu;
  out[ "sig_one_muon_above_thresh" ] = sig_one_muon_above_thresh;
  out[ "sig_one_proton_above_thresh" ] = sig_one_proton_above_thresh;
  out[ "sig_no_pions" ] = sig_no_pions;
  out[ "sig_no_heavy_mesons" ] = sig_no_heavy_mesons;
  out[ "mc_n_threshold_muon" ] = sig_mc_n_threshold_muon;
  out[ "mc_n_threshold_proton" ] = sig_mc_n_threshold_proton;
  out[ "mc_n_threshold_pion0" ] = sig_mc_n_threshold_pion0;
  out[ "mc_n_threshold_pionpm" ] = sig_mc_n_threshold_pionpm;
  out[ "mc_n_heaviermeson" ] = sig_mc_n_heaviermeson;
  out[ "true_muon_index" ] = true_muon_index;
  out[ "true_proton_index" ] = true_proton_index;

  out[ "True_pT" ] = true_pT;
  out[ "True_pTx" ] = true_pTx;
  out[ "True_pTy" ] = true_pTy;
  out[ "True_pL" ] = true_pL;
  out[ "True_pn" ] = true_pn;
  out[ "True_pn_perp" ] = true_pn_perp;
  out[ "True_pn_perpx" ] = true_pn_perpx;
  out[ "True_pn_perpy" ] = true_pn_perpy;
  out[ "True_pn_par" ] = true_pn_par;
  out[ "True_delta_alphaT" ] = true_delta_alphaT;
  out[ "True_delta_alpha3D_q" ] = true_delta_alpha3D_q;
  out[ "True_delta_alpha3D_mu" ] = true_delta_alpha3D_mu;
  out[ "True_delta_phiT" ] = true_delta_phiT;
  out[ "True_delta_phi3D" ] = true_delta_phi3D;
  out[ "True_Ecal" ] = true_Ecal;
  out[ "True_E_QE" ] = true_E_QE;
  out[ "True_Q2" ] = true_Q2;
  out[ "True_A" ] = true_A;
  out[ "True_E_miss" ] = true_E_miss;
  out[ "True_k_miss" ] = true_k_miss;
  out[ "True_p_miss" ] = true_p_miss;
  out[ "True_p_miss_minus" ] = true_p_miss_minus;

  // Categorize by interaction mode and whether the event fulfills the
  // signal definition
  std::string categ;
  if ( is_signal ) categ += "S ";
  else categ += "B ";

  if ( mc_interaction_type == QE_INTERACTION ) categ += "QE";
  else if ( mc_interaction_type == MEC_INTERACTION ) categ += "MEC";
  else if ( mc_interaction_type == RES_INTERACTION ) categ += "RES";
  else categ += "DIS";

  return categ;
}

bool CC1mu1p0pi::is_selected( AnalysisEvent& ev ) {

  // Get access to variables from the PeLEE ntuple
  const auto& in = ev.in();

  int nslice, num_pf_particles;
  std::vector< int >* pfp_reco_pdg;
  std::vector< unsigned int >* gen_vec;
  std::vector< float > *pfp_track_score, *track_llr_pid_score;
  in.at( "nslice" ) >> nslice;
  in.at( "n_pfps" ) >> num_pf_particles;
  in.at( "pfp_generation_v" ) >> gen_vec;
  in.at( "trk_score_v" ) >> pfp_track_score;
  in.at( "trk_llr_pid_score_v" ) >> track_llr_pid_score;
  in.at( "pfpdg" ) >> pfp_reco_pdg;

  // Requirement for exactly one neutrino slice
  bool sel_nslice_eq_1 = ( nslice == 1 );

  // Requirement for exactly 2 tracks and 0 showers for a CC1p0pi selection
  int reco_shower_count = 0;
  int reco_track_count = 0;
  std::vector< int > candidate_indices;

  for ( int p = 0; p < num_pf_particles; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = gen_vec->at( p );
    if ( generation != 2u ) continue;

    float tscore = pfp_track_score->at( p );
    if ( tscore <= TRACK_SCORE_CUT ) {
      ++reco_shower_count;
    } else {
      ++reco_track_count;
      candidate_indices.push_back( p );
    }

  }

  bool sel_nshower_eq_0 = ( reco_shower_count == 0 );
  bool sel_ntrack_eq_2 = ( reco_track_count == 2 );

  // Identify candidate muon & proton
  // Muon = the one with the highest LLR PID Score
  // Proton = the one with the lowest LLR PID Score

  int mu_idx = BOGUS_INDEX;
  int p_idx = BOGUS_INDEX;

  // Default values of these cuts
  bool sel_muon_candidate_tracklike = false;
  bool sel_proton_candidate_tracklike = false;

  if ( sel_ntrack_eq_2 ) {
    float first_pid_score = track_llr_pid_score->at( candidate_indices.at(0) );
    float second_pid_score = track_llr_pid_score->at( candidate_indices.at(1) );

    if ( first_pid_score > second_pid_score ) {
      mu_idx = candidate_indices.at(0);
      p_idx = candidate_indices.at(1);
    } else {
      mu_idx = candidate_indices.at(1);
      p_idx = candidate_indices.at(0);
    }


    sel_muon_candidate_tracklike = ( pfp_reco_pdg->at( mu_idx ) == MUON );
    sel_proton_candidate_tracklike = ( pfp_reco_pdg->at( p_idx ) == MUON );
  }

  // Neutrino vertex in FV?
  float nu_vx, nu_vy, nu_vz;
  in.at( "reco_nu_vtx_sce_x" ) >> nu_vx;
  in.at( "reco_nu_vtx_sce_y" ) >> nu_vy;
  in.at( "reco_nu_vtx_sce_z" ) >> nu_vz;

  bool sel_nu_vertex_contained = this->get_fv().is_inside( nu_vx,
    nu_vy, nu_vz );

  // Containment check on the muon and proton
  std::vector< float > *track_start_x, *track_start_y, *track_start_z,
    *track_end_x, *track_end_y, *track_end_z, *muon_mom_range, *muon_mom_mcs,
    *track_KE_proton;

  in.at( "trk_sce_start_x_v" ) >> track_start_x;
  in.at( "trk_sce_start_y_v" ) >> track_start_y;
  in.at( "trk_sce_start_z_v" ) >> track_start_z;

  in.at( "trk_sce_end_x_v" ) >> track_end_x;
  in.at( "trk_sce_end_y_v" ) >> track_end_y;
  in.at( "trk_sce_end_z_v" ) >> track_end_z;

  in.at( "trk_range_muon_mom_v" ) >> muon_mom_range;
  in.at( "trk_mcs_muon_mom_v" ) >> muon_mom_mcs;
  in.at( "trk_energy_proton_v" ) >> track_KE_proton;

  bool sel_muon_candidate_above_p_thresh = false,
    sel_proton_candidate_above_p_thresh = false, sel_mu_contained = false,
    sel_p_contained = false;

  double mu_mom = BOGUS, p_mom = BOGUS;

  if ( sel_ntrack_eq_2 ) {

    // Muon variables
    double mu_trk_start_x = track_start_x->at( mu_idx );
    double mu_trk_start_y = track_start_y->at( mu_idx );
    double mu_trk_start_z = track_start_z->at( mu_idx );

    double mu_trk_end_x = track_end_x->at( mu_idx );
    double mu_trk_end_y = track_end_y->at( mu_idx );
    double mu_trk_end_z = track_end_z->at( mu_idx );

    bool mu_start_contained = this->get_fv().is_inside(
      mu_trk_start_x, mu_trk_start_y, mu_trk_start_z );

    bool mu_end_contained = this->get_fv().is_inside(
      mu_trk_end_x, mu_trk_end_y, mu_trk_end_z );

    double mu_mom = muon_mom_range->at( mu_idx ); // GeV/c

    // If exiting muon, switch to MCS and recalculate the momentum
    if ( !mu_end_contained ) {
      mu_mom = muon_mom_mcs->at( mu_idx );
    }

    // Proton variables
    double p_trk_start_x = track_start_x->at( p_idx );
    double p_trk_start_y = track_start_y->at( p_idx );
    double p_trk_start_z = track_start_z->at( p_idx );

    double p_trk_end_x = track_end_x->at( p_idx );
    double p_trk_end_y = track_end_y->at( p_idx );
    double p_trk_end_z = track_end_z->at( p_idx );

    bool p_start_contained = this->get_fv().is_inside( p_trk_start_x,
      p_trk_start_y, p_trk_start_z );

    bool p_end_contained = this->get_fv().is_inside( p_trk_end_x, p_trk_end_y,
      p_trk_end_z );

    double p_KE = track_KE_proton->at( p_idx ); // GeV
    p_mom = real_sqrt( p_KE*p_KE + 2.*PROTON_MASS*p_KE );

    // Redefinition if p_p < 0.5 where biases have been observed
    if ( p_mom < 0.5 ) {
      p_mom = ( 1.-0.01*( 29.354172*p_mom - 14.674918 ) )* p_mom;
    }

    // Comes from /uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs
    // /ubana/ubana/myClasses/Constants.hh:L1762
    sel_muon_candidate_above_p_thresh = ( mu_mom >= MUON_P_MIN_MOM_CUT );

    // Comes from /uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/
    // ubana/ubana/myClasses/Constants.hh:L1775
    sel_proton_candidate_above_p_thresh = ( p_mom >= 0.3 );

    sel_mu_contained = ( mu_start_contained && mu_end_contained );
    sel_p_contained = ( p_start_contained && p_end_contained );

  }

  // Check that the muon momentum estimators agree within 25%

  bool sel_muon_momentum_quality = false;
  if ( mu_idx != BOGUS_INDEX ) {
    double p_mcs = muon_mom_mcs->at( mu_idx );
    double p_range = muon_mom_range->at( mu_idx );

    double frac_diff = std::abs( p_mcs - p_range ) / p_range;
    sel_muon_momentum_quality = ( frac_diff <= 0.25 );
  }

  // Check for flipped tracks
  bool sel_no_flipped_tracks = false;

  if ( mu_idx != BOGUS_INDEX && p_idx != BOGUS_INDEX ) {
    TVector3 x3_vtx( nu_vx, nu_vy, nu_vz );

    TVector3 x3_mu_start( track_start_x->at( mu_idx ),
      track_start_y->at( mu_idx ), track_start_z->at( mu_idx ) );

    TVector3 x3_mu_end( track_end_x->at( mu_idx ), track_end_y->at( mu_idx ),
      track_end_z->at( mu_idx ) );

    TVector3 x3_p_start( track_start_x->at( p_idx ),
      track_start_y->at( p_idx ), track_start_z->at( p_idx ) );

    TVector3 x3_p_end( track_end_x->at( p_idx ), track_end_y->at( p_idx ),
      track_end_z->at( p_idx ) );

    double vtx_to_mu_start = ( x3_vtx - x3_mu_start ).Mag();
    double vtx_to_mu_end = ( x3_vtx - x3_mu_end ).Mag();
    double vtx_to_p_start = ( x3_vtx - x3_p_start ).Mag();
    double vtx_to_p_end = ( x3_vtx - x3_p_end ).Mag();

    double mu_start_to_p_start = ( x3_mu_start - x3_p_start ).Mag();
    double mu_end_to_p_end = ( x3_mu_end - x3_p_end ).Mag();

    sel_no_flipped_tracks = !( ( vtx_to_mu_start > vtx_to_mu_end )
      || ( vtx_to_p_start > vtx_to_p_end )
      || ( mu_start_to_p_start > mu_end_to_p_end ) );
  }

  // Check Proton Candidate's LLH to be a proton
  bool sel_proton_cand_passed_llr_cut = false;
  if ( mu_idx != BOGUS_INDEX && p_idx != BOGUS_INDEX ) {
    double llr_score = track_llr_pid_score->at( p_idx );
    sel_proton_cand_passed_llr_cut = ( llr_score < 0.05 );
  }

  // Apply kinematic cuts
  std::vector< float > *track_theta, *track_phi;
  in.at( "trk_theta_v" ) >> track_theta;
  in.at( "trk_phi_v" ) >> track_phi;

  bool sel_muon_momentum_in_range = false, sel_muon_costheta_in_range = false,
    sel_muon_phi_in_range = false, sel_proton_momentum_in_range = false,
    sel_proton_costheta_in_range = false, sel_proton_phi_in_range = false;

  TVector3 p3_mu;
  if ( mu_idx != BOGUS_INDEX ) {
    // This is essentially duplicate of sel_muon_candidate_above_p_thresh cut
    double mu_momentum = muon_mom_range->at( mu_idx );
    if ( ! ( mu_momentum < 0.1 || mu_momentum > 1.2 ) ) {
      sel_muon_momentum_in_range = true;
    }

    double mu_cth = std::cos( track_theta->at( mu_idx ) );
    if ( !( mu_cth < -1. || mu_cth > 1. ) ) {
      sel_muon_costheta_in_range = true;
    }

    double mu_phi = track_phi->at( mu_idx ) * 180. / M_PI;
    if ( !( mu_phi < -180. || mu_phi > 180. ) ) {
      sel_muon_phi_in_range = true;
    }

    p3_mu = TVector3( 0., 0., 1. );
    p3_mu.SetMag( mu_mom );
    p3_mu.SetTheta( std::acos( mu_cth ) );
    p3_mu.SetPhi( mu_phi);
  }

  TVector3 p3_p;
  if ( p_idx != BOGUS_INDEX ) {

    if ( ! ( p_mom < 0.3 || p_mom > 1. ) ) {
      sel_proton_momentum_in_range = true;
    }

    double p_cth = std::cos( track_theta->at( p_idx ) );

    if ( !( p_cth < -1. || p_cth > 1. ) ) {
      sel_proton_costheta_in_range = true;
    }

    double p_phi = track_phi->at( p_idx ) * 180. / M_PI;
    if ( !( p_phi < -180. || p_phi > 180. ) ) {
      sel_proton_phi_in_range = true;
    }

    p3_p = TVector3( 0., 0., 1. );
    p3_p.SetMag( p_mom );
    p3_p.SetTheta( std::acos( p_cth ) );
    p3_p.SetPhi( p_phi);

  }

  // Compute some observables
  double reco_pT = BOGUS, reco_pTx = BOGUS, reco_pTy = BOGUS, reco_pL = BOGUS,
    reco_pn = BOGUS, reco_pn_perp = BOGUS, reco_pn_perpx = BOGUS,
    reco_pn_perpy = BOGUS, reco_pn_par = BOGUS, reco_delta_alphaT = BOGUS,
    reco_delta_alpha3D_q = BOGUS, reco_delta_alpha3D_mu = BOGUS,
    reco_delta_phiT = BOGUS, reco_delta_phi3D = BOGUS, reco_Ecal = BOGUS,
    reco_E_QE = BOGUS, reco_Q2 = BOGUS, reco_A = BOGUS, reco_E_miss = BOGUS,
    reco_k_miss = BOGUS, reco_p_miss = BOGUS, reco_p_miss_minus = BOGUS,
    backTrack_pT = BOGUS, backTrack_pTx = BOGUS, backTrack_pTy = BOGUS,
    backTrack_pL = BOGUS, backTrack_pn = BOGUS, backTrack_pn_perp = BOGUS,
    backTrack_pn_perpx = BOGUS, backTrack_pn_perpy = BOGUS,
    backTrack_pn_par = BOGUS, backTrack_delta_alphaT = BOGUS,
    backTrack_delta_alpha3D_q = BOGUS, backTrack_delta_alpha3D_mu = BOGUS,
    backTrack_delta_phiT = BOGUS, backTrack_delta_phi3D = BOGUS,
    backTrack_Ecal = BOGUS, backTrack_E_QE = BOGUS, backTrack_Q2 = BOGUS,
    backTrack_A = BOGUS, backTrack_E_miss = BOGUS, backTrack_k_miss = BOGUS,
    backTrack_p_miss = BOGUS, backTrack_p_miss_minus = BOGUS;

  if ( mu_idx != BOGUS_INDEX && p_idx != BOGUS_INDEX) {

    KICalculator ki_calc( p3_mu, p3_p, ki_calc_type_ );

    reco_pT = ki_calc.pT();
    reco_pTx = ki_calc.pTx();
    reco_pTy = ki_calc.pTy();
    reco_pL = ki_calc.pL();
    reco_pn = ki_calc.pn();
    reco_pn_perp = ki_calc.pn_perp();
    reco_pn_perpx = ki_calc.pn_perpx();
    reco_pn_perpy = ki_calc.pn_perpy();
    reco_pn_par = ki_calc.pn_par();
    reco_delta_alphaT = ki_calc.delta_alphaT();
    reco_delta_alpha3D_q = ki_calc.delta_alpha3D_q();
    reco_delta_alpha3D_mu = ki_calc.delta_alpha3D_mu();
    reco_delta_phiT = ki_calc.delta_phiT();
    reco_delta_phi3D = ki_calc.delta_phi3D();
    reco_Ecal = ki_calc.Ecal();
    reco_E_QE = ki_calc.E_QE();
    reco_Q2 = ki_calc.Q2();
    reco_A = ki_calc.A();
    reco_E_miss = ki_calc.E_miss();
    reco_k_miss = ki_calc.k_miss();
    reco_p_miss = ki_calc.p_miss();
    reco_p_miss_minus = ki_calc.p_miss_minus();

    // Now deal with the backtracking
    std::vector< float > *bt_px, *bt_py, *bt_pz;
    in.at( "backtracked_px" ) >> bt_px;
    in.at( "backtracked_py" ) >> bt_py;
    in.at( "backtracked_pz" ) >> bt_pz;

    TVector3 p3_mu_bt( bt_px->at( mu_idx ),
      bt_py->at( mu_idx ), bt_pz->at( mu_idx ) );

    TVector3 p3_p_bt( bt_px->at( p_idx ), bt_py->at( p_idx ),
      bt_pz->at( p_idx ) );

    KICalculator bt_ki_calc( p3_mu_bt, p3_p_bt, ki_calc_type_ );

    backTrack_pT = bt_ki_calc.pT();
    backTrack_pTx = bt_ki_calc.pTx();
    backTrack_pTy = bt_ki_calc.pTy();
    backTrack_pL = bt_ki_calc.pL();
    backTrack_pn = bt_ki_calc.pn();
    backTrack_pn_perp = bt_ki_calc.pn_perp();
    backTrack_pn_perpx = bt_ki_calc.pn_perpx();
    backTrack_pn_perpy = bt_ki_calc.pn_perpy();
    backTrack_pn_par = bt_ki_calc.pn_par();
    backTrack_delta_alphaT = bt_ki_calc.delta_alphaT();
    backTrack_delta_alpha3D_q = bt_ki_calc.delta_alpha3D_q();
    backTrack_delta_alpha3D_mu = bt_ki_calc.delta_alpha3D_mu();
    backTrack_delta_phiT = bt_ki_calc.delta_phiT();
    backTrack_delta_phi3D = bt_ki_calc.delta_phi3D();
    backTrack_Ecal = bt_ki_calc.Ecal();
    backTrack_E_QE = bt_ki_calc.E_QE();
    backTrack_Q2 = bt_ki_calc.Q2();
    backTrack_A = bt_ki_calc.A();
    backTrack_E_miss = bt_ki_calc.E_miss();
    backTrack_k_miss = bt_ki_calc.k_miss();
    backTrack_p_miss = bt_ki_calc.p_miss();
    backTrack_p_miss_minus = bt_ki_calc.p_miss_minus();
  }

  auto& out = ev.out();

  out[ "nslice_eq_1" ] = sel_nslice_eq_1;
  out[ "nshower_eq_0" ] = sel_nshower_eq_0;
  out[ "ntrack_eq_2" ] = sel_ntrack_eq_2;
  out[ "muon_candidate_tracklike" ] = sel_muon_candidate_tracklike;
  out[ "proton_candidate_tracklike" ] = sel_proton_candidate_tracklike;
  out[ "nu_vertex_contained" ] = sel_nu_vertex_contained;

  out[ "muon_candidate_above_p_thresh" ] = sel_muon_candidate_above_p_thresh;

  out[ "proton_candidate_above_p_thresh" ]
    = sel_proton_candidate_above_p_thresh;

  out[ "muon_candidate_contained" ] = sel_mu_contained;
  out[ "proton_candidate_contained" ] = sel_p_contained;
  out[ "sel_muon_momentum_quality" ] = sel_muon_momentum_quality;
  out[ "sel_no_flipped_tracks" ] = sel_no_flipped_tracks;
  out[ "sel_proton_cand_passed_llr_cut" ] = sel_proton_cand_passed_llr_cut;
  out[ "sel_muon_momentum_in_range" ] = sel_muon_momentum_in_range;
  out[ "sel_muon_costheta_in_range" ] = sel_muon_costheta_in_range;
  out[ "sel_muon_phi_in_range" ] = sel_muon_phi_in_range;
  out[ "sel_proton_momentum_in_range" ] = sel_proton_momentum_in_range;
  out[ "sel_proton_costheta_in_range" ] = sel_proton_costheta_in_range;
  out[ "sel_proton_phi_in_range" ] = sel_proton_phi_in_range;
  out[ "CandidateMuonIndex" ] = mu_idx;
  out[ "CandidateProtonIndex" ] = p_idx;

  out[ "Reco_pT" ] = reco_pT;
  out[ "Reco_pTx" ] = reco_pTx;
  out[ "Reco_pTy" ] = reco_pTy;
  out[ "Reco_pL" ] = reco_pL;
  out[ "Reco_pn" ] = reco_pn;
  out[ "Reco_pn_perp" ] = reco_pn_perp;
  out[ "Reco_pn_perpx" ] = reco_pn_perpx;
  out[ "Reco_pn_perpy" ] = reco_pn_perpy;
  out[ "Reco_pn_par" ] = reco_pn_par;
  out[ "Reco_delta_alphaT" ] = reco_delta_alphaT;
  out[ "Reco_delta_alpha3D_q" ] = reco_delta_alpha3D_q;
  out[ "Reco_delta_alpha3D_mu" ] = reco_delta_alpha3D_mu;
  out[ "Reco_delta_phiT" ] = reco_delta_phiT;
  out[ "Reco_delta_phi3D" ] = reco_delta_phi3D;
  out[ "Reco_Ecal" ] = reco_Ecal;
  out[ "Reco_E_QE" ] = reco_E_QE;
  out[ "Reco_Q2" ] = reco_Q2;
  out[ "Reco_A" ] = reco_A;
  out[ "Reco_E_miss" ] = reco_E_miss;
  out[ "Reco_k_miss" ] = reco_k_miss;
  out[ "Reco_p_miss" ] = reco_p_miss;
  out[ "Reco_p_miss_minus" ] = reco_p_miss_minus;

  out[ "BackTrack_pT" ] = backTrack_pT;
  out[ "BackTrack_pTx" ] = backTrack_pTx;
  out[ "BackTrack_pTy" ] = backTrack_pTy;
  out[ "BackTrack_pL" ] = backTrack_pL;
  out[ "BackTrack_pn" ] = backTrack_pn;
  out[ "BackTrack_pn_perp" ] = backTrack_pn_perp;
  out[ "BackTrack_pn_perpx" ] = backTrack_pn_perpx;
  out[ "BackTrack_pn_perpy" ] = backTrack_pn_perpy;
  out[ "BackTrack_pn_par" ] = backTrack_pn_par;
  out[ "BackTrack_delta_alphaT" ] = backTrack_delta_alphaT;
  out[ "BackTrack_delta_alpha3D_q" ] = backTrack_delta_alpha3D_q;
  out[ "BackTrack_delta_alpha3D_mu" ] = backTrack_delta_alpha3D_mu;
  out[ "BackTrack_delta_phiT" ] = backTrack_delta_phiT;
  out[ "BackTrack_delta_phi3D" ] = backTrack_delta_phi3D;
  out[ "BackTrack_Ecal" ] = backTrack_Ecal;
  out[ "BackTrack_E_QE" ] = backTrack_E_QE;
  out[ "BackTrack_Q2" ] = backTrack_Q2;
  out[ "BackTrack_A" ] = backTrack_A;
  out[ "BackTrack_E_miss" ] = backTrack_E_miss;
  out[ "BackTrack_k_miss" ] = backTrack_k_miss;
  out[ "BackTrack_p_miss" ] = backTrack_p_miss;
  out[ "BackTrack_p_miss_minus" ] = backTrack_p_miss_minus;

  // Does everything pass selection?
  bool passed = sel_nslice_eq_1 && sel_nshower_eq_0 && sel_ntrack_eq_2
    && sel_muon_candidate_tracklike && sel_proton_candidate_tracklike
    && sel_nu_vertex_contained && sel_muon_candidate_above_p_thresh
    && sel_proton_candidate_above_p_thresh && sel_mu_contained
    && sel_p_contained && sel_muon_momentum_quality
    && sel_no_flipped_tracks && sel_proton_cand_passed_llr_cut
    && sel_muon_momentum_in_range && sel_muon_costheta_in_range
    && sel_muon_phi_in_range && sel_proton_momentum_in_range
    && sel_proton_costheta_in_range && sel_proton_phi_in_range;
  return passed;
}
