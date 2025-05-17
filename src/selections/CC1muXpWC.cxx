// XSecAnalyzer includes
#include "XSecAnalyzer/Constants.hh"
#include "XSecAnalyzer/Selections/CC1muXpWC.hh"

namespace {
  constexpr float NUMU_SCORE_CUT = 0.9f;
  constexpr float PROTON_ENERGY_THRESHOLD = 0.973f; // GeV
  constexpr int PRIMARY_PARTICLE = 0;
}

CC1muXpWC::CC1muXpWC() : SelectionBase( "CC1muXpWC" ) {
  this->define_fv( 3., 253., -113., 113., 3., 1034. ); // cm
}

bool CC1muXpWC::is_selected( AnalysisEvent& ev ) {

  // Get access to the input WireCell TTrees
  const auto& bdt_vars = ev.in( "wcpselection/T_BDTvars" );
  const auto& pf_eval = ev.in( "wcpselection/T_PFeval" );
  const auto& eval = ev.in( "wcpselection/T_eval" );

  // Check whether the space-charge-corrected position of the reconstructed
  // neutrino vertex lies within the fiducial volume
  float nu_vx, nu_vy, nu_vz;
  pf_eval.at( "reco_nuvtxX" ) >> nu_vx;
  pf_eval.at( "reco_nuvtxY" ) >> nu_vy;
  pf_eval.at( "reco_nuvtxZ" ) >> nu_vz;

  // Is the reco neutrino vertex in the fiducial volume
  bool sel_reco_vertex_in_FV = this->get_fv().is_inside( nu_vx, nu_vy, nu_vz );

  static const std::string inclusive_cuts( "numu_score > "
    + std::to_string( NUMU_SCORE_CUT ) + " && numu_cc_flag >= 0" );

  // Does the event pass the numu CC inclusive selection?
  bool sel_is_numu_CC = bdt_vars.formula( inclusive_cuts );

  // Is the event fully (true) or partially (false) contained?
  bool sel_is_FC = eval.formula( "match_isFC" );

  static const std::string np_cuts( "reco_mother == 0 && reco_pdg == "
    + std::to_string( PROTON ) + " && reco_startMomentum[][ 3 ] > "
    + std::to_string( PROTON_ENERGY_THRESHOLD ) );

  // Classify by 0p/Np by checking for a proton candidate track
  bool sel_is_Np = pf_eval.formula( np_cuts ).use_or();

  // Apply the "tight" version of the selection by explicitly requiring a muon
  bool sel_has_muon = pf_eval.formula( "reco_muonMomentum[ 3 ]  > 0." );

  // Get access to the output tree
  auto& out = ev.out();
  out[ "sel_reco_vertex_in_FV" ] = sel_reco_vertex_in_FV;
  out[ "sel_is_numu_CC" ] = sel_is_numu_CC;
  out[ "sel_is_FC" ] = sel_is_FC;
  out[ "sel_is_Np" ] = sel_is_Np;
  out[ "sel_has_muon" ] = sel_has_muon;
  out[ "sel_CC1muXp_tight" ] = sel_is_numu_CC && sel_has_muon;

  return sel_is_numu_CC;
}

std::string CC1muXpWC::categorize_event( AnalysisEvent& ev ) {

  static const std::string signal_cuts( "truth_vtxInside && truth_nuPdg == "
    + std::to_string( MUON_NEUTRINO ) + " && truth_isCC" );

  static const std::string np_cuts( "truth_mother == 0 && truth_pdg == "
    + std::to_string( PROTON ) + " && truth_startMomentum[][ 3 ] > "
    + std::to_string( PROTON_ENERGY_THRESHOLD ) );

  const auto& eval = ev.in( "wcpselection/T_eval" );
  const auto& pf_eval = ev.in( "wcpselection/T_PFeval" );

  bool is_signal = eval.formula( signal_cuts );
  bool is_Np = pf_eval.formula( np_cuts ).use_or();

  auto& out = ev.out();
  out[ "sig_is_Np" ] = is_Np;

  if ( is_signal ) {
    if ( is_Np ) return "CC Np";
    else return "CC 0p";
  }
  return "Background";
}
