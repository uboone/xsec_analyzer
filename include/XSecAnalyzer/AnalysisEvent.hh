#pragma once

// STV analysis includes
#include "TreeUtils.hh"
#include "FiducialVolume.hh"
#include "Constants.hh"

#include <vector>
#include <map>

#include "TVector3.h"

class AnalysisEvent{
public:
  AnalysisEvent() {}
  ~AnalysisEvent() {}

  /* EVENT LEVEL RECO VARIABLES */
  //================================================================================================================

  // Event scores needed for numu CC selection
  float topological_score_ = BOGUS;
  float cosmic_impact_parameter_ = BOGUS;

  // Backtracked purity and completeness of hits (MC only)
  float nu_completeness_from_pfp_ = BOGUS;
  float nu_purity_from_pfp_ = BOGUS;

  // Reco PDG code of the neutrino candidate
  int nu_pdg_ = BOGUS_INT;

  // Number of neutrino slices identified by the SliceID. Allowed values
  // are zero or one.
  int nslice_ = BOGUS_INT;

  // Reco neutrino vertex coordinates (cm). Space charge corrections have
  // been applied for these.
  float nu_vx_ = BOGUS;
  float nu_vy_ = BOGUS;
  float nu_vz_ = BOGUS;

  // Reconstructed object counts
  int num_pf_particles_ = BOGUS_INT;
  int num_tracks_ = BOGUS_INT;
  int num_showers_ = BOGUS_INT;

  //================================================================================================================
  /* EVENT LEVEL TRUE VARIABLES */

  // True neutrino PDG code
  int mc_nu_pdg_ = BOGUS_INT;

  // True neutrino vertex coordinates (cm)
  float mc_nu_vx_ = BOGUS;
  float mc_nu_vy_ = BOGUS;
  float mc_nu_vz_ = BOGUS;

  float mc_nu_sce_vx_ = BOGUS;
  float mc_nu_sce_vy_ = BOGUS;
  float mc_nu_sce_vz_ = BOGUS;

  // True neutrino 4-momentum
  float mc_nu_energy_ = BOGUS;

  // Whether the event is CC (0) or NC (1)
  int mc_nu_ccnc_ = false;

  // Interaction mode (QE, MEC, etc.)
  int mc_nu_interaction_type_ = BOGUS_INT;

  //================================================================================================================
  /* PARTICLE LEVEL RECO VARIABLES */

  // BDT scores
  MyPointer< std::vector<float> > proton_BDT_score_;
  MyPointer< std::vector<float> > muon_BDT_score_;

  // pion/golden pion handles
  MyPointer< std::vector<int> > pfp_n_descendents_;
  MyPointer< std::vector<int> > trk_end_spacepoints_;
  MyPointer< std::vector<float> > trk_avg_deflection_stdev_;
  
  // Calo features, 3-plane Bragg
  /*
  MyPointer< std::vector<float> > bragg_mip_uvw_; // TODO: output only 

  MyPointer< std::vector<float> > bragg_p_fwd_uvw_; // TODO: output only 
  MyPointer< std::vector<float> > bragg_p_bwd_uvw_; // TODO: output only 
  MyPointer< std::vector<float> > bragg_p_to_MIP_; // TODO: output only 
  MyPointer< std::vector<float> > bragg_p_fwd_2_bwd_; // TODO: output only 

  MyPointer< std::vector<float> > bragg_pion_fwd_uvw_; // TODO: output only 
  MyPointer< std::vector<float> > bragg_pion_bwd_uvw_; // TODO: output only 
  MyPointer< std::vector<float> > bragg_pion_to_MIP_; // TODO: output only 
  MyPointer< std::vector<float> > bragg_pion_fwd_2_bwd_; // TODO: output only 

  MyPointer< std::vector<float> > bragg_mu_fwd_uvw_; // TODO: output only 
  MyPointer< std::vector<float> > bragg_mu_bwd_uvw_; // TODO: output only 

  MyPointer< std::vector<float> > truncated_mean_dEdx_; // TODO: output only 

  // Calo feature single plane Bragg

  MyPointer< std::vector<float> > trk_trunk_dEdx_w_;
  MyPointer< std::vector<float> > trk_trunk_dEdx_u_;
  MyPointer< std::vector<float> > trk_trunk_dEdx_v_;
  */
  /*
  MyPointer< std::vector<float> > bragg_p_fwd_w_; // TODO: output only 
  MyPointer< std::vector<float> > bragg_p_bwd_w_; // TODO: output only 

  MyPointer< std::vector<float> > bragg_mu_fwd_w_; // TODO: output only 
  MyPointer< std::vector<float> > bragg_mu_bwd_w_; // TODO: output only 

  MyPointer< std::vector<float> > bragg_pion_fwd_w_; // TODO: output only 
  MyPointer< std::vector<float> > bragg_pion_bwd_w_; // TODO: output only 

  MyPointer< std::vector<float> > bragg_p_fwd_v_; // TODO: output only 
  MyPointer< std::vector<float> > bragg_p_bwd_v_; // TODO: output only 

  MyPointer< std::vector<float> > bragg_mu_fwd_v_; // TODO: output only 
  MyPointer< std::vector<float> > bragg_mu_bwd_v_; // TODO: output only 

  MyPointer< std::vector<float> > bragg_pion_fwd_v_; // TODO: output only 
  MyPointer< std::vector<float> > bragg_pion_bwd_v_; // TODO: output only 

  MyPointer< std::vector<float> > bragg_p_fwd_u_; // TODO: output only 
  MyPointer< std::vector<float> > bragg_p_bwd_u_; // TODO: output only 

  MyPointer< std::vector<float> > bragg_mu_fwd_u_; // TODO: output only 
  MyPointer< std::vector<float> > bragg_mu_bwd_u_; // TODO: output only 

  MyPointer< std::vector<float> > bragg_pion_fwd_u_; // TODO: output only 
  MyPointer< std::vector<float> > bragg_pion_bwd_u_; // TODO: output only 

  
  MyPointer< std::vector<bool> > trk_bragg_p_fwd_preferred_w_;
  MyPointer< std::vector<float> > trk_bragg_p_w_;
  MyPointer< std::vector<float> > trk_bragg_p_alt_dir_w_;

  MyPointer< std::vector<bool> > trk_bragg_pion_fwd_preferred_w_;
  MyPointer< std::vector<float> > trk_bragg_pion_w_;
  MyPointer< std::vector<float> > trk_bragg_pion_alt_dir_w_;

  MyPointer< std::vector<bool> > trk_bragg_mu_fwd_preferred_w_;
  MyPointer< std::vector<float> > trk_bragg_mu_w_;
  MyPointer< std::vector<float> > trk_bragg_mu_alt_dir_w_;

  MyPointer< std::vector<bool> > trk_bragg_p_fwd_preferred_u_;
  MyPointer< std::vector<float> > trk_bragg_p_u_;
  MyPointer< std::vector<float> > trk_bragg_p_alt_dir_u_;

  MyPointer< std::vector<bool> > trk_bragg_pion_fwd_preferred_u_;
  MyPointer< std::vector<float> > trk_bragg_pion_u_;
  MyPointer< std::vector<float> > trk_bragg_pion_alt_dir_u_;

  MyPointer< std::vector<bool> > trk_bragg_mu_fwd_preferred_u_;
  MyPointer< std::vector<float> > trk_bragg_mu_u_;
  MyPointer< std::vector<float> > trk_bragg_mu_alt_dir_u_;

  MyPointer< std::vector<bool> > trk_bragg_p_fwd_preferred_v_;
  MyPointer< std::vector<float> > trk_bragg_p_v_;
  MyPointer< std::vector<float> > trk_bragg_p_alt_dir_v_;

  MyPointer< std::vector<bool> > trk_bragg_pion_fwd_preferred_v_;
  MyPointer< std::vector<float> > trk_bragg_pion_v_;
  MyPointer< std::vector<float> > trk_bragg_pion_alt_dir_v_;

  MyPointer< std::vector<bool> > trk_bragg_mu_fwd_preferred_v_;
  MyPointer< std::vector<float> > trk_bragg_mu_v_;
  MyPointer< std::vector<float> > trk_bragg_mu_alt_dir_v_;

  MyPointer< std::vector<float> > bragg_mip_w_;
  MyPointer< std::vector<float> > bragg_mip_v_;
  MyPointer< std::vector<float> > bragg_mip_u_;
  */

  // PFParticle properties
  MyPointer< std::vector<unsigned int> > pfp_generation_;
  MyPointer< std::vector<unsigned int> > pfp_trk_daughters_count_;
  MyPointer< std::vector<unsigned int> > pfp_shr_daughters_count_;

  MyPointer< std::vector<float> > pfp_track_score_;

  // Reco PDG code assigned by Pandora
  MyPointer< std::vector<int> > pfp_reco_pdg_;

  // Total number of wire plane hits associated with each PFParticle
  MyPointer< std::vector<int> > pfp_hits_;

  // Number of hits on the three individual planes
  // (Y is the collection plane)
  MyPointer< std::vector<int> > pfp_hitsU_;
  MyPointer< std::vector<int> > pfp_hitsV_;
  MyPointer< std::vector<int> > pfp_hitsY_;

  // Shower properties
  MyPointer< std::vector<unsigned long> > shower_pfp_id_;
  MyPointer< std::vector<float> > shower_startx_;
  MyPointer< std::vector<float> > shower_starty_;
  MyPointer< std::vector<float> > shower_startz_;
  MyPointer< std::vector<float> > shower_start_distance_;

  // Track properties
  MyPointer< std::vector<unsigned long> > track_pfp_id_;
  MyPointer< std::vector<float> > track_length_;
  MyPointer< std::vector<float> > track_startx_;
  MyPointer< std::vector<float> > track_starty_;
  MyPointer< std::vector<float> > track_startz_;
  MyPointer< std::vector<float> > track_start_distance_;
  MyPointer< std::vector<float> > track_endx_;
  MyPointer< std::vector<float> > track_endy_;
  MyPointer< std::vector<float> > track_endz_;
  MyPointer< std::vector<float> > track_dirx_;
  MyPointer< std::vector<float> > track_diry_;
  MyPointer< std::vector<float> > track_dirz_;
  MyPointer< std::vector<float> > track_theta_;
  MyPointer< std::vector<float> > track_phi_;

  // Proton *kinetic* energy using range-based momentum reconstruction
  MyPointer< std::vector<float> > track_kinetic_energy_p_;

  MyPointer< std::vector<float> > track_range_mom_mu_;
  MyPointer< std::vector<float> > track_mcs_mom_mu_;
  MyPointer< std::vector<float> > track_chi2_proton_;

  // Log-likelihood ratio particle ID information

  // Product of muon/proton log-likelihood ratios from all wire three planes
  MyPointer< std::vector<float> > track_llr_pid_;

  // Individual wire plane muon/proton log-likelihood ratios
  MyPointer< std::vector<float> > track_llr_pid_U_;
  MyPointer< std::vector<float> > track_llr_pid_V_;
  MyPointer< std::vector<float> > track_llr_pid_Y_;

  // Rescaled overall PID score (all three planes) that lies
  // on the interval [-1, 1]
  MyPointer< std::vector<float> > track_llr_pid_score_;

  // True PDG code found using the backtracker
  MyPointer< std::vector<int> > pfp_true_pdg_;
  //================================================================================================================
  /* PARTICLE LEVEL TRUE VARIABLES */

  // pion/golden pion handles
  MyPointer< std::vector<float> > mc_end_p_;
  MyPointer< std::vector<int> > mc_n_inelastic_;
  MyPointer< std::vector<int> > mc_n_elastic_;

  // True 4-momentum components found using the backtracker
  MyPointer< std::vector<float> > pfp_true_E_;
  MyPointer< std::vector<float> > pfp_true_px_;
  MyPointer< std::vector<float> > pfp_true_py_;
  MyPointer< std::vector<float> > pfp_true_pz_;

  // Final-state particle PDG codes and energies (post-FSIs)
  MyPointer< std::vector<int> > mc_nu_daughter_pdg_;
  MyPointer< std::vector<float> > mc_nu_daughter_energy_;
  MyPointer< std::vector<float> > mc_nu_daughter_px_;
  MyPointer< std::vector<float> > mc_nu_daughter_py_;
  MyPointer< std::vector<float> > mc_nu_daughter_pz_;

  //================================================================================================================
  /* OTHER */

  // General systematic weights
  MyPointer< std::map< std::string, std::vector<double> > > mc_weights_map_;
  // Map of pointers used to set output branch addresses for the elements
  // of the weights map. Hacky, but it works.
  // TODO: revisit this to make something more elegant
  std::map< std::string, std::vector<double>* > mc_weights_ptr_map_;

  // GENIE weights
  float spline_weight_ = DEFAULT_WEIGHT;
  float tuned_cv_weight_ = DEFAULT_WEIGHT;

  // Signal definition requirements
  bool is_mc_ = false;

};
