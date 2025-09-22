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

  // True PDG code found using the backtracker
  MyPointer< std::vector<int> > pfp_true_pdg_;

  // True 4-momentum components found using the backtracker
  MyPointer< std::vector<float> > pfp_true_E_;
  MyPointer< std::vector<float> > pfp_true_px_;
  MyPointer< std::vector<float> > pfp_true_py_;
  MyPointer< std::vector<float> > pfp_true_pz_;

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

  // Final-state particle PDG codes and energies (post-FSIs)
  MyPointer< std::vector<int> > mc_nu_daughter_pdg_;
  MyPointer< std::vector<float> > mc_nu_daughter_energy_;
  MyPointer< std::vector<float> > mc_nu_daughter_px_;
  MyPointer< std::vector<float> > mc_nu_daughter_py_;
  MyPointer< std::vector<float> > mc_nu_daughter_pz_;

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

  //Begin Burke Additions
  //Metadata info
  MyPointer< std::vector<float> > slice_topo_score_v_;
  int event_;
  int run_;
  int subrun_;

  //All MC Variables
  MyPointer< std::vector<int> >  all_mc_trkid_;
  MyPointer< std::vector<int> >  all_mc_pdg_;
  MyPointer< std::vector<int> >  all_mc_mother_;
  MyPointer< std::vector<float> >  all_mc_vx_;
  MyPointer< std::vector<float> >  all_mc_vy_;
  MyPointer< std::vector<float> >  all_mc_vz_;
  MyPointer< std::vector<float> >  all_mc_endx_;
  MyPointer< std::vector<float> >  all_mc_endy_;
  MyPointer< std::vector<float> >  all_mc_endz_;
  MyPointer< std::vector<float> >  all_mc_px_;
  MyPointer< std::vector<float> >  all_mc_py_;
  MyPointer< std::vector<float> >  all_mc_pz_;
  MyPointer< std::vector<double> >  all_mc_E_;
  MyPointer< std::vector<std::__cxx11::string > > all_mc_process_;
  MyPointer< std::vector<std::__cxx11::string > > all_mc_end_process_;
  //MyPointer< std::vector<float> >  all_mc_distance_;

  //All PFP reco variables
  unsigned int slice_id_ = BOGUS_INDEX;
  //MyPointer< std::vector<int> > nonprim_pfp_parent_v_;
  //MyPointer< std::vector<int> > nonprim_pfp_ID_v_;
  MyPointer< std::vector<int> >  nonprim_slc_id_v_;
  MyPointer< std::vector<int> >  backtracked_tid_;
  MyPointer< std::vector<int> >  nonprim_backtracked_tid_;
  //int n_nonprim_pfps_ = BOGUS_INDEX;
  MyPointer< std::vector<float> >  nonprim_trk_llr_pid_score_v_;
  MyPointer< std::vector<float> >  nonprim_trk_len_v_;
  MyPointer< std::vector<float> >  nonprim_trk_distance_v_;
  MyPointer< std::vector<float> >  nonprim_trk_score_v_;
  //MyPointer< std::vector<float> >  nonprim_backtracked_purity_;

  MyPointer< std::vector<unsigned int> > nonprim_pfp_generation_;
  MyPointer< std::vector<float> >  nonprim_trk_sce_start_x_v_;
  MyPointer< std::vector<float> >  nonprim_trk_sce_start_y_v_;
  MyPointer< std::vector<float> >  nonprim_trk_sce_start_z_v_;
  MyPointer< std::vector<float> >  nonprim_trk_sce_end_x_v_;
  MyPointer< std::vector<float> >  nonprim_trk_sce_end_y_v_;
  MyPointer< std::vector<float> >  nonprim_trk_sce_end_z_v_;
  //MyPointer< std::vector<int> >  nonprim_trk_nhits_v_v_;
  //MyPointer< std::vector<int> >  nonprim_trk_nhits_u_v_;
  //MyPointer< std::vector<int> >  nonprim_trk_nhits_y_v_;
  //MyPointer< std::vector<float> >  nonprim_trk_charge_v_;
  //MyPointer<std::vector<float> >  pfp_direction_v;
  //MyPointer<std::vector<float> >  pfp_proximity_v;
  //MyPointer<std::vector<int> >  secondary_proton_trk_candidate_indices;

  int secondary_proton_candidate_idx_ = BOGUS_INDEX;
  //End Burke Additions


  //================================================================================================================
  // ** Reconstructed observables **
};
