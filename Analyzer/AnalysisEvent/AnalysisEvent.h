#ifndef __ANALYSIS_EVENT_H__
#define __ANALYSIS_EVENT_H__

// STV analysis includes
#include "EventCategory.hh"
#include "TreeUtils.hh"
#include "FiducialVolume.hh"
#include "Constants.h"

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
  
  EventCategory category_ = kUnknown;
  
  //================================================================================================================
  // ** Reconstructed observables **
  
  // 3-momenta
  MyPointer< TVector3 > p3_mu_;
  MyPointer< TVector3 > p3_lead_p_;
  
  // Reconstructed 3-momenta for all proton candidates,
  // ordered from highest to lowest by magnitude
  MyPointer< std::vector<TVector3> > p3_p_vec_;
  
  // Reco STVs and other variables of interest
  float delta_pT_ = BOGUS;
  float delta_phiT_ = BOGUS;
  float delta_alphaT_ = BOGUS;
  float delta_pL_ = BOGUS;
  float pn_ = BOGUS;
  float delta_pTx_ = BOGUS;
  float delta_pTy_ = BOGUS;
  float theta_mu_p_ = BOGUS;
  
  // ** MC truth observables **
  // These are loaded for signal events whenever we have MC information
  // to use
  
  // 3-momenta
  MyPointer< TVector3 > mc_p3_mu_;
  MyPointer< TVector3 > mc_p3_lead_p_;
  
  // True 3-momenta for all true MC protons, ordered from highest to lowest
  // by magnitude
  MyPointer< std::vector<TVector3> > mc_p3_p_vec_;
  
  // MC truth STVs and other variables of interest
  float mc_delta_pT_ = BOGUS;
  float mc_delta_phiT_ = BOGUS;
  float mc_delta_alphaT_ = BOGUS;
  float mc_delta_pL_ = BOGUS;
  float mc_pn_ = BOGUS;
  float mc_delta_pTx_ = BOGUS;
  float mc_delta_pTy_ = BOGUS;
  float mc_theta_mu_p_ = BOGUS;
};

#endif
