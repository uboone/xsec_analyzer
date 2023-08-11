// Analysis macro for use in the CCNp0pi single transverse variable analysis
// Designed for use with the PeLEE group's "searchingfornues" ntuples
//
// Updated 22 April 2023
// Steven Gardiner <gardiner@fnal.gov>

// Standard library includes
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

// ROOT includes
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"

//Dan's includes
#include "AnalysisEvent.h"
#include "Constants.h"
#include "Functions.h"

#include "SelectionBase.h"
#include "CC1mu1p0pi.h"
#include "CC1mu2p0pi.h"
#include "CC1muNp0pi.h"

void SetBranchAddress(TTree& etree, std::string BranchName, void* Variable) {
  etree.SetBranchAddress(BranchName.c_str(),Variable);
}

// Helper function to set branch addresses for reading information
// from the Event TTree
void set_event_branch_addresses(TTree& etree, AnalysisEvent& ev)
{
  // Reco PDG code of primary PFParticle in slice (i.e., the neutrino
  // candidate)
  SetBranchAddress(etree, "slpdg", &ev.nu_pdg_ );

  // Number of neutrino slices identified by the SliceID. Allowed values
  // are zero or one.
  //SetBranchAddress(etree,"nslice",&ev.nslice_);
  SetBranchAddress(etree, "nslice", &ev.nslice_ );

  // Topological score
  SetBranchAddress(etree, "topological_score", &ev.topological_score_ );
  SetBranchAddress(etree, "CosmicIP", &ev.cosmic_impact_parameter_ );
  //SetBranchAddress(etree, "CosmicIPAll3D", &ev.CosmicIPAll3D_ );

  // Reconstructed neutrino vertex position (with corrections for
  // space charge applied)
  SetBranchAddress(etree, "reco_nu_vtx_sce_x", &ev.nu_vx_ );
  SetBranchAddress(etree, "reco_nu_vtx_sce_y", &ev.nu_vy_ );
  SetBranchAddress(etree, "reco_nu_vtx_sce_z", &ev.nu_vz_ );

  // Reconstructed object counts
  SetBranchAddress(etree, "n_pfps", &ev.num_pf_particles_ );
  SetBranchAddress(etree, "n_tracks", &ev.num_tracks_ );
  SetBranchAddress(etree, "n_showers", &ev.num_showers_ );

  // PFParticle properties
  set_object_input_branch_address( etree, "pfp_generation_v",
    ev.pfp_generation_ );

  set_object_input_branch_address( etree, "pfp_trk_daughters_v",
    ev.pfp_trk_daughters_count_ );

  set_object_input_branch_address( etree, "pfp_shr_daughters_v",
    ev.pfp_shr_daughters_count_ );

  set_object_input_branch_address( etree, "trk_score_v", ev.pfp_track_score_ );
  set_object_input_branch_address( etree, "pfpdg", ev.pfp_reco_pdg_ );
  set_object_input_branch_address( etree, "pfnhits", ev.pfp_hits_ );
  set_object_input_branch_address( etree, "pfnplanehits_U", ev.pfp_hitsU_ );
  set_object_input_branch_address( etree, "pfnplanehits_V", ev.pfp_hitsV_ );
  set_object_input_branch_address( etree, "pfnplanehits_Y", ev.pfp_hitsY_ );

  // Backtracked PFParticle properties
  set_object_input_branch_address( etree, "backtracked_pdg", ev.pfp_true_pdg_ );
  set_object_input_branch_address( etree, "backtracked_e", ev.pfp_true_E_ );
  set_object_input_branch_address( etree, "backtracked_px", ev.pfp_true_px_ );
  set_object_input_branch_address( etree, "backtracked_py", ev.pfp_true_py_ );
  set_object_input_branch_address( etree, "backtracked_pz", ev.pfp_true_pz_ );

  // Shower properties
  // These are excluded from some ntuples to ensure blindness for the LEE
  // analyses. We will skip them when not available.
  bool has_shower_branches = ( etree.GetBranch("shr_pfp_id_v") != nullptr );
  if ( has_shower_branches ) {
    set_object_input_branch_address( etree, "shr_pfp_id_v", ev.shower_pfp_id_ );
    set_object_input_branch_address( etree, "shr_start_x_v", ev.shower_startx_ );
    set_object_input_branch_address( etree, "shr_start_y_v", ev.shower_starty_ );
    set_object_input_branch_address( etree, "shr_start_z_v", ev.shower_startz_ );
    // Shower start distance from reco neutrino vertex (pre-calculated for
    // convenience)
    set_object_input_branch_address( etree, "shr_dist_v",
      ev.shower_start_distance_ );
  }
  else {
    // When the shower information is not available, delete the owned vectors
    // to signal that the associated branches should not be written to the
    // output TTree
    ev.shower_pfp_id_.reset( nullptr );
    ev.shower_startx_.reset( nullptr );
    ev.shower_starty_.reset( nullptr );
    ev.shower_startz_.reset( nullptr );
    ev.shower_start_distance_.reset( nullptr );
  }

  // Track properties
  set_object_input_branch_address( etree, "trk_pfp_id_v", ev.track_pfp_id_ );
  set_object_input_branch_address( etree, "trk_len_v", ev.track_length_ );
  set_object_input_branch_address( etree, "trk_sce_start_x_v", ev.track_startx_ );
  set_object_input_branch_address( etree, "trk_sce_start_y_v", ev.track_starty_ );
  set_object_input_branch_address( etree, "trk_sce_start_z_v", ev.track_startz_ );

  // Track start distance from reco neutrino vertex (pre-calculated for
  // convenience)
  set_object_input_branch_address( etree, "trk_distance_v",
    ev.track_start_distance_ );

  set_object_input_branch_address( etree, "trk_sce_end_x_v", ev.track_endx_ );
  set_object_input_branch_address( etree, "trk_sce_end_y_v", ev.track_endy_ );
  set_object_input_branch_address( etree, "trk_sce_end_z_v", ev.track_endz_ );

  set_object_input_branch_address( etree, "trk_dir_x_v", ev.track_dirx_ );
  set_object_input_branch_address( etree, "trk_dir_y_v", ev.track_diry_ );
  set_object_input_branch_address( etree, "trk_dir_z_v", ev.track_dirz_ );

  set_object_input_branch_address( etree, "trk_energy_proton_v",
    ev.track_kinetic_energy_p_ );

  set_object_input_branch_address( etree, "trk_range_muon_mom_v",
    ev.track_range_mom_mu_ );

  set_object_input_branch_address( etree, "trk_mcs_muon_mom_v",
    ev.track_mcs_mom_mu_ );

  // Some ntuples exclude the old proton chi^2 PID score. Only include it
  // in the output if this branch is available.
  bool has_chipr = ( etree.GetBranch("trk_pid_chipr_v") != nullptr );
  if ( has_chipr ) {
    set_object_input_branch_address( etree, "trk_pid_chipr_v",
      ev.track_chi2_proton_ );
  }
  else {
    ev.track_chi2_proton_.reset( nullptr );
  }

  // Log-likelihood-based particle ID information
  set_object_input_branch_address( etree, "trk_llr_pid_v", ev.track_llr_pid_ );

  set_object_input_branch_address( etree, "trk_llr_pid_u_v",
    ev.track_llr_pid_U_ );

  set_object_input_branch_address( etree, "trk_llr_pid_v_v",
    ev.track_llr_pid_V_ );

  set_object_input_branch_address( etree, "trk_llr_pid_y_v",
    ev.track_llr_pid_Y_ );

  set_object_input_branch_address( etree, "trk_llr_pid_score_v",
    ev.track_llr_pid_score_ );

  // MC truth information for the neutrino
  SetBranchAddress(etree, "nu_pdg", &ev.mc_nu_pdg_ );
  SetBranchAddress(etree, "true_nu_vtx_x", &ev.mc_nu_vx_ );
  SetBranchAddress(etree, "true_nu_vtx_y", &ev.mc_nu_vy_ );
  SetBranchAddress(etree, "true_nu_vtx_z", &ev.mc_nu_vz_ );
  SetBranchAddress(etree, "nu_e", &ev.mc_nu_energy_ );
  SetBranchAddress(etree, "ccnc", &ev.mc_nu_ccnc_ );
  SetBranchAddress(etree, "interaction", &ev.mc_nu_interaction_type_ );

  // MC truth information for the final-state primary particles
  set_object_input_branch_address( etree, "mc_pdg", ev.mc_nu_daughter_pdg_ );
  set_object_input_branch_address( etree, "mc_E", ev.mc_nu_daughter_energy_ );
  set_object_input_branch_address( etree, "mc_px", ev.mc_nu_daughter_px_ );
  set_object_input_branch_address( etree, "mc_py", ev.mc_nu_daughter_py_ );
  set_object_input_branch_address( etree, "mc_pz", ev.mc_nu_daughter_pz_ );

  // GENIE and other systematic variation weights
  bool has_genie_mc_weights = ( etree.GetBranch("weightSpline") != nullptr );
  if ( has_genie_mc_weights ) {
    SetBranchAddress(etree, "weightSpline", &ev.spline_weight_ );
    SetBranchAddress(etree, "weightTune", &ev.tuned_cv_weight_ );
  }

  bool has_weight_map = ( etree.GetBranch("weights") != nullptr );
  if ( has_weight_map ) {
    set_object_input_branch_address( etree, "weights", ev.mc_weights_map_ );
  }
  else {
    ev.mc_weights_map_.reset( nullptr );
  }

  // Purity and completeness of the backtracked hits in the neutrino slice
  bool has_pfp_backtracked_purity = ( etree.GetBranch("nu_purity_from_pfp")
    != nullptr );
  if ( has_pfp_backtracked_purity ) {

    SetBranchAddress(etree, "nu_completeness_from_pfp",
      &ev.nu_completeness_from_pfp_ );

    SetBranchAddress(etree, "nu_purity_from_pfp", &ev.nu_purity_from_pfp_ );

  }

}

// Helper function to set branch addresses for the output TTree
void set_event_output_branch_addresses(TTree& out_tree, AnalysisEvent& ev,
  bool create = false)
{
  // Signal definition flags
  set_output_branch_address( out_tree, "is_mc", &ev.is_mc_, create, "is_mc/O" );

  set_output_branch_address( out_tree, "mc_neutrino_is_numu",
    &ev.mc_neutrino_is_numu_, create, "mc_neutrino_is_numu/O" );

  set_output_branch_address( out_tree, "mc_vertex_in_FV",
    &ev.mc_vertex_in_FV_, create, "mc_vertex_in_FV/O" );

  set_output_branch_address( out_tree, "mc_muon_in_mom_range",
    &ev.mc_muon_in_mom_range_, create, "mc_muon_in_mom_range/O" );

  set_output_branch_address( out_tree, "mc_lead_p_in_mom_range",
    &ev.mc_lead_p_in_mom_range_, create, "mc_lead_p_in_mom_range/O" );

  set_output_branch_address( out_tree, "mc_no_fs_pi0",
    &ev.mc_no_fs_pi0_, create, "mc_no_fs_pi0/O" );

  set_output_branch_address( out_tree, "mc_no_charged_pi_above_threshold",
    &ev.mc_no_charged_pi_above_threshold_, create,
    "mc_no_charged_pi_above_threshold/O" );

  set_output_branch_address( out_tree, "mc_no_fs_mesons",
    &ev.mc_no_fs_mesons_, create, "mc_no_fs_mesons/O" );

  set_output_branch_address( out_tree, "mc_is_signal",
    &ev.mc_is_signal_, create, "mc_is_signal/O" );

  // MC event category
  set_output_branch_address( out_tree, "category",
    &ev.category_, create, "category/I" );

  // Event weights
  set_output_branch_address( out_tree, "spline_weight",
    &ev.spline_weight_, create, "spline_weight/F" );

  set_output_branch_address( out_tree, "tuned_cv_weight",
    &ev.tuned_cv_weight_, create, "tuned_cv_weight/F" );

  // If MC weights are available, prepare to store them in the output TTree
  if ( ev.mc_weights_map_ ) {

    // Make separate branches for the various sets of systematic variation
    // weights in the map
    for ( auto& pair : *ev.mc_weights_map_ ) {

      // Prepend "weight_" to the name of the vector of weights in the map
      std::string weight_branch_name = "weight_" + pair.first;

      // Store a pointer to the vector of weights (needed to set the branch
      // address properly) in the temporary map of pointers
      ev.mc_weights_ptr_map_[ weight_branch_name ] = &pair.second;

      // Set the branch address for this vector of weights
      set_object_output_branch_address< std::vector<double> >( out_tree,
        weight_branch_name, ev.mc_weights_ptr_map_.at(weight_branch_name),
        create );
    }
  }

  // Backtracked neutrino purity and completeness
  set_output_branch_address( out_tree, "nu_completeness_from_pfp",
    &ev.nu_completeness_from_pfp_, create, "nu_completeness_from_pfp/F" );

  set_output_branch_address( out_tree, "nu_purity_from_pfp",
    &ev.nu_purity_from_pfp_, create, "nu_purity_from_pfp/F" );

  // Number of neutrino slices identified by the SliceID
  set_output_branch_address( out_tree, "nslice", &ev.nslice_, create,
    "nslice/I" );

  /*
  // Reco 3-momenta (muon, leading proton)
  set_object_output_branch_address< TVector3 >( out_tree,
    "p3_mu", ev.p3_mu_, create );

  set_object_output_branch_address< TVector3 >( out_tree,
    "p3_lead_p", ev.p3_lead_p_, create );

  // Reco 3-momenta (all proton candidates, ordered from highest to lowest
  // magnitude)
  set_object_output_branch_address< std::vector<TVector3> >( out_tree,
    "p3_p_vec", ev.p3_p_vec_, create );

  // True 3-momenta (muon, leading proton)
  set_object_output_branch_address< TVector3 >( out_tree,
    "mc_p3_mu", ev.mc_p3_mu_, create );

  set_object_output_branch_address< TVector3 >( out_tree,
    "mc_p3_lead_p", ev.mc_p3_lead_p_, create );

  // True 3-momenta (all protons, ordered from highest to lowest magnitude)
  set_object_output_branch_address< std::vector<TVector3> >( out_tree,
    "mc_p3_p_vec", ev.mc_p3_p_vec_, create );

  // Reco STVs
  set_output_branch_address( out_tree, "delta_pT",
    &ev.delta_pT_, create, "delta_pT/F" );

  set_output_branch_address( out_tree, "delta_phiT",
    &ev.delta_phiT_, create, "delta_phiT/F" );

  set_output_branch_address( out_tree, "delta_alphaT",
    &ev.delta_alphaT_, create, "delta_alphaT/F" );

  set_output_branch_address( out_tree, "delta_pL",
    &ev.delta_pL_, create, "delta_pL/F" );

  set_output_branch_address( out_tree, "pn",
    &ev.pn_, create, "pn/F" );

  set_output_branch_address( out_tree, "delta_pTx",
    &ev.delta_pTx_, create, "delta_pTx/F" );

  set_output_branch_address( out_tree, "delta_pTy",
    &ev.delta_pTy_, create, "delta_pTy/F" );

  set_output_branch_address( out_tree, "theta_mu_p",
    &ev.theta_mu_p_, create, "theta_mu_p/F" );

  // MC STVs (only filled for signal events)
  set_output_branch_address( out_tree, "mc_delta_pT",
    &ev.mc_delta_pT_, create, "mc_delta_pT/F" );

  set_output_branch_address( out_tree, "mc_delta_phiT",
    &ev.mc_delta_phiT_, create, "mc_delta_phiT/F" );

  set_output_branch_address( out_tree, "mc_delta_alphaT",
    &ev.mc_delta_alphaT_, create, "mc_delta_alphaT/F" );

  set_output_branch_address( out_tree, "mc_delta_pL",
    &ev.mc_delta_pL_, create, "mc_delta_pL/F" );

  set_output_branch_address( out_tree, "mc_pn",
    &ev.mc_pn_, create, "mc_pn/F" );

  set_output_branch_address( out_tree, "mc_delta_pTx",
    &ev.mc_delta_pTx_, create, "mc_delta_pTx/F" );

  set_output_branch_address( out_tree, "mc_delta_pTy",
    &ev.mc_delta_pTy_, create, "mc_delta_pTy/F" );

  set_output_branch_address( out_tree, "mc_theta_mu_p",
    &ev.mc_theta_mu_p_, create, "mc_theta_mu_p/F" );
  */
  
  // *** Branches copied directly from the input ***

  // Cosmic rejection parameters for numu CC inclusive selection
  set_output_branch_address( out_tree, "topological_score",
    &ev.topological_score_, create, "topological_score/F" );

  set_output_branch_address( out_tree, "CosmicIP",
    &ev.cosmic_impact_parameter_, create, "CosmicIP/F" );

  // Reconstructed neutrino vertex position
  set_output_branch_address( out_tree, "reco_nu_vtx_sce_x",
    &ev.nu_vx_, create, "reco_nu_vtx_sce_x/F" );

  set_output_branch_address( out_tree, "reco_nu_vtx_sce_y",
    &ev.nu_vy_, create, "reco_nu_vtx_sce_y/F" );

  set_output_branch_address( out_tree, "reco_nu_vtx_sce_z",
    &ev.nu_vz_, create, "reco_nu_vtx_sce_z/F" );

  // MC truth information for the neutrino
  set_output_branch_address( out_tree, "mc_nu_pdg", &ev.mc_nu_pdg_,
    create, "mc_nu_pdg/I" );

  set_output_branch_address( out_tree, "mc_nu_vtx_x", &ev.mc_nu_vx_,
    create, "mc_nu_vtx_x/F" );

  set_output_branch_address( out_tree, "mc_nu_vtx_y", &ev.mc_nu_vy_,
    create, "mc_nu_vtx_y/F" );

  set_output_branch_address( out_tree, "mc_nu_vtx_z", &ev.mc_nu_vz_,
    create, "mc_nu_vtx_z/F" );

  set_output_branch_address( out_tree, "mc_nu_energy", &ev.mc_nu_energy_,
    create, "mc_nu_energy/F" );

  set_output_branch_address( out_tree, "mc_ccnc", &ev.mc_nu_ccnc_,
    create, "mc_ccnc/I" );

  set_output_branch_address( out_tree, "mc_interaction",
    &ev.mc_nu_interaction_type_, create, "mc_interaction/I" );

  // PFParticle properties
  set_object_output_branch_address< std::vector<unsigned int> >( out_tree,
    "pfp_generation_v", ev.pfp_generation_, create );

  set_object_output_branch_address< std::vector<unsigned int> >( out_tree,
    "pfp_trk_daughters_v", ev.pfp_trk_daughters_count_, create );

  set_object_output_branch_address< std::vector<unsigned int> >( out_tree,
    "pfp_shr_daughters_v", ev.pfp_shr_daughters_count_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_score_v", ev.pfp_track_score_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "pfpdg", ev.pfp_reco_pdg_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "pfnhits", ev.pfp_hits_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "pfnplanehits_U", ev.pfp_hitsU_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "pfnplanehits_V", ev.pfp_hitsV_, create );

  set_object_output_branch_address< std::vector<int> >( out_tree,
    "pfnplanehits_Y", ev.pfp_hitsY_, create );

  // Backtracked PFParticle properties
  set_object_output_branch_address< std::vector<int> >( out_tree,
    "backtracked_pdg", ev.pfp_true_pdg_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "backtracked_e", ev.pfp_true_E_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "backtracked_px", ev.pfp_true_px_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "backtracked_py", ev.pfp_true_py_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "backtracked_pz", ev.pfp_true_pz_, create );

  // Shower properties
  // For some ntuples, reconstructed shower information is excluded.
  // In such cases, skip writing these branches to the output TTree.
  if ( ev.shower_startx_ ) {
    set_object_output_branch_address< std::vector<float> >( out_tree,
      "shr_start_x_v", ev.shower_startx_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "shr_start_y_v", ev.shower_starty_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "shr_start_z_v", ev.shower_startz_, create );

    // Shower start distance from reco neutrino vertex (pre-calculated for
    // convenience)
    set_object_output_branch_address< std::vector<float> >( out_tree,
      "shr_dist_v", ev.shower_start_distance_, create );
  }

  // Track properties
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_len_v", ev.track_length_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_start_x_v", ev.track_startx_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_start_y_v", ev.track_starty_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_start_z_v", ev.track_startz_, create );

  // Track start distance from reco neutrino vertex (pre-calculated for
  // convenience)
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_distance_v", ev.track_start_distance_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_end_x_v", ev.track_endx_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_end_y_v", ev.track_endy_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_sce_end_z_v", ev.track_endz_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_dir_x_v", ev.track_dirx_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_dir_y_v", ev.track_diry_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_dir_z_v", ev.track_dirz_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_energy_proton_v", ev.track_kinetic_energy_p_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_range_muon_mom_v", ev.track_range_mom_mu_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_mcs_muon_mom_v", ev.track_mcs_mom_mu_, create );

  // Some ntuples exclude the old chi^2 proton PID score. Only include it in
  // the output if it is available.
  if ( ev.track_chi2_proton_ ) {
    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_pid_chipr_v", ev.track_chi2_proton_, create );
  }

  // Log-likelihood-based particle ID information
  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_llr_pid_v", ev.track_llr_pid_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_llr_pid_u_v", ev.track_llr_pid_U_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_llr_pid_v_v", ev.track_llr_pid_V_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_llr_pid_y_v", ev.track_llr_pid_Y_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree,
    "trk_llr_pid_score_v", ev.track_llr_pid_score_, create );

  // MC truth information for the final-state primary particles
  set_object_output_branch_address< std::vector<int> >( out_tree, "mc_pdg",
    ev.mc_nu_daughter_pdg_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree, "mc_E",
    ev.mc_nu_daughter_energy_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree, "mc_px",
    ev.mc_nu_daughter_px_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree, "mc_py",
    ev.mc_nu_daughter_py_, create );

  set_object_output_branch_address< std::vector<float> >( out_tree, "mc_pz",
    ev.mc_nu_daughter_pz_, create );
}

void analyze(const std::vector<std::string>& in_file_names,
  const std::string& output_filename)
{
  // Get the TTrees containing the event ntuples and subrun POT information
  // Use TChain objects for simplicity in manipulating multiple files
  TChain events_ch( "nuselection/NeutrinoSelectionFilter" );
  TChain subruns_ch( "nuselection/SubRun" );

  for ( const auto& f_name : in_file_names ) {
    events_ch.Add( f_name.c_str() );
    subruns_ch.Add( f_name.c_str() );
  }
  
  // OUTPUT TTREE
  // Make an output TTree for plotting (one entry per event)
  TFile* out_file = new TFile( output_filename.c_str(), "recreate" );
  out_file->cd();
  TTree* out_tree = new TTree( "stv_tree", "STV analysis tree" );

  // Get the total POT from the subruns TTree. Save it in the output
  // TFile as a TParameter<float>. Real data doesn't have this TTree,
  // so check that it exists first.
  float pot;
  float summed_pot = 0.;
  bool has_pot_branch = ( subruns_ch.GetBranch("pot") != nullptr );
  if ( has_pot_branch ) {
    subruns_ch.SetBranchAddress( "pot", &pot );
    for ( int se = 0; se < subruns_ch.GetEntries(); ++se ) {
      subruns_ch.GetEntry( se );
      summed_pot += pot;
    }
  }

  TParameter<float>* summed_pot_param = new TParameter<float>( "summed_pot",
    summed_pot );

  summed_pot_param->Write();

  std::vector<SelectionBase*> Selections;

  CC1mu1p0pi* CC1mu1p0piObj = new CC1mu1p0pi();
  Selections.push_back((SelectionBase*)CC1mu1p0piObj);
  
  CC1mu2p0pi* CC1mu2p0piObj = new CC1mu2p0pi();
  Selections.push_back((SelectionBase*)CC1mu2p0piObj);
  
  CC1muNp0pi* CC1muNp0piObj = new CC1muNp0pi();
  Selections.push_back((SelectionBase*)CC1muNp0piObj);

  for (size_t i=0;i<Selections.size();i++) {
    Selections[i]->SetupTree(out_tree);
  }

  // EVENT LOOP
  // TChains can potentially be really big (and spread out over multiple
  // files). When that's the case, calling TChain::GetEntries() can be very
  // slow. I get around this by using a while loop instead of a for loop.
  bool created_output_branches = false;
  long events_entry = 0;

  while ( true ) {

    //if ( events_entry > 10000) break;
    
    if ( events_entry % 1000 == 0 ) {
      std::cout << "Processing event #" << events_entry << '\n';
    }
    
    // Create a new AnalysisEvent object. This will reset all analysis
    // variables for the current event.
    AnalysisEvent cur_event;

    // Set branch addresses for the member variables that will be read
    // directly from the Event TTree.
    set_event_branch_addresses( events_ch, cur_event );

    // TChain::LoadTree() returns the entry number that should be used with
    // the current TTree object, which (together with the TBranch objects
    // that it owns) doesn't know about the other TTrees in the TChain.
    // If the return value is negative, there was an I/O error, or we've
    // attempted to read past the end of the TChain.
    int local_entry = events_ch.LoadTree( events_entry );

    // If we've reached the end of the TChain (or encountered an I/O error),
    // then terminate the event loop
    if ( local_entry < 0 ) break;

    // Load all of the branches for which we've called
    // TChain::SetBranchAddress() above
    events_ch.GetEntry( events_entry );

    // Set the output TTree branch addresses, creating the branches if needed
    // (during the first event loop iteration)
    bool create_them = false;
    if ( !created_output_branches ) {
      create_them = true;
      created_output_branches = true;
    }
    set_event_output_branch_addresses(*out_tree, cur_event, create_them );

    cur_event.categorize_event();
    
    for (size_t i=0;i<Selections.size();i++) {
      Selections[i]->ApplySelection(&(cur_event));
    }
    
    // We're done. Save the results and move on to the next event.
    out_tree->Fill();
    ++events_entry;
  }
  
  for (size_t i=0;i<Selections.size();i++) {
    Selections[i]->Print();
  }
  
  out_tree->Write();
  out_file->Close();
  delete out_file;
}

void analyzer(const std::string& in_file_name,
 const std::string& output_filename)
{
  std::vector<std::string> in_files = { in_file_name };
  analyze( in_files, output_filename );
}

int main( int argc, char* argv[] ) {

  if ( argc != 3 ) {
    std::cout << "Usage: analyzer INPUT_PELEE_NTUPLE_FILE OUTPUT_FILE\n";
    return 1;
  }

  std::string input_file_name( argv[1] );
  std::string output_file_name( argv[2] );

  analyzer( input_file_name, output_file_name );

  return 0;
}
