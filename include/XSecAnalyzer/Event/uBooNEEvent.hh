/*
Description: Implementation of AnalysisEvent class for MicroBooNE
Date: 2024-10-16
Authors: Steven Gardiner
*/

#include "AnalysisEvent.hh"

class uBooNEEvent : public AnalysisEvent {
public:
  uBooNEEvent() {}
  ~uBooNEEvent() {}

  // Helper function to set branch addresses for reading information
  // from the Event TTree
  void set_event_branch_addresses(TTree& etree) override
  {
    std::cout<< "Setting branch addresses for uBooNEEvent" << std::endl;
    // Reco PDG code of primary PFParticle in slice (i.e., the neutrino
    // candidate)
    SetBranchAddress(etree, "slpdg", &this->nu_pdg_ );

    // Number of neutrino slices identified by the SliceID. Allowed values
    // are zero or one.
    //SetBranchAddress(etree,"nslice",&this->nslice_);
    SetBranchAddress(etree, "nslice", &this->nslice_ );

    // Topological score
    SetBranchAddress(etree, "topological_score", &this->topological_score_ );
    SetBranchAddress(etree, "CosmicIP", &this->cosmic_impact_parameter_ );
    //SetBranchAddress(etree, "CosmicIPAll3D", &this->CosmicIPAll3D_ );

    // Reconstructed neutrino vertex position (with corrections for
    // space charge applied)
    SetBranchAddress(etree, "reco_nu_vtx_sce_x", &this->nu_vx_ );
    SetBranchAddress(etree, "reco_nu_vtx_sce_y", &this->nu_vy_ );
    SetBranchAddress(etree, "reco_nu_vtx_sce_z", &this->nu_vz_ );

    // Reconstructed object counts
    SetBranchAddress(etree, "n_pfps", &this->num_pf_particles_ );
    SetBranchAddress(etree, "n_tracks", &this->num_tracks_ );
    SetBranchAddress(etree, "n_showers", &this->num_showers_ );

    // PFParticle properties
    set_object_input_branch_address( etree, "pfp_generation_v",
        this->pfp_generation_ );

    set_object_input_branch_address( etree, "pfp_trk_daughters_v",
        this->pfp_trk_daughters_count_ );

    set_object_input_branch_address( etree, "pfp_shr_daughters_v",
        this->pfp_shr_daughters_count_ );

    set_object_input_branch_address( etree, "trk_score_v", this->pfp_track_score_ );
    set_object_input_branch_address( etree, "pfpdg", this->pfp_reco_pdg_ );
    set_object_input_branch_address( etree, "pfnhits", this->pfp_hits_ );
    set_object_input_branch_address( etree, "pfnplanehits_U", this->pfp_hitsU_ );
    set_object_input_branch_address( etree, "pfnplanehits_V", this->pfp_hitsV_ );
    set_object_input_branch_address( etree, "pfnplanehits_Y", this->pfp_hitsY_ );

    // Backtracked PFParticle properties
    set_object_input_branch_address( etree, "backtracked_pdg", this->pfp_true_pdg_ );
    set_object_input_branch_address( etree, "backtracked_e", this->pfp_true_E_ );
    set_object_input_branch_address( etree, "backtracked_px", this->pfp_true_px_ );
    set_object_input_branch_address( etree, "backtracked_py", this->pfp_true_py_ );
    set_object_input_branch_address( etree, "backtracked_pz", this->pfp_true_pz_ );

    // Shower properties
    // These are excluded from some ntuples to ensure blindness for the LEE
    // analyses. We will skip them when not available.
    bool has_shower_branches = ( etree.GetBranch("shr_pfp_id_v") != nullptr );
    if ( has_shower_branches ) {
        set_object_input_branch_address( etree, "shr_pfp_id_v", this->shower_pfp_id_ );
        set_object_input_branch_address( etree, "shr_start_x_v", this->shower_startx_ );
        set_object_input_branch_address( etree, "shr_start_y_v", this->shower_starty_ );
        set_object_input_branch_address( etree, "shr_start_z_v", this->shower_startz_ );
        // Shower start distance from reco neutrino vertex (pre-calculated for
        // convenience)
        set_object_input_branch_address( etree, "shr_dist_v",
        this->shower_start_distance_ );
    }
    else {
        // When the shower information is not available, delete the owned vectors
        // to signal that the associated branches should not be written to the
        // output TTree
        this->shower_pfp_id_.reset( nullptr );
        this->shower_startx_.reset( nullptr );
        this->shower_starty_.reset( nullptr );
        this->shower_startz_.reset( nullptr );
        this->shower_start_distance_.reset( nullptr );
    }

    // Track properties
    set_object_input_branch_address( etree, "trk_pfp_id_v", this->track_pfp_id_ );
    set_object_input_branch_address( etree, "trk_len_v", this->track_length_ );
    set_object_input_branch_address( etree, "trk_sce_start_x_v", this->track_startx_ );
    set_object_input_branch_address( etree, "trk_sce_start_y_v", this->track_starty_ );
    set_object_input_branch_address( etree, "trk_sce_start_z_v", this->track_startz_ );

    // Track start distance from reco neutrino vertex (pre-calculated for
    // convenience)
    set_object_input_branch_address( etree, "trk_distance_v",
        this->track_start_distance_ );

    set_object_input_branch_address( etree, "trk_sce_end_x_v", this->track_endx_ );
    set_object_input_branch_address( etree, "trk_sce_end_y_v", this->track_endy_ );
    set_object_input_branch_address( etree, "trk_sce_end_z_v", this->track_endz_ );

    set_object_input_branch_address( etree, "trk_dir_x_v", this->track_dirx_ );
    set_object_input_branch_address( etree, "trk_dir_y_v", this->track_diry_ );
    set_object_input_branch_address( etree, "trk_dir_z_v", this->track_dirz_ );

    set_object_input_branch_address( etree, "trk_theta_v", this->track_theta_ );
    set_object_input_branch_address( etree, "trk_phi_v", this->track_phi_ );

    set_object_input_branch_address( etree, "trk_energy_proton_v",
        this->track_kinetic_energy_p_ );

    set_object_input_branch_address( etree, "trk_range_muon_mom_v",
        this->track_range_mom_mu_ );

    set_object_input_branch_address( etree, "trk_mcs_muon_mom_v",
        this->track_mcs_mom_mu_ );

    // Some ntuples exclude the old proton chi^2 PID score. Only include it
    // in the output if this branch is available.
    bool has_chipr = ( etree.GetBranch("trk_pid_chipr_v") != nullptr );
    if ( has_chipr ) {
        set_object_input_branch_address( etree, "trk_pid_chipr_v",
        this->track_chi2_proton_ );
    }
    else {
        this->track_chi2_proton_.reset( nullptr );
    }

    // Log-likelihood-based particle ID information
    set_object_input_branch_address( etree, "trk_llr_pid_v", this->track_llr_pid_ );

    set_object_input_branch_address( etree, "trk_llr_pid_u_v",
        this->track_llr_pid_U_ );

    set_object_input_branch_address( etree, "trk_llr_pid_v_v",
        this->track_llr_pid_V_ );

    set_object_input_branch_address( etree, "trk_llr_pid_y_v",
        this->track_llr_pid_Y_ );

    set_object_input_branch_address( etree, "trk_llr_pid_score_v",
        this->track_llr_pid_score_ );

    // MC truth information for the neutrino
    SetBranchAddress(etree, "nu_pdg", &this->mc_nu_pdg_ );
    SetBranchAddress(etree, "true_nu_vtx_x", &this->mc_nu_vx_ );
    SetBranchAddress(etree, "true_nu_vtx_y", &this->mc_nu_vy_ );
    SetBranchAddress(etree, "true_nu_vtx_z", &this->mc_nu_vz_ );
    SetBranchAddress(etree, "nu_e", &this->mc_nu_energy_ );
    SetBranchAddress(etree, "ccnc", &this->mc_nu_ccnc_ );
    SetBranchAddress(etree, "interaction", &this->mc_nu_interaction_type_ );

    //=============================================
    //DB Added to match Samantha's Signal defintion

    SetBranchAddress(etree, "true_nu_vtx_sce_x", &this->mc_nu_sce_vx_ );
    SetBranchAddress(etree, "true_nu_vtx_sce_y", &this->mc_nu_sce_vy_ );
    SetBranchAddress(etree, "true_nu_vtx_sce_z", &this->mc_nu_sce_vz_ );

    //=============================================

    // MC truth information for the final-state primary particles
    set_object_input_branch_address( etree, "mc_pdg", this->mc_nu_daughter_pdg_ );
    set_object_input_branch_address( etree, "mc_E", this->mc_nu_daughter_energy_ );
    set_object_input_branch_address( etree, "mc_px", this->mc_nu_daughter_px_ );
    set_object_input_branch_address( etree, "mc_py", this->mc_nu_daughter_py_ );
    set_object_input_branch_address( etree, "mc_pz", this->mc_nu_daughter_pz_ );

    // GENIE and other systematic variation weights
    bool has_genie_mc_weights = ( etree.GetBranch("weightSpline") != nullptr );
    if ( has_genie_mc_weights ) {
        SetBranchAddress(etree, "weightSpline", &this->spline_weight_ );
        SetBranchAddress(etree, "weightTune", &this->tuned_cv_weight_ );
    }

    bool has_weight_map = ( etree.GetBranch("weights") != nullptr );
    if ( has_weight_map ) {
        set_object_input_branch_address( etree, "weights", this->mc_weights_map_ );
    }
    else {
        this->mc_weights_map_.reset( nullptr );
    }

    // Purity and completeness of the backtracked hits in the neutrino slice
    bool has_pfp_backtracked_purity = ( etree.GetBranch("nu_purity_from_pfp")
        != nullptr );
    if ( has_pfp_backtracked_purity ) {

        SetBranchAddress(etree, "nu_completeness_from_pfp",
        &this->nu_completeness_from_pfp_ );

        SetBranchAddress(etree, "nu_purity_from_pfp", &this->nu_purity_from_pfp_ );

    }
  }

  // Helper function to set branch addresses for the output TTree
  void set_event_output_branch_addresses(TTree& out_tree, bool create = false) override
  {
    std::cout << "create " << create << std::endl;
    // Signal definition flags
    set_output_branch_address( out_tree, "is_mc", &this->is_mc_, create, "is_mc/O" );
    // Event weights
    set_output_branch_address( out_tree, "spline_weight",
      &this->spline_weight_, create, "spline_weight/F" );

    std::cout<< "(-1) setting output branch addresses for uBooNEEvent" << std::endl;
    set_output_branch_address( out_tree, "tuned_cv_weight",
      &this->tuned_cv_weight_, create, "tuned_cv_weight/F" );
    std::cout << "(0) setting output branch addresses for uBooNEEvent" << std::endl;

    // If MC weights are available, prepare to store them in the output TTree
    if (!this->mc_weights_map_.empty() ) {
      
      // Make separate branches for the various sets of systematic variation
      // weights in the map
      for ( auto& pair : *this->mc_weights_map_ ) {

        // Prepend "weight_" to the name of the vector of weights in the map
        std::string weight_branch_name = "weight_" + pair.first;

        // Store a pointer to the vector of weights (needed to set the branch
        // address properly) in the temporary map of pointers
        this->mc_weights_ptr_map_[ weight_branch_name ] = &pair.second;

        // Set the branch address for &this vector of weights
        set_object_output_branch_address< std::vector<double> >( out_tree,
          weight_branch_name, this->mc_weights_ptr_map_.at(weight_branch_name),
          create );
      }
    }
    std::cout << "(1) setting output branch addresses for uBooNEEvent" << std::endl;

    // Backtracked neutrino purity and completeness
    set_output_branch_address( out_tree, "nu_completeness_from_pfp",
      &this->nu_completeness_from_pfp_, create, "nu_completeness_from_pfp/F" );

    set_output_branch_address( out_tree, "nu_purity_from_pfp",
      &this->nu_purity_from_pfp_, create, "nu_purity_from_pfp/F" );

    // Number of neutrino slices identified by the SliceID
    set_output_branch_address( out_tree, "nslice", &this->nslice_, create,
      "nslice/I" );

    // *** Branches copied directly from the input ***

    // Cosmic rejection parameters for numu CC inclusive selection
    set_output_branch_address( out_tree, "topological_score",
      &this->topological_score_, create, "topological_score/F" );

    set_output_branch_address( out_tree, "CosmicIP",
      &this->cosmic_impact_parameter_, create, "CosmicIP/F" );

    // Reconstructed neutrino vertex position
    set_output_branch_address( out_tree, "reco_nu_vtx_sce_x",
      &this->nu_vx_, create, "reco_nu_vtx_sce_x/F" );

    set_output_branch_address( out_tree, "reco_nu_vtx_sce_y",
      &this->nu_vy_, create, "reco_nu_vtx_sce_y/F" );

    set_output_branch_address( out_tree, "reco_nu_vtx_sce_z",
      &this->nu_vz_, create, "reco_nu_vtx_sce_z/F" );

    // MC truth information for the neutrino
    set_output_branch_address( out_tree, "mc_nu_pdg", &this->mc_nu_pdg_,
      create, "mc_nu_pdg/I" );

    set_output_branch_address( out_tree, "mc_nu_vtx_x", &this->mc_nu_vx_,
      create, "mc_nu_vtx_x/F" );

    set_output_branch_address( out_tree, "mc_nu_vtx_y", &this->mc_nu_vy_,
      create, "mc_nu_vtx_y/F" );

    set_output_branch_address( out_tree, "mc_nu_vtx_z", &this->mc_nu_vz_,
      create, "mc_nu_vtx_z/F" );

    set_output_branch_address( out_tree, "mc_nu_energy", &this->mc_nu_energy_,
      create, "mc_nu_energy/F" );

    set_output_branch_address( out_tree, "mc_ccnc", &this->mc_nu_ccnc_,
      create, "mc_ccnc/I" );

    set_output_branch_address( out_tree, "mc_interaction",
      &this->mc_nu_interaction_type_, create, "mc_interaction/I" );

    // PFParticle properties
    set_object_output_branch_address< std::vector<unsigned int> >( out_tree,
      "pfp_generation_v", this->pfp_generation_, create );

    set_object_output_branch_address< std::vector<unsigned int> >( out_tree,
      "pfp_trk_daughters_v", this->pfp_trk_daughters_count_, create );

    set_object_output_branch_address< std::vector<unsigned int> >( out_tree,
      "pfp_shr_daughters_v", this->pfp_shr_daughters_count_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_score_v", this->pfp_track_score_, create );

    set_object_output_branch_address< std::vector<int> >( out_tree,
      "pfpdg", this->pfp_reco_pdg_, create );

    set_object_output_branch_address< std::vector<int> >( out_tree,
      "pfnhits", this->pfp_hits_, create );

    set_object_output_branch_address< std::vector<int> >( out_tree,
      "pfnplanehits_U", this->pfp_hitsU_, create );

    set_object_output_branch_address< std::vector<int> >( out_tree,
      "pfnplanehits_V", this->pfp_hitsV_, create );

    set_object_output_branch_address< std::vector<int> >( out_tree,
      "pfnplanehits_Y", this->pfp_hitsY_, create );

    // Backtracked PFParticle properties
    set_object_output_branch_address< std::vector<int> >( out_tree,
      "backtracked_pdg", this->pfp_true_pdg_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "backtracked_e", this->pfp_true_E_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "backtracked_px", this->pfp_true_px_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "backtracked_py", this->pfp_true_py_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "backtracked_pz", this->pfp_true_pz_, create );

    // Shower properties
    // For some ntuples, reconstructed shower information is excluded.
    // In such cases, skip writing these branches to the output TTree.
    if ( this->shower_startx_ ) {
      set_object_output_branch_address< std::vector<float> >( out_tree,
        "shr_start_x_v", this->shower_startx_, create );

      set_object_output_branch_address< std::vector<float> >( out_tree,
        "shr_start_y_v", this->shower_starty_, create );

      set_object_output_branch_address< std::vector<float> >( out_tree,
        "shr_start_z_v", this->shower_startz_, create );

      // Shower start distance from reco neutrino vertex (pre-calculated for
      // convenience)
      set_object_output_branch_address< std::vector<float> >( out_tree,
        "shr_dist_v", this->shower_start_distance_, create );
    }

    // Track properties
    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_len_v", this->track_length_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_sce_start_x_v", this->track_startx_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_sce_start_y_v", this->track_starty_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_sce_start_z_v", this->track_startz_, create );

    // Track start distance from reco neutrino vertex (pre-calculated for
    // convenience)
    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_distance_v", this->track_start_distance_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_sce_end_x_v", this->track_endx_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_sce_end_y_v", this->track_endy_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_sce_end_z_v", this->track_endz_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_dir_x_v", this->track_dirx_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_dir_y_v", this->track_diry_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_dir_z_v", this->track_dirz_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_energy_proton_v", this->track_kinetic_energy_p_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_range_muon_mom_v", this->track_range_mom_mu_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_mcs_muon_mom_v", this->track_mcs_mom_mu_, create );

    // Some ntuples exclude the old chi^2 proton PID score. Only include it in
    // the output if it is available.
    if ( this->track_chi2_proton_ ) {
      set_object_output_branch_address< std::vector<float> >( out_tree,
        "trk_pid_chipr_v", this->track_chi2_proton_, create );
    }

    // Log-likelihood-based particle ID information
    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_llr_pid_v", this->track_llr_pid_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_llr_pid_u_v", this->track_llr_pid_U_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_llr_pid_v_v", this->track_llr_pid_V_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_llr_pid_y_v", this->track_llr_pid_Y_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree,
      "trk_llr_pid_score_v", this->track_llr_pid_score_, create );

    // MC truth information for the final-state primary particles
    set_object_output_branch_address< std::vector<int> >( out_tree, "mc_pdg",
      this->mc_nu_daughter_pdg_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree, "mc_E",
      this->mc_nu_daughter_energy_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree, "mc_px",
      this->mc_nu_daughter_px_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree, "mc_py",
      this->mc_nu_daughter_py_, create );

    set_object_output_branch_address< std::vector<float> >( out_tree, "mc_pz",
      this->mc_nu_daughter_pz_, create );
  }

};