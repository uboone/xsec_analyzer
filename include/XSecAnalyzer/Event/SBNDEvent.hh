/*
Description: Implementation of AnalysisEvent class for SBND
Date: 2024-10-16
Authors: Brinden Carlson (bcarlson1@ufl.edu)
*/

#include "AnalysisEvent.hh"
#include "XSecAnalyzer/SBND/Constants.hh"

class SBNDEvent : public AnalysisEvent {
public:
  SBNDEvent() {
    reserve_event_branch_space();
  }
  ~SBNDEvent() {}

  // Branches

  // Slice properties, overwrite
  MyPointer< std::vector<float> > topological_score_;

  // Shower properties
  MyPointer< std::vector<float> > shower_dirx_;
  MyPointer< std::vector<float> > shower_diry_;
  MyPointer< std::vector<float> > shower_dirz_;
  // -other properties (still important)
  MyPointer< std::vector<float> > shower_openangle_;
  MyPointer< std::vector<float> > shower_energy_;
  MyPointer< std::vector<float> > shower_dEdx_;

  // Track properties

  // -range momentum
  MyPointer< std::vector<float> > track_range_mom_pi_;
  MyPointer< std::vector<float> > track_range_mom_p_;
  // -MCS momentum
  MyPointer< std::vector<float> > track_mcs_bwd_mom_pi_;
  MyPointer< std::vector<float> > track_mcs_bwd_mom_p_;
  MyPointer< std::vector<float> > track_mcs_bwd_mom_mu_;
  MyPointer< std::vector<float> > track_mcs_fwd_mom_pi_;
  MyPointer< std::vector<float> > track_mcs_fwd_mom_p_;
  MyPointer< std::vector<float> > track_mcs_fwd_mom_mu_;
  // -chi2 proton
  MyPointer< std::vector<float> > track_chi2_u_proton_;
  MyPointer< std::vector<float> > track_chi2_v_proton_;
  MyPointer< std::vector<float> > track_chi2_y_proton_;
  // -chi2 muon
  MyPointer< std::vector<float> > track_chi2_u_muon_;
  MyPointer< std::vector<float> > track_chi2_v_muon_;
  MyPointer< std::vector<float> > track_chi2_y_muon_;
  // -chi2 pion
  MyPointer< std::vector<float> > track_chi2_u_pion_;
  MyPointer< std::vector<float> > track_chi2_v_pion_;
  MyPointer< std::vector<float> > track_chi2_y_pion_;
  // -other properties (still important)
  MyPointer< std::vector<float> > track_bestplane_;

  // Final-state primary particle properties
  MyPointer< std::vector<unsigned long> > mc_nu_daughter_startx_;
  MyPointer< std::vector<unsigned long> > mc_nu_daughter_starty_;
  MyPointer< std::vector<unsigned long> > mc_nu_daughter_startz_;

  //Redefine some branches to have the proper type

  //Helper function to reserve space for the vectors
  void reserve_event_branch_space(){

    //PFP branches
    this->pfp_generation_->reserve(PFP_SIZE);
    this->pfp_trk_daughters_count_->reserve(PFP_SIZE);
    this->pfp_track_score_->reserve(PFP_SIZE);
    this->pfp_hitsU_->reserve(PFP_SIZE);
    this->pfp_hitsV_->reserve(PFP_SIZE);
    this->pfp_hitsY_->reserve(PFP_SIZE);
    this->pfp_true_pdg_->reserve(PFP_SIZE);
    this->pfp_true_E_->reserve(PFP_SIZE);
    this->pfp_true_px_->reserve(PFP_SIZE);
    this->pfp_true_py_->reserve(PFP_SIZE);
    this->pfp_true_pz_->reserve(PFP_SIZE);
    this->shower_pfp_id_->reserve(PFP_SIZE);
    this->shower_startx_->reserve(PFP_SIZE);
    this->shower_starty_->reserve(PFP_SIZE);
    this->shower_startz_->reserve(PFP_SIZE);
    this->shower_dirx_->reserve(PFP_SIZE);
    this->shower_diry_->reserve(PFP_SIZE);
    this->shower_dirz_->reserve(PFP_SIZE);
    this->shower_openangle_->reserve(PFP_SIZE);
    this->shower_start_distance_->reserve(PFP_SIZE);
    this->shower_energy_->reserve(PFP_SIZE);
    this->shower_dEdx_->reserve(PFP_SIZE);
    this->track_pfp_id_->reserve(PFP_SIZE);
    this->track_length_->reserve(PFP_SIZE);
    this->track_startx_->reserve(PFP_SIZE);
    this->track_starty_->reserve(PFP_SIZE);
    this->track_startz_->reserve(PFP_SIZE);
    this->track_start_distance_->reserve(PFP_SIZE);
    this->track_endx_->reserve(PFP_SIZE);
    this->track_endy_->reserve(PFP_SIZE);
    this->track_endz_->reserve(PFP_SIZE);
    this->track_dirx_->reserve(PFP_SIZE);
    this->track_diry_->reserve(PFP_SIZE);
    this->track_dirz_->reserve(PFP_SIZE);
    this->track_theta_->reserve(PFP_SIZE);
    this->track_phi_->reserve(PFP_SIZE);
    this->track_range_mom_mu_->reserve(PFP_SIZE);
    this->track_range_mom_p_->reserve(PFP_SIZE);
    this->track_range_mom_pi_->reserve(PFP_SIZE);
    this->track_mcs_bwd_mom_mu_->reserve(PFP_SIZE);
    this->track_mcs_bwd_mom_p_->reserve(PFP_SIZE);
    this->track_mcs_bwd_mom_pi_->reserve(PFP_SIZE);
    this->track_mcs_fwd_mom_mu_->reserve(PFP_SIZE);
    this->track_mcs_fwd_mom_p_->reserve(PFP_SIZE);
    this->track_mcs_fwd_mom_pi_->reserve(PFP_SIZE);
    this->track_chi2_u_proton_->reserve(PFP_SIZE);
    this->track_chi2_v_proton_->reserve(PFP_SIZE);
    this->track_chi2_y_proton_->reserve(PFP_SIZE);
    this->track_chi2_u_muon_->reserve(PFP_SIZE);
    this->track_chi2_v_muon_->reserve(PFP_SIZE);
    this->track_chi2_y_muon_->reserve(PFP_SIZE);
    this->track_chi2_u_pion_->reserve(PFP_SIZE);
    this->track_chi2_v_pion_->reserve(PFP_SIZE);
    this->track_chi2_y_pion_->reserve(PFP_SIZE);
    this->track_bestplane_->reserve(PFP_SIZE);

    //MC branches
    this->mc_nu_daughter_pdg_->reserve(EVENT_SIZE);
    
    //Slice branches
    // this->nu_pdg_->reserve(SLC_SIZE);
    this->topological_score_->reserve(SLC_SIZE);
    // this->nu_vx_->reserve(SLC_SIZE);
    // this->nu_vy_->reserve(SLC_SIZE);
    // this->nu_vz_->reserve(SLC_SIZE);
    // this->num_pf_particles_->reserve(SLC_SIZE);
  }

  // Helper function to set branch addresses for reading information
  // from the Event TTree
  void set_event_branch_addresses(TTree &etree) override
  {
    // Number of neutrino slices identified by the SliceID.
    SetBranchAddress(etree, "rec.nslc", &this->nslice_);

    // Is mc
    //SetBranchAddress(etree, "rec.hdr.ismc", &this->is_mc_);

    // Reco PDG code of primary PFParticle in slice (i.e., the neutrino
    // candidate)
    //SetBranchAddress(etree, "rec.slc.nu_pdg", &this->nu_pdg_); //error

    // Topological score
    SetBranchAddress(etree, "rec.slc.nu_score", this->topological_score_->data());

    // Reconstructed neutrino vertex position (with corrections for
    // space charge applied)
    // SetBranchAddress(etree, "rec.slc.vertex.x", &this->nu_vx_);
    // SetBranchAddress(etree, "rec.slc.vertex.y", &this->nu_vy_);
    // SetBranchAddress(etree, "rec.slc.vertex.z", &this->nu_vz_);

    // // Reconstructed object counts
    // SetBranchAddress(etree, "rec.slc.reco.npfp", &this->num_pf_particles_);
    
    // PFParticle properties
    //SetBranchAddress( etree, "rec.slc.reco.pfp.parent_is_primary", this->pfp_generation_->data() );
    //SetBranchAddress( etree, "rec.slc.reco.pfp.ndaughters", this->pfp_trk_daughters_count_->data() );
    //SetBranchAddress( etree, "rec.slc.reco.pfp.trackScore", this->pfp_track_score_->data() );
    //SetBranchAddress( etree, "rec.slc.reco.pfp.pdg", this->pfp_reco_pdg_->data() );
    //SetBranchAddress( etree, "rec.slc.reco.pfp.nhits", this->pfp_hits_->data() );
    
    //SetBranchAddress(etree, "rec.slc.reco.pfp.trk.calo.0.nhit", this->pfp_hitsU_->data()); // U plane
    //SetBranchAddress( etree, "rec.slc.reco.pfp.trk.calo.1.nhit", this->pfp_hitsV_->data() ); // V plane
    //SetBranchAddress( etree, "rec.slc.reco.pfp.trk.calo.2.nhit", this->pfp_hitsY_->data() ); // Y plane

    // Backtracked (MC truth) PFParticle properties
    // SetBranchAddress( etree, "rec.slc.reco.pfp.trk.truth.p.pdg", this->pfp_true_pdg_->data() );
    // SetBranchAddress( etree, "rec.slc.reco.pfp.trk.truth.p.genE", this->pfp_true_E_->data() );
    // SetBranchAddress( etree, "rec.slc.reco.pfp.trk.truth.p.genp.x", this->pfp_true_px_->data() );
    // SetBranchAddress( etree, "rec.slc.reco.pfp.trk.truth.p.genp.y", this->pfp_true_py_->data() );
    // SetBranchAddress( etree, "rec.slc.reco.pfp.trk.truth.p.genp.z", this->pfp_true_pz_->data() );

    // Showers (save for all PFParticles)
    // SetBranchAddress( etree, "rec.slc.reco.pfp.id", this->shower_pfp_id_->data() );
    // -start point
    SetBranchAddress( etree, "rec.slc.reco.pfp.shw.start.x", this->shower_startx_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.shw.start.y", this->shower_starty_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.shw.start.z", this->shower_startz_->data() );
    // -direction
    SetBranchAddress( etree, "rec.slc.reco.pfp.shw.dir.x", this->shower_dirx_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.shw.dir.y", this->shower_diry_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.shw.dir.z", this->shower_dirz_->data() );
    // -other properties (still important)
    SetBranchAddress( etree, "rec.slc.reco.pfp.shw.open_angle", this->shower_openangle_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.shw.len", this->shower_start_distance_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.shw.bestplane_energy", this->shower_energy_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.shw.bestplane_dEdx", this->shower_dEdx_->data() );

    // Tracks (save for all PFParticles)
    SetBranchAddress( etree, "rec.slc.reco.pfp.id", this->track_pfp_id_->data() );
    // -start point
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.start.x", this->track_startx_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.start.y", this->track_starty_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.start.z", this->track_startz_->data() );
    // -direction
    std::vector<float>* my_temp_ptr = this->track_dirx_.get_bare_ptr();
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.dir.x", &my_temp_ptr );
    //SetBranchAddress( etree, "rec.slc.reco.pfp.trk.dir.x", &(this->track_dirx_->data()) );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.dir.y", this->track_diry_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.dir.z", this->track_dirz_->data() );
    // -end point
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.end.x", this->track_endx_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.end.y", this->track_endy_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.end.z", this->track_endz_->data() );
    // -energy and momentum
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.rangeP.p_muon", this->track_range_mom_mu_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.rangeP.p_proton", this->track_range_mom_p_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.rangeP.p_pion", this->track_range_mom_pi_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.mcsP.bwdP_muon", this->track_mcs_bwd_mom_mu_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.mcsP.bwdP_proton", this->track_mcs_bwd_mom_p_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.mcsP.bwdP_pion", this->track_mcs_bwd_mom_pi_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.mcsP.fwdP_muon", this->track_mcs_fwd_mom_mu_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.mcsP.fwdP_proton", this->track_mcs_fwd_mom_p_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.mcsP.fwdP_pion", this->track_mcs_fwd_mom_pi_->data() );

    //-chi2 PIDs
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.chi2pid.0.chi2_proton", this->track_chi2_u_proton_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.chi2pid.1.chi2_proton", this->track_chi2_v_proton_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.chi2pid.2.chi2_proton", this->track_chi2_y_proton_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.chi2pid.0.chi2_muon", this->track_chi2_u_muon_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.chi2pid.1.chi2_muon", this->track_chi2_v_muon_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.chi2pid.2.chi2_muon", this->track_chi2_y_muon_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.chi2pid.0.chi2_pion", this->track_chi2_u_pion_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.chi2pid.1.chi2_pion", this->track_chi2_v_pion_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.chi2pid.2.chi2_pion", this->track_chi2_y_pion_->data() );

    // // -other properties (still important)
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.len", this->track_length_->data() );
    //SetBranchAddress( etree, "rec.slc.reco.pfp.trk.theta", this->track_theta_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.phi", this->track_phi_->data() );
    SetBranchAddress( etree, "rec.slc.reco.pfp.trk.bestplane", this->track_bestplane_->data() );

    // // MC truth information for the neutrino
    SetBranchAddress(etree, "rec.mc.nu.pdg", &this->mc_nu_pdg_);
    SetBranchAddress(etree, "rec.mc.nu.E", &this->mc_nu_energy_);
    SetBranchAddress(etree, "rec.mc.nu.isnc", &this->mc_nu_ccnc_); // 0=CC, 1=NC
    //SetBranchAddress(etree, "rec.mc.nu.genie_inttype", this->mc_nu_interaction_type_->data()); //error
    SetBranchAddress(etree, "rec.mc.nu.vtx.x", &this->mc_nu_vx_);
    SetBranchAddress(etree, "rec.mc.nu.vtx.y", &this->mc_nu_vy_);
    SetBranchAddress(etree, "rec.mc.nu.vtx.z", &this->mc_nu_vz_);
    // Error ends here

    // MC truth information for the final-state primary particles
    // SetBranchAddress( etree, "rec.mc.nu.prim.pdg", this->mc_nu_daughter_pdg_->data() );
    // SetBranchAddress( etree, "rec.mc.nu.prim.genE", this->mc_nu_daughter_energy_->data() );
    // // -Momenta
    // SetBranchAddress( etree, "rec.mc.nu.prim.genp.x", this->mc_nu_daughter_px_->data() );
    // SetBranchAddress( etree, "rec.mc.nu.prim.genp.y", this->mc_nu_daughter_py_->data() );
    // SetBranchAddress( etree, "rec.mc.nu.prim.genp.z", this->mc_nu_daughter_pz_->data() );
    // // -Start
    // SetBranchAddress( etree, "rec.mc.nu.prim.gen.x", this->mc_nu_daughter_startx_->data() );
    // SetBranchAddress( etree, "rec.mc.nu.prim.gen.y", this->mc_nu_daughter_starty_->data() );
    // SetBranchAddress( etree, "rec.mc.nu.prim.gen.z", this->mc_nu_daughter_startz_->data() );

    // GENIE and other systematic variation weights - todo
    bool has_genie_mc_weights = ( etree.GetBranch("globalTree") != nullptr );
    std::cout << "has_genie_mc_weights " << has_genie_mc_weights << std::endl;





  }

  // Helper function to set branch addresses for the output TTree
  void set_event_output_branch_addresses(TTree &out_tree, bool create = false) override
  {
    std::cout << "create " << create << std::endl;
    // Signal definition flags
    set_output_branch_address(out_tree, "is_mc", &this->is_mc_, create, "is_mc/O");

    // Event weights

    // Reco PDG code of primary PFParticle in slice (i.e., the neutrino
    // candidate)
    //set_output_branch_address(out_tree, "nu_pdg", &this->nu_pdg_, create, "nu_pdg/I"); //error

    // Number of neutrino slices identified by the SliceID. Allowed values
    // are zero or one.
    set_output_branch_address(out_tree, "nslice", &this->nslice_, create, "nslice/I");

    // Topological score
    set_output_branch_address(out_tree, "topological_score", &this->topological_score_, create, "topological_score/F");

    // Reconstructed neutrino vertex position (with corrections for
    // space charge applied)
    //set_output_branch_address(out_tree, "nu_vx", &this->nu_vx_, create, "nu_vx/F");
    //set_output_branch_address(out_tree, "nu_vy", &this->nu_vy_, create, "nu_vy/F");
    //set_output_branch_address(out_tree, "nu_vz", &this->nu_vz_, create, "nu_vz/F");

    // Reconstructed object counts
    //set_output_branch_address(out_tree, "num_pf_particles", &this->num_pf_particles_, create, "num_pf_particles/I");

    // PFParticle properties
    //set_output_branch_address(out_tree, "pfp_generation", this->pfp_generation_->data(), create, "pfp_generation/I");
    //set_output_branch_address(out_tree, "pfp_trk_daughters_count", this->pfp_trk_daughters_count_->data(), create, "pfp_trk_daughters_count/I");
    //set_output_branch_address(out_tree, "pfp_track_score", this->pfp_track_score_->data(), create, "pfp_track_score/F");
    //set_output_branch_address(out_tree, "pfp_hitsU", this->pfp_hitsU_->data(), create, "pfp_hitsU/I");
    //set_output_branch_address(out_tree, "pfp_hitsV", this->pfp_hitsV_->data(), create, "pfp_hitsV/I");
    //set_output_branch_address(out_tree, "pfp_hitsY", this->pfp_hitsY_->data(), create, "pfp_hitsY/I");
    
    // Backtracked (MC truth) PFParticle properties
    // set_output_branch_address(out_tree, "pfp_true_pdg", this->pfp_true_pdg_->data(), create, "pfp_true_pdg/I");
    // set_output_branch_address(out_tree, "pfp_true_E", this->pfp_true_E_->data(), create, "pfp_true_E/F");
    // set_output_branch_address(out_tree, "pfp_true_px", this->pfp_true_px_->data(), create, "pfp_true_px/F");
    // set_output_branch_address(out_tree, "pfp_true_py", this->pfp_true_py_->data(), create, "pfp_true_py/F");
    // set_output_branch_address(out_tree, "pfp_true_pz", this->pfp_true_pz_->data(), create, "pfp_true_pz/F");

    // Shower properties (save for all PFParticles)
    set_output_branch_address( out_tree, "shower_pfp_id", this->shower_pfp_id_->data(), create, "shower_pfp_id/I" );
    // -direction
    set_output_branch_address( out_tree, "shower_dirx", this->shower_dirx_->data(), create, "shower_dirx/F" );
    set_output_branch_address( out_tree, "shower_diry", this->shower_diry_->data(), create, "shower_diry/F" );
    set_output_branch_address( out_tree, "shower_dirz", this->shower_dirz_->data(), create, "shower_dirz/F" );
    // -other properties (still important)
    set_output_branch_address( out_tree, "shower_openangle", this->shower_openangle_->data(), create, "shower_openangle/F" );
    set_output_branch_address( out_tree, "shower_energy", this->shower_energy_->data(), create, "shower_energy/F" );
    set_output_branch_address( out_tree, "shower_dEdx", this->shower_dEdx_->data(), create, "shower_dEdx/F" );


    // Tracks (save for all PFParticles)
    set_output_branch_address( out_tree, "track_pfp_id", this->track_pfp_id_->data(), create, "track_pfp_id/I" );
    // -start point
    set_output_branch_address( out_tree, "track_startx", this->track_startx_->data(), create, "track_startx/F" );
    set_output_branch_address( out_tree, "track_starty", this->track_starty_->data(), create, "track_starty/F" );
    set_output_branch_address( out_tree, "track_startz", this->track_startz_->data(), create, "track_startz/F" );
    // -direction
    set_output_branch_address( out_tree, "track_dirx", this->track_dirx_->data(), create, "track_dirx/F" );
    set_output_branch_address( out_tree, "track_diry", this->track_diry_->data(), create, "track_diry/F" );
    set_output_branch_address( out_tree, "track_dirz", this->track_dirz_->data(), create, "track_dirz/F" );
    // -end point
    set_output_branch_address( out_tree, "track_endx", this->track_endx_->data(), create, "track_endx/F" );
    set_output_branch_address( out_tree, "track_endy", this->track_endy_->data(), create, "track_endy/F" );
    set_output_branch_address( out_tree, "track_endz", this->track_endz_->data(), create, "track_endz/F" );
    // -energy and momentum
    set_output_branch_address( out_tree, "track_range_mom_mu", this->track_range_mom_mu_->data(), create, "track_range_mom_mu/F" );
    set_output_branch_address( out_tree, "track_range_mom_p", this->track_range_mom_p_->data(), create, "track_range_mom_p/F" );
    set_output_branch_address( out_tree, "track_range_mom_pi", this->track_range_mom_pi_->data(), create, "track_range_mom_pi/F" );
    set_output_branch_address( out_tree, "track_mcs_bwd_mom_mu", this->track_mcs_bwd_mom_mu_->data(), create, "track_mcs_bwd_mom_mu/F" );
    set_output_branch_address( out_tree, "track_mcs_bwd_mom_p", this->track_mcs_bwd_mom_p_->data(), create, "track_mcs_bwd_mom_p/F" );
    set_output_branch_address( out_tree, "track_mcs_bwd_mom_pi", this->track_mcs_bwd_mom_pi_->data(), create, "track_mcs_bwd_mom_pi/F" );
    set_output_branch_address( out_tree, "track_mcs_fwd_mom_mu", this->track_mcs_fwd_mom_mu_->data(), create, "track_mcs_fwd_mom_mu/F" );
    set_output_branch_address( out_tree, "track_mcs_fwd_mom_p", this->track_mcs_fwd_mom_p_->data(), create, "track_mcs_fwd_mom_p/F" );
    set_output_branch_address( out_tree, "track_mcs_fwd_mom_pi", this->track_mcs_fwd_mom_pi_->data(), create, "track_mcs_fwd_mom_pi/F" );
    // -chi2 PIDs
    set_output_branch_address( out_tree, "track_chi2_u_proton", this->track_chi2_u_proton_->data(), create, "track_chi2_u_proton/F" );
    set_output_branch_address( out_tree, "track_chi2_v_proton", this->track_chi2_v_proton_->data(), create, "track_chi2_v_proton/F" );
    set_output_branch_address( out_tree, "track_chi2_y_proton", this->track_chi2_y_proton_->data(), create, "track_chi2_y_proton/F" );
    set_output_branch_address( out_tree, "track_chi2_u_muon", this->track_chi2_u_muon_->data(), create, "track_chi2_u_muon/F" );
    set_output_branch_address( out_tree, "track_chi2_v_muon", this->track_chi2_v_muon_->data(), create, "track_chi2_v_muon/F" );
    set_output_branch_address( out_tree, "track_chi2_y_muon", this->track_chi2_y_muon_->data(), create, "track_chi2_y_muon/F" );
    set_output_branch_address( out_tree, "track_chi2_u_pion", this->track_chi2_u_pion_->data(), create, "track_chi2_u_pion/F" );
    set_output_branch_address( out_tree, "track_chi2_v_pion", this->track_chi2_v_pion_->data(), create, "track_chi2_v_pion/F" );
    set_output_branch_address( out_tree, "track_chi2_y_pion", this->track_chi2_y_pion_->data(), create, "track_chi2_y_pion/F" );
    // -other properties (still important)
    set_output_branch_address( out_tree, "track_length", this->track_length_->data(), create, "track_length/F" );
    //set_output_branch_address( out_tree, "track_theta", this->track_theta_->data(), create, "track_theta/F" );
    set_output_branch_address( out_tree, "track_phi", this->track_phi_->data(), create, "track_phi/F" );
    set_output_branch_address( out_tree, "track_bestplane", this->track_bestplane_->data(), create, "track_bestplane/I" );


    // MC truth information for the final-state primary particles
    // set_output_branch_address( out_tree, "mc_nu_daughter_pdg", this->mc_nu_daughter_pdg_->data(), create, "mc_nu_daughter_pdg/I" );
    // set_output_branch_address( out_tree, "mc_nu_daughter_energy", this->mc_nu_daughter_energy_->data(), create, "mc_nu_daughter_energy/F" );
    // set_output_branch_address( out_tree, "mc_nu_daughter_px", this->mc_nu_daughter_px_->data(), create, "mc_nu_daughter_px/F" );
    // set_output_branch_address( out_tree, "mc_nu_daughter_py", this->mc_nu_daughter_py_->data(), create, "mc_nu_daughter_py/F" );
    // set_output_branch_address( out_tree, "mc_nu_daughter_pz", this->mc_nu_daughter_pz_->data(), create, "mc_nu_daughter_pz/F" );
    // set_output_branch_address( out_tree, "mc_nu_daughter_startx", this->mc_nu_daughter_startx_->data(), create, "mc_nu_daughter_startx/F" );
    // set_output_branch_address( out_tree, "mc_nu_daughter_starty", this->mc_nu_daughter_starty_->data(), create, "mc_nu_daughter_starty/F" );
    // set_output_branch_address( out_tree, "mc_nu_daughter_startz", this->mc_nu_daughter_startz_->data(), create, "mc_nu_daughter_startz/F" );
  }
};


