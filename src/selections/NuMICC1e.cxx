// XSecAnalyzer includes
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/TreeUtils.hh"

#include "XSecAnalyzer/Selections/EventCategoriesNuMICC1e.hh"
#include "XSecAnalyzer/Selections/NuMICC1e.hh"

NuMICC1e::NuMICC1e() : SelectionBase( "NuMICC1e" ) {}

void NuMICC1e::define_constants() { 
    // FV definitions as in PeLEE analysis
    // x_min, x_max, y_min, y_max, z_min, z_max
    this->define_true_FV( 10., 246., -105., 105., 10., 1026. );
    this->define_reco_FV( 10., 246., -105., 105., 10., 1026. );
}

void NuMICC1e::compute_reco_observables( AnalysisEvent* Event ) {
    // Evaluate the reconstructed kinematic variables of interest for the xsec
    // measurement
  
    // Set the reco energy of the electron candidate if we found one
    if ( !sel_pass_shower_identification_ ) return;

    // set reco energy, applying shower energy correction factor
    reco_electron_energy_ = Event->shr_energy_cali_ / 0.83;
}

void NuMICC1e::compute_true_observables( AnalysisEvent* Event ) {
    // Evaluate the true kinematic variables of interest
  
    // Check if there is a true final-state electron in this event
    bool has_true_electron = ( sig_isNuE_ && Event->mc_nu_ccnc_ == CHARGED_CURRENT );
  
    // If there isn't one, we don't need to do anything else
    if ( !has_true_electron ) return;

    // set true electron energy
    mc_electron_energy_ = Event->mc_elec_e_;
  
}

int NuMICC1e::categorize_event( AnalysisEvent* Event ) {
    // Identify the event category of the selected event
  
    // Real data has a bogus true neutrino PDG code that is not one of the
    // allowed values (±12, ±14, ±16)
    // is this EXT only? can change name
    int abs_mc_nu_pdg = std::abs( Event->mc_nu_pdg_ );
    Event->is_mc_ = ( abs_mc_nu_pdg == ELECTRON_NEUTRINO ||
      abs_mc_nu_pdg == MUON_NEUTRINO || abs_mc_nu_pdg == TAU_NEUTRINO );
    if ( !Event->is_mc_ ) {
      return kEXT;
    }
  
    // All events outside of the true fiducial volume should be categorized
    // as "out of fiducial volume"
    bool mcVertexInFV = point_inside_FV( this->true_FV(),
      Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_ );
    if ( !mcVertexInFV ) {
      return kOOFV;
    }
  
    // neutral current
    if ( Event->mc_nu_ccnc_ == NEUTRAL_CURRENT ) {
        if(Event->mc_npi0_ > 0) return kNCPi0;
        else return kNCOther;
    }

    // CC muon (anti)neutrinos
    if ( std::abs(Event->mc_nu_pdg_) == MUON_NEUTRINO ) {
        if(Event->mc_npi0_ > 0) return kNuMuCCPi0;
        else return kNuMuCCOther;
    }
  
    // CC electron (anti)neutrinos
    if (std::abs(Event->mc_nu_pdg_) == ELECTRON_NEUTRINO) {
        // signal events
        if ( this->is_event_mc_signal() ) return kNuECC;
        // non-signal nues
        else return kNuECCOther;
    }
  
    // We shouldn't ever get here, but return "unknown" just in case
    std::cout << "Warning: Unknown event! Check the categorization logic." << std::endl;
    return kUnknown;
}
  
bool NuMICC1e::define_signal( AnalysisEvent* Event ) {
    // Determine whether or not the input event matches the signal definition
    // for this selection. This determination should be done only with MC truth
    // information.
  
    // Require signal events to be inside the true fiducial volume
    sig_inFV_ = point_inside_FV( this->true_FV(), Event->mc_nu_vx_,
      Event->mc_nu_vy_, Event->mc_nu_vz_ );
  
    // Require an incident electron (anti)neutrino above threshold
    sig_isNuE_ = ( std::abs(Event->mc_nu_pdg_) == ELECTRON_NEUTRINO );

    // Require charged-current interaction
    sig_isCC_ = Event->mc_nu_ccnc_ == CHARGED_CURRENT;
  
    // Require a final-state electron above threshold
    if (Event->mc_nelec_ > 0) sig_has_fs_electron_ = true; // 30 MeV threshold
  
    sig_is_signal_ = sig_inFV_ && sig_isNuE_ && sig_isCC_ && sig_has_fs_electron_;
    
    return sig_is_signal_;
}

bool NuMICC1e::selection(AnalysisEvent* Event) {
    // Apply the selection cuts on reco variables only to determine whether this
    // event is a candidate signal event or not.

    // reset variables
    sel_pass_preselection_ = false;
    sel_pass_cosmic_rejection_ = false;
    sel_pass_shower_identification_ = false;
    sel_pass_electron_identification_ = false;
    sel_nu_e_cc_ = false;

    // pre-selection
    // neutrino slice
    if (Event->nslice_ != 1) return false;
    // vertex inside FV
    if (!point_inside_FV( this->reco_FV(), Event->nu_vx_, Event->nu_vy_, Event->nu_vz_ )) return false;
    // at least one shower
    if (Event->num_showers_ < 1) return false;
    // contained fraction
    if (Event->contained_fraction_ < 0.85) return false;

    sel_pass_preselection_ = true;

    // cosmic rejection
    // topological score
    if (Event->topological_score_ < 0.2) return false;
    // cosmic impact parameter
    if (Event->cosmic_impact_parameter_ < 10) return false;

    sel_pass_cosmic_rejection_ = true;

    // shower identification
    // valid shower ID
    if (Event->shr_id_ == 0) return false; // default value, not filled
    // shower pfp generation
    if (Event->pfp_generation_->at(Event->shr_id_-1) != 2) return false;
    // shower energy (to match signal definition)
    if ( (Event->shr_energy_cali_ / 0.83) < 0.03) return false;
    // shower score
    if (Event->shr_score_ > 0.15) return false;
    // shower hits ratio
    if (Event->hits_ratio_ < 0.5) return false;
  
    sel_pass_shower_identification_ = true;

    // electron identification
    // moliere average angle
    if (Event->shrmoliereavg_ > 7) return false;
    // shower distance and dE/dx
    if (Event->num_tracks_ > 0) {
        // track present, 2D distance-dE/dx cut
        if (Event->shr_tkfit_gap10_dedx_Y_ >= 0 && Event->shr_tkfit_gap10_dedx_Y_ < 1.75) {
            if (Event->shr_distance_ > 3) return false;
        }
        else if (Event->shr_tkfit_gap10_dedx_Y_ >= 1.75 && Event->shr_tkfit_gap10_dedx_Y_ < 2.5) {
            if (Event->shr_distance_ > 12) return false;
        }
        else if (Event->shr_tkfit_gap10_dedx_Y_ >= 2.5 && Event->shr_tkfit_gap10_dedx_Y_ < 3.5) {
            if (Event->shr_distance_ > 3) return false;
        }
        else if (Event->shr_tkfit_gap10_dedx_Y_ >= 3.5 && Event->shr_tkfit_gap10_dedx_Y_ < 4.7) {
            return false;
        }
        else if (Event->shr_tkfit_gap10_dedx_Y_ >= 4.7) {
            if (Event->shr_distance_ > 3) return false;
        }
        else return false;
    }
    else {
        // no track, 1D dE/dx cut
        if (Event->shr_tkfit_gap10_dedx_Y_ < 1.7) return false;
        if (Event->shr_tkfit_gap10_dedx_Y_ > 2.7 && Event->shr_tkfit_gap10_dedx_Y_ < 5.5) return false;
    }
   
    // if gets to here, it has passed the selection
    sel_nu_e_cc_ = true;
  
    return sel_nu_e_cc_;
  }


  void NuMICC1e::define_output_branches() {
    // Save any additional variables to the output TTree
    this->set_branch( &sig_isNuE_, "mc_is_nue" );
    this->set_branch( &sig_isCC_, "mc_is_CC" );
    this->set_branch( &sig_inFV_, "mc_vertex_in_FV" );
    this->set_branch( &sig_has_fs_electron_, "mc_has_fs_electron" );
    this->set_branch( &sig_is_signal_, "mc_is_signal" );

    this->set_branch( &sel_pass_preselection_, "sel_pass_preselection");
    this->set_branch( &sel_pass_cosmic_rejection_, "sel_pass_cosmic_rejection");
    this->set_branch( &sel_pass_shower_identification_, "sel_pass_shower_identification");
    this->set_branch( &sel_nu_e_cc_, "sel_nu_e_cc");
    
    this->set_branch( &reco_electron_energy_, "reco_electron_energy" );
    this->set_branch( &mc_electron_energy_, "mc_electron_energy" );
  }
  
  void NuMICC1e::define_category_map() {
    // Use NuMI CC 1e event categories
    categ_map_ = NuMICC1e_MAP;
  }

  void NuMICC1e::reset() {

    sig_isNuE_ = false;
    sig_isCC_ = false;
    sig_inFV_ = false;
    sig_has_fs_electron_ = false;
    sig_is_signal_ = false;
  
    sel_pass_preselection_ = false;
    sel_pass_cosmic_rejection_ = false;
    sel_pass_shower_identification_ = false;
    sel_pass_electron_identification_ = false;
    sel_nu_e_cc_ = false;

    reco_electron_energy_ = BOGUS;
    mc_electron_energy_ = BOGUS;
}