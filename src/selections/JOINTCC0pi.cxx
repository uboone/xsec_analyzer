// XSecAnalyzer includes
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/TreeUtils.hh"

#include "XSecAnalyzer/Selections/JOINTCC0pi.hh"
#include "XSecAnalyzer/Selections/EventCategoriesJOINTCC0pi.hh"

// XGBoost
//#include <xgboost/c_api.h>

// XGBoost model
//BoosterHandle* booster;

JOINTCC0pi::JOINTCC0pi() : SelectionBase("JOINTCC0pi") {
  //CalcType = kOpt1;
  booster = new BoosterHandle;
  XGBoosterCreate(NULL, 0, booster);
  assert(XGBoosterLoadModel(*booster, "/exp/uboone/app/users/cnguyen/stv-analysis-II/xsec_analyzer/lib/cc0pi_fstrack_pid_softprob.json") == 0);
  std::cout<<"Got XGBoosterCreate  "<<std::endl;
  std::cout<<"Finished with Constuctor  "<<std::endl;
  
  xgb_pid_vec_.reset( new std::vector<int>() );
  
}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void JOINTCC0pi::DefineConstants() {
  DefineTrueFV(10.,246.,-105.,105.,10.,1026.);
  DefineRecoFV(10.,246.,-105.,105.,10.,1026.);
  std::cout<<"Finished ::DefineConstants "<< std::endl;
}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void JOINTCC0pi::ComputeRecoObservables(AnalysisEvent* Event){

 if(verbosal > 0) std::cout<< "FInished truth observables"<< std::endl;
  // In cases where we failed to find a muon candidate, check whether there are
  // at least two generation == 2 PFParticles. If there are, then compute the
  // usual observables using the longest track as the muon candidate and the
  // second-longest track as the leading proton candidate. This will enable
  // sideband studies of NC backgrounds in the STV phase space.
  if ( !sel_has_muon_candidate_ ) {

    float max_trk_len = LOW_FLOAT;
    int max_trk_idx = BOGUS_INDEX;

    float next_to_max_trk_len = LOW_FLOAT;
    int next_to_max_trk_idx = BOGUS_INDEX;
    
    if(verbosal > 0) std::cout<< "Stating  num_pf_particles_ loop "<< std::endl;

    for ( int p = 0; p < num_pf_particles_; ++p ) {

      // Only include direct neutrino daughters (generation == 2)
      unsigned int generation = Event->pfp_generation_->at( p );
      if ( generation != 2u ) continue;

      float trk_len = Event->track_length_->at( p );

      if ( trk_len > next_to_max_trk_len ) {

        next_to_max_trk_len = trk_len;
        next_to_max_trk_idx = p;

        if ( next_to_max_trk_len > max_trk_len ) {

          next_to_max_trk_len = max_trk_len;
          next_to_max_trk_idx = max_trk_idx;

          max_trk_len = trk_len;
          max_trk_idx = p;
        }
      }
    }

    // If we found at least two usable PFParticles, then assign the indices to
    // be used below
    if ( max_trk_idx != BOGUS_INDEX && next_to_max_trk_idx != BOGUS_INDEX ) {
      muon_candidate_idx_ = max_trk_idx;
      lead_p_candidate_idx_ = next_to_max_trk_idx;
    }
  }

  // Abbreviate some of the calculations below by using these handy
  // references to the muon and leading proton 3-momenta
      if(verbosal > 0)std::cout<< "Seting p3mu   "<< std::endl;
  auto& p3mu = *p3_mu_;
      if(verbosal > 0) std::cout<< "Seting p3mu_rangle   "<< std::endl;
  auto& p3mu_range = *p3_mu_range_;
    if(verbosal > 0)  std::cout<< "Seting p3mu_mcs   "<< std::endl;
  auto& p3mu_mcs = *p3_mu_mcs_;
  
  auto& p3p = *p3_lead_p_;

  // Set the reco 3-momentum of the muon candidate if we found one
  bool muon = muon_candidate_idx_ != BOGUS_INDEX;
  if ( muon ) {
    float mu_dirx = Event->track_dirx_->at( muon_candidate_idx_ );
    float mu_diry = Event->track_diry_->at( muon_candidate_idx_ );
    float mu_dirz = Event->track_dirz_->at( muon_candidate_idx_ );

    p3mu = TVector3( mu_dirx, mu_diry, mu_dirz );
    p3mu_mcs = TVector3( mu_dirx, mu_diry, mu_dirz );
    p3mu_range= TVector3( mu_dirx, mu_diry, mu_dirz );
    // The selection flag indicating whether the muon candidate is contained
    // was already set when the selection was applied. Use it to choose the
    // best momentum estimator to use.
    float muon_mom = LOW_FLOAT;
    float muon_muon_range = Event->track_range_mom_mu_->at( muon_candidate_idx_ );
    float muon_muon_mcs = Event->track_mcs_mom_mu_->at( muon_candidate_idx_ );
    /// Set Muon Tracklen
    muontrklen_ = Event->track_length_->at( muon_candidate_idx_ );
    muondistancetovectex_ =  Event->track_start_distance_->at( muon_candidate_idx_ );
    muon_llr_pid_score_ =  Event->track_llr_pid_score_->at( muon_candidate_idx_ );
    
    muon_tkscore_ = Event->pfp_track_score_->at( muon_candidate_idx_ );
    
    muon_trkstart_x_ = Event->track_startx_->at(muon_candidate_idx_); 
    muon_trkstart_y_ = Event->track_starty_->at(muon_candidate_idx_); 
    muon_trkstart_z_ = Event->track_startz_->at(muon_candidate_idx_); 
        
    muon_trkend_x_ = Event->track_endx_->at(muon_candidate_idx_);
    muon_trkend_y_ = Event->track_endx_->at(muon_candidate_idx_);
    muon_trkend_z_ = Event->track_endx_->at(muon_candidate_idx_);
    
    muon_trkchi2muon_ = Event->track_chi2_muon_->at(muon_candidate_idx_);
    
    
     if(verbosal > 0)std::cout<< "Finished getting float muon_mom  float muon_muon_range loat muon_muon_mcs   "<< std::endl;
    
    
    
   ////No Corrections////// 
   /*if ( sel_muon_contained_ ) {
       muon_mom = track_range_mom_mu_->at( muon_candidate_idx_ );
    }
    else{
       muon_mom = track_mcs_mom_mu_->at( muon_candidate_idx_ );
    }*/

    if ( sel_muon_contained_ ) {
    // Unsure how to if I should apply it since not all root files have these branches 
    if (Event->track_length_->size() == 1 && Event->trk_bragg_mu_fwd_preferred_v_->at(muon_candidate_idx_) == 0 && 6. < Event->track_chi2_muon_->at(muon_candidate_idx_)){
            std::cout<<"Rejecting flipped tracks"<<std::endl; 
             muon_mom = LOW_FLOAT;
            }
    else{
       muon_mom = Event->track_range_mom_mu_->at( muon_candidate_idx_ );
       }
     }
     
    else{
    if(p3mu.CosTheta() > -0.9){
    
    
          if (Event->track_mcs_mom_mu_->at( muon_candidate_idx_ ) < 0.11){
              muon_mom = LOW_FLOAT;
          }
          else{
              muon_mom = Event->track_mcs_mom_mu_->at( muon_candidate_idx_ ) - 0.0361*Event->track_mcs_mom_mu_->at( muon_candidate_idx_ ) + 0.04;
          }
    }
    else{ 
    muon_mom = LOW_FLOAT;
    }
          
    }
    
    if(verbosal > 0)std::cout<<"  muon_mom = "<<  muon_mom << std::endl;
    


    p3mu = p3mu.Unit() * muon_mom;    
    p3mu_mcs = p3mu.Unit() * muon_muon_mcs;
    p3mu_range = p3mu.Unit() * muon_muon_range;

    if(verbosal > 0) std::cout<<" stating distance_FV_surface_" << std::endl;
    distance_FV_surface_ = this->reco_distance_to_FV_Surface(Event);
    
    if(verbosal > 0)std::cout<<"  Finsished ::distance_FV_surface_ =" << distance_FV_surface_<<std::endl;


  }
  // Set the reco 3-momentum of the leading proton candidate if we found one
  bool lead_p = lead_p_candidate_idx_ != BOGUS_INDEX;
  if ( lead_p ) {

    float p_dirx = Event->track_dirx_->at( lead_p_candidate_idx_ );
    float p_diry = Event->track_diry_->at( lead_p_candidate_idx_ );
    float p_dirz = Event->track_dirz_->at( lead_p_candidate_idx_ );
    float KEp = Event->track_kinetic_energy_p_->at( lead_p_candidate_idx_ );
    float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

    p3p = TVector3( p_dirx, p_diry, p_dirz );
    p3p = p3p.Unit() * p_mom;
  }

  // Reset the vector of reconstructed proton candidate 3-momenta
  p3_p_vec_->clear();

  // Set the reco 3-momenta of all proton candidates (i.e., all generation == 2
  // tracks except the muon candidate) assuming we found both a muon candidate
  // and at least one proton candidate.
  if ( muon && lead_p ) {
    for ( int p = 0; p < Event->num_pf_particles_; ++p ) {
      // Skip the muon candidate
      if ( p == muon_candidate_idx_ ) continue;

      // Only include direct neutrino daughters (generation == 2)
      unsigned int generation = Event->pfp_generation_->at( p );
      if ( generation != 2u ) continue;

      float p_dirx = Event->track_dirx_->at( p );
      float p_diry = Event->track_diry_->at( p );
      float p_dirz = Event->track_dirz_->at( p );
      float KEp = Event->track_kinetic_energy_p_->at( p );
      float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

      TVector3 p3_temp( p_dirx, p_diry, p_dirz );
      p3_temp = p3_temp.Unit() * p_mom;

      p3_p_vec_->push_back( p3_temp );
    }

    // TODO: reduce code duplication by just getting the leading proton
    // 3-momentum from this sorted vector
    // Sort the reco proton 3-momenta in order from highest to lowest magnitude
    std::sort( p3_p_vec_->begin(), p3_p_vec_->end(), [](const TVector3& a,
      const TVector3& b) -> bool { return a.Mag() > b.Mag(); } );
  }

  // Compute reco STVs if we have both a muon candidate
  // and a leading proton candidate in the event
  if ( muon && lead_p ) {
    this->compute_stvs( p3mu, p3p, delta_pT_, delta_phiT_,
      delta_alphaT_, delta_pL_, pn_, delta_pTx_,delta_pTy_ );

    theta_mu_p_ = std::acos( p3mu.Dot(p3p) / p3mu.Mag() / p3p.Mag() );
  }


}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void JOINTCC0pi::ComputeTrueObservables(AnalysisEvent* Event){

 // If this is not an MC event, then just return without doing anything
  if ( !Event->is_mc_ ) return;

  size_t num_mc_daughters = Event->mc_nu_daughter_pdg_->size();

  // Set the true 3-momentum of the final-state muon if there is one
  bool true_muon = ( mc_neutrino_is_numu_ && Event->mc_nu_ccnc_ == CHARGED_CURRENT );
  if ( true_muon ) {
    // Loop over the MC neutrino daughters, find the muon, and get its
    // true 3-momentum. Note that we assume there is only one muon in
    // this loop.
    bool found_muon = false;
    for ( size_t d = 0u; d < num_mc_daughters; ++d ) {
      int pdg = Event->mc_nu_daughter_pdg_->at( d );
      if ( pdg == MUON ) {
        found_muon = true;
        float px = Event->mc_nu_daughter_px_->at( d );
        float py = Event->mc_nu_daughter_py_->at( d );
        float pz = Event->mc_nu_daughter_pz_->at( d );
        *mc_p3_mu_ = TVector3( px, py, pz );
        mc_muontrklen_ = distanceBetweentwopoint(Event->mc_nu_vx_,
        Event->mc_nu_vy_,
        Event->mc_nu_vz_,
        Event->mc_nu_daughter_endx_->at( d ),
        Event->mc_nu_daughter_endy_->at( d ),
        Event->mc_nu_daughter_endz_->at( d ));
        
        break;
      }
    }

    if ( !found_muon ) {
      std::cout << "WARNING: Missing muon in MC signal event!\n";
      return;
    }
  }

  // Reset the vector of true MC proton 3-momenta
  mc_p3_p_vec_->clear();

  // Set the true 3-momentum of the leading proton (if there is one)
  float max_mom = LOW_FLOAT;
  for ( int p = 0; p < num_mc_daughters; ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    if ( pdg == PROTON )
    {
      float px = Event->mc_nu_daughter_px_->at( p );
      float py = Event->mc_nu_daughter_py_->at( p );
      float pz = Event->mc_nu_daughter_pz_->at( p );
      TVector3 temp_p3 = TVector3( px, py, pz );

      mc_p3_p_vec_->push_back( temp_p3 );

      float mom = temp_p3.Mag();
      if ( mom > max_mom ) {
        max_mom = mom;
        *mc_p3_lead_p_ = temp_p3;
      }
    }
  }

  // TODO: reduce code duplication by just getting the leading proton
  // 3-momentum from this sorted vector
  // Sort the true proton 3-momenta in order from highest to lowest magnitude
  std::sort( mc_p3_p_vec_->begin(), mc_p3_p_vec_->end(), [](const TVector3& a,
    const TVector3& b) -> bool { return a.Mag() > b.Mag(); } );

  mc_distance_FV_surface_ = this->mc_distance_to_FV_Surface(Event);


  // If the event contains a leading proton, then set the 3-momentum
  // accordingly
  bool true_lead_p = ( max_mom != LOW_FLOAT );
  if ( !true_lead_p && mc_is_signal_ ) {
    // If it doesn't for a signal event, then something is wrong.
    std::cout << "WARNING: Missing leading proton in MC signal event!\n";
    return;
  }

  // Compute true STVs if the event contains both a muon and a leading
  // proton
  if ( true_muon && true_lead_p ) {
    this->compute_stvs( *mc_p3_mu_, *mc_p3_lead_p_, mc_delta_pT_, mc_delta_phiT_,
      mc_delta_alphaT_, mc_delta_pL_, mc_pn_, mc_delta_pTx_, mc_delta_pTy_ );

    mc_theta_mu_p_ = std::acos( mc_p3_mu_->Dot(*mc_p3_lead_p_)
      / mc_p3_mu_->Mag() / mc_p3_lead_p_->Mag() );
  }

 return; 

}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
int JOINTCC0pi::CategorizeEvent(AnalysisEvent* Event){
 
  // Real data has a bogus true neutrino PDG code that is not one of the
  // allowed values (±12, ±14, ±16)
  int abs_mc_nu_pdg = std::abs( Event->mc_nu_pdg_ );
  Event->is_mc_ = ( abs_mc_nu_pdg == ELECTRON_NEUTRINO
    || abs_mc_nu_pdg == MUON_NEUTRINO || abs_mc_nu_pdg == TAU_NEUTRINO );
  if ( !Event->is_mc_ ) return kUnknown;
  
  mc_neutrino_is_numu_ = ( Event->mc_nu_pdg_ == MUON_NEUTRINO );

  if ( !mc_vertex_in_FV_ ) {
    mc_is_signal_ = false;
    mc_is_cc0pi_signal_ = false;
    mc_is_cc0pi_wc_signal_ = false;
    return kOOFV;
  }
  else if ( Event->mc_nu_ccnc_ == NEUTRAL_CURRENT ) {
    mc_is_signal_ = false;
    mc_is_cc0pi_signal_ = false;
    mc_is_cc0pi_wc_signal_ = false;
    return kNC;
  }
  else if ( !mc_neutrino_is_numu_ ) {
    mc_is_signal_ = false;
    mc_is_cc0pi_signal_ = false;
    mc_is_cc0pi_wc_signal_ = false;
    if ( Event->mc_nu_pdg_ == ELECTRON_NEUTRINO
      && Event->mc_nu_ccnc_ == CHARGED_CURRENT ) return kNuECC;
    else return kOther;
  }


  // Sort signal by interaction mode
  if ( mc_is_cc0pi_signal_ ) {
    if ( Event->mc_nu_interaction_type_ == 0 ) return kSignalCCQE; // QE
    else if ( Event->mc_nu_interaction_type_ == 10 ) return kSignalCCMEC; // MEC
    else if ( Event->mc_nu_interaction_type_ == 1 ) return kSignalCCRES; // RES
    //else if ( Event->mc_nu_interaction_type_ == 2 ) // DIS
    //else if ( Event->mc_nu_interaction_type_ == 3 ) // COH
    else return kSignalOther;
  }
  else if ( !mc_no_fs_pi0_ || !mc_no_charged_pi_above_threshold_ ) {
    return kNuMuCCNpi;
  }
  /*else if ( !mc_lead_p_in_mom_range_ ) {
    return kNuMuCC0pi0p;
  }*/
  else return kNuMuCCOther;

}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
bool JOINTCC0pi::DefineSignal(AnalysisEvent* Event){

  // Set flags to their default values here
  mc_muon_in_mom_range_ = false;
  mc_muon_in_wc_mom_range_ = false;
  mc_lead_p_in_mom_range_ = false;
  mc_no_fs_pi0_ = true;
  mc_no_charged_pi_above_threshold_ = true;
  mc_no_charged_pi_above_wc_threshold_ = true;
  mc_no_fs_mesons_ = true;
  mc_vertex_in_FV_ = this->mc_vertex_inside_FV(Event);
  
  double lead_p_mom = LOW_FLOAT;

  for ( size_t p = 0u; p < Event->mc_nu_daughter_pdg_->size(); ++p ) {
    int pdg = Event->mc_nu_daughter_pdg_->at( p );
    float energy = Event->mc_nu_daughter_energy_->at( p );

    // Do the general check for (anti)mesons first before considering
    // any individual PDG codes
    if ( is_meson_or_antimeson(pdg) ) {
      mc_no_fs_mesons_ = false;
    }

    // Check that the muon has a momentum within the allowed range
    if ( pdg == MUON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(MUON_MASS, 2) );
      if ( mom >= MUON_P_MIN_MOM_CUT_jointcc0pi && mom <= MUON_P_MAX_MOM_CUT_jointcc0pi ) {
        mc_muon_in_mom_range_ = true;
        
      }
      if ( mom >= MUON_P_MIN_WC_MOM_CUT_jointcc0pi && mom <= MUON_P_MAX_MOM_CUT_jointcc0pi ) {
        mc_muon_in_wc_mom_range_ = true;
      
      }
    }
    else if ( pdg == PROTON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PROTON_MASS, 2) );
      if ( mom > lead_p_mom ) lead_p_mom = mom;
      mc_num_protons_++;
    }
    else if ( pdg == NEUTRON ) {
      mc_num_neutrons_++;
    }
    else if ( pdg == PI_ZERO ) {
      mc_no_fs_pi0_ = false;
    }
    else if ( std::abs(pdg) == PI_PLUS ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PI_PLUS_MASS, 2) );
      if ( mom > CHARGED_PI_MOM_CUT_jointcc0pi ) {
        mc_no_charged_pi_above_threshold_ = false;
        mc_num_charged_pions_++;
      }
      if ( mom > CHARGED_PI_WC_MOM_CUT_jointcc0pi ) {
        mc_no_charged_pi_above_wc_threshold_ = false;
        mc_num_wc_charged_pions_++;
      }
    }
  }

  // Check that the leading proton has a momentum within the allowed range
  if ( lead_p_mom >= LEAD_P_MIN_MOM_CUT && lead_p_mom <= LEAD_P_MAX_MOM_CUT ) {
    mc_lead_p_in_mom_range_ = true;
  }

  mc_is_signal_ = mc_vertex_in_FV_ && mc_neutrino_is_numu_
    && mc_muon_in_mom_range_ && mc_lead_p_in_mom_range_
    && mc_no_fs_mesons_;

  mc_is_cc0pi_signal_ = mc_vertex_in_FV_ && mc_neutrino_is_numu_
    && mc_muon_in_mom_range_ && mc_no_fs_pi0_
    && mc_no_charged_pi_above_threshold_;

  mc_is_cc0pi_wc_signal_ = mc_vertex_in_FV_ && mc_neutrino_is_numu_
    && mc_muon_in_wc_mom_range_ && mc_no_fs_pi0_
    && mc_no_charged_pi_above_wc_threshold_;

  return mc_is_cc0pi_signal_;
}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
bool JOINTCC0pi::Selection(AnalysisEvent* Event){
  
  int reco_shower_count = 0;
   if(verbosal > 0) std::cout<< "Starting loop "<<std::endl;
  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    float tscore = Event->pfp_track_score_->at( p );
    if ( tscore <= TRACK_SCORE_CUT ) ++reco_shower_count;
  }
  
 if(verbosal > 0) std::cout<< "Finished loop "<<std::endl;
  // Check the shower cut
  sel_no_reco_showers_ = ( reco_shower_count == 0 );
  
  /////////////////////////////
  // Apply Numu CC Selection
  ///////////////////////////
 this->apply_numu_CC_selection(Event);

 if(verbosal > 0) std::cout<< "Finished apply_numu_CC_selection"<<std::endl;
 
  // Set flags that default to true here
  sel_passed_proton_pid_cut_ = true;
  sel_protons_contained_ = true;

  // Set flags that default to false here
  sel_muon_passed_mom_cuts_ = false;
  sel_muon_passed_wc_mom_cuts_ = false;
  sel_muon_contained_ = false;
  sel_has_pion_candidate_ = false;
  sel_has_pion_wc_candidate_ = false;

  sel_num_proton_candidates_ = 0;
  sel_num_pion_candidates_ = 0;
  sel_num_pion_wc_candidates_ = 0;

  if(verbosal > 0)std::cout<< "Starting loop "<<std::endl;

  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {
      if(verbosal > 0) std::cout<<"p  = "<< p<< std::endl;  

     // Only worry about direct neutrino daughters (PFParticles considered
     // daughters of the reconstructed neutrino)
     unsigned int generation = Event->pfp_generation_->at( p );
    
     if(verbosal > 0)std::cout<<"generation  = "<< generation<< std::endl;  
     if ( generation != 2u ) continue;

     // Check that we can find a muon candidate in the event. If more than
     // one is found, also fail the cut.
    if ( p == muon_candidate_idx_ ) {
     if(verbosal > 0) std::cout<<" muon_candidate_idx_  "<< std::endl; 
      // Check whether the muon candidate is contained. Use the same
      // containment volume as the protons. TODO: revisit this as needed.
      float endx = Event->track_endx_->at( p );
      float endy = Event->track_endy_->at( p );
      float endz = Event->track_endz_->at( p );
      bool end_contained = this->in_proton_containment_vol( endx, endy, endz );
      if ( end_contained ) sel_muon_contained_ = true;
      
      auto& p3mu = *p3_mu_; 
      float mu_dirx = Event->track_dirx_->at( muon_candidate_idx_ );
      float mu_diry = Event->track_diry_->at( muon_candidate_idx_ );
      float mu_dirz = Event->track_dirz_->at( muon_candidate_idx_ );
       p3mu = TVector3( mu_dirx, mu_diry, mu_dirz );
       double angle = p3mu.CosTheta();
      float muon_mom = LOW_FLOAT;
      float range_muon_mom = Event->track_range_mom_mu_->at( p );
      float mcs_muon_mom = Event->track_mcs_mom_mu_->at( p );
    
      ///////No Corrections////////
      /*if ( sel_muon_contained_ ){
         muon_mom = range_muon_mom;  
      }
      else {
         muon_mom = mcs_muon_mom; 
      }*/
      
      if ( sel_muon_contained_ ){
            
           if (Event->track_length_->size() == 1 && Event->trk_bragg_mu_fwd_preferred_v_->at(muon_candidate_idx_) == 0 && 6. < Event->track_chi2_muon_->at(muon_candidate_idx_)){
            std::cout<<"Rejecting flipped tracks"<<std::endl; 
            muon_mom = LOW_FLOAT;
           }
            
     else{  
          muon_mom = range_muon_mom;
          } 
      } 
      else{
          if(p3mu.CosTheta() > -0.9){
      
                if (mcs_muon_mom < 0.11){
                    muon_mom = LOW_FLOAT;
                }
                else{
                    muon_mom = mcs_muon_mom - 0.0361*mcs_muon_mom + 0.04;
                    }   
          }
          else{
                muon_mom = LOW_FLOAT;
              }
        }
      
      if ( muon_mom >= MUON_P_MIN_MOM_CUT_jointcc0pi && muon_mom <= MUON_P_MAX_MOM_CUT_jointcc0pi ) {
        sel_muon_passed_mom_cuts_ = true;
      }

      if ( muon_mom >= MUON_P_MIN_WC_MOM_CUT_jointcc0pi && muon_mom <= MUON_P_MAX_MOM_CUT_jointcc0pi ) {
        sel_muon_passed_wc_mom_cuts_ = true;
      }

      // Apply muon candidate quality cut by comparing MCS and range-based
      // momentum estimators. Default to failing the cut.
      sel_muon_quality_ok_ = false;

      double frac_diff_range_mcs = std::abs( range_muon_mom - mcs_muon_mom );
      if ( range_muon_mom > 0. ) {
        frac_diff_range_mcs /= range_muon_mom;
        if ( frac_diff_range_mcs < MUON_MOM_QUALITY_CUT ) {
          sel_muon_quality_ok_ = true;
        }
      }
    }
    // Pion candidate
    else if (xgb_pid_vec_->at(p) == (int) BDTClass::kMu ||
             xgb_pid_vec_->at(p) == (int) BDTClass::kPi) {

      float track_length = Event->track_length_->at( p );
      if ( track_length <= 0. ) continue;

      float endx = Event->track_endx_->at( p );
      float endy = Event->track_endy_->at( p );
      float endz = Event->track_endz_->at( p );
      bool end_contained = this->in_proton_containment_vol( endx, endy, endz );

      float pi_mom;
      if ( end_contained ) {
        pi_mom = Event->track_range_mom_mu_->at(p);
      }
      else {
        pi_mom = Event->track_mcs_mom_mu_->at(p);
      }

      if (pi_mom > CHARGED_PI_MOM_CUT_jointcc0pi) {
        sel_num_pion_candidates_++;
        sel_has_pion_candidate_ = true;
      }

      if (pi_mom > CHARGED_PI_WC_MOM_CUT_jointcc0pi) {
        sel_num_pion_wc_candidates_++;
        sel_has_pion_wc_candidate_ = true;
      }
    }

    else if (xgb_pid_vec_->at(p) == (int) BDTClass::kP) {
      float track_length = Event->track_length_->at( p );
      if ( track_length <= 0. ) continue;

      sel_has_p_candidate_ = true;
      sel_passed_proton_pid_cut_ = true;
      sel_num_proton_candidates_++;

      // Check whether the current proton candidate fails the containment cut
      float endx = Event->track_endx_->at( p );
      float endy = Event->track_endy_->at( p );
      float endz = Event->track_endz_->at( p );
      bool end_contained = this->in_proton_containment_vol( endx, endy, endz );
      if ( !end_contained ) sel_protons_contained_ = false;
    }

    else if (xgb_pid_vec_->at(p) == (int) BDTClass::kOther) {}
    else if (xgb_pid_vec_->at(p) == (int) BDTClass::kInvalid) {}
    else {
      // hmm
    }
  }
  if(verbosal > 0) std::cout<< "Finished loop:line1945 "<<std::endl;

  sel_CC0pi_ = sel_nu_mu_cc_ && sel_no_reco_showers_
    && sel_muon_passed_mom_cuts_ && !sel_has_pion_candidate_
    && sel_n_bdt_other_ == 0 && sel_n_bdt_invalid_ == 0;

  sel_CC0pi_wc_ = sel_nu_mu_cc_ && sel_no_reco_showers_
    && sel_muon_passed_wc_mom_cuts_ && !sel_has_pion_wc_candidate_
    && sel_n_bdt_other_ == 0 && sel_n_bdt_invalid_ == 0;
  
  // Don't bother to apply the cuts that involve the leading
  // proton candidate if we don't have one
  if ( !sel_has_p_candidate_ ) {
    sel_CCNp0pi_ = false;
  }
  else {
    // All that remains is to apply the leading proton candidate cuts. We could
    // search for it above, but doing it here makes the code more readable (with
    // likely negligible impact on performance)
    this->find_lead_p_candidate(Event);

    // Check the range-based reco momentum for the leading proton candidate
    float lead_p_KE = Event->track_kinetic_energy_p_->at( lead_p_candidate_idx_ );
    float range_mom_lead_p = real_sqrt( lead_p_KE*lead_p_KE
      + 2.*PROTON_MASS*lead_p_KE );
    if ( range_mom_lead_p >= LEAD_P_MIN_MOM_CUT
      && range_mom_lead_p <= LEAD_P_MAX_MOM_CUT )
    {
      sel_lead_p_passed_mom_cuts_ = true;
    }

    // All right, we've applied all selection cuts. Set the flag that indicates
    // whether all were passed (and thus the event is selected as a CCNp0pi
    // candidate)
    sel_CCNp0pi_ = sel_nu_mu_cc_ && sel_no_reco_showers_
      && sel_muon_passed_mom_cuts_ && sel_muon_contained_ && sel_muon_quality_ok_
      && sel_has_p_candidate_ && sel_passed_proton_pid_cut_
      && sel_protons_contained_ && sel_lead_p_passed_mom_cuts_;
  }

 return sel_CC0pi_; 
}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void JOINTCC0pi::DefineOutputBranches(){

  // Signal definition flags
    
  SetBranch(&mc_neutrino_is_numu_,"mc_neutrino_is_numu",kBool);  
  SetBranch(&mc_vertex_in_FV_,"mc_vertex_in_FV",kBool);
  SetBranch(&mc_muon_in_mom_range_,"mc_muon_in_mom_range",kBool);
  SetBranch(&mc_lead_p_in_mom_range_,"mc_lead_p_in_mom_range",kBool);
  SetBranch(&mc_no_fs_pi0_,"mc_no_fs_pi0",kBool);
  SetBranch(&mc_no_charged_pi_above_threshold_,"mc_no_charged_pi_above_threshold",kBool);
  
  SetBranch(&mc_num_protons_,"mc_num_protons",kInteger);
  SetBranch(&mc_num_neutrons_,"mc_num_neutrons",kInteger);
  SetBranch(&mc_num_charged_pions_,"mc_num_charged_pions",kInteger);
  SetBranch(&mc_num_wc_charged_pions_,"mc_num_wc_charged_pions",kInteger);
 
 
  SetBranch(&mc_no_fs_mesons_,"mc_no_fs_mesons",kBool);
  SetBranch(&mc_is_signal_,"mc_is_signal",kBool);
  SetBranch(&mc_is_cc0pi_signal_,"mc_is_cc0pi_signal",kBool);
  SetBranch(&mc_is_cc0pi_wc_signal_,"mc_is_cc0pi_wc_signal",kBool);

  // MC event category
   SetBranch(&sel_nu_mu_cc_,"sel_nu_mu_cc",kBool);
   SetBranch(&sel_reco_vertex_in_FV_,"sel_reco_vertex_in_FV",kBool);
   SetBranch(&sel_topo_cut_passed_,"sel_topo_cut_passed",kBool);
   SetBranch(&sel_cosmic_ip_cut_passed_,"sel_cosmic_ip_cut_passed",kBool);
   SetBranch(&sel_pfp_starts_in_PCV_,"sel_pfp_starts_in_PCV",kBool);
   SetBranch(&sel_no_reco_showers_,"sel_no_reco_showers",kBool);
   SetBranch(&sel_has_muon_candidate_,"sel_has_muon_candidate",kBool);
   SetBranch(&sel_muon_contained_,"sel_muon_contained",kBool);
   SetBranch(&sel_muon_passed_mom_cuts_,"sel_muon_passed_mom_cuts",kBool);
   SetBranch(&sel_muon_quality_ok_,"sel_muon_quality_ok",kBool);
   SetBranch(&sel_has_p_candidate_,"sel_has_p_candidate",kBool);
   SetBranch(&sel_passed_proton_pid_cut_,"sel_passed_proton_pid_cut",kBool);
   SetBranch(&sel_protons_contained_,"sel_protons_contained",kBool);
   SetBranch(&sel_lead_p_passed_mom_cuts_,"sel_lead_p_passed_mom_cuts",kBool);
   SetBranch(&sel_CCNp0pi_,"sel_CCNp0pi",kBool);
   SetBranch(&sel_presel_,"sel_presel",kBool);  
   SetBranch(&sel_has_pion_candidate_,"sel_has_pion_candidate",kBool);
   SetBranch(&sel_CC0pi_,"sel_CC0pi",kBool);
   SetBranch(&sel_CC0pi_wc_,"sel_CC0pi_wc",kBool);
  
  // Still on the fence if I should make the WC conditions as a different selection
  SetBranch(&sel_muon_passed_wc_mom_cuts_,"sel_muon_passed_wc_mom_cuts",kBool);
  SetBranch(&sel_CCNp0pi_wc_,"sel_CCNp0pi_wc",kBool);
  SetBranch(&sel_has_pion_wc_candidate_,"sel_has_pion_wc_candidate",kBool);
  SetBranch(&mc_muon_in_wc_mom_range_,"mc_muon_in_wc_mom_range",kBool);
  SetBranch(&mc_no_charged_pi_above_wc_threshold_,"mc_no_charged_pi_above_wc_threshold",kBool);


  // XGBoost outputs

  SetBranch(&sel_n_bdt_other_,"sel_n_bdt_other",kInteger);
  SetBranch(&sel_n_bdt_muon_,"sel_n_bdt_muon",kInteger);
  SetBranch(&sel_n_bdt_pion_,"sel_n_bdt_pion",kInteger);
  SetBranch(&sel_n_bdt_proton_,"sel_n_bdt_proton",kInteger);
  SetBranch(&sel_n_bdt_invalid_,"sel_n_bdt_invalid",kInteger);


  // CC0pi
  SetBranch(&sel_num_proton_candidates_,"sel_num_proton_candidates",kInteger);
  SetBranch(&sel_num_pion_candidates_,"sel_num_pion_candidates",kInteger);
  SetBranch(&sel_num_pion_wc_candidates_,"sel_num_pion_wc_candidates",kInteger);
  SetBranch(&muon_candidate_idx_,"muon_candidate_idx",kInteger);
  SetBranch(&lead_p_candidate_idx_,"lead_p_candidate_idx",kInteger);

  // Index for the leading proton candidate in the vectors of PFParticles

   SetBranch(&distance_FV_surface_,"distance_FV_surface",kFloat);
   SetBranch(&muontrklen_,"muontrklen",kFloat);
   SetBranch(&muondistancetovectex_,"muondistancetovectex",kFloat);
   SetBranch(&muon_llr_pid_score_,"muon_llr_pid_score",kFloat);
   SetBranch(&muon_tkscore_,"muon_tkscore_",kFloat);

   SetBranch(&muon_trkend_x_,"muon_trkend_x",kFloat);
   SetBranch(&muon_trkend_y_,"muon_trkend_y",kFloat);
   SetBranch(&muon_trkend_z_,"muon_trkend_z",kFloat);

   SetBranch(&muon_trkstart_x_,"muon_trkstart_x",kFloat);
   SetBranch(&muon_trkstart_y_,"muon_trkstart_y",kFloat);
   SetBranch(&muon_trkstart_z_,"muon_trkstart_z",kFloat);
   SetBranch(&muon_trkchi2muon_,"muon_trkchi2muon",kFloat);  



  // Reco 3-momenta (muon, leading proton)

     SetBranch(&p3_mu_,"p3_mu",kTVector);
     SetBranch(&p3_lead_p_,"p3_lead_p",kTVector);
     SetBranch(&p3_lead_pi_,"p3_lead_pi",kTVector);

  //set_output_branch_address( out_tree, "p3_mag_corr_mu",
  //  &ev.p3_mag_corr_mu_, create, "p3_mag_corr_mu/F" );

  //set_output_branch_address( out_tree, "p3_costheta_corr_mu",
  //  &ev.p3_costheta_corr_mu_, create, "p3_costheta_corr_mu/F" );



  // Reco 3-momenta (all proton candidates, ordered from highest to lowest
  // magnitude)
  
    SetBranch(&p3_p_vec_,"p3_p_vec",kTVector);
    SetBranch(&mc_p3_mu_,"mc_p3_mu",kTVector);
    
    SetBranch(&mc_p3_lead_p_,"mc_p3_lead_p",kTVector);
    SetBranch(&mc_p3_lead_pi_,"mc_p3_lead_pi",kTVector);
    SetBranch(&mc_p3_p_vec_,"mc_p3_p_vec",kTVector);
  // True 3-momenta (muon, leading proton)

    SetBranch(&mc_muontrklen_,"mc_muontrklen",kFloat);


  // True 3-momenta (all protons, ordered from highest to lowest magnitude)

 
  // XGBoost scores
  //set_object_output_branch_address< std::vector<int> >( out_tree,
  //  "xgb_pid_vec", &xgb_pid_vec_, create );
   SetBranch(&xgb_pid_vec_,"xgb_pid_vec",kVectorInteger);

  //set_object_output_branch_address< std::vector<std::vector<float > > >( out_tree,
    //"xgb_score_vec", ev.xgb_score_vec_, create );
   
   SetBranch(&xgb_score_vec_,"xgb_score_vec",kVectorVectorFloat);
  // Reco STVs
   SetBranch(&delta_pT_,"delta_pT",kFloat);
   SetBranch(&delta_phiT_,"delta_phiT",kFloat);
   SetBranch(&delta_alphaT_,"delta_alphaT",kFloat);
   SetBranch(&delta_pL_,"delta_pL",kFloat);
   SetBranch(&pn_,"pn",kFloat);
   SetBranch(&delta_pTx_,"delta_pTx",kFloat);
   
   SetBranch(&delta_pTy_,"delta_pTy",kFloat);
   SetBranch(&theta_mu_p_,"theta_mu_p",kFloat);
     // MC STVs (only filled for signal events)
   SetBranch(&mc_delta_pT_,"mc_delta_pT",kFloat);
   SetBranch(&mc_delta_phiT_,"mc_delta_phiT",kFloat);
   SetBranch(&mc_delta_alphaT_,"mc_delta_alphaT",kFloat);
   SetBranch(&mc_delta_pL_,"mc_delta_pL",kFloat);
   SetBranch(&mc_pn_,"mc_pn",kFloat);
   SetBranch(&mc_delta_pTx_,"mc_delta_pTx",kFloat);
   SetBranch(&mc_delta_pTy_,"mc_delta_pTy",kFloat);
   SetBranch(&mc_theta_mu_p_,"mc_theta_mu_p",kFloat);
   SetBranch(&mc_distance_FV_surface_,"mc_distance_FV_surface",kFloat);
  // *** Branches copied directly from the input ***

  // Cosmic rejection parameters for numu CC inclusive selection
 
 // MC truth information for the neutrino
  






  // Log-likelihood-based particle ID information
    
  // MC truth information for the final-state primary particles



}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void JOINTCC0pi::DefineCategoryMap(){
  // Use the shared category map for 1p/2p/Np/Xp
  categ_map_ = JOINTCC0Pi_MAP;
  std::cout<<"Finished :: DefineCategoryMap "<< std::endl;
}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void JOINTCC0pi::classify_tracks(AnalysisEvent* Event){
  xgb_pid_vec_->clear();
  xgb_pid_vec_->resize(Event->num_pf_particles_, (int) BDTClass::kInvalid);

  xgb_score_vec_->clear();
  xgb_score_vec_->resize(Event->num_pf_particles_);

  sel_n_bdt_other_ = 0;
  sel_n_bdt_muon_ = 0;
  sel_n_bdt_pion_ = 0;
  sel_n_bdt_proton_ = 0;
  sel_n_bdt_invalid_ = 0;

  // Only make predictions for events passing the BDT pre-selection cuts
  if (!sel_presel_) return;

  // Classify all FS objects
  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {
    // Only consider generation 2 tracks
    float track_score =  Event->pfp_track_score_->at( p );
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2 || track_score <= TRACK_SCORE_CUT_jointcc0pi ) continue;

    // Observables vector
    float pmcs = Event->track_mcs_mom_mu_->at(p);
    float prange = Event->track_range_mom_mu_->at(p);
    float trk_relative_dif_mcs_range = (pmcs - prange) / prange;

    float vx = Event->track_endx_->at(p);
    float vy = Event->track_endx_->at(p);
    float vz = Event->track_endx_->at(p);
    bool trk_contained = point_inside_FV(this->ReturnRecoFV(), vx, vy, vz );

    unsigned int num_daughters = Event->pfp_trk_daughters_count_->at(p) +
                                 Event->pfp_shr_daughters_count_->at(p);

    float fs_v[] = {
      Event->track_start_distance_->at(p),
      Event->track_llr_pid_score_->at(p),
      Event->track_chi2_proton_->at(p),
      Event->track_kinetic_energy_p_->at(p),
      (float) trk_contained,
      (float) num_daughters,
      trk_relative_dif_mcs_range
    };

    // Check for inf/nans
    bool fs_v_ok = true;
    for (size_t i=0; i<7; i++) {
      if (std::isnan(fs_v[i]) || std::isinf(fs_v[i])) {
        fs_v_ok = false;
        break;
      }
    }

    if (fs_v_ok) {
      int r;  // XGBoost return codes

      // Create input DMatrix
      DMatrixHandle dmat;
      r = XGDMatrixCreateFromMat(fs_v, 1, 7, -1, &dmat);
      if (r != 0) std::cout << XGBGetLastError() << std::endl;
      assert(r == 0);

      // Load trained model from JSON
      const char* c_json_config = "{\"type\": 0,\"training\": false,\"iteration_begin\": 0,\"iteration_end\": 0,\"strict_shape\": false}";
      const unsigned long int* out_shape;
      unsigned long int out_dim;  
      const float* out_result;

      r = XGBoosterPredictFromDMatrix(*booster, dmat, c_json_config,
                                      &out_shape, &out_dim, &out_result);
      if (r != 0) std::cout << XGBGetLastError() << std::endl;
      assert(r == 0);

      XGDMatrixFree(dmat);

      int idx_max = (int) BDTClass::kInvalid;
      float v_max = -1;

      //for (int k=0; k<out_dim; k++) {
      for (int k=0; k < 4; k++) {
          xgb_score_vec_->at(p).push_back(out_result[k]);
          if (out_result[k] > v_max) {
             v_max = out_result[k];
             idx_max = k;
          }
       }
      xgb_pid_vec_->at(p) = idx_max;

      switch ((BDTClass) idx_max) {
        case BDTClass::kOther: sel_n_bdt_other_++; break;
        case BDTClass::kMu: sel_n_bdt_muon_++; break;
        case BDTClass::kPi: sel_n_bdt_pion_++; break;
        case BDTClass::kP: sel_n_bdt_proton_++; break;
        case BDTClass::kInvalid: sel_n_bdt_invalid_++; break;
      }
    }
  }
}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void JOINTCC0pi::find_muon_candidate(AnalysisEvent* Event){

  std::vector<int> muon_candidate_indices;
  std::vector<int> muon_pid_scores;

  if (sel_presel_) {
    for ( int p = 0; p < Event->num_pf_particles_; ++p ) {
      // Only direct neutrino daughters (generation == 2) will be considered as
      // possible muon candidates
      unsigned int generation = Event->pfp_generation_->at( p );
      if ( generation != 2u ) continue;

      if (xgb_pid_vec_->at(p) == (int) BDTClass::kMu) {
        muon_candidate_indices.push_back( p );
        muon_pid_scores.push_back( xgb_score_vec_->at(p)[1] );
      }
    }
  }

  size_t num_candidates = muon_candidate_indices.size();
  if ( num_candidates > 0u ) sel_has_muon_candidate_ = true;

  if ( num_candidates == 1u ) {
    muon_candidate_idx_ = muon_candidate_indices.front();
  }
  else if ( num_candidates > 1u ) {
    // In the case of multiple muon candidates, choose the one with the highest
    // PID score (most muon-like) as the one to use
    float highest_score = LOW_FLOAT;
    int chosen_index = BOGUS_INDEX;
    for ( size_t c = 0; c < num_candidates; ++c ) {
      float score = muon_pid_scores.at( c );
      if ( highest_score < score ) {
        highest_score = score;
        chosen_index = muon_candidate_indices.at( c );
      }
    }
    muon_candidate_idx_ = chosen_index;
  }
  else {
    muon_candidate_idx_ = BOGUS_INDEX;
  }
}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void JOINTCC0pi::find_lead_p_candidate(AnalysisEvent* Event){
  float lead_p_track_length = LOW_FLOAT;
  size_t lead_p_index = 0u;
  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    // Skip the muon candidate reco track (this function assumes that it has
    // already been found)
    if ( p == muon_candidate_idx_ ) continue;

    // Skip PFParticles that are shower-like (track scores near 0)
    float track_score = Event->pfp_track_score_->at( p );
    if ( track_score <= TRACK_SCORE_CUT_jointcc0pi ) continue;

    // All non-muon-candidate reco tracks are considered proton candidates
    float track_length = Event->track_length_->at( p );
    if ( track_length <= 0. ) continue;

    if ( track_length > lead_p_track_length ) {
      lead_p_track_length = track_length;
      lead_p_index = p;
    }
  }

  // If the leading proton track length changed from its initial
  // value, then we found one. Set the index appropriately.
  if ( lead_p_track_length != LOW_FLOAT ) lead_p_candidate_idx_ = lead_p_index;
  // Otherwise, set the index to BOGUS_INDEX
  else lead_p_candidate_idx_ = BOGUS_INDEX;
}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void JOINTCC0pi::compute_stvs( const TVector3& p3mu, 
   const TVector3& p3p, float& delta_pT,
  float& delta_phiT, float& delta_alphaT, 
  float& delta_pL, float& pn,
  float& delta_pTx, float& delta_pTy )
{
  delta_pT = (p3mu + p3p).Perp();

  delta_phiT = std::acos( (-p3mu.X()*p3p.X() - p3mu.Y()*p3p.Y())
    / (p3mu.XYvector().Mod() * p3p.XYvector().Mod()) );

  TVector2 delta_pT_vec = (p3mu + p3p).XYvector();
  delta_alphaT = std::acos( (-p3mu.X()*delta_pT_vec.X()
    - p3mu.Y()*delta_pT_vec.Y())
    / (p3mu.XYvector().Mod() * delta_pT_vec.Mod()) );

  float Emu = std::sqrt(std::pow(MUON_MASS, 2) + p3mu.Mag2());
  float Ep = std::sqrt(std::pow(PROTON_MASS, 2) + p3p.Mag2());
  float R = TARGET_MASS + p3mu.Z() + p3p.Z() - Emu - Ep;

  // Estimated mass of the final remnant nucleus (CCQE assumption)
  float mf = TARGET_MASS - NEUTRON_MASS + BINDING_ENERGY;
  delta_pL = 0.5*R - (std::pow(mf, 2) + std::pow(delta_pT, 2)) / (2.*R);

  pn = std::sqrt( std::pow(delta_pL, 2) + std::pow(delta_pT, 2) );

  // Components of the 2D delta_pT vector (see arXiv:1910.08658)

  // We assume that the neutrino travels along the +z direction (also done
  // in the other expressions above)
  TVector3 zUnit( 0., 0., 1. );

  // Defines the x direction for the components of the delta_pT vector
  TVector2 xTUnit = zUnit.Cross( p3mu ).XYvector().Unit();

  delta_pTx = xTUnit.X()*delta_pT_vec.X() + xTUnit.Y()*delta_pT_vec.Y();

  // Defines the y direction for the components of the delta_T vector
  TVector2 yTUnit = ( -p3mu ).XYvector().Unit();

  delta_pTy = yTUnit.X()*delta_pT_vec.X() + yTUnit.Y()*delta_pT_vec.Y();
}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
float JOINTCC0pi::point_distance_to_FV( float x, float y, float z )
{
    double FV_X_MIN =   21.5;
    double FV_X_MAX =  234.85;
    
    double FV_Y_MIN = -95.0;
    double FV_Y_MAX =  95.0;
    
    double FV_Z_MIN =   21.5;
    double FV_Z_MAX =  966.8;
    // Distances to the surfaces along each axis
    float dist_to_x_min = x - FV_X_MIN;
    float dist_to_x_max = FV_X_MAX - x;
    float dist_to_y_min = y - FV_Y_MIN;
    float dist_to_y_max = FV_Y_MAX - y;
    float dist_to_z_min = z - FV_Z_MIN;
    float dist_to_z_max = FV_Z_MAX - z;

    // Find the minimum distance to any surface
    float min_distance = std::min({
        dist_to_x_min, dist_to_x_max,
        dist_to_y_min, dist_to_y_max,
        dist_to_z_min, dist_to_z_max
    });

    return min_distance;
    
}
/////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void JOINTCC0pi::apply_numu_CC_selection(AnalysisEvent* Event)
{

  sel_reco_vertex_in_FV_ = this->reco_vertex_inside_FV(Event);
  sel_topo_cut_passed_ = Event->topological_score_ > TOPO_SCORE_CUT_jointcc0pi;
  sel_cosmic_ip_cut_passed_ = Event->cosmic_impact_parameter_ > COSMIC_IP_CUT_jointcc0pi;

  // Apply the containment cut to the starting positions of all
  // reconstructed tracks and showers. Pass this cut by default.
  sel_pfp_starts_in_PCV_ = true;

  // Loop over each PFParticle in the event
  for ( int p = 0; p < Event->num_pf_particles_; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = Event->pfp_generation_->at( p );
    if ( generation != 2u ) continue;

    // Use the track reconstruction results to get the start point for
    // every PFParticle for the purpose of verifying containment. We could
    // in principle differentiate between tracks and showers here, but
    // (1) we cut out all showers later on in the selection anyway, and
    // (2) the blinded PeLEE data ntuples do not include shower information.
    // We therefore apply the track reconstruction here unconditionally.
    float x = Event->track_startx_->at( p );
    float y = Event->track_starty_->at( p );
    float z = Event->track_startz_->at( p );

    // Verify that the start of the PFParticle lies within the containment
    // volume.
    // TODO: revisit which containment volume to use for PFParticle start
    // positions. See https://stackoverflow.com/a/2488507 for an explanation
    // of the use of &= here. Don't worry, it's type-safe since both operands
    // are bool.
    sel_pfp_starts_in_PCV_ &= this->in_proton_containment_vol( x, y, z );
  }

  sel_presel_ = (Event->nslice_ == 1 &&
                 sel_reco_vertex_in_FV_ &&
                 sel_pfp_starts_in_PCV_ &&
                 sel_topo_cut_passed_ &&
                 sel_no_reco_showers_);

  // Apply BDT to classify final state tracks
  this->classify_tracks(Event);

  // Sets the sel_has_muon_candidate_ flag as appropriate. The threshold check
  // is handled later.
  this->find_muon_candidate(Event);

  sel_nu_mu_cc_ = sel_presel_ && sel_has_muon_candidate_;
}
/////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////

bool JOINTCC0pi::reco_vertex_inside_FV(AnalysisEvent* Event) {
      return point_inside_FV( this->ReturnRecoFV(),Event->nu_vx_, Event->nu_vy_, Event->nu_vz_ );
}
/////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
float JOINTCC0pi::reco_distance_to_FV_Surface(AnalysisEvent* Event){
   return point_distance_to_FV( Event->nu_vx_, Event->nu_vy_, Event->nu_vz_  );
}
/////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
bool JOINTCC0pi::mc_vertex_inside_FV(AnalysisEvent* Event) {
      return point_inside_FV( this->ReturnTrueFV(),Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_ );
}
/////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
float JOINTCC0pi::mc_distance_to_FV_Surface(AnalysisEvent* Event){
       return point_distance_to_FV( Event->mc_nu_vx_, Event->mc_nu_vy_, Event->mc_nu_vz_  );
}
/////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
bool JOINTCC0pi::in_proton_containment_vol( float x, float y, float z ) {
      bool x_inside_PCV = ( PCV_X_MIN < x ) && ( x < PCV_X_MAX );
      bool y_inside_PCV = ( PCV_Y_MIN < y ) && ( y < PCV_Y_MAX );
      bool z_inside_PCV = ( PCV_Z_MIN < z ) && ( z < PCV_Z_MAX );
      return ( x_inside_PCV && y_inside_PCV && z_inside_PCV );
}
/////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
float JOINTCC0pi::distanceBetweentwopoint(float x1, float y1, float z1, float x2, float y2, float z2) {
    // Calculate the differences in each coordinate
    float dx = x2 - x1;
    float dy = y2 - y1;
    float dz = z2 - z1;

    // Calculate the Euclidean distance
    return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
}