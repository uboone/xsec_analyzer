#ifndef __JOINTCC0pi_h__
#define __JOINTCC0pi_h__


#include "XSecAnalyzer/Selections/SelectionBase.hh"

// XGBoost
#include <xgboost/c_api.h>

#define DEFAULT_NUM 20

class JOINTCC0pi : public SelectionBase {
 public:
  /////Constructor
  JOINTCC0pi();
  
  ///destructor 
   ~JOINTCC0pi() {
        XGBoosterFree(*booster);
        std::cout << "Destructor: Memory deallocated." << std::endl;
    }
    
  /////////////
  //BDT Model// 
  BoosterHandle* booster;
  
  int CategorizeEvent(AnalysisEvent* Event);
  bool Selection(AnalysisEvent* Event);
  bool DefineSignal(AnalysisEvent* Event);
  void ComputeRecoObservables(AnalysisEvent* Event);
  void ComputeTrueObservables(AnalysisEvent* Event);
  void DefineOutputBranches();
  void DefineConstants();
  void DefineCategoryMap();
  void apply_numu_CC_selection(AnalysisEvent* Event);
  void find_muon_candidate(AnalysisEvent* Event);
  void find_lead_p_candidate(AnalysisEvent* Event);
  void classify_tracks(AnalysisEvent* Event);
  bool in_proton_containment_vol( float x, float y, float z );
  float distanceBetweentwopoint(float x1, float y1, 
                                float z1, float x2,
                                float y2, float z2);
  void compute_stvs( const TVector3& p3mu, 
                    const TVector3& p3p, float& delta_pT,
                     float& delta_phiT, float& delta_alphaT, 
                     float& delta_pL, float& pn,
                     float& delta_pTx, float& delta_pTy );
  bool reco_vertex_inside_FV(AnalysisEvent* Event);
  bool mc_vertex_inside_FV(AnalysisEvent* Event);
  float reco_distance_to_FV_Surface(AnalysisEvent* Event);
  float mc_distance_to_FV_Surface(AnalysisEvent* Event);
 float point_distance_to_FV( float x, float y, float z );
 void SetVerbosal(int input ){verbosal = input;}
  //////////////////////////////////////////////////////
  // Print Error Messages > 0
  int verbosal = 1;
  ///////////////////////////////////////////////////////
  // cc0pi Cut Values 
  //////////////////////////////////////////////////////
    const float MUON_P_MIN_MOM_CUT_jointcc0pi =.100; // GeV/c
    const float MUON_P_MIN_WC_MOM_CUT_jointcc0pi =.121; // GeV/c
    const float MUON_P_MAX_MOM_CUT_jointcc0pi = 2.000; // GeV/c
    const float CHARGED_PI_MOM_CUT_jointcc0pi =.07; // GeV/c
    const float CHARGED_PI_WC_MOM_CUT_jointcc0pi =.161; // GeV/c
    const float MUON_MOM_QUALITY_CUT_jointcc0pi= .25; // fractional difference
   
    const float TOPO_SCORE_CUT_jointcc0pi =  .15;
    const float COSMIC_IP_CUT_jointcc0pi = 25.; // cm
    const float TRACK_SCORE_CUT_jointcc0pi  = .5;
// REMOVE ME
//constexpr float MUON_TRACK_SCORE_CUT .8;
//constexpr float MUON_VTX_DISTANCE_CUT = 4.; // cm
//constexpr float MUON_LENGTH_CUT = 10.; // cm
//constexpr float MUON_PID_CUT .2;

 
  
  float nu_vx_;
  float nu_vy_;
  float nu_vz_;
  
  int num_pf_particles_ ;
  float track_chi2_muon_; 
  

  


////////////////////////////////////////////////////////////
///// Only need to declear varibles that are exclusive to the Joint CC0pi 
////////////////////////////////////////////////////////////


    // Signal definition requirements
    bool mc_neutrino_is_numu_ ;
    bool mc_vertex_in_FV_ ;
    bool mc_muon_in_mom_range_ ;
    bool mc_lead_p_in_mom_range_ ;
    bool mc_no_fs_mesons_ ;
    // Intersection of all of these requirements
    bool mc_is_signal_ ;

    // Extra flags for looking specifically at final-state pions
    bool mc_no_fs_pi0_ ;
    bool mc_no_charged_pi_above_threshold_ ;

    bool mc_is_cc0pi_signal_ ;

    int mc_num_protons_ ;
    int mc_num_neutrons_ ;
    int mc_num_charged_pions_ ;
    int mc_num_wc_charged_pions_ ;

    // Water Cherenkov cuts
    bool mc_muon_in_wc_mom_range_ ;
    bool mc_no_charged_pi_above_wc_threshold_ ;
    bool mc_is_cc0pi_wc_signal_ ;


    // **** Reco selection requirements ****

    // Whether the event passed the numu CC selection (a subset of the cuts
    // used for the full analysis)
    bool sel_nu_mu_cc_ ;

    // Whether the reconstructed neutrino vertex lies within the fiducial
    // volume
    bool sel_reco_vertex_in_FV_;
    // Whether the event passed the topological score cut
    bool sel_topo_cut_passed_;
    // Whether the event passed the cosmic impact parameter cut
    bool sel_cosmic_ip_cut_passed_ ;
    // Whether the start points for all PFParticles lie within the
    // proton containment volume
    bool sel_pfp_starts_in_PCV_ ;

    // True if a generation == 2 muon candidate was identified
    bool sel_has_muon_candidate_ ;

    // Whether the end point of the muon candidate track is contained
    // in the "containment volume"
    bool sel_muon_contained_ ;

    // Whether the muon candidate has MCS- and range-based reco momenta
    // that agree within a given tolerance
    bool sel_muon_quality_ok_ ;

    // Whether the muon candidate has a reco momentum above threshold
    bool sel_muon_passed_mom_cuts_ ;
    bool sel_muon_passed_wc_mom_cuts_ ;

    // False if at least one generation == 2 shower was reconstructed
    bool sel_no_reco_showers_ ;

    // Whether at least one generation == 2 reco track exists that is not the
    // muon candidate
    bool sel_has_p_candidate_ ;

    // Whether all proton candidates (i.e., all tracks which are not the muon
    // candidate) pass the proton PID cut or not
    bool sel_passed_proton_pid_cut_ ;

    // Whether all proton candidates have track end coordinates that lie within
    // the "containment volume"
    bool sel_protons_contained_ ;

    // Whether the leading proton candidate has a range-based reco momentum
    // above LEAD_P_MIN_MOM_CUT and below LEAD_P_MAX_MOM_CUT
    bool sel_lead_p_passed_mom_cuts_ ;

    // Intersection of all of the above requirements
    bool sel_CCNp0pi_ ;
    bool sel_CCNp0pi_wc_ ;
    bool sel_presel_ ;

    int sel_n_bdt_other_ ;
    int sel_n_bdt_muon_ ;
    int sel_n_bdt_pion_ ;
    int sel_n_bdt_proton_ ;
    int sel_n_bdt_invalid_ ;

    int sel_num_proton_candidates_ ;
    int sel_num_pion_candidates_ ;
    int sel_num_pion_wc_candidates_ ;
    bool sel_has_pion_candidate_ ;
    bool sel_has_pion_wc_candidate_ ;
    bool sel_CC0pi_ ;
    bool sel_CC0pi_wc_ ;

    // Muon and leading proton candidate indices (BOGUS_INDEX if not present)
    // in the reco track arrays
    int muon_candidate_idx_ ;
    int lead_p_candidate_idx_ ;

    // ** Reconstructed observables **

    // 3-momenta
    MyPointer< TVector3 > p3_mu_;
    MyPointer< TVector3 > p3_mu_range_;
    MyPointer< TVector3 > p3_mu_mcs_;
    //float p3_mag_corr_mu_ ;
    //float p3_costheta_corr_mu_ ;

    MyPointer< TVector3 > p3_lead_p_;
    MyPointer< TVector3 > p3_lead_pi_;

    MyPointer< std::vector<bool> > trk_bragg_mu_fwd_preferred_v_;
    // Reconstructed 3-momenta for all proton candidates,
    // ordered from highest to lowest by magnitude
    MyPointer< std::vector<TVector3> > p3_p_vec_;

    // XGBoost class predictions
    MyPointer< std::vector<std::vector<float> > > xgb_score_vec_;
    MyPointer< std::vector<int> > xgb_pid_vec_;

    // Reco STVs and other variables of interest
    float delta_pT_;
    float delta_phiT_;
    float delta_alphaT_;
    float delta_pL_;
    float pn_;
    float delta_pTx_;
    float delta_pTy_;
    float theta_mu_p_;

    float muontrklen_;
    float distance_FV_surface_;
    float muondistancetovectex_;
    float muon_llr_pid_score_;
    
    float muon_trkstart_x_;
    float muon_trkstart_y_;
    float muon_trkstart_z_;
    float muon_trkend_x_;
    float muon_trkend_y_;
    float muon_trkend_z_;
    float muon_tkscore_; 
    float muon_trkchi2muon_;
    
    
    
    
    // ** MC truth observables **
    // These are loaded for signal events whenever we have MC information
    // to use

    // 3-momenta
    MyPointer< TVector3 > mc_p3_mu_;
    MyPointer< TVector3 > mc_p3_lead_p_;
    MyPointer< TVector3 > mc_p3_lead_pi_;
    // True 3-momenta for all true MC protons, ordered from highest to lowest
    // by magnitude
    MyPointer< std::vector<TVector3> > mc_p3_p_vec_;

    // MC truth STVs and other variables of interest
    float mc_delta_pT_ ;
    float mc_delta_phiT_ ;
    float mc_delta_alphaT_ ;
    float mc_delta_pL_ ;
    float mc_pn_;
    float mc_delta_pTx_ ;
    float mc_delta_pTy_ ;
    float mc_theta_mu_p_ ;
    float mc_distance_FV_surface_ ;
    float mc_muontrklen_;

  STVCalcType CalcType;
  
  private:
};

#endif
