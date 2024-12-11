#include "TVector3.h"

#include "/exp/uboone/app/users/kwresilo/CCNp1pi/FiducialVolume.hh"
#include "/exp/uboone/app/users/kwresilo/CCNp1pi/TreeUtils.hh"


struct TempEvent {

  TempEvent() {}

  void set_branch_addresses( TTree& etree );

  /*
  float delta_pT_;
  float delta_pL_;
  float delta_pL2_;
  float pn_;
  float pn2_;
  float theta_mu_p_;
  float theta_mu_cpi_;
  float delta_alpha3D_;
  float delta_phi3d_had_;
  float delta_phi3d_mu_;

  float mc_delta_pT_;
  float mc_delta_pL_;
  float mc_delta_pL2_;
  float mc_pn_;
  float mc_pn2_;
  float mc_theta_mu_p_;
  float mc_theta_mu_cpi_;
  float mc_delta_alpha3D_;
  float mc_delta_phi3d_had_;
  float mc_delta_phi3d_mu_;

  MyPointer< std::vector<int> > pfp_generation_;
  MyPointer< std::vector<int> > pfp_true_pdg_;
  MyPointer< std::vector<float> > pfp_true_px_;
  MyPointer< std::vector<float> > pfp_true_py_;
  MyPointer< std::vector<float> > pfp_true_pz_;
  MyPointer< std::vector<float> > pfp_track_score_;
  MyPointer< std::vector<float> > track_range_mom_mu_;
  MyPointer< std::vector<float> > track_startx_;
  MyPointer< std::vector<float> > track_starty_;
  MyPointer< std::vector<float> > track_startz_;
  MyPointer< std::vector<float> > track_endx_;
  MyPointer< std::vector<float> > track_endy_;
  MyPointer< std::vector<float> > track_endz_;
  MyPointer< std::vector<float> > track_dirx_;
  MyPointer< std::vector<float> > track_diry_;
  MyPointer< std::vector<float> > track_dirz_;
  MyPointer< std::vector<float> > track_length_;

  */
  MyPointer< std::vector<int> > mc_pdg_;

  int pion_candidate_idx_;
  int muon_candidate_pid_idx_;
  int lead_p_candidate_idx_;

  Bool_t mc_is_signal_;
  Bool_t mc_pi_stopping_;
  Bool_t mc_golden_;

  Bool_t sel_CCNp1pi_;
  Bool_t sel_has_pion_candidate_;
  Bool_t sel_has_muon_candidate_;
  Bool_t sel_has_p_candidate_;

  TVector3 *mc_p3_cpi_ = new TVector3();
  TVector3 *p3_cpi_ = new TVector3();
  TVector3 *mc_p3_mu_ = new TVector3();
  TVector3 *p3_mu_ = new TVector3();
  TVector3 *mc_p3_lead_p_ = new TVector3();
  TVector3 *p3_lead_p_ = new TVector3();

  double gki_proton_KE_, gki_Ecal_, gki_Q_, gki_Pt_, gki_Pl_;
  double gki_PtMuon_, gki_PtProton_, gki_PtPion_, gki_PlMuon_, gki_PlProton_, gki_PlPion_;
  double gki_Pn_, gki_DeltaAlpha3D_, gki_DeltaAlpha3DMu_, gki_DeltaPhi3D_;
  double gki_DeltaPhi3D_pion_, gki_DeltaPhi3D_proton_;

  double gki_Total_KE_, gki_Total_Ecal_, gki_Total_Q_, gki_Total_Pt_, gki_Total_Pl_;
  double gki_Total_PtMuon_, gki_Total_PtProton_, gki_Total_PtPion_;
  double gki_Total_PlMuon_, gki_Total_PlProton_, gki_Total_PlPion_;
  double gki_Total_Pn_, gki_Total_DeltaAlpha3D_, gki_Total_DeltaAlpha3DMu_;
  double gki_Total_DeltaPhi3D_, gki_Total_DeltaPhi3D_pion_, gki_Total_DeltaPhi3D_proton_;

  double mc_gki_proton_KE_, mc_gki_Ecal_, mc_gki_Q_, mc_gki_Pt_, mc_gki_Pl_;
  double mc_gki_PtMuon_, mc_gki_PtProton_, mc_gki_PtPion_, mc_gki_PlMuon_, mc_gki_PlProton_;
  double mc_gki_PlPion_, mc_gki_Pn_, mc_gki_DeltaAlpha3D_, mc_gki_DeltaAlpha3DMu_;
  double mc_gki_DeltaPhi3D_, mc_gki_DeltaPhi3D_pion_, mc_gki_DeltaPhi3D_proton_;

  double mc_gki_Total_KE_, mc_gki_Total_Ecal_, mc_gki_Total_Q_, mc_gki_Total_Pt_, mc_gki_Total_Pl_;
  double mc_gki_Total_PtMuon_, mc_gki_Total_PtProton_, mc_gki_Total_PtPion_;
  double mc_gki_Total_PlMuon_, mc_gki_Total_PlProton_, mc_gki_Total_PlPion_;
  double mc_gki_Total_Pn_, mc_gki_Total_DeltaAlpha3D_, mc_gki_Total_DeltaAlpha3DMu_;
  double mc_gki_Total_DeltaPhi3D_, mc_gki_Total_DeltaPhi3D_pion_, mc_gki_Total_DeltaPhi3D_proton_
};

void TempEvent::set_branch_addresses( TTree& etree ) {
  set_object_input_branch_address( etree, "mc_pdg", mc_pdg_ );

  etree.SetBranchAddress("pion_candidate_idx", &pion_candidate_idx_);
  etree.SetBranchAddress("muon_candidate_pid_idx", &muon_candidate_pid_idx_);
  etree.SetBranchAddress("lead_p_candidate_idx", &lead_p_candidate_idx_);

  etree.SetBranchAddress("MC_Signal", &mc_is_signal_);
  etree.SetBranchAddress("true_pi_stopping", &mc_pi_stopping_);
  etree.SetBranchAddress("true_pi_golden", &mc_golden_);

  etree.SetBranchAddress("Selected", &sel_CCNp1pi_);
  etree.SetBranchAddress("has_pion_candidate", &sel_has_pion_candidate_);
  etree.SetBranchAddress("has_muon_candidate", &sel_has_muon_candidate_);
  etree.SetBranchAddress("has_p_candidate", &sel_has_p_candidate_);

  etree.SetBranchAddress("true_p3_cpi", &mc_p3_cpi_);
  etree.SetBranchAddress("reco_p3_cpi", &p3_cpi_);
  etree.SetBranchAddress("true_p3_mu", &mc_p3_mu_);
  etree.SetBranchAddress("reco_p3_mu", &p3_mu_);
  etree.SetBranchAddress("true_p3_lead_p", &mc_p3_lead_p_);
  etree.SetBranchAddress("reco_p3_lead_p", &p3_lead_p_);

}

void GoldenPionDependency(){

  TChain stv_tree( "stv_tree" );
  stv_tree.Add("/exp/uboone/data/users/kwresilo/CCNp1pi/data/GKI_v1/xsec-ana-nu_overlay_run1.roottest.root");
  //stv_tree.Add("/exp/uboone/data/users/kwresilo/CCNp1pi/data/GKI_v1/xsec-ana-nu_overlay_run2.roottest.root");
  //stv_tree.Add("/exp/uboone/data/users/kwresilo/CCNp1pi/data/GKI_v1/xsec-ana-nu_overlay_run3.roottest.root"); 

  TH1D *h1 = new TH1D("pi momentum", "", 41, -1., 1.);  
  TH1D *h2 = new TH1D("pi momentum golden", "", 41, -1., 1.);

  
  long entry = 0;
  while (true) {
    TempEvent cur_event;
      cur_event.set_branch_addresses( stv_tree );
      //if ( entry > 200000) break; 
      
      int local_entry = stv_tree.LoadTree( entry );
      if ( local_entry < 0 ) break;
      
      if (entry % 10000 == 0) std::cout << "Event " << entry <<std::endl; 
      stv_tree.GetEntry( entry );
      ++entry;

      if ( !(cur_event.mc_is_signal_ && cur_event.sel_has_pion_candidate_ && cur_event.sel_has_p_candidate_ && cur_event.sel_has_muon_candidate_) ) continue;
      else {
        std::cout << "Event " << entry << " is a signal event with pion, muon and proton candidate" << std::endl;
      }
  }
}