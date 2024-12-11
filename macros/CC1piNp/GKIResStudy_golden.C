#include "TVector3.h"

#include "/exp/uboone/app/users/kwresilo/CCNp1pi/FiducialVolume.hh"
#include "/exp/uboone/app/users/kwresilo/CCNp1pi/TreeUtils.hh"

constexpr int MUON_PDG = 13;
constexpr int PION_PDG = 211;

// Boundaries of the containment volume in cm. Chosen to leave a 10-cm border
// between all edges and the active volume boundary
double VOL_X_MIN =   21.5;
double VOL_X_MAX =  234.85;

double VOL_Y_MIN = -95.0;
double VOL_Y_MAX =  95.0;

double VOL_Z_MIN =   21.5;
double VOL_Z_MAX =  966.8;

const int FONT_STYLE = 132;
const double TEXT_SIZE  = 0.06;

const double PION_MOM_MIN = 0.05;
const double PION_MOM_MAX = 0.7;

const double MUON_MOM_MIN = 0.15;
const double MUON_MOM_MAX = 1.2;

const double PROTON_MOM_MIN = 0.3;
const double PROTON_MOM_MAX = 1.;

bool in_FV( double x, double y, double z ) {
  bool x_inside_V = ( VOL_X_MIN < x ) && ( x < VOL_X_MAX );
  bool y_inside_V = ( VOL_Y_MIN < y ) && ( y < VOL_Y_MAX );
  bool z_inside_V = ( VOL_Z_MIN < z ) && ( z < VOL_Z_MAX );
  return ( x_inside_V && y_inside_V && z_inside_V );
}

// Set up temporary event memebrs
// Useful in this way since variables all relate to same event
struct TempEvent {

  TempEvent() {}

  void set_branch_addresses( TTree& etree );

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

  MyPointer< std::vector<int> > mc_pdg_;
  MyPointer< std::vector<int> > mc_n_inelastic_;
  MyPointer< std::vector<int> > mc_n_elastic_;
  MyPointer< std::vector<float> > mc_end_p_;

  int pion_candidate_idx_;
  int muon_candidate_pid_idx_;
  int lead_p_candidate_idx_;

  Bool_t mc_is_signal_;

  Bool_t sel_CCNp1pi_;
  Bool_t sel_has_pion_candidate_;

  Bool_t sel_has_muon_candidate_;
  Bool_t sel_muon_contained_;
  Bool_t sel_muon_quality_ok_;
  Bool_t sel_muon_passed_mom_cuts_;
  Bool_t sel_has_p_candidate_;

  TVector3 *mc_p3_cpi_ = new TVector3();
  TVector3 *p3_cpi_ = new TVector3();
  TVector3 *mc_p3_mu_ = new TVector3();
  TVector3 *p3_mu_ = new TVector3();
  TVector3 *mc_p3_lead_p_ = new TVector3();
  TVector3 *p3_lead_p_ = new TVector3();
};

void TempEvent::set_branch_addresses( TTree& etree ) {

  etree.SetBranchAddress("delta_pT", &delta_pT_);
  etree.SetBranchAddress("delta_pL", &delta_pL_);
  etree.SetBranchAddress("delta_pL2", &delta_pL2_);
  etree.SetBranchAddress("pn", &pn_);
  etree.SetBranchAddress("pn2", &pn2_);
  etree.SetBranchAddress("theta_mu_p", &theta_mu_p_);
  etree.SetBranchAddress("theta_mu_cpi", &theta_mu_cpi_);
  etree.SetBranchAddress("delta_alpha3D", &delta_alpha3D_);
  etree.SetBranchAddress("delta_phi3d_had", &delta_phi3d_had_);
  etree.SetBranchAddress("delta_phi3d_mu", &delta_phi3d_mu_);

  etree.SetBranchAddress("mc_delta_pT", &mc_delta_pT_);
  etree.SetBranchAddress("mc_delta_pL", &mc_delta_pL_);
  etree.SetBranchAddress("mc_delta_pL2", &mc_delta_pL2_);
  etree.SetBranchAddress("mc_pn", &mc_pn_);
  etree.SetBranchAddress("mc_pn2", &mc_pn2_);
  etree.SetBranchAddress("mc_theta_mu_p", &mc_theta_mu_p_);
  etree.SetBranchAddress("mc_theta_mu_cpi", &mc_theta_mu_cpi_);
  etree.SetBranchAddress("mc_delta_alpha3D", &mc_delta_alpha3D_);
  etree.SetBranchAddress("mc_delta_phi3d_had", &mc_delta_phi3d_had_);
  etree.SetBranchAddress("mc_delta_phi3d_mu", &mc_delta_phi3d_mu_);

  set_object_input_branch_address( etree,"pfp_generation_v", pfp_generation_ );
  set_object_input_branch_address( etree, "backtracked_pdg", pfp_true_pdg_ );
  set_object_input_branch_address( etree, "backtracked_px", pfp_true_px_ );
  set_object_input_branch_address( etree, "backtracked_py", pfp_true_py_ );
  set_object_input_branch_address( etree, "backtracked_pz", pfp_true_pz_ );
  set_object_input_branch_address( etree, "trk_range_muon_mom_v",
    track_range_mom_mu_ );
  set_object_input_branch_address( etree, "trk_score_v", pfp_track_score_ );

  set_object_input_branch_address( etree, "trk_sce_start_x_v", track_startx_ );
  set_object_input_branch_address( etree, "trk_sce_start_y_v", track_starty_ );
  set_object_input_branch_address( etree, "trk_sce_start_z_v", track_startz_ );

  set_object_input_branch_address( etree, "trk_sce_end_x_v", track_endx_ );
  set_object_input_branch_address( etree, "trk_sce_end_y_v", track_endy_ );
  set_object_input_branch_address( etree, "trk_sce_end_z_v", track_endz_ );

  set_object_input_branch_address( etree, "trk_dir_x_v", track_dirx_ );
  set_object_input_branch_address( etree, "trk_dir_y_v", track_diry_ );
  set_object_input_branch_address( etree, "trk_dir_z_v", track_dirz_ );

  set_object_input_branch_address( etree, "trk_len_v", track_length_ );

  set_object_input_branch_address( etree, "mc_pdg", mc_pdg_ );
  set_object_input_branch_address( etree, "mc_n_elastic", mc_n_elastic_ );
  set_object_input_branch_address( etree, "mc_n_inelastic", mc_n_inelastic_ );
  set_object_input_branch_address( etree, "mc_end_p", mc_end_p_ );

 
  etree.SetBranchAddress("pion_candidate_idx", &pion_candidate_idx_);
  etree.SetBranchAddress("muon_candidate_pid_idx", &muon_candidate_pid_idx_); 
  etree.SetBranchAddress("lead_p_candidate_idx", &lead_p_candidate_idx_ );
  etree.SetBranchAddress("mc_is_signal", &mc_is_signal_);

  etree.SetBranchAddress("sel_CCNp1pi", &sel_CCNp1pi_);
  etree.SetBranchAddress("sel_has_pion_candidate", &sel_has_pion_candidate_);
 
  etree.SetBranchAddress("sel_has_muon_candidate", &sel_has_muon_candidate_);
  etree.SetBranchAddress("sel_muon_contained", &sel_muon_contained_);
  etree.SetBranchAddress("sel_muon_quality_ok", &sel_muon_quality_ok_);
  etree.SetBranchAddress("sel_muon_passed_mom_cuts", &sel_muon_passed_mom_cuts_); 
  etree.SetBranchAddress("sel_has_p_candidate", &sel_has_p_candidate_);
  etree.SetBranchAddress("mc_p3_cpi", &mc_p3_cpi_);
  etree.SetBranchAddress("p3_cpi", &p3_cpi_);
  etree.SetBranchAddress("mc_p3_mu", &mc_p3_mu_);
  etree.SetBranchAddress("p3_mu", &p3_mu_);
  etree.SetBranchAddress("mc_p3_lead_p", &mc_p3_lead_p_);
  etree.SetBranchAddress("p3_lead_p", &p3_lead_p_);
}

int get_best_matched_truth_index_pdg(TempEvent& event, int pdg){

  int truth_index = -1;

  int n_mc_truth_particles = event.mc_pdg_->size();
  for (int i = 0; i<n_mc_truth_particles; i++){
    if( std::abs(event.mc_pdg_->at(i)) == pdg) truth_index = i;
  }
  return truth_index;
}


void GKIResStudy_golden(){

  TChain stv_tree( "stv_tree" );
  stv_tree.Add("/exp/uboone/app/users/kwresilo/CCNp1pi/gki_study2.root");
  //stv_tree.Add( "/exp/uboone/data/users/kwresilo/CCNp1pi/data/stv_ntuples_pass8/stv-nu_overlay_run1.roottest.root");
  //stv_tree.Add("/exp/uboone/data/users/kwresilo/CCNp1pi/data/stv_ntuples_pass8/stv-nu_overlay_run2.roottest.root");
  //stv_tree.Add("/exp/uboone/data/users/kwresilo/CCNp1pi/data/stv_ntuples_pass8/stv-nu_overlay_run3.roottest.root");

  //TH2D *h1 = new TH2D("pn resolution", "", 40, 0., 1.0,  40, 0., 1.0);
  /*TH2D *h2 = new TH2D("pn2 resolution", "", 40, 0., 1.0,  40, 0., 1.0); 
  TH2D *h3 = new TH2D("alpha res", "", 40, 0., 3.14,  40,  0., 3.14);
  TH2D *h4 = new TH2D("phi_had res", "", 40, 0., 3.14,  40,  0., 3.1); 
  TH2D *h6 = new TH2D("phi_mu res", "", 40, 0., 3.14,  40,  0., 3.1);
  */
  //TH2D *h7 = new TH2D("pt res", "", 40,  0., 1.0,  40, 0., 1.0); 
  //TH2D *h8 = new TH2D("pl res", "", 40,  -1.0, 0.,  40, -1.0 , 0.);
  //TH2D *h9 = new TH2D("pl2 res", "", 40,  -1.0, 0.,  40, -1.0, 0.);

  //TH1D *h21 = new TH1D("pn fractional error before", "", 41, -1., 1.);
  TH1D *h22 = new TH1D("pn2 fractional error", "", 41, -1., 1.);
  TH1D *h23 = new TH1D("alpha fractional error", "", 41, -1., 1.);
  TH1D *h24 = new TH1D("phi_had fractional error", "", 41, -1., 1.);
  TH1D *h26 = new TH1D("phi_mu fractional error", "", 41, -1., 1.);

  TH1D *h221 = new TH1D("pn2 fractional error stop", "", 41, -1., 1.);
  TH1D *h231 = new TH1D("alpha fractional error stop", "", 41, -1., 1.);
  TH1D *h241 = new TH1D("phi_had fractional error stop", "", 41, -1., 1.);
  TH1D *h261 = new TH1D("phi_mu fractional error stop", "", 41, -1., 1.);

  TH1D *h222 = new TH1D("pn2 fractional error golden", "", 41, -1., 1.);
  TH1D *h232 = new TH1D("alpha fractional error golden", "", 41, -1., 1.);
  TH1D *h242 = new TH1D("phi_had fractional error golden", "", 41, -1., 1.);
  TH1D *h262 = new TH1D("phi_mu fractional error golden", "", 41, -1., 1.);
  //TH1D *h27 = new TH1D("pt fractional error ",  "", 41, -1., 1.);
  //TH1D *h28 = new TH1D("pl fractional error ",  "", 41, -1., 1.);
  //TH1D *h29 = new TH1D("pl2 fractional error ",  "", 41, -1., 1.);
 
  h22->SetStats(0);
  h23->SetStats(0);
  h24->SetStats(0);
  h26->SetStats(0);
  
  h221->SetStats(0);
  h231->SetStats(0);
  h241->SetStats(0);
  h261->SetStats(0);

  h222->SetStats(0);
  h232->SetStats(0);
  h242->SetStats(0);
  h262->SetStats(0);

  long entry = 0;
  while (true) {
    
    TempEvent cur_event;
    cur_event.set_branch_addresses( stv_tree );
    //if ( entry > 200000) break; 
    
    int local_entry = stv_tree.LoadTree( entry );
    if ( local_entry < 0 ) break;
    
    if (entry % 1000 == 0) std::cout << "Event " << entry <<std::endl; 
    stv_tree.GetEntry( entry );
    ++entry;

    size_t num_pfparticles = cur_event.pfp_generation_->size();
    
    // Skip to next event if not signal event that has pion, muon and proton candidate
    if ( !(cur_event.mc_is_signal_ && cur_event.sel_has_pion_candidate_ && cur_event.sel_has_p_candidate_ && cur_event.sel_has_muon_candidate_) ) continue;
    h22->Fill( (cur_event.pn_ - cur_event.mc_pn_) / cur_event.mc_pn_ );

    h23->Fill( (cur_event.delta_alpha3D_ - cur_event.mc_delta_alpha3D_ ) / cur_event.mc_delta_alpha3D_ );

    h24->Fill( (cur_event.delta_phi3d_had_- cur_event.mc_delta_phi3d_had_) / cur_event.mc_delta_phi3d_had_);

    h26->Fill((cur_event.delta_phi3d_mu_ - cur_event.mc_delta_phi3d_mu_) / cur_event.mc_delta_phi3d_mu_ );

    int truth_idx = get_best_matched_truth_index_pdg(cur_event, PION_PDG);
    if (truth_idx != -1){
      bool is_stopping = ( cur_event.mc_end_p_->at(truth_idx) <= std::numeric_limits<float>::epsilon() );
      int n_scatters = cur_event.mc_n_elastic_->at(truth_idx) + cur_event.mc_n_inelastic_->at(truth_idx);
      if (is_stopping) {
        h221->Fill( (cur_event.pn_ - cur_event.mc_pn_) / cur_event.mc_pn_ );
        h231->Fill( (cur_event.delta_alpha3D_ - cur_event.mc_delta_alpha3D_ ) / cur_event.mc_delta_alpha3D_ );
        h241->Fill( (cur_event.delta_phi3d_had_- cur_event.mc_delta_phi3d_had_) / cur_event.mc_delta_phi3d_had_);
        h261->Fill((cur_event.delta_phi3d_mu_ - cur_event.mc_delta_phi3d_mu_) / cur_event.mc_delta_phi3d_mu_ );    
        if (n_scatters == 0){
          h222->Fill( (cur_event.pn_ - cur_event.mc_pn_) / cur_event.mc_pn_ );
          h232->Fill( (cur_event.delta_alpha3D_ - cur_event.mc_delta_alpha3D_ ) / cur_event.mc_delta_alpha3D_ );
          h242->Fill( (cur_event.delta_phi3d_had_- cur_event.mc_delta_phi3d_had_) / cur_event.mc_delta_phi3d_had_);
          h262->Fill((cur_event.delta_phi3d_mu_ - cur_event.mc_delta_phi3d_mu_) / cur_event.mc_delta_phi3d_mu_ );
        }
      }
    }
  }


  h22->GetXaxis()->SetTitleFont( FONT_STYLE );
  h22->GetYaxis()->SetTitleFont( FONT_STYLE );
  h22->GetXaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h22->GetYaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h22->GetXaxis()->SetLabelFont( FONT_STYLE );
  h22->GetYaxis()->SetLabelFont( FONT_STYLE );
  h22->GetXaxis()->CenterTitle();
  h22->GetYaxis()->CenterTitle();
  h22->GetXaxis()->SetTitle("p_{n} Fractional Error");
  h22->GetYaxis()->SetTitle("Events");
  h22->SetLineColor(kAzure+2);

  h23->GetXaxis()->SetTitleFont( FONT_STYLE );
  h23->GetYaxis()->SetTitleFont( FONT_STYLE );
  h23->GetXaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h23->GetYaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h23->GetXaxis()->SetLabelFont( FONT_STYLE );
  h23->GetYaxis()->SetLabelFont( FONT_STYLE );
  h23->GetXaxis()->CenterTitle();
  h23->GetYaxis()->CenterTitle();
  h23->GetXaxis()->SetTitle("#alpha_{3D} Fractional Error");
  h23->GetYaxis()->SetTitle("Events");
  h23->SetLineColor(kAzure+2);

  h24->GetXaxis()->SetTitleFont( FONT_STYLE );
  h24->GetYaxis()->SetTitleFont( FONT_STYLE );
  h24->GetXaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h24->GetYaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h24->GetXaxis()->SetLabelFont( FONT_STYLE );
  h24->GetYaxis()->SetLabelFont( FONT_STYLE );
  h24->GetXaxis()->CenterTitle();
  h24->GetYaxis()->CenterTitle();
  h24->GetXaxis()->SetTitle("#Phi_{3D}^{had} Fractional Error");
  h24->GetYaxis()->SetTitle("Events");
  h24->SetLineColor(kAzure+2);

  h26->GetXaxis()->SetTitleFont( FONT_STYLE );
  h26->GetYaxis()->SetTitleFont( FONT_STYLE );
  h26->GetXaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h26->GetYaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h26->GetXaxis()->SetLabelFont( FONT_STYLE );
  h26->GetYaxis()->SetLabelFont( FONT_STYLE );
  h26->GetXaxis()->CenterTitle();
  h26->GetYaxis()->CenterTitle();
  h26->GetXaxis()->SetTitle("#Phi_{3D}^{#mu} Fractional Error");
  h26->GetYaxis()->SetTitle("Events");
  h26->SetLineColor(kAzure+2);
  
  h221->GetXaxis()->SetTitleFont( FONT_STYLE );
  h221->GetYaxis()->SetTitleFont( FONT_STYLE );
  h221->GetXaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h221->GetYaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h221->GetXaxis()->SetLabelFont( FONT_STYLE );
  h221->GetYaxis()->SetLabelFont( FONT_STYLE );
  h221->GetXaxis()->CenterTitle();
  h221->GetYaxis()->CenterTitle();
  h221->GetXaxis()->SetTitle("p_{n} Fractional Error");
  h221->GetYaxis()->SetTitle("Events");
  h221->SetLineColor(kOrange-3);

  h231->GetXaxis()->SetTitleFont( FONT_STYLE );
  h231->GetYaxis()->SetTitleFont( FONT_STYLE );
  h231->GetXaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h231->GetYaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h231->GetXaxis()->SetLabelFont( FONT_STYLE );
  h231->GetYaxis()->SetLabelFont( FONT_STYLE );
  h231->GetXaxis()->CenterTitle();
  h231->GetYaxis()->CenterTitle();
  h231->GetXaxis()->SetTitle("#alpha_{3D} Fractional Error");
  h231->GetYaxis()->SetTitle("Events");
  h231->SetLineColor(kOrange-3);

  h241->GetXaxis()->SetTitleFont( FONT_STYLE );
  h241->GetYaxis()->SetTitleFont( FONT_STYLE );
  h241->GetXaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h241->GetYaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h241->GetXaxis()->SetLabelFont( FONT_STYLE );
  h241->GetYaxis()->SetLabelFont( FONT_STYLE );
  h241->GetXaxis()->CenterTitle();
  h241->GetYaxis()->CenterTitle();
  h241->GetXaxis()->SetTitle("#Phi_{3D}^{had} Fractional Error");
  h241->GetYaxis()->SetTitle("Events");
  h241->SetLineColor(kOrange-3);

  h261->GetXaxis()->SetTitleFont( FONT_STYLE );
  h261->GetYaxis()->SetTitleFont( FONT_STYLE );
  h261->GetXaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h261->GetYaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h261->GetXaxis()->SetLabelFont( FONT_STYLE );
  h261->GetYaxis()->SetLabelFont( FONT_STYLE );
  h261->GetXaxis()->CenterTitle();
  h261->GetYaxis()->CenterTitle();
  h261->GetXaxis()->SetTitle("#Phi_{3D}^{#mu} Fractional Error");
  h261->GetYaxis()->SetTitle("Events");
  h261->SetLineColor(kOrange-3);  

  h222->GetXaxis()->SetTitleFont( FONT_STYLE );
  h222->GetYaxis()->SetTitleFont( FONT_STYLE );
  h222->GetXaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h222->GetYaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h222->GetXaxis()->SetLabelFont( FONT_STYLE );
  h222->GetYaxis()->SetLabelFont( FONT_STYLE );
  h222->GetXaxis()->CenterTitle();
  h222->GetYaxis()->CenterTitle();
  h222->GetXaxis()->SetTitle("p_{n} Fractional Error");
  h222->GetYaxis()->SetTitle("Events");
  h222->SetLineColor(kGreen+2);

  h232->GetXaxis()->SetTitleFont( FONT_STYLE );
  h232->GetYaxis()->SetTitleFont( FONT_STYLE );
  h232->GetXaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h232->GetYaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h232->GetXaxis()->SetLabelFont( FONT_STYLE );
  h232->GetYaxis()->SetLabelFont( FONT_STYLE );
  h232->GetXaxis()->CenterTitle();
  h232->GetYaxis()->CenterTitle();
  h232->GetXaxis()->SetTitle("#alpha_{3D} Fractional Error");
  h232->GetYaxis()->SetTitle("Events");
  h232->SetLineColor(kGreen+2);

  h242->GetXaxis()->SetTitleFont( FONT_STYLE );
  h242->GetYaxis()->SetTitleFont( FONT_STYLE );
  h242->GetXaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h242->GetYaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h242->GetXaxis()->SetLabelFont( FONT_STYLE );
  h242->GetYaxis()->SetLabelFont( FONT_STYLE );
  h242->GetXaxis()->CenterTitle();
  h242->GetYaxis()->CenterTitle();
  h242->GetXaxis()->SetTitle("#Phi_{3D}^{had} Fractional Error");
  h242->GetYaxis()->SetTitle("Events");
  h242->SetLineColor(kGreen+2);

  h262->GetXaxis()->SetTitleFont( FONT_STYLE );
  h262->GetYaxis()->SetTitleFont( FONT_STYLE );
  h262->GetXaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h262->GetYaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h262->GetXaxis()->SetLabelFont( FONT_STYLE );
  h262->GetYaxis()->SetLabelFont( FONT_STYLE );
  h262->GetXaxis()->CenterTitle();
  h262->GetYaxis()->CenterTitle();
  h262->GetXaxis()->SetTitle("#Phi_{3D}^{#mu} Fractional Error");
  h262->GetYaxis()->SetTitle("Events");
  h262->SetLineColor(kGreen+2); 
 
  TCanvas *c22 = new TCanvas("c22", "c22", 900, 600);
  h22->Draw();
  h221->Draw("same");
  h222->Draw("same");
  TLegend *l1 = new TLegend(0.6, 0.7, 0.88, 0.88);
  l1->AddEntry(h22, "Contained Pions");
  l1->AddEntry(h221, "AND Stopping");
  l1->AddEntry(h222, "AND No scatters");
  l1->SetTextSize(0.04);
  l1->SetTextFont(FONT_STYLE);
  l1->SetBorderSize( 0 );
  l1->Draw();
  c22->SaveAs("/exp/uboone/app/users/kwresilo/CCNp1pi/plots/stv_plots/GKIRes/pn_frac_err_golden1.pdf");


  TCanvas *c23 = new TCanvas("c23", "c23", 900, 600);
  h23->Draw();
  h231->Draw("same");
  h232->Draw("same");
  TLegend *l2 = new TLegend(0.6, 0.7, 0.88, 0.88);
  l2->AddEntry(h23, "Contained Pions");
  l2->AddEntry(h231, "AND Stopping");
  l2->AddEntry(h232, "AND No scatters");
  l2->SetTextSize(0.04);
  l2->SetTextFont(FONT_STYLE);
  l2->SetBorderSize( 0 );
  l2->Draw();
  c23->SaveAs("/exp/uboone/app/users/kwresilo/CCNp1pi/plots/stv_plots/GKIRes/alpha_frac_err_golden1.pdf");

  TCanvas *c24 = new TCanvas("c24", "c24", 900, 600);
  h24->Draw();
  h241->Draw("same");
  h242->Draw("same");
  TLegend *l3 = new TLegend(0.6, 0.7, 0.88, 0.88);
  l3->AddEntry(h24, "Contained Pions");
  l3->AddEntry(h241, "AND Stopping");
  l3->AddEntry(h242, "AND No scatters");
  l3->SetTextSize(0.04);
  l3->SetTextFont(FONT_STYLE);
  l3->SetBorderSize( 0 );
  l3->Draw();
  c24->SaveAs("/exp/uboone/app/users/kwresilo/CCNp1pi/plots/stv_plots/GKIRes/phi_had_frac_err_golden1.pdf");
 

  TCanvas *c26 = new TCanvas("c26", "c26", 900, 600);
  h26->Draw();
  h261->Draw("same");
  h262->Draw("same");
  TLegend *l4 = new TLegend(0.6, 0.7, 0.88, 0.88);
  l4->AddEntry(h26, "Contained Pions");
  l4->AddEntry(h261, "AND Stopping");
  l4->AddEntry(h262, "AND No scatters");
  l4->SetTextSize(0.04);
  l4->SetTextFont(FONT_STYLE);
  l4->SetBorderSize( 0 );
  l4->Draw();
  c26->SaveAs("/exp/uboone/app/users/kwresilo/CCNp1pi/plots/stv_plots/GKIRes/phi_mu_frac_err_golden1.pdf");

}




