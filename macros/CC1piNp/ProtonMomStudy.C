#include "TVector3.h"

#include "/exp/uboone/app/users/kwresilo/CCNp1pi/FiducialVolume.hh"
#include "/exp/uboone/app/users/kwresilo/CCNp1pi/TreeUtils.hh"


constexpr int MUON_PDG = 13;
constexpr int PION_PDG = 211;
constexpr int PROTON_PDG = 2212;
constexpr float PROTON_MASS = 0.93827208;

const int FONT_STYLE = 132;
const double TEXT_SIZE  = 0.06;

struct TempEvent {

  TempEvent() {}

  void set_branch_addresses( TTree& etree );

  MyPointer< std::vector<int> > pfp_generation_;
  MyPointer< std::vector<int> > pfp_true_pdg_;
  MyPointer< std::vector<float> > track_kinetic_energy_p_; 
  
  Bool_t mc_is_signal_;
  Bool_t sel_nu_mu_cc_;
  Bool_t mc_vertex_in_FV_;
  Bool_t mc_neutrino_is_numu_;
  Bool_t mc_muon_in_mom_range_;
  Bool_t mc_no_fs_pi0_;
  Bool_t mc_1cpi_;
  Bool_t mc_no_fs_mesons_;
  
  TVector3 *mc_p3_lead_p_ = new TVector3();
};


void TempEvent::set_branch_addresses( TTree& etree ) {

  set_object_input_branch_address( etree,"pfp_generation_v", pfp_generation_ );
  set_object_input_branch_address( etree, "backtracked_pdg", pfp_true_pdg_ );
  set_object_input_branch_address( etree, "trk_energy_proton_v", track_kinetic_energy_p_);

  etree.SetBranchAddress("MC_Signal", &mc_is_signal_);
  etree.SetBranchAddress("sel_nu_mu_cc", &sel_nu_mu_cc_);
  etree.SetBranchAddress("mc_vertex_in_FV", &mc_vertex_in_FV_);
  etree.SetBranchAddress("mc_neutrino_is_numu", &mc_neutrino_is_numu_);
  etree.SetBranchAddress("mc_muon_in_mom_range", &mc_muon_in_mom_range_);
  etree.SetBranchAddress("mc_no_fs_pi0", &mc_no_fs_pi0_);
  etree.SetBranchAddress("mc_1cpi", &mc_1cpi_);
  etree.SetBranchAddress("mc_no_fs_mesons", &mc_no_fs_mesons_); 
  etree.SetBranchAddress("mc_p3_lead_p", &mc_p3_lead_p_);
}


void ProtonMomRes(){

  TChain stv_tree( "stv_tree" );
  stv_tree.Add( "/exp/uboone/data/users/kwresilo/CCNp1pi/data/stv_ntuples_pass6/stv-nu_overlay_run1.roottest.root");
  stv_tree.Add("/exp/uboone/data/users/kwresilo/CCNp1pi/data/stv_ntuples_pass6/stv-nu_overlay_run2.roottest.root");
  stv_tree.Add("/exp/uboone/data/users/kwresilo/CCNp1pi/data/stv_ntuples_pass6/stv-nu_overlay_run3.roottest.root");

  TH1D *h1 = new TH1D("proton fractional error", "", 41, -1., 1.);
  TH1D *h2 = new TH1D("proton fractional error2", "", 41, -1., 1.);
  TH1D *h3 = new TH1D("proton fractional error3", "", 41, -1, 1.);
  TH1D *h4 = new TH1D("proton fractional error3", "", 41, -1, 1.);
  h1->SetStats(0);
  h2->SetStats(0);
  h3->SetStats(0);
  h4->SetStats(0);

  long entry = 0;
  while (true) {

    TempEvent cur_event;
    cur_event.set_branch_addresses( stv_tree );

    int local_entry = stv_tree.LoadTree( entry );
    if ( local_entry < 0 ) break;

    //if (entry > 200000) break;

    if (entry % 1000 == 0) std::cout << "Event " << entry <<std::endl;
    stv_tree.GetEntry( entry );
    ++entry;

    size_t num_pfparticles = cur_event.pfp_generation_->size();
    //std::cout << cur_event.pfp_generation_->size() << std::endl;
    //if (cur_event.mc_is_signal_) std::cout << "signal event found" << std::endl;
    //if (cur_event.sel_nu_mu_cc_ == 1) std::cout << "cc sel event found"  << std::endl;
 
    if( !(cur_event.mc_is_signal) ) continue;

    //std::cout << "signal event found" << std::endl;   
    for (size_t p = 0u; p < num_pfparticles; ++p ) {
      if ( cur_event.pfp_generation_->at( p ) != 2) continue;

      if ( std::abs( cur_event.pfp_true_pdg_->at( p ) )  == PROTON_PDG) {
	//std::cout << "proton found" << std::endl;
        float mc_mom = cur_event.mc_p3_lead_p_->Mag();
        float reco_ke = cur_event.track_kinetic_energy_p_->at( p );	
	      float reco_mom = cur_event.p3_lead_p_->Mag();
        if (mc_mom >= 0.2 && mc_mom < 0.25) h1->Fill( (reco_mom - mc_mom) / mc_mom );
        if (mc_mom >= 0.25 && mc_mom < 0.3) h2->Fill( (reco_mom - mc_mom) / mc_mom );
        if (mc_mom >= 0.3 && mc_mom < 0.35) h3->Fill( (reco_mom - mc_mom) / mc_mom );
        if (mc_mom >= 0.35 && mc_mom < 0.4) h4->Fill( (reco_mom - mc_mom) / mc_mom );
      } //proton loop
    } //pfparticle
  } //event


  h1->GetXaxis()->SetTitleFont( FONT_STYLE );
  h1->GetYaxis()->SetTitleFont( FONT_STYLE );
  h1->GetXaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h1->GetYaxis()->SetTitleSize( TEXT_SIZE - 0.015);
  h1->GetXaxis()->SetLabelFont( FONT_STYLE );
  h1->GetYaxis()->SetLabelFont( FONT_STYLE );
  h1->GetZaxis()->SetLabelFont( FONT_STYLE );
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->SetTitle("Momentum Fractional Error");
  h1->GetYaxis()->SetTitle("# Protons");

  h1->SetLineColor(kBlack);
  h2->SetLineColor(kRed);
  h3->SetLineColor(kBlue);
  h4->SetLineColor(kOrange);
  h1->Scale(1. / h1->Integral() );
  h2->Scale(1. / h2->Integral() );
  h3->Scale(1. / h3->Integral() );
  h4->Scale(1. / h4->Integral() );
  h1->Sumw2(kFALSE);
  h2->Sumw2(kFALSE);
  h3->Sumw2(kFALSE);
  h4->Sumw2(kFALSE);

  float max1 = h1->GetBinContent( h1->GetMaximumBin() );
  h1->SetMaximum(max1 + max1*0.2);
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1100, 800);
  c1->SetLeftMargin(0.15);

  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  h4->Draw("same");

  TLegend *l1 = new TLegend(0.7, 0.7, 0.9, 0.9);
  l1->AddEntry(h1, "p_{p, true} #in [0.20, 0.25] GeV");
  l1->AddEntry(h2, "p_{p, true} #in [0.25, 0.30] GeV");
  l1->AddEntry(h3, "p_{p, true} #in [0.30, 0.35] GeV");
  l1->AddEntry(h4, "p_{p, true} #in [0.35, 0.40] GeV");
  l1->Draw("same");
  c1->SaveAs("/exp/uboone/app/users/kwresilo/xsec_analyzer/plots/proton_study/proton_low_momentum_study.pdf");

}


