#include <vector>
#include <iostream>
#include <iomanip>

TTree* GetTree(TString FileName, TString TreeName) {
  TFile* File = new TFile(FileName);
  if (!File) {
    std::cout << "Did not find file:" << FileName << std::endl;
    throw;
  }

  TTree* Tree = (TTree*)File->Get(TreeName);
  if (!Tree) {
    std::cout << "Did not find TTree:" << TreeName << " in File:" << FileName << std::endl;
    throw;
  }

  return Tree;
}

void TreeComparer() {
  TString NewTreeName = "stv_tree";
  TString NewTreeShortName = "New Code";
  TString NewFileName = "/uboone/data/users/barrow/XSecAnalyzer_Validation/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";
  
  TString OldTreeName = "stv_tree";
  TString OldTreeShortName = "Old Code";
  TString OldFileName = "/uboone/data/users/barrow/XSecAnalyzer_Orig/stv-prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root";

  TTree* NewTree = GetTree(NewFileName,NewTreeName);
  TTree* OldTree = GetTree(OldFileName,OldTreeName);

  TObjArray* NewTreeBranches = (TObjArray*)NewTree->GetListOfBranches();
  int NewTreeNBranches = NewTreeBranches->GetEntries();
  std::vector<bool> NewTreeBranchFound(NewTreeNBranches,false);
  std::vector<TString> NewTreeBranchMatch(NewTreeNBranches,"");

  TObjArray* OldTreeBranches = (TObjArray*)OldTree->GetListOfBranches();
  int OldTreeNBranches = OldTreeBranches->GetEntries();
  std::vector<bool> OldTreeBranchFound(OldTreeNBranches,false);
  std::vector<TString> OldTreeBranchMatch(OldTreeNBranches,"");

  std::vector<TString> NewBranchesToSkipDrawing;
  NewBranchesToSkipDrawing.push_back("weight_All_UBGenie");
  NewBranchesToSkipDrawing.push_back("weight_flux_all");
  NewBranchesToSkipDrawing.push_back("weight_reint_all");
  
  std::cout << "Number of branches in NewTree:" << NewTreeNBranches << std::endl;
  std::cout << "Number of branches in OldTree:" << OldTreeNBranches << std::endl;

  //====================================================================================================================================
  //First check for identical matching Branch Names
  
  //NewTree
  
  for (int iBr=0;iBr<NewTreeNBranches;iBr++) {
    TString BranchToFind = NewTreeBranches->At(iBr)->GetName();
    for (int jBr=0;jBr<OldTreeNBranches;jBr++) {
      TString BranchName = OldTreeBranches->At(jBr)->GetName();

      if (BranchToFind == BranchName) {
	NewTreeBranchFound[iBr] = true;
	NewTreeBranchMatch[iBr] = BranchName;
	break;
      }
    }
  }

  //OldTree
  
  for (int iBr=0;iBr<OldTreeNBranches;iBr++) {
    TString BranchToFind = OldTreeBranches->At(iBr)->GetName();
    for (int jBr=0;jBr<NewTreeNBranches;jBr++) {
      TString BranchName = NewTreeBranches->At(jBr)->GetName();
      
      if (BranchToFind == BranchName) {
        OldTreeBranchFound[iBr] = true;
	OldTreeBranchMatch[iBr]	= BranchName;
        break;
      }
    }
  }

  //====================================================================================================================================

  std::vector<TString> Map_NewBranchNames;
  std::vector<TString> Map_OldBranchNames;

  Map_NewBranchNames.push_back("CC1muNp0pi_Selected");
  Map_OldBranchNames.push_back("sel_CCNp0pi");

  Map_NewBranchNames.push_back("CC1muNp0pi_MC_Signal");
  Map_OldBranchNames.push_back("mc_is_signal");

  Map_NewBranchNames.push_back("CC1muNp0pi_EventCategory");
  Map_OldBranchNames.push_back("category");

  Map_NewBranchNames.push_back("CC1muNp0pi_reco_vertex_in_FV");
  Map_OldBranchNames.push_back("sel_reco_vertex_in_FV");

  Map_NewBranchNames.push_back("CC1muNp0pi_pfp_starts_in_PCV");
  Map_OldBranchNames.push_back("sel_pfp_starts_in_PCV");

  Map_NewBranchNames.push_back("CC1muNp0pi_has_muon_candidate");
  Map_OldBranchNames.push_back("sel_has_muon_candidate");

  Map_NewBranchNames.push_back("CC1muNp0pi_topo_cut_passed");
  Map_OldBranchNames.push_back("sel_topo_cut_passed");

  Map_NewBranchNames.push_back("CC1muNp0pi_nu_mu_cc");
  Map_OldBranchNames.push_back("sel_nu_mu_cc");

  Map_NewBranchNames.push_back("CC1muNp0pi_muon_contained");
  Map_OldBranchNames.push_back("sel_muon_contained");

  Map_NewBranchNames.push_back("CC1muNp0pi_muon_passed_mom_cuts");
  Map_OldBranchNames.push_back("sel_muon_passed_mom_cuts");

  Map_NewBranchNames.push_back("CC1muNp0pi_no_reco_showers");
  Map_OldBranchNames.push_back("sel_no_reco_showers");

  Map_NewBranchNames.push_back("CC1muNp0pi_has_p_candidate");
  Map_OldBranchNames.push_back("sel_has_p_candidate");

  Map_NewBranchNames.push_back("CC1muNp0pi_muon_quality_ok");
  Map_OldBranchNames.push_back("sel_muon_quality_ok");

  Map_NewBranchNames.push_back("CC1muNp0pi_protons_contained");
  Map_OldBranchNames.push_back("sel_protons_contained");

  Map_NewBranchNames.push_back("CC1muNp0pi_passed_proton_pid_cut");
  Map_OldBranchNames.push_back("sel_passed_proton_pid_cut");

  Map_NewBranchNames.push_back("CC1muNp0pi_lead_p_passed_mom_cuts");
  Map_OldBranchNames.push_back("sel_lead_p_passed_mom_cuts");

  Map_NewBranchNames.push_back("CC1muNp0pi_lead_p_candidate_idx");
  Map_OldBranchNames.push_back("lead_p_candidate_idx");

  Map_NewBranchNames.push_back("CC1muNp0pi_muon_candidate_idx");
  Map_OldBranchNames.push_back("muon_candidate_idx");

  Map_NewBranchNames.push_back("CC1muNp0pi_mc_is_numu");
  Map_OldBranchNames.push_back("mc_neutrino_is_numu");

  Map_NewBranchNames.push_back("CC1muNp0pi_mc_vertex_in_FV");
  Map_OldBranchNames.push_back("mc_vertex_in_FV");

  Map_NewBranchNames.push_back("CC1muNp0pi_mc_lead_p_in_range");
  Map_OldBranchNames.push_back("mc_lead_p_in_mom_range");

  Map_NewBranchNames.push_back("CC1muNp0pi_mc_muon_in_mom_range");
  Map_OldBranchNames.push_back("mc_muon_in_mom_range");

  Map_NewBranchNames.push_back("CC1muNp0pi_mc_no_pi0s");
  Map_OldBranchNames.push_back("mc_no_fs_pi0");

  Map_NewBranchNames.push_back("CC1muNp0pi_mc_no_FS_mesons");
  Map_OldBranchNames.push_back("mc_no_fs_mesons");

  Map_NewBranchNames.push_back("CC1muNp0pi_mc_no_charged_pions_above_thres");
  Map_OldBranchNames.push_back("mc_no_charged_pi_above_threshold");

  Map_NewBranchNames.push_back("CC1muNp0pi_cosmic_ip_cut_passed");
  Map_OldBranchNames.push_back("sel_cosmic_ip_cut_passed");

  Map_NewBranchNames.push_back("CC1muNp0pi_reco_delta_pT");
  Map_OldBranchNames.push_back("delta_pT");

  Map_NewBranchNames.push_back("CC1muNp0pi_reco_delta_phiT");
  Map_OldBranchNames.push_back("delta_phiT");

  Map_NewBranchNames.push_back("CC1muNp0pi_reco_delta_pTx");
  Map_OldBranchNames.push_back("delta_pTx");

  Map_NewBranchNames.push_back("CC1muNp0pi_reco_delta_pTy");
  Map_OldBranchNames.push_back("delta_pTy");

  Map_NewBranchNames.push_back("CC1muNp0pi_reco_theta_mu_p");
  Map_OldBranchNames.push_back("theta_mu_p");

  Map_NewBranchNames.push_back("CC1muNp0pi_reco_pn");
  Map_OldBranchNames.push_back("pn");

  Map_NewBranchNames.push_back("CC1muNp0pi_reco_delta_pL");
  Map_OldBranchNames.push_back("delta_pL");

  Map_NewBranchNames.push_back("CC1muNp0pi_reco_delta_alphaT");
  Map_OldBranchNames.push_back("delta_alphaT");

  Map_NewBranchNames.push_back("CC1muNp0pi_reco_p3_lead_p");
  Map_OldBranchNames.push_back("p3_lead_p");

  Map_NewBranchNames.push_back("CC1muNp0pi_reco_p3_mu");
  Map_OldBranchNames.push_back("p3_mu");

  Map_NewBranchNames.push_back("CC1muNp0pi_reco_p3_p_vec");
  Map_OldBranchNames.push_back("p3_p_vec");

  Map_NewBranchNames.push_back("CC1muNp0pi_true_delta_pT");
  Map_OldBranchNames.push_back("mc_delta_pT");

  Map_NewBranchNames.push_back("CC1muNp0pi_true_delta_phiT");
  Map_OldBranchNames.push_back("mc_delta_phiT");

  Map_NewBranchNames.push_back("CC1muNp0pi_true_delta_pTx");
  Map_OldBranchNames.push_back("mc_delta_pTx");

  Map_NewBranchNames.push_back("CC1muNp0pi_true_delta_pTy");
  Map_OldBranchNames.push_back("mc_delta_pTy");

  Map_NewBranchNames.push_back("CC1muNp0pi_true_theta_mu_p");
  Map_OldBranchNames.push_back("mc_theta_mu_p");

  Map_NewBranchNames.push_back("CC1muNp0pi_true_pn");
  Map_OldBranchNames.push_back("mc_pn");

  Map_NewBranchNames.push_back("CC1muNp0pi_true_delta_pL");
  Map_OldBranchNames.push_back("mc_delta_pL");

  Map_NewBranchNames.push_back("CC1muNp0pi_true_delta_alphaT");
  Map_OldBranchNames.push_back("mc_delta_alphaT");

  Map_NewBranchNames.push_back("CC1muNp0pi_true_p3_lead_p");
  Map_OldBranchNames.push_back("mc_p3_lead_p");

  Map_NewBranchNames.push_back("CC1muNp0pi_true_p3_mu");
  Map_OldBranchNames.push_back("mc_p3_mu");

  Map_NewBranchNames.push_back("CC1muNp0pi_true_p3_p_vec");
  Map_OldBranchNames.push_back("mc_p3_p_vec");
  
  for (int iBr=0;iBr<Map_NewBranchNames.size();iBr++) {
    TString BranchOfInterest = Map_NewBranchNames[iBr];
    TString BranchToFind = Map_OldBranchNames[iBr];
    
    for (int jBr=0;jBr<NewTreeNBranches;jBr++) {
      TString NewBranchName = NewTreeBranches->At(jBr)->GetName();
      if (BranchOfInterest == NewBranchName) {

	for (int kBr=0;kBr<OldTreeNBranches;kBr++) {
	  TString OldBranchName = OldTreeBranches->At(kBr)->GetName();
	  if (OldBranchName == BranchToFind) {
	    NewTreeBranchFound[jBr] = true;
	    NewTreeBranchMatch[jBr] = BranchToFind;
	    break;
	  }
	}
	break;
      }
    }    
  }

  for (int iBr=0;iBr<Map_OldBranchNames.size();iBr++) {
    TString BranchOfInterest = Map_OldBranchNames[iBr];
    TString BranchToFind = Map_NewBranchNames[iBr];

    for (int jBr=0;jBr<OldTreeNBranches;jBr++) {
      TString OldBranchName = OldTreeBranches->At(jBr)->GetName();
      if (BranchOfInterest == OldBranchName) {

        for (int kBr=0;kBr<NewTreeNBranches;kBr++) {
          TString NewBranchName = NewTreeBranches->At(kBr)->GetName();
          if (NewBranchName == BranchToFind) {
            OldTreeBranchFound[jBr] = true;
            OldTreeBranchMatch[jBr] = BranchToFind;
            break;
          }
        }
        break;
      }
    }
  }
  
  //====================================================================================================================================

  std::cout << "\n\nNewTreeBranches:" << std::endl;
  for (int iBr=0;iBr<NewTreeNBranches;iBr++) {
    TString BranchName = NewTreeBranches->At(iBr)->GetName();
    std::cout << "NewTree | BranchName: " << std::setw(50) << BranchName << " | Found: " << NewTreeBranchFound[iBr] << " (" << std::setw(50) << NewTreeBranchMatch[iBr] << ")" << std::endl; 
  }
  std::cout << "\n\nOldTreeBranches:" << std::endl;
  for (int iBr=0;iBr<OldTreeNBranches;iBr++) {
    TString BranchName = OldTreeBranches->At(iBr)->GetName();
    std::cout << "OldTree | BranchName: " << std::setw(50) << BranchName << " | Found: " << OldTreeBranchFound[iBr] << " (" << std::setw(50) << OldTreeBranchMatch[iBr] << ")" << std::endl; 
  }

  //====================================================================================================================================

  gErrorIgnoreLevel = kFatal;
  
  TH1D* NewTreeHist;
  TH1D* OldTreeHist;

  TCanvas* Canv = new TCanvas("Canv","");
  TString OutputName = "Comparison.pdf";
  Canv->Print(OutputName+"[");

  bool IsFirstPrint = true;
  
  for (int iBr=0;iBr<NewTreeNBranches;iBr++) {
    TString NewBranchName = NewTreeBranches->At(iBr)->GetName();
    TString OldBranchName = NewTreeBranchMatch[iBr];
    
    bool Skip = false;
    for (int kBr=0;kBr<NewBranchesToSkipDrawing.size();kBr++) {
      if (NewBranchName == NewBranchesToSkipDrawing[kBr]) {
	Skip = true;
	break;
      }
    }
    if (Skip) {continue;}

    if (NewTreeBranchFound[iBr]) {
      TString NewHistName = Form("hNew_%i",iBr);
      TString OldHistName = Form("hOld_%i",iBr);

      if (NewBranchName.Contains("p3")) {
	NewBranchName += ".fX";
	OldBranchName += ".fX";
      }
      
      //std::cout << "Drawing: " << std::setw(50) << NewBranchName << " | " << std::setw(50) << OldBranchName << std::endl;
      TString Cut;

      Cut = "( "+NewBranchName+" != 9999) && ( "+NewBranchName+" != 9999.)";
      NewTree->Draw(NewBranchName+">>"+NewHistName,Cut,"goff");
      
      //std::cout << Cut << std::endl;
      Cut = "( "+OldBranchName+" != 9999) && ( "+OldBranchName+" != 9999.)";
      OldTree->Draw(OldBranchName+">>"+OldHistName,Cut,"goff");

      NewTreeHist = (TH1D*)gDirectory->Get(NewHistName);
      NewTreeHist->SetLineColor(kRed);
      
      OldTreeHist = (TH1D*)gDirectory->Get(OldHistName);
      OldTreeHist->SetLineColor(kBlack);
      OldTreeHist->SetLineStyle(kDashed);

      NewTreeHist->SetStats(false);
      OldTreeHist->SetStats(false);
      NewTreeHist->GetXaxis()->SetTitle(NewBranchName);
      OldTreeHist->GetXaxis()->SetTitle(NewBranchName);

      if (IsFirstPrint) {
	Canv->cd();
	TLegend* Text = new TLegend(0.1,0.1,0.9,0.9);
	Text->AddEntry(NewTreeHist,NewTreeShortName,"l");
	Text->AddEntry(OldTreeHist,OldTreeShortName,"l");
	Text->Draw();
	Canv->Print(OutputName);
	IsFirstPrint = false;
      }
      
      NewTreeHist->Draw();
      OldTreeHist->Draw("SAME");
      Canv->Print(OutputName);

      int nBins_New = NewTreeHist->GetNbinsX();
      int nBins_Old = OldTreeHist->GetNbinsX();
      if (nBins_New != nBins_Old) {
	std::cout << "Different number of bins - " << std::setw(50) << NewBranchName << "(" << nBins_New << ")" << " | " << std::setw(50) << OldBranchName << "(" << nBins_Old << ")" << std::endl;
      } else {
	bool IsIdentical = true;
	for (int iBin=1;iBin<nBins_New+1;iBin++) {
	  if (NewTreeHist->GetBinContent(iBin) != OldTreeHist->GetBinContent(iBin)) {
	    IsIdentical = false;
	    break;
	  }
	}
	if (!IsIdentical) {
	  std::cout << "Different bin contents   - " << std::setw(50) << NewBranchName << " | " << std::setw(50) << OldBranchName << std::endl;

	  for (int iBin=1;iBin<nBins_New+1;iBin++) {
	    if (NewTreeHist->GetBinContent(iBin) != OldTreeHist->GetBinContent(iBin)) {
	      std::cout << "\t" << iBin << " " << NewTreeHist->GetBinContent(iBin) << " " << OldTreeHist->GetBinContent(iBin) << std::endl;
	    }
	  }

	}
      }
    }
  }

  Canv->Print(OutputName+"]");
  
}
