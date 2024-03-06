#include <string>
#include <fstream>
#include <iomanip>
#include <vector>

std::vector< std::vector<std::string> > ReadConfig(std::string ConfigFileName) {
  std::vector< std::vector<std::string> > SelectionCuts;

  ifstream input(ConfigFileName);
  if (!input.is_open()) {
    std::cout << "Did not find Config File named:" << ConfigFileName << std::endl;
    throw;
  }

  std::string tmp;
  for(std::string line; getline( input, line );) {
    std::vector<std::string> Selection;

    stringstream ss(line);
    while (getline(ss, tmp, ' ')) {
      Selection.push_back(tmp);
    }

    if (Selection.size()!=0) {
      SelectionCuts.push_back(Selection);
    }
  }

  input.close();

  return SelectionCuts;
}

void CalcEfficiency(TTree* Tree, std::vector<std::string> SelectionCuts) {
  std::vector<TH1D*> Hists(SelectionCuts.size());
  std::vector<float> Efficiency(SelectionCuts.size(),0);
  
  int HistCounter = 0;
  TString HistName;
  TString Cut = SelectionCuts[0]+"_MC_Signal==1";

  HistName = Form("Hist_Eff_%s_%i",SelectionCuts[0].c_str(),HistCounter);
  HistCounter += 1;

  Tree->Draw(TString(SelectionCuts[0])+"_Selected>>"+HistName,Cut,"goff");
  Hists[0] = (TH1D*)gDirectory->Get(HistName);

  for (int i=1;i<SelectionCuts.size();i++) {
    Cut += " && "+SelectionCuts[i]+"==1";

    HistName = Form("Hist_%s_%i",SelectionCuts[0].c_str(),HistCounter);
    HistCounter += 1;

    Tree->Draw(TString(SelectionCuts[0])+"_Selected>>"+HistName,Cut,"goff");
    Hists[i] = (TH1D*)gDirectory->Get(HistName);
  }

  for (int i=0;i<Hists.size();i++) {
    Efficiency[i] = 100.0*float(Hists[i]->Integral())/float(Hists[0]->Integral());
  }
  
  for (int i=0;i<Hists.size();i++) {
    if (i==0) {
      std::cout << "Sample: " << std::setw(30) << SelectionCuts[0] << " | Cut: " << std::setw(60) << "No Cut" << " | #Signal: " << std::setw(10) << Hists[i]->Integral() << " | Efficiency (%):" << std::setw(10) << Efficiency[i] << std::endl;
    } else {
      std::cout << "Sample: " << std::setw(30) << SelectionCuts[0] << " | Cut: " << std::setw(60) << SelectionCuts[i] << " | #Signal: " << std::setw(10) << Hists[i]->Integral() << " | Efficiency (%):" << std::setw(10) << Efficiency[i] << std::endl;
    }
  }

  std::cout << "\n\n" << std::endl;
}

void CalcPurity(TTree* Tree, std::vector<std::string> SelectionCuts) {
  std::vector<TH1D*> Hists_Selected(SelectionCuts.size());
  std::vector<TH1D*> Hists_SelectedSignal(SelectionCuts.size());
  std::vector<float> Purity(SelectionCuts.size());
  
  int HistCounter = 0;
  TString HistName;

  TString Cut_Selected, Cut_SelectedSignal;

  Cut_Selected = "";
  HistName = Form("Hist_Pur_Sel_%s_%i",SelectionCuts[0].c_str(),HistCounter);
  Tree->Draw(TString(SelectionCuts[0])+"_Selected>>"+HistName,Cut_Selected,"goff");
  Hists_Selected[0] = (TH1D*)gDirectory->Get(HistName);

  Cut_SelectedSignal = TString(SelectionCuts[0])+"_MC_Signal==1";
  HistName = Form("Hist_Pur_SelSig_%s_%i",SelectionCuts[0].c_str(),HistCounter);
  Tree->Draw(TString(SelectionCuts[0])+"_Selected>>"+HistName,Cut_SelectedSignal,"goff");
  Hists_SelectedSignal[0] = (TH1D*)gDirectory->Get(HistName);

  Purity[0] = 100.0*float(Hists_SelectedSignal[0]->Integral())/float(Hists_Selected[0]->Integral());
  HistCounter += 1;

  for (int i=1;i<SelectionCuts.size();i++) {
    if (i == 1) {
      Cut_Selected += SelectionCuts[i]+"==1";
    } else {
      Cut_Selected += " && "+SelectionCuts[i]+"==1";
    }
    HistName = Form("Hist_Pur_Sel_%s_%i",SelectionCuts[i].c_str(),HistCounter);
    Tree->Draw(TString(SelectionCuts[0])+"_Selected>>"+HistName,Cut_Selected,"goff");
    Hists_Selected[i] = (TH1D*)gDirectory->Get(HistName);
    
    Cut_SelectedSignal += " && "+SelectionCuts[i]+"==1";
    HistName = Form("Hist_Pur_SelSig_%s_%i",SelectionCuts[i].c_str(),HistCounter);
    Tree->Draw(TString(SelectionCuts[0])+"_Selected>>"+HistName,Cut_SelectedSignal,"goff");
    Hists_SelectedSignal[i] = (TH1D*)gDirectory->Get(HistName);
    
    Purity[i] = 100.0*float(Hists_SelectedSignal[i]->Integral())/float(Hists_Selected[i]->Integral());
    HistCounter += 1;
  }

  for (int i=0;i<Hists_Selected.size();i++) {
    if (i == 0) {
      std::cout << "Sample: " << std::setw(30) << SelectionCuts[0] << " | Cut: " << std::setw(60) << "No Cut" << " | #Selected: " << std::setw(15) << Hists_Selected[i]->Integral() << " | #(Selected && Signal): " << std::setw(15) << Hists_SelectedSignal[i]->Integral() << " | Purity (%):" << std::setw(10) << Purity[i] << std::endl;
    } else {
      std::cout << "Sample: " << std::setw(30) << SelectionCuts[0] << " | Cut: " << std::setw(60) << SelectionCuts[i] << " | #Selected: " << std::setw(15) << Hists_Selected[i]->Integral() << " | #(Selected && Signal): " << std::setw(15) << Hists_SelectedSignal[i]->Integral() << " | Purity (%):" << std::setw(10) << Purity[i] << std::endl;
    }
  }

  std::cout << "\n\n" << std::endl;
}

void PlotEfficiency() {
  std::string ConfigFileName = "Selections.txt";
  std::vector< std::vector<std::string> > SelectionCuts = ReadConfig(ConfigFileName);

  TString FileName = "/Users/barrowd/Desktop/CorrelatedXsecAnalysis/OutputFiles/prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2__STV-ANALYSIS-PROCESSED__.root";
  TFile* File = new TFile(FileName);
  if (!File) {std::cout << "File:" << FileName << " not found!" << std::endl; throw;}
  TTree* Tree = (TTree*)File->Get("stv_tree");
  if (!Tree) {std::cout << "Did not find stv_tree in File:" << FileName << std::endl; throw;}

  int nSelections = SelectionCuts.size();

  for (int i=0;i<nSelections;i++) {
    std::cout << "SelectionName:" << SelectionCuts[i][0] << std::endl;
    for (int j=1;j<SelectionCuts[i].size();j++) {
      std::cout << "\t Cut " << j << " = " << SelectionCuts[i][j] << std::endl;
    }
  }
  std::cout << "\n\n" << std::endl;

  for (int i=0;i<nSelections;i++) {
    CalcEfficiency(Tree,SelectionCuts[i]);
    CalcPurity(Tree,SelectionCuts[i]);
  }
}
