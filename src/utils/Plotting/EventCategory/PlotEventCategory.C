#include <iostream>
#include <iomanip>
#include <vector>

#include "TTree.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TString.h"
#include "TCanvas.h"

#include "../../Utils/EventCategory.hh"

void PlotEventCategory() {
  TString FileName = "/Users/barrowd/Desktop/CorrelatedXsecAnalysis/OutputFiles/prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2__STV-ANALYSIS-PROCESSED__.root";
  TString TreeName = "stv_tree";

  TFile* File = new TFile(FileName);
  if (!File) {std::cerr << "Did not find a File named:" << FileName << std::endl;}
  TTree* Tree = (TTree*)File->Get(TreeName);
  if (!Tree) {std::cerr << "Did not find a TTree named:" << TreeName << " in File named:" << FileName << std::endl;}

  TObjArray* Branches = Tree->GetListOfBranches();
  int nBranches = Branches->GetEntries();

  std::vector<TString> BranchesOfInterest;

  for (int iBranch=0;iBranch<nBranches;iBranch++) {
    TString BranchName = ((TBranch*)(Branches->At(iBranch)))->GetName();
    if (BranchName.Contains("_EventCategory")) {
      BranchesOfInterest.push_back(BranchName);
    }
  }

  int nBranchesOfInterest = BranchesOfInterest.size();
  std::cout << "Branches which will be drawn:" << std::endl;
  for (int iBranch=0;iBranch<nBranchesOfInterest;iBranch++) {
    std::cout << "\t" << BranchesOfInterest[iBranch] << std::endl;
  }

  const EventCategoryInterpreter& eci = EventCategoryInterpreter::Instance();
  int nExpectedCategories = eci.get_number_categories();

  TH1D* Hist;
  TCanvas* Canv = new TCanvas("Canv","");
  Canv->SetBottomMargin(0.4);
  Canv->SetLogy(true);

  TString OutputName = "EventCategory.pdf";
  Canv->Print(OutputName+"[");

  for (int iBranch=0;iBranch<nBranchesOfInterest;iBranch++) {
    TString HistName = "Hist_"+BranchesOfInterest[iBranch];
    TString SampleName = TString(BranchesOfInterest[iBranch]).ReplaceAll("_EventCategory","");

    Tree->Draw(BranchesOfInterest[iBranch]+">>"+HistName+Form("(%i,%4.2f,%4.2f)",nExpectedCategories,-0.5,nExpectedCategories-0.5),"","goff");
    Hist = (TH1D*)gDirectory->Get(HistName);

    std::cout << SampleName << "  ===========================================================================================================" << std::endl;
    for (int xBin=1;xBin<=Hist->GetNbinsX();xBin++) {
      std::cout << "Bin:" << std::setw(40) << (eci.label((EventCategory)(xBin-1))).c_str() << " | Content:" << Hist->GetBinContent(xBin) << std::endl;
    }

    for (int xBin=1;xBin<=Hist->GetNbinsX();xBin++) {
      Hist->GetXaxis()->SetBinLabel(xBin,(eci.label((EventCategory)(xBin-1))).c_str());
    }
    Hist->LabelsOption("v");
    Hist->SetStats(0);
    Hist->Draw("TEXT");

    //Canv->Print(OutputName);
  }

  Canv->Print(OutputName+"]");
}
