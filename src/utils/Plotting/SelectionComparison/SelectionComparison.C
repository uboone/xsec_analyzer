#include <string>
#include <fstream>
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

void SelectionComparison() {
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

  //Hists[Selection][EveryOtherSelection][CutInOtherSelection]
  std::vector< std::vector< std::vector< TH1D* > > > Hists;
  std::vector< std::vector< std::vector< std::string > > > HistDescript;

  Hists.resize(nSelections);
  HistDescript.resize(nSelections);
  for (int iSel=0;iSel<nSelections;iSel++) {
    Hists[iSel].resize(nSelections-1);
    HistDescript[iSel].resize(nSelections-1);
    int Counter = 0;
    for (int jSel=0;jSel<nSelections;jSel++) {
      if (iSel == jSel) continue;
      Hists[iSel][Counter].resize(SelectionCuts[jSel].size()-1);
      HistDescript[iSel][Counter].resize(SelectionCuts[jSel].size()-1);
      Counter += 1;
    }
  }

  int HistCounter = 0;
  for (int iSel=0;iSel<nSelections;iSel++) {
    int Counter = 0;
    for	(int jSel=0;jSel<nSelections;jSel++) {
      if (iSel == jSel) continue;
      for (int kCut=1;kCut<SelectionCuts[jSel].size();kCut++) {
	TString HistName = Form("Hist_%i",HistCounter);
	HistCounter += 1;
	
	TString DrawString = TString(SelectionCuts[iSel][0].c_str())+"_Selected>>"+HistName;
	TString DrawOpt = SelectionCuts[iSel][0]+"_Selected==1 && "+SelectionCuts[jSel][kCut]+"!=1";

	HistDescript[iSel][Counter][kCut-1] = DrawOpt;
	
	Tree->Draw(DrawString,DrawOpt,"goff");
	Hists[iSel][Counter][kCut-1] = (TH1D*)gDirectory->Get(HistName);
      }
      Counter += 1;
    }
  }

  for (int iSel=0;iSel<Hists.size();iSel++) {
    for (int jSel=0;jSel<Hists[iSel].size();jSel++) {
      for (int kSel=0;kSel<Hists[iSel][jSel].size();kSel++) {
	std::cout << HistDescript[iSel][jSel][kSel] << " : " << Hists[iSel][jSel][kSel]->Integral() << std::endl;
	//std::cout << HistDescript[iSel][jSel][kSel] << std::endl;
      }
      std::cout << "" << std::endl;
    }
    std::cout << "\n\n" << std::endl;
  }

}
