
// Univmake expects all MC files to have identical systematic weights structure.
// For the original processing of NuMI Dirt MC the weights were not included.
// This script adds a fake set weights with identical structure to the NuMI 
// Dirt MC files to make them compatible.
// This script should be run after addBeamlineGeometryWeights.cpp and before ProcessNTuples.C.


#include <iostream>
#include <string>
#include <map>
#include <vector>

#include "TFile.h"
#include "TTree.h"

int main(int argc, char *argv[]) {
	
	// parse arguments
	if ( argc != 4 ) {
    	std::cout << "Usage: AddFakeWeights INPUT_FILE INPUT_FULL_WEIGHTS_FILE OUTPUT_FILE" << std::endl;
    	return 1;
  	}

  	std::string input_filename( argv[1] );
  	std::string input_full_weights_filename( argv[2] );
  	std::string output_filename( argv[3] );

    std::cout << "Adding fake weights to file: " << input_filename << ", using structure from: " << input_full_weights_filename << std::endl;

    // load file containing real weights structure
	TFile *f1 = new TFile(input_full_weights_filename.c_str());
	if(!f1->IsOpen()) {std::cout << "Could not open input file containing weights structure." << std::endl; exit(1);}

    // load tree
	TTree *t1 = (TTree*)f1->Get("nuselection/NeutrinoSelectionFilter");

	// set addresses to weights map
	std::map<std::string, std::vector<double>> *weights = nullptr;
	t1->SetBranchAddress("weights", &weights);

	// get first entry to populate weights
	t1->GetEntry(0);
	
	// number of weights
	std::cout <<  "Number of weight pairs in map: " << weights->size() << std::endl;

    // create new set of fake weights to populate
	std::map<std::string, std::vector<double>> fake_weights;

	// loop through weights
	std::cout << std::endl;
	for (auto pair : *weights) {
		std::cout << "Name: " << pair.first << ", Size: " << (pair.second).size() << std::endl;
		// reset all weights to 1
		std::fill((pair.second).begin(), (pair.second).end(), 1);
		// add to new vector
		fake_weights.insert({pair.first, pair.second});
	}
	std::cout << std::endl;

	f1->Close();
    delete f1;    

    // open second file weights are to be added to
	TFile *f2 = new TFile(input_filename.c_str());
	if(!f2->IsOpen()) {std::cout << "Could not open input file." << std::endl; exit(1);}
	TTree *t2_nu = (TTree*)f2->Get("nuselection/NeutrinoSelectionFilter");
    TTree *t2_pot = (TTree*)f2->Get("nuselection/SubRun");

	// open third file to write out, cloning tree from second file
	TFile *f3 = new TFile(output_filename.c_str(),"recreate");
	f3->cd();
    f3->mkdir("nuselection");
    f3->cd("/nuselection/");

    // create new tree, disabling partially filled weights branch
    t2_nu->SetBranchStatus("weights",0);
    TTree *t3_nu = t2_nu->CloneTree(0);
    TTree *t3_pot = t2_pot->CloneTree();

    // create new branch weights
    std::map<std::string, std::vector<double>> new_weights;
    t3_nu->Branch( "weights", "std::map<std::string, std::vector<double>>", &new_weights );

    // loop over events adding fake weights to new branch
    int n_entries = t2_nu->GetEntries();

    for (int e = 0; e < n_entries; e++) {

		// get current entry
	  	t2_nu->GetEntry(e);    	

	    if ( (e != 0) && (n_entries >= 10) &&  (e % (n_entries/10) == 0) ) {
	      std::cout << Form("%i0%% Completed...\n", e / (n_entries/10));
	    }

	    new_weights = fake_weights;
	    t3_nu->Fill();
	}

	// write output
    t3_nu->Write("NeutrinoSelectionFilter");
    t3_pot->Write("SubRun");

    f2->Close();
    f3->Close();

    delete f2;
    delete f3;
    
}