// This script adds beamline geometry weights to NuMI  MC files
// It should be run before ProcessNTuples.C
// The beamline geometry weights are stored in a separate ROOT file, required as input
// Note - path might need fixing

#include <string>
#include <iostream>
#include <sstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TRotation.h"

float GetNuMIAngle(double px, double py, double pz, std::string direction); 

std::vector<float> getWeights(float nu_e, float nu_angle, std::vector<TH2F> &h_weights);

float checkWeight(float weight);

int main(int argc, char *argv[]) {

    // parse arguments
    if ( argc != 4 ) {
        std::cout << "Usage: AddBeamlineGeometryWeights INPUT_FILE HORN_CURRENT_MODE OUTPUT_FILE" << std::endl;
        return 1;
        }

        std::string input_filename( argv[1] );
        std::string horn_current_mode( argv[2] );
        std::string output_filename( argv[3] );

    std::cout << "Adding beamline geometry weights to file: " << input_filename << std::endl;
    if (horn_current_mode == "FHC" || horn_current_mode == "RHC") {
        std::cout << "Horn current mode: " << horn_current_mode << std::endl;
    }
    else {
        std::cout << "Error: invalid horn current mode. Valid modes: FHC, RHC" << std::endl;
        exit(1);
    }

    // load beamline variation histograms for chosen horn current mode
    TFile *beamlineVariationsFile = new TFile("src/app/NuMI/NuMI_Geometry_Weights_Histograms.root");

    if (!beamlineVariationsFile || beamlineVariationsFile->IsZombie()) {
        std::cerr << "Error: Could not open file 'NuMI_Geometry_Weights_Histograms.root'. Please check the file path and ensure the file exists." << std::endl;
        exit(1);
    }

    std::vector<TH2F> h_nue;
	std::vector<TH2F> h_nuebar;
	std::vector<TH2F> h_numu;
	std::vector<TH2F> h_numubar;

    for (int i = 1; i <= 20; i++) {
        // construct names
        std::stringstream name_nue_ss; name_nue_ss << "EnergyTheta2D/ratio_run" << i <<  "_" << horn_current_mode << "_nue_CV_AV_TPC_2D";
        std::stringstream name_nuebar_ss; name_nuebar_ss << "EnergyTheta2D/ratio_run" << i <<  "_" << horn_current_mode << "_nuebar_CV_AV_TPC_2D";
        std::stringstream name_numu_ss; name_numu_ss << "EnergyTheta2D/ratio_run" << i <<  "_" << horn_current_mode << "_numu_CV_AV_TPC_2D";
        std::stringstream name_numubar_ss; name_numubar_ss << "EnergyTheta2D/ratio_run" << i <<  "_" << horn_current_mode << "_numubar_CV_AV_TPC_2D";
        std::string name_nue = name_nue_ss.str();
        std::string name_nuebar = name_nuebar_ss.str();
        std::string name_numu = name_numu_ss.str();
        std::string name_numubar = name_numubar_ss.str();		

        // load histograms
        h_nue.push_back(*(TH2F*)beamlineVariationsFile->Get(name_nue.c_str()));
        h_nuebar.push_back(*(TH2F*)beamlineVariationsFile->Get(name_nuebar.c_str()));
        h_numu.push_back(*(TH2F*)beamlineVariationsFile->Get(name_numu.c_str()));
        h_numubar.push_back(*(TH2F*)beamlineVariationsFile->Get(name_numubar.c_str()));  
    }

	beamlineVariationsFile->Close();
    delete beamlineVariationsFile;

    // load input file
    // input file
    TFile *f = new TFile(input_filename.c_str());
    if(!f->IsOpen()) {std::cout << "Could not open input file." << std::endl; exit(1);}
    
    TTree *t_nu = (TTree*)f->Get("nuselection/NeutrinoSelectionFilter");
    TTree *t_pot = (TTree*)f->Get("nuselection/SubRun");

    // set input branch addresses
    // neutrino truth information
    int nu_pdg;
    float nu_e;
    float true_nu_px;
    float true_nu_py;
    float true_nu_pz;

    t_nu->SetBranchAddress("nu_pdg", &nu_pdg);
    t_nu->SetBranchAddress("nu_e", &nu_e);
    t_nu->SetBranchAddress("true_nu_px", &true_nu_px);
    t_nu->SetBranchAddress("true_nu_py", &true_nu_py);
    t_nu->SetBranchAddress("true_nu_pz", &true_nu_pz);

    // set addresses to weights map
	std::map<std::string, std::vector<double>> *mc_weights_map = nullptr;
	t_nu->SetBranchAddress("weights", &mc_weights_map);

    // check whether beamline geometry weights already exist
    // get first entry to populate weights
	t_nu->GetEntry(0);
    /*
	std::cout << std::endl;
	for (auto pair : *mc_weights_map) {
		std::cout << "Name: " << pair.first << ", Size: " << (pair.second).size() << std::endl;
	}
	std::cout << std::endl;
    */
    auto Horn_2kA_index = mc_weights_map->find("Horn_2kA");
    if (Horn_2kA_index != mc_weights_map->end()) {
        std::cout << "Beamline geometry weights already present." << std::endl;
        f->Close();
        delete f;
        exit(0);
    } 

    // open ouptut file, cloning tree from input file
	TFile *f_out = new TFile(output_filename.c_str(),"recreate");
	f_out->cd();
    f_out->mkdir("nuselection");
    f_out->cd("/nuselection/");

    // create new tree
    TTree *t_out_nu = t_nu->CloneTree(0);
    TTree *t_out_pot = t_pot->CloneTree();

    // event loop
    int nEntries = t_nu->GetEntries();

    std::cout << "Number entries: " << nEntries << std::endl;
    
    for (int iEntry = 0; iEntry < nEntries; iEntry++) {
    //for (int iEntry = 0; iEntry < 1; iEntry++) {

        t_nu->GetEntry(iEntry);

        if ( (iEntry != 0) && (nEntries >= 10) &&  (iEntry % (nEntries/10) == 0) ) {
            std::cout << Form("%i0%% Completed...\n", iEntry / (nEntries/10));
        }
    
        // calculate neutrino direction, relative to NuMI beam direction
        float nu_angle = GetNuMIAngle(true_nu_px, true_nu_py, true_nu_pz, "beam");

        // calculate weights
        std::vector<float> weights;
        if (nu_pdg == 12) weights = getWeights(nu_e, nu_angle, h_nue);
        else if (nu_pdg == -12) weights = getWeights(nu_e, nu_angle, h_nuebar);
        else if (nu_pdg == 14) weights = getWeights(nu_e, nu_angle, h_numu);
        else if (nu_pdg == -14) weights = getWeights(nu_e, nu_angle, h_numubar);
        else {
            std::cout << "Error: cannot get beamline variation weights" << std::endl;
            exit(1);
        }

        // sanity check
        if (weights.size() != 20) {
            std::cout << "Error: missing expected beamline variation weights" << std::endl;
            exit(1);
        }

        // populate weights
        std::vector<double> _weight_Horn_2kA {checkWeight(weights[0]), checkWeight(weights[1])};
        std::vector<double> _weight_Horn1_x_3mm {checkWeight(weights[2]), checkWeight(weights[3])};
        std::vector<double> _weight_Horn1_y_3mm {checkWeight(weights[4]), checkWeight(weights[5])};
        std::vector<double> _weight_Beam_spot_1_1mm {checkWeight(weights[6])};
        std::vector<double> _weight_Beam_spot_1_5mm {checkWeight(weights[7])};
        std::vector<double> _weight_Horn2_x_3mm {checkWeight(weights[8]), checkWeight(weights[9])};
        std::vector<double> _weight_Horn2_y_3mm {checkWeight(weights[10]), checkWeight(weights[11])};
        std::vector<double> _weight_Horns_0mm_water {checkWeight(weights[12])};;
        std::vector<double> _weight_Horns_2mm_water {checkWeight(weights[13])};;
        std::vector<double> _weight_Beam_shift_x_1mm {checkWeight(weights[14]), checkWeight(weights[15])};
        std::vector<double> _weight_Beam_shift_y_1mm {checkWeight(weights[16]), checkWeight(weights[17])};
        std::vector<double> _weight_Target_z_7mm {checkWeight(weights[18]), checkWeight(weights[19])};

        // add to map
        mc_weights_map->insert( std::pair<std::string, std::vector<double>> ("Horn_2kA", _weight_Horn_2kA) );
        mc_weights_map->insert( std::pair<std::string, std::vector<double>> ("Horn1_x_3mm", _weight_Horn1_x_3mm) );
        mc_weights_map->insert( std::pair<std::string, std::vector<double>> ("Horn1_y_3mm", _weight_Horn1_y_3mm) );
        mc_weights_map->insert( std::pair<std::string, std::vector<double>> ("Beam_spot_1_1mm", _weight_Beam_spot_1_1mm) );
        mc_weights_map->insert( std::pair<std::string, std::vector<double>> ("Beam_spot_1_5mm", _weight_Beam_spot_1_5mm) );
        mc_weights_map->insert( std::pair<std::string, std::vector<double>> ("Horn2_x_3mm", _weight_Horn2_x_3mm) );
        mc_weights_map->insert( std::pair<std::string, std::vector<double>> ("Horn2_y_3mm", _weight_Horn2_y_3mm) );
        mc_weights_map->insert( std::pair<std::string, std::vector<double>> ("Horns_0mm_water", _weight_Horns_0mm_water) );
        mc_weights_map->insert( std::pair<std::string, std::vector<double>> ("Horns_2mm_water", _weight_Horns_2mm_water) );
        mc_weights_map->insert( std::pair<std::string, std::vector<double>> ("Beam_shift_x_1mm", _weight_Beam_shift_x_1mm) );
        mc_weights_map->insert( std::pair<std::string, std::vector<double>> ("Beam_shift_y_1mm", _weight_Beam_shift_y_1mm) );
        mc_weights_map->insert( std::pair<std::string, std::vector<double>> ("Target_z_7mm", _weight_Target_z_7mm) );

        // fill new branches
        t_out_nu->Fill();
    }

   // write output
   t_out_nu->Write("NeutrinoSelectionFilter");
   t_out_pot->Write("SubRun");

   f->Close();
   f_out->Close();

   delete f;
   delete f_out;
}


float GetNuMIAngle(double px, double py, double pz, std::string direction) {

    // Variables
    TRotation RotDet2Beam;             // Rotations
    TVector3  detxyz, BeamCoords;      // Translations
    std::vector<double> rotmatrix;     // Inputs

    // input detector coordinates to translate
    detxyz = {px, py, pz};     

    // From beam to detector rotation matrix
    rotmatrix = {
        0.92103853804025681562, 0.022713504803924120662, 0.38880857519374290021,
        4.6254001262154668408e-05, 0.99829162468141474651, -0.058427989452906302359,
        -0.38947144863934973769, 0.053832413938664107345, 0.91946400794392302291 };

    // Return the TRotation
    TVector3 newX, newY, newZ;
    newX = TVector3(rotmatrix[0], rotmatrix[1], rotmatrix[2]);
    newY = TVector3(rotmatrix[3], rotmatrix[4], rotmatrix[5]);
    newZ = TVector3(rotmatrix[6], rotmatrix[7], rotmatrix[8]);

    RotDet2Beam.RotateAxes(newX, newY, newZ); // Return the TRotation now det to beam

    // Rotate to beam coords
    BeamCoords = RotDet2Beam * detxyz;

    TVector3 beamdir = {0 , 0 , 1};;
    
    // Get the angle wrt to the beam
    if (direction == "beam") beamdir = {0 , 0 , 1};
    
    // Get the angle wrt to the target to detector direction
    else if (direction == "target") {
        beamdir = {5502, 7259, 67270};
        beamdir = beamdir.Unit(); // Get the direction
    }
    else {
        std::cout << "Warning unknown angle type specified, you should check this" << std::endl;
    }
    
    double angle = BeamCoords.Angle(beamdir) * 180 / 3.1415926;

    return angle;
}

std::vector<float> getWeights(float nu_e, float nu_angle, std::vector<TH2F> &h_weights) {

    std::vector<float> weights; weights.reserve(20);

    for (int i = 0; i < h_weights.size(); i++) {
        int binx = h_weights[i].GetXaxis()->FindBin(nu_e);
        int biny = h_weights[i].GetYaxis()->FindBin(nu_angle);
        float weight = h_weights[i].GetBinContent(binx,biny);
        weights.push_back(weight);
    }

    return weights;
}

float checkWeight(float weight) {

	// infinite weight
	if (std::isinf(weight)) weight = 1.0;
    
    // nan weight
    else if (std::isnan(weight)) weight = 1.0;

	// overly large weight
	else if (weight > 30.0) weight = 1.0;

	// negative weight
	else if (weight < 0.0) weight = 1.0;

	// approximately zero weight
	else if (weight < 1e-4) weight = 0.0;

	return weight;
}