// This script adds beamline geometry weights to NuMI  MC files
// It should be run before ProcessNTuples.C
// The beamline geometry weights are stored in a separate ROOT file, required as input

#include <string>
#include <iostream>
#include <sstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TRotation.h"


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

void addBeamlineGeometryWeights(std::string input_filename, std::string horn_current_mode) {

    std::cout << "Adding beamline geometry weights to file: " << input_filename << std::endl;
    if (horn_current_mode == "FHC" || horn_current_mode == "RHC") {
        std::cout << "Horn current mode: " << horn_current_mode << std::endl;
    }
    else {
        std::cout << "Error: invalid horn current mode. Valid modes: FHC, RHC" << std::endl;
        exit(1);
    }

    // load beamline variation histograms for chosen horn current mode
    TFile *beamlineVariationsFile = new TFile("/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/NuMIBeamLineWeights/NuMI_Geometry_Weights_Histograms.root");

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
    TFile *f = new TFile(input_filename.c_str(), "UPDATE");
    TTree *t_nu = (TTree*)f->Get("nuselection/NeutrinoSelectionFilter");

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

    // check whether branches already exist
    auto branches = t_nu->GetListOfBranches();
    if (branches->FindObject("_weight_Horn_2kA")) {
        std::cout << "Beamline geometry weights already present." << std::endl;
        f->Close();
        delete f;
        exit(0);
    }

    // set output branch addresses
    std::vector<float> _weight_Horn_2kA;
    std::vector<float> _weight_Horn1_x_3mm;
    std::vector<float> _weight_Horn1_y_3mm;
    std::vector<float> _weight_Beam_spot_1_1mm;
    std::vector<float> _weight_Beam_spot_1_5mm;
    std::vector<float> _weight_Horn2_x_3mm;
    std::vector<float> _weight_Horn2_y_3mm;
    std::vector<float> _weight_Horns_0mm_water;
    std::vector<float> _weight_Horns_2mm_water;
    std::vector<float> _weight_Beam_shift_x_1mm;
    std::vector<float> _weight_Beam_shift_y_1mm;
    std::vector<float> _weight_Target_z_7mm;

    auto branch_weight_Horn_2kA = t_nu->Branch("_weight_Horn_2kA", &_weight_Horn_2kA);
    auto branch_weight_Horn1_x_3mm = t_nu->Branch("_weight_Horn1_x_3mm", &_weight_Horn1_x_3mm);
    auto branch_weight_Horn1_y_3mm = t_nu->Branch("_weight_Horn1_y_3mm", &_weight_Horn1_y_3mm);
    auto branch_weight_Beam_spot_1_1mm = t_nu->Branch("_weight_Beam_spot_1_1mm", &_weight_Beam_spot_1_1mm);
    auto branch_weight_Beam_spot_1_5mm = t_nu->Branch("_weight_Beam_spot_1_5mm", &_weight_Beam_spot_1_5mm);
    auto branch_weight_Horn2_x_3mm = t_nu->Branch("_weight_Horn2_x_3mm", &_weight_Horn2_x_3mm);
    auto branch_weight_Horn2_y_3mm = t_nu->Branch("_weight_Horn2_y_3mm", &_weight_Horn2_y_3mm);
    auto branch_weight_Horns_0mm_water = t_nu->Branch("_weight_Horns_0mm_water", &_weight_Horns_0mm_water);
    auto branch_weight_Horns_2mm_water = t_nu->Branch("_weight_Horns_2mm_water", &_weight_Horns_2mm_water);
    auto branch_weight_Beam_shift_x_1mm = t_nu->Branch("_weight_Beam_shift_x_1mm", &_weight_Beam_shift_x_1mm);
    auto branch_weight_Beam_shift_y_1mm = t_nu->Branch("_weight_Beam_shift_y_1mm", &_weight_Beam_shift_y_1mm);
    auto branch_weight_Target_z_7mm = t_nu->Branch("_weight_Target_z_7mm", &_weight_Target_z_7mm);

    // event loop
    int nEntries = t_nu->GetEntries();

    std::cout << "Number entries: " << nEntries << std::endl;
    
    for (int iEntry = 0; iEntry < nEntries; iEntry++) {

        t_nu->GetEntry(iEntry);

        if ( (iEntry != 0) && (nEntries >= 10) &&  (iEntry % (nEntries/10) == 0) ) {
            std::cout << Form("%i0%% Completed...\n", iEntry / (nEntries/10));
        }

        // clear weights from current entry
        _weight_Horn_2kA.clear();
        _weight_Horn1_x_3mm.clear();
        _weight_Horn1_y_3mm.clear();
        _weight_Beam_spot_1_1mm.clear();
        _weight_Beam_spot_1_5mm.clear();
        _weight_Horn2_x_3mm.clear();
        _weight_Horn2_y_3mm.clear();
        _weight_Horns_0mm_water.clear();
        _weight_Horns_2mm_water.clear();
        _weight_Beam_shift_x_1mm.clear();
        _weight_Beam_shift_y_1mm.clear();
        _weight_Target_z_7mm.clear();
    
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
        _weight_Horn_2kA.push_back(checkWeight(weights[0]));
        _weight_Horn_2kA.push_back(checkWeight(weights[1]));
        _weight_Horn1_x_3mm.push_back(checkWeight(weights[2]));
        _weight_Horn1_x_3mm.push_back(checkWeight(weights[3]));
        _weight_Horn1_y_3mm.push_back(checkWeight(weights[4]));
        _weight_Horn1_y_3mm.push_back(checkWeight(weights[5]));
        _weight_Beam_spot_1_1mm.push_back(checkWeight(weights[6]));
        _weight_Beam_spot_1_5mm.push_back(checkWeight(weights[7]));
        _weight_Horn2_x_3mm.push_back(checkWeight(weights[8]));
        _weight_Horn2_x_3mm.push_back(checkWeight(weights[9]));
        _weight_Horn2_y_3mm.push_back(checkWeight(weights[10]));
        _weight_Horn2_y_3mm.push_back(checkWeight(weights[11]));
        _weight_Horns_0mm_water.push_back(checkWeight(weights[12]));
        _weight_Horns_2mm_water.push_back(checkWeight(weights[13]));
        _weight_Beam_shift_x_1mm.push_back(checkWeight(weights[14]));
        _weight_Beam_shift_x_1mm.push_back(checkWeight(weights[15]));
        _weight_Beam_shift_y_1mm.push_back(checkWeight(weights[16]));
        _weight_Beam_shift_y_1mm.push_back(checkWeight(weights[17]));
        _weight_Target_z_7mm.push_back(checkWeight(weights[18]));
        _weight_Target_z_7mm.push_back(checkWeight(weights[19]));

        // fill new branches
        branch_weight_Horn_2kA->Fill();
        branch_weight_Horn1_x_3mm->Fill();
        branch_weight_Horn1_y_3mm->Fill();
        branch_weight_Beam_spot_1_1mm->Fill();
        branch_weight_Beam_spot_1_5mm->Fill();
        branch_weight_Horn2_x_3mm->Fill();
        branch_weight_Horn2_y_3mm->Fill();
        branch_weight_Horns_0mm_water->Fill();
        branch_weight_Horns_2mm_water->Fill();
        branch_weight_Beam_shift_x_1mm->Fill();
        branch_weight_Beam_shift_y_1mm->Fill();
        branch_weight_Target_z_7mm->Fill();
    }

    // write file
    f->cd("/nuselection/");
    t_nu->Write("", TObject::kOverwrite);
    f->Close();
    delete f;    
}