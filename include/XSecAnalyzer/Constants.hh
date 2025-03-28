#pragma once
#include <cmath>

// Boundaries of the neutrino vertex fiducial volume (cm)
// This is handled the same way for reco and in MC
constexpr double FV_X_MIN =   21.5;
constexpr double FV_X_MAX =  234.85;

constexpr double FV_Y_MIN = -95.0;
constexpr double FV_Y_MAX =  95.0;

constexpr double FV_Z_MIN =   21.5;
constexpr double FV_Z_MAX =  966.8;

// Pion momentum calibration from range-based momentum
constexpr float A = 1.6353;
constexpr float B = 0.0022385;
constexpr float C = 1.5668;
constexpr float D = 0.013327;

// constants needed for multi-plane Bragg calculations
const float FLOAT_MIN_VALUE = -340282346638528859811704183484516925440.000000f;
const auto SIN2_ANGLE_THRESHOLD = 0.175f;
const auto PI_BY_3 = std::acos(0.5f);
const auto LENGTH_FRACTION = 0.333333333f;
const auto N_HITS_TO_SKIP = 3;

// A few helpful dummy constants
constexpr float BOGUS = 9999.;
constexpr int BOGUS_INT = 9999;
constexpr int BOGUS_INDEX = -1;
constexpr float LOW_FLOAT = -1e30;
constexpr float DEFAULT_WEIGHT = 1.;

// Integer representation of CC versus NC for the ccnc branch
constexpr int CHARGED_CURRENT = 0;
constexpr int NEUTRAL_CURRENT = 1;

// Useful PDG codes
constexpr int ELECTRON_NEUTRINO = 12;
constexpr int MUON = 13;
constexpr int MUON_NEUTRINO = 14;
constexpr int TAU_NEUTRINO = 16;
constexpr int PROTON = 2212;
constexpr int PI_ZERO = 111;
constexpr int PI_PLUS = 211;
constexpr int NEUTRON = 2112;

// Values of parameters to use in analysis cuts
constexpr float DEFAULT_PROTON_PID_CUT = 0.2;
constexpr float LEAD_P_MIN_MOM_CUT = 0.300; // GeV/c
constexpr float LEAD_P_MAX_MOM_CUT = 1.; // GeV/c
constexpr float MUON_P_MIN_MOM_CUT = 0.150; // GeV/c
constexpr float MUON_P_MAX_MOM_CUT = 1.200; // GeV/c
constexpr float CHARGED_PI_MOM_CUT = 0.05; // GeV/c
constexpr float MUON_MOM_QUALITY_CUT = 0.25; // fractional difference
constexpr float PROTON_MIN_MOM_CUT = 0.3; //GeV/c
constexpr float PROTON_MAX_MOM_CUT = 1.0; //GeV/c
constexpr float TOPO_SCORE_CUT = 0.1;
constexpr float COSMIC_IP_CUT = 10.; // cm
constexpr float MUON_TRACK_SCORE_CUT = 0.8;
constexpr float MUON_VTX_DISTANCE_CUT = 4.; // cm
constexpr float MUON_LENGTH_CUT = 10.; // cm
constexpr float MUON_PID_CUT = 0.2;
constexpr float TRACK_SCORE_CUT = 0.5;
constexpr float PROTON_BDT_CUT = 0.01;
constexpr float MUON_BDT_CUT = 0.0;
constexpr float TOPO_SCORE_CUT_2 = 0.67; // TODO: check this value/remove as unusued
constexpr float PFP_DISTANCE_CUT = 9.5;  //cm 

// Boundaries of the proton containment volume (used in reco only) in cm
constexpr double PCV_X_MIN =   10.;
constexpr double PCV_X_MAX =  246.35;
constexpr double PCV_Y_MIN = -106.5;
constexpr double PCV_Y_MAX =  106.5;
constexpr double PCV_Z_MIN =   10.;
constexpr double PCV_Z_MAX = 1026.8;

// Mass values from GENIE v3.0.6
constexpr double TARGET_MASS = 37.215526; // 40Ar, GeV
constexpr double NEUTRON_MASS = 0.93956541; // GeV
constexpr double PROTON_MASS = 0.93827208; // GeV
constexpr double MUON_MASS = 0.10565837; // GeV
constexpr double PI_PLUS_MASS = 0.13957000; // GeV

// This binding energy value is used in GENIE v3.0.6
//constexpr double BINDING_ENERGY = 0.0295; // 40Ar, GeV

// This value is the shell-occupancy-weighted mean of the $E_{\alpha}$ values
// listed for 40Ar in Table II of arXiv:1609.03530. MINERvA uses an identical
// procedure for 12C to obtain the binding energy value of 27.13 MeV, which is
// adopted in their STV analysis described in arXiv:1910.08658.
//constexpr double BINDING_ENERGY = 0.02478; // 40Ar, GeV

constexpr double BINDING_ENERGY = 0.34381; // 40Ar, GeV
constexpr double MEAN_EXCITATION_ENERGY = 0.0309; 

enum VarType{kString, kDouble, kFloat, kInteger, kBool, kTVector, kSTDVector};

enum STVCalcType{kOpt1,kOpt2,kOpt3,kOpt4,nOpts};
