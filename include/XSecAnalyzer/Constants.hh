#pragma once

// Enable NuMI mode
constexpr bool useNuMI = false;

// Boundaries of the neutrino vertex fiducial volume (cm)
// This is handled the same way for reco and in MC
constexpr double FV_X_MIN =   21.5;
constexpr double FV_X_MAX =  234.85;

constexpr double FV_Y_MIN = -95.0;
constexpr double FV_Y_MAX =  95.0;

constexpr double FV_Z_MIN =   21.5;
constexpr double FV_Z_MAX =  966.8;

// PeLEE FV
/*
constexpr double FV_X_MIN =   10.0;
constexpr double FV_X_MAX =  246.0;

constexpr double FV_Y_MIN = -101.0;
constexpr double FV_Y_MAX =  101.0;

constexpr double FV_Z_MIN =   10.0;
constexpr double FV_Z_MAX =  986.0;
*/

// A few helpful dummy constants
constexpr float BOGUS = 9999.;
constexpr int BOGUS_INT = 9999;
constexpr int BOGUS_INDEX = -1;
constexpr float LOW_FLOAT = -1e30;
constexpr float DEFAULT_WEIGHT = 1.;

// Integer representation of CC versus NC used in the "ccnc" branch
// from the PeLEE ntuples
constexpr int CHARGED_CURRENT = 0;
constexpr int NEUTRAL_CURRENT = 1;

// Integer labels for primary interaction modes used in the "interaction"
// branch from the PeLEE ntuples
constexpr int QE_INTERACTION = 0;
constexpr int MEC_INTERACTION = 10;
constexpr int RES_INTERACTION = 1;

// Useful PDG codes
constexpr int ELECTRON_NEUTRINO = 12;
constexpr int MUON = 13;
constexpr int MUON_NEUTRINO = 14;
constexpr int TAU_NEUTRINO = 16;
constexpr int PROTON = 2212;
constexpr int PI_ZERO = 111;
constexpr int PI_PLUS = 211;

// Values of parameters to use in analysis cuts
constexpr float DEFAULT_PROTON_PID_CUT = 0.2;
constexpr float LEAD_P_MIN_MOM_CUT = 0.250; // GeV/c
constexpr float LEAD_P_MAX_MOM_CUT = 1.; // GeV/c
constexpr float MUON_P_MIN_MOM_CUT = 0.100; // GeV/c
constexpr float MUON_P_MAX_MOM_CUT = 1.200; // GeV/c
constexpr float CHARGED_PI_MOM_CUT = 0.; // GeV/c
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

// Mass values from GENIE v3.0.6
constexpr double TARGET_MASS = 37.215526; // 40Ar, GeV
constexpr double NEUTRON_MASS = 0.93956541; // GeV
constexpr double PROTON_MASS = 0.93827208; // GeV
constexpr double MUON_MASS = 0.10565837; // GeV
constexpr double PI_PLUS_MASS = 0.13957000; // GeV
