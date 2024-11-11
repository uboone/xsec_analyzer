#pragma once

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

#include "TMath.h"
#include <TVector3.h>
#include <TLorentzVector.h>

#include "Constants.hh"

class STVTools {

private:

  double fkMiss;
  double fEMiss;
  double fPMissMinus;
  double fPMiss;
  double fPt;
  double fPL;
  double fPn;
  double fDeltaAlphaT;
  double fDeltaAlpha3Dq;
  double fDeltaAlpha3DMu;
  double fDeltaPhiT;
  double fDeltaPhi3D;
  double fECal;
  double fECalMB;
  double fEQE;
  double fQ2;
  double fA;
  double fPtx;
  double fPty;
  double fPnPerp;
  double fPnPerpx;
  double fPnPerpy;
  double fPnPar;

public:

  // Default constructor
  STVTools() {};
  void CalculateSTVs( TVector3 MuonVector, TVector3 ProtonVector,
    double MuonEnergy, double ProtonEnergy, STVCalcType CalcOption = kOpt1 );

  // Default destructor
  ~STVTools() {}

  inline double ReturnkMiss() {return fkMiss;}
  inline double ReturnEMiss() {return fEMiss;}
  inline double ReturnPMissMinus() {return fPMissMinus;}
  inline double ReturnPMiss() {return fPMiss;}
  inline double ReturnPt() {return fPt;}
  inline double ReturnPL() {return fPL;}
  inline double ReturnPn() {return fPn;}
  inline double ReturnDeltaAlphaT() {return fDeltaAlphaT;}
  inline double ReturnDeltaAlpha3Dq() {return fDeltaAlpha3Dq;}
  inline double ReturnDeltaAlpha3DMu() {return fDeltaAlpha3DMu;}
  inline double ReturnDeltaPhiT() {return fDeltaPhiT;}
  inline double ReturnDeltaPhi3D() {return fDeltaPhi3D;}
  inline double ReturnECal() {return fECal;}
  inline double ReturnECalMB() {return fECalMB;}
  inline double ReturnEQE() {return fEQE;}
  inline double ReturnQ2() {return fQ2;}
  inline double ReturnA() {return fA;}
  inline double ReturnPtx() {return fPtx;}
  inline double ReturnPty() {return fPty;}
  inline double ReturnPnPerp() {return fPnPerp;}
  inline double ReturnPnPerpx() {return fPnPerpx;}
  inline double ReturnPnPerpy() {return fPnPerpy;}
  inline double ReturnPnPar() {return fPnPar;}
};
