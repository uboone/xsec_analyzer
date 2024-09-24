#ifndef STV_TOOLS_H
#define STV_TOOLS_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

#include "TMath.h"
#include <TVector3.h>
#include <TLorentzVector.h>

#include "Constants.hh"

class STV_Tools {
  
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
  STV_Tools() {};
  void CalculateSTVs(TVector3 MuonVector,TVector3 ProtonVector, double MuonEnergy, double ProtonEnergy, STVCalcType CalcOption=kOpt1);
  
  // Default destructor
  ~STV_Tools(){}
  
  double ReturnkMiss() {return fkMiss;}
  double ReturnEMiss() {return fEMiss;}
  double ReturnPMissMinus() {return fPMissMinus;}
  double ReturnPMiss() {return fPMiss;}
  double ReturnPt() {return fPt;}
  double ReturnPL() {return fPL;}
  double ReturnPn() {return fPn;}
  double ReturnDeltaAlphaT() {return fDeltaAlphaT;}
  double ReturnDeltaAlpha3Dq() {return fDeltaAlpha3Dq;}
  double ReturnDeltaAlpha3DMu() {return fDeltaAlpha3DMu;}
  double ReturnDeltaPhiT() {return fDeltaPhiT;}
  double ReturnDeltaPhi3D() {return fDeltaPhi3D;}
  double ReturnECal() {return fECal;}
  double ReturnECalMB() {return fECalMB;}
  double ReturnEQE() {return fEQE;}
  double ReturnQ2() {return fQ2;}
  double ReturnA() {return fA;}
  double ReturnPtx() {return fPtx;}
  double ReturnPty() {return fPty;}
  double ReturnPnPerp() {return fPnPerp;}
  double ReturnPnPerpx() {return fPnPerpx;}
  double ReturnPnPerpy() {return fPnPerpy;}
  double ReturnPnPar() {return fPnPar;}
};

#endif
