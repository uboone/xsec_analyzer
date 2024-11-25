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

class GKITools {

private:

  double fLeadProtonKE;
  double fEcalMB;
  double fQ;
  double fPt;
  double fPl;
  double fPtMuon;
  double fPtLeadProton;
  double fPtPion;
  double fPlMuon;
  double fPlLeadProton;
  double fPlPion;
  double fPn;
  double fDeltaAlpha3D;
  double fDeltaAlpha3DMu;
  double fDeltaPhi3D;
  double fDeltaPhi3D_pion;
  double fDeltaPhi3D_proton;
  double fDeltaPhi3D_muon;

  double fTotalProtonKE;
  double fEcalMBTotal;
  double fQTotal;
  double fPtTotal;
  double fPlTotal;
  double fPtMuonTotal;
  double fPtProtonTotal;
  double fPtPionTotal;
  double fPlMuonTotal;
  double fPlProtonTotal;
  double fPlPionTotal;
  double fPnTotal;
  double fDeltaAlpha3DTotal;
  double fDeltaAlpha3DMuTotal;
  double fDeltaPhi3DTotal;
  double fDeltaPhi3D_pionTotal;
  double fDeltaPhi3D_protonTotal;
  double fDeltaPhi3D_muonTotal;

public:

  // Default constructor
  GKITools() {};
  void CalculateGKIs( TVector3 MuonVector, TVector3 ProtonVector, TVector3 PionVector,
    double MuonEnergy, double ProtonEnergy, double PionEnergy, const std::vector<TVector3>& ProtonMomenta);

  void ComputeObservables(const TVector3 MuonVector, const TVector3 ProtonVector, const TVector3 PionVector,
               double MuonEnergy, double ProtonKE, double PionEnergy,
               double& EcalMB, double& Q, double& Pt, double& Pl,
               double& PtMuon, double& PtProton, double& PtPion,
               double& PlMuon, double& PlProton, double& PlPion,
               double& Pn, double& DeltaAlpha3D, double& DeltaAlpha3DMu,
               double& DeltaPhi3D, double& DeltaPhi3D_pion, double& DeltaPhi3D_proton, double& DeltaPhi3D_muon);
  // Default destructor
  ~GKITools() {}

  inline double ReturnLeadProtonKE() {return fLeadProtonKE;}
  inline double ReturnEcalMB() {return fEcalMB;}
  inline double ReturnQ() {return fQ;}
  inline double ReturnPt() {return fPt;}
  inline double ReturnPl() {return fPl;}
  inline double ReturnPtMuon() {return fPtMuon;}
  inline double ReturnPtProton() {return fPtLeadProton;}
  inline double ReturnPtPion() {return fPtPion;}
  inline double ReturnPlMuon() {return fPlMuon;}
  inline double ReturnPlProton() {return fPlLeadProton;}
  inline double ReturnPlPion() {return fPlPion;}
  inline double ReturnPn() {return fPn;}
  inline double ReturnDeltaAlpha3D() {return fDeltaAlpha3D;}
  inline double ReturnDeltaAlpha3DMu() {return fDeltaAlpha3DMu;}
  inline double ReturnDeltaPhi3D() {return fDeltaPhi3D;}
  inline double ReturnDeltaPhi3DPion() {return fDeltaPhi3D_pion;}
  inline double ReturnDeltaPhi3DProton() {return fDeltaPhi3D_proton;}
  inline double ReturnDeltaPhi3DMuon() {return fDeltaPhi3D_muon;}

  inline double ReturnProtonKETotal() {return fTotalProtonKE;}
  inline double ReturnEcalMBTotal() {return fEcalMBTotal;}
  inline double ReturnQTotal() {return fQTotal;}
  inline double ReturnPtTotal() {return fPtTotal;}
  inline double ReturnPlTotal() {return fPlTotal;}
  inline double ReturnPtMuonTotal() {return fPtMuonTotal;}
  inline double ReturnPtProtonTotal() {return fPtProtonTotal;}
  inline double ReturnPtPionTotal() {return fPtPionTotal;}
  inline double ReturnPlMuonTotal() {return fPlMuonTotal;}
  inline double ReturnPlProtonTotal() {return fPlProtonTotal;}
  inline double ReturnPlPionTotal() {return fPlPionTotal;}
  inline double ReturnPnTotal() {return fPnTotal;}
  inline double ReturnDeltaAlpha3DTotal() {return fDeltaAlpha3DTotal;}
  inline double ReturnDeltaAlpha3DMuTotal() {return fDeltaAlpha3DMuTotal;}
  inline double ReturnDeltaPhi3DTotal() {return fDeltaPhi3DTotal;}
  inline double ReturnDeltaPhi3DPionTotal() {return fDeltaPhi3D_pionTotal;}
  inline double ReturnDeltaPhi3DProtonTotal() {return fDeltaPhi3D_protonTotal;}
  inline double ReturnDeltaPhi3DMuonTotal() {return fDeltaPhi3D_muonTotal;}
};;
