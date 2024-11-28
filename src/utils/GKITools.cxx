//Class created by Afroditi Papadopoulou (apapadop@mit.edu, apapadopoulou@anl.gov)

// _________________________________________________________________________________________________________________________________________________

#include "XSecAnalyzer/GKITools.hh"

// __________________________________________________________________________________________________________________________________________________
void GKITools::CalculateGKIs(TVector3 MuonVector, TVector3 LeadProtonVector, TVector3 PionVector,
  double MuonEnergy, double LeadProtonEnergy, double PionEnergy, const std::vector<TVector3>& ProtonMomenta){

  fLeadProtonKE = LeadProtonEnergy - PROTON_MASS;
  fTotalProtonKE = 0.0;
  TVector3 totalProtonMomentum;
  totalProtonMomentum.SetXYZ(0.0, 0.0, 0.0);

  for (const auto& protonMomentum : ProtonMomenta)
  {
    double protonEnergy = sqrt(protonMomentum.Mag2() + PROTON_MASS * PROTON_MASS);
    double protonKE = protonEnergy - PROTON_MASS;
    fTotalProtonKE += protonKE;
    totalProtonMomentum += protonMomentum;
  }
  
  // Leading Proton
  ComputeObservables(MuonVector, LeadProtonVector, PionVector, MuonEnergy, fLeadProtonKE, PionEnergy, 
    fEcalMB, fQ, fPt, fPl, fPtMuon, fPtLeadProton, fPtPion, fPlMuon, fPlLeadProton, fPlPion, fPn,
    fDeltaAlpha3D, fDeltaAlpha3DMu, fDeltaPhi3D, fDeltaPhi3D_pion, fDeltaPhi3D_proton, fDeltaPhi3D_muon);

  ComputeObservables(MuonVector, totalProtonMomentum, PionVector, MuonEnergy, fTotalProtonKE, PionEnergy, 
    fEcalMBTotal, fQTotal, fPtTotal, fPlTotal, fPtMuonTotal, fPtProtonTotal, fPtPionTotal, fPlMuonTotal, fPlProtonTotal, fPlPionTotal, fPnTotal,
    fDeltaAlpha3DTotal, fDeltaAlpha3DMuTotal, fDeltaPhi3DTotal, fDeltaPhi3D_pionTotal, fDeltaPhi3D_protonTotal, fDeltaPhi3D_muonTotal);

}

void GKITools::ComputeObservables(const TVector3 MuonVector, const TVector3 ProtonVector, const TVector3 PionVector,
  double MuonEnergy, double ProtonKE, double PionEnergy, double &ECalMB, double &Q, double &Pt, double &Pl, double &PtMuon, double &PtProton, double &PtPion, double &PlMuon, double &PlLeadProton, double &PlPion, double &Pn, 
  double &DeltaAlpha3D, double &DeltaAlpha3DMu, double &DeltaPhi3D, double &DeltaPhi3D_pion, double &DeltaPhi3D_proton, double &DeltaPhi3D_muon)
{
  double MeanExcitationEnergy = 0.0309; // GeV

  // E_cal estimator
  ECalMB = MuonEnergy + ProtonKE + PionEnergy + MeanExcitationEnergy; 

  // Momentm transfer
  TLorentzVector nuLorentzVectorMB(0.,0.,ECalMB,ECalMB);
  TLorentzVector MuonLorentzVector(MuonVector,MuonEnergy);
  TLorentzVector qLorentzVectorMB = nuLorentzVectorMB - MuonLorentzVector;

  TVector3 qVector = qLorentzVectorMB.Vect();
  Q = qVector.Mag();

  // Transverse particle momenta
  TVector3 MuonVectorTrans;
  MuonVectorTrans.SetXYZ(MuonVector.X(),MuonVector.Y(),0.);
  TVector3 ProtonVectorTrans;
  ProtonVectorTrans.SetXYZ(ProtonVector.X(), ProtonVector.Y(),0.);
  TVector3 PionVectorTrans;
  PionVectorTrans.SetXYZ(PionVector.X(),PionVector.Y(),0.);

  PtMuon = MuonVectorTrans.Mag();
  PtProton = ProtonVectorTrans.Mag();
  PtPion = PionVectorTrans.Mag();

  // Total missing transverse momentum
  TVector3 PtVector = MuonVectorTrans + ProtonVectorTrans + PionVectorTrans;
  Pt = PtVector.Mag();

  // Longitudinal particle momenta
  PlMuon = MuonVector.Z();
  PlLeadProton = ProtonVector.Z();
  PlPion = PionVector.Z();

  // Total missing longitudinal momentum
  Pl = MuonVector.Z() + ProtonVector.Z() + PionVector.Z() - ECalMB;

  // Total missing momentum
  Pn = TMath::Sqrt( Pt * Pt + Pl * Pl );
  
  // GKI Calculations
  TVector3 PnVector(PtVector.X(),PtVector.Y(),Pl);

  // Angle between q and Pn
  DeltaAlpha3D = TMath::ACos( (qVector * PnVector) / ( Q * Pn ) ) * 180./TMath::Pi();
  if (DeltaAlpha3D > 180.) { DeltaAlpha3D -= 180.; }
  if (DeltaAlpha3D < 0.) { DeltaAlpha3D += 180.; }

  // Angle between muon and Pn
  DeltaAlpha3DMu = TMath::ACos( (MuonVector * PnVector) / ( MuonVector.Mag() * Pn ) ) * 180./TMath::Pi();
  if (fDeltaAlpha3DMu > 180.) { fDeltaAlpha3DMu -= 180.; }
  if (fDeltaAlpha3DMu < 0.) { fDeltaAlpha3DMu += 180.; }

  // Opening angle between q and hadronic system
  TVector3 HadronVector = ProtonVector + PionVector;
  fDeltaPhi3D = TMath::ACos( (qVector * HadronVector) / ( Q * HadronVector.Mag() ) ) * 180./TMath::Pi();
  if (fDeltaPhi3D > 180.) { fDeltaPhi3D -= 180.; }
  if (fDeltaPhi3D < 0.) { fDeltaPhi3D += 180.; }

  // Opening angle between q and proton
  fDeltaPhi3D_proton = TMath::ACos( (qVector * ProtonVector) / ( Q * ProtonVector.Mag() ) ) * 180./TMath::Pi();
  if (fDeltaPhi3D_proton > 180.) { fDeltaPhi3D_proton -= 180.; }
  if (fDeltaPhi3D_proton < 0.) { fDeltaPhi3D_proton += 180.; }

  // Opening angle between q and pion
  fDeltaPhi3D_pion = TMath::ACos( (qVector * PionVector) / ( Q * PionVector.Mag() ) ) * 180./TMath::Pi();
  if (fDeltaPhi3D_pion > 180.) { fDeltaPhi3D_pion -= 180.; }
  if (fDeltaPhi3D_pion < 0.) { fDeltaPhi3D_pion += 180.; }

  fDeltaPhi3D_muon = TMath::ACos( (qVector * PionVector) / ( Q * MuonVector.Mag() ) ) * 180./TMath::Pi();
  if (fDeltaPhi3D_pion > 180.) { fDeltaPhi3D_pion -= 180.; }
  if (fDeltaPhi3D_pion < 0.) { fDeltaPhi3D_pion += 180.; }

}
