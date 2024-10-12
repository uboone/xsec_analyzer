//Class created by Afroditi Papadopoulou (apapadop@mit.edu, apapadopoulou@anl.gov)

// _________________________________________________________________________________________________________________________________________________

#include "XSecAnalyzer/STVTools.hh"

// __________________________________________________________________________________________________________________________________________________

void STVTools::CalculateSTVs( TVector3 MuonVector, TVector3 ProtonVector,
  double MuonEnergy, double ProtonEnergy, STVCalcType CalcOpt )
{
  double BindingEnergy_GeV;

  switch (CalcOpt) {
  case kOpt1:
  case kOpt4:
    // This value is the shell-occupancy-weighted mean of the $E_{\alpha}$ values
    // listed for 40Ar in Table II of arXiv:1609.03530. MINERvA uses an identical
    // procedure for 12C to obtain the binding energy value of 27.13 MeV, which is
    // adopted in their STV analysis described in arXiv:1910.08658
    BindingEnergy_GeV = 0.02478;
    break;
  case kOpt2:
    // For the calculation of the excitation energies
    // https://doi.org/10.1140/epjc/s10052-019-6750-3
    BindingEnergy_GeV = 0.0309;
    break;
  case kOpt3:
    // https://github.com/afropapp13/myClasses/blob/de02b1da3d7629146a9a3f2244e3f9227d4059c8/STV_Tools.cxx#L20
    BindingEnergy_GeV = 0.04;
  default:
    std::cerr << "STVCalcOpt not defined:" << CalcOpt << std::endl;
    throw;
  }

  double DeltaM2 = TMath::Power(NEUTRON_MASS,2.) - TMath::Power(PROTON_MASS,2.);

  TVector3 MuonVectorTrans;
  MuonVectorTrans.SetXYZ(MuonVector.X(),MuonVector.Y(),0.);
  double MuonVectorTransMag = MuonVectorTrans.Mag();

  TVector3 MuonVectorLong;
  MuonVectorLong.SetXYZ(0.,0.,MuonVector.Z());
  double MuonVectorLongMag = MuonVectorLong.Mag();

  TLorentzVector MuonLorentzVector(MuonVector,MuonEnergy);

  TVector3 ProtonVectorTrans;
  ProtonVectorTrans.SetXYZ(ProtonVector.X(),ProtonVector.Y(),0.);
  double ProtonVectorTransMag = ProtonVectorTrans.Mag();

  TVector3 ProtonVectorLong;
  ProtonVectorLong.SetXYZ(0.,0.,ProtonVector.Z());
  double ProtonVectorLongMag = ProtonVectorLong.Mag();

  TLorentzVector ProtonLorentzVector(ProtonVector,ProtonEnergy);
  double ProtonKE = ProtonEnergy - PROTON_MASS;

  TVector3 PtVector = MuonVectorTrans + ProtonVectorTrans;

  fPt = PtVector.Mag();
  TVector2 fPt_2DVec = (MuonVector + ProtonVector).XYvector();

  fDeltaAlphaT = TMath::ACos( (- MuonVectorTrans * PtVector) / ( MuonVectorTransMag * fPt ) ) * 180./TMath::Pi();
  if (fDeltaAlphaT > 180.) { fDeltaAlphaT -= 180.; }
  if (fDeltaAlphaT < 0.) { fDeltaAlphaT += 180.; }

  fDeltaPhiT = TMath::ACos( (- MuonVectorTrans * ProtonVectorTrans) / ( MuonVectorTransMag * ProtonVectorTransMag ) ) * 180./TMath::Pi();
  if (fDeltaPhiT > 180.) { fDeltaPhiT -= 180.; }
  if (fDeltaPhiT < 0.) { fDeltaPhiT += 180.; }

  // -------------------------------------------------------------------------------------------------------------------------
  // Calorimetric Energy Reconstruction

  fECal = MuonEnergy + ProtonKE + BindingEnergy_GeV; // GeV

  // QE Energy Reconstruction

  double EQENum = 2 * (NEUTRON_MASS - BindingEnergy_GeV) * MuonEnergy - (BindingEnergy_GeV*BindingEnergy_GeV - 2 * NEUTRON_MASS *BindingEnergy_GeV + MUON_MASS * MUON_MASS + DeltaM2);
  double EQEDen = 2 * ( NEUTRON_MASS - BindingEnergy_GeV - MuonEnergy + MuonVector.Mag() * MuonVector.CosTheta() );
  fEQE = EQENum / EQEDen;

  // Reconstructed Q2

  TLorentzVector nuLorentzVector(0.,0.,fECal,fECal);
  TLorentzVector qLorentzVector = nuLorentzVector - MuonLorentzVector;
  fQ2 = - qLorentzVector.Mag2();

  // https://journals.aps.org/prd/pdf/10.1103/PhysRevD.101.092001

  // Just the magnitudes
  // fPtx = fPt * TMath::Sin(fDeltaAlphaT * TMath::Pi() / 180.);
  // fPty = fPt * TMath::Cos(fDeltaAlphaT * TMath::Pi() / 180.);

  TVector3 zUnit(0.,0.,1.);

  // Small differences found compared to Stepehen's code when MuonVectorTransMag ~ 0. Moved to using Stepehen's code below
  //fPtx = ( zUnit.Cross(MuonVectorTrans) ).Dot(PtVector) / MuonVectorTransMag;
  //fPty = - (MuonVectorTrans).Dot(PtVector) / MuonVectorTransMag;

  TVector2 xTUnit = zUnit.Cross(MuonVector).XYvector().Unit();
  fPtx = xTUnit.X()*fPt_2DVec.X() + xTUnit.Y()*fPt_2DVec.Y();

  TVector2 yTUnit = (-MuonVector).XYvector().Unit();
  fPty= yTUnit.X()*fPt_2DVec.X() + yTUnit.Y()*fPt_2DVec.Y();

  // -------------------------------------------------------------------------------------------------------------------------

  // JLab Light Cone Variables

  TLorentzVector MissLorentzVector = MuonLorentzVector + ProtonLorentzVector - nuLorentzVector;

  fEMiss = TMath::Abs(MissLorentzVector.E());
  fPMiss = (MissLorentzVector.Vect()).Mag();

  // fPMissMinus = fEMiss - MissLorentzVector.Z();

  // Suggestion from Jackson to avoid Ecal assumption
  fPMissMinus = (MuonEnergy - MuonVector.Z()) + (ProtonEnergy - ProtonVector.Z());

  double fkMissNum = ( TMath::Power(fPt,2.) + TMath::Power(PROTON_MASS,2.) );
  double fkMissDen = ( fPMissMinus * (2*PROTON_MASS - fPMissMinus) );

  double fkMiss2 = TMath::Power(PROTON_MASS,2.) * fkMissNum / fkMissDen - TMath::Power(PROTON_MASS,2.); // Jackson's GlueX note

  fkMiss = sqrt(fkMiss2);

  fA = fPMissMinus / PROTON_MASS;

  // -------------------------------------------------------------------------------------------------------------------------

  // Minerva longitudinal & total variables

  // For the calculation of the masses
  //https://journals.aps.org/prc/pdf/10.1103/PhysRevC.95.065501

  double MA = 22 * NEUTRON_MASS + 18 * PROTON_MASS - 0.34381; // GeV

  // For the calculation of the excitation energies
  // https://doi.org/10.1140/epjc/s10052-019-6750-3

  double MAPrime = MA - NEUTRON_MASS + BindingEnergy_GeV; // GeV, constant obtained from table 7

  // For the calculation of p_n, back to the Minerva PRL
  // https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.121.022504

  // Original code from Afro's calculation framework was found to have a bug - now fixed
  //double R = MA + (MuonVectorLong + ProtonVectorLong).Mag() - MuonEnergy - ProtonEnergy; // Equation 8
  double R = MA + MuonVectorLong.Z() + ProtonVectorLong.Z() - MuonEnergy - ProtonEnergy; // Equation 8

  // -------------------------------------------------------------------------------------------------------------------------

  // Beyond the transverse variables
  // Based on Andy F's xsec meeting presentation
  // https://microboone-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=38090&filename=BeyondTransverseVariables_xsec_2022_06_14.pdf&version=1

  fECalMB = MuonEnergy + ProtonKE + BindingEnergy_GeV; // GeV, after discussion with Andy F who got the numbers from Jan S
  TLorentzVector nuLorentzVectorMB(0.,0.,fECalMB,fECalMB);
  TLorentzVector qLorentzVectorMB = nuLorentzVectorMB - MuonLorentzVector;

  // Equation 7
  // Abandoned this expression on Mar 6 2023
  switch (CalcOpt) {
  case kOpt1:
  case kOpt2:
  case kOpt3:
    fPL = 0.5 * R - (MAPrime * MAPrime + fPt * fPt) / (2 * R);
    break;
  case kOpt4:
    fPL = MuonVector.Z() + ProtonVector.Z() - fECalMB;
    break;
  default:
    std::cerr << "STVCalcOpt not defined:" << CalcOpt << std::endl;
    throw;
  }
  TVector3 PnVector(PtVector.X(),PtVector.Y(),fPL);

  TVector3 qVector = qLorentzVectorMB.Vect();
  TVector3 qTVector(qVector.X(), qVector.X(), 0.);
  TVector3 qVectorUnit = qVector.Unit();
  TVector3 qTVectorUnit = qTVector.Unit();

  fPn = TMath::Sqrt( fPt * fPt + fPL * fPL );

  double qMag = qVector.Mag();
  fDeltaAlpha3Dq = TMath::ACos( (qVector * PnVector) / ( qMag * fPn ) ) * 180./TMath::Pi();
  if (fDeltaAlpha3Dq > 180.) { fDeltaAlpha3Dq -= 180.; }
  if (fDeltaAlpha3Dq < 0.) { fDeltaAlpha3Dq += 180.; }
  //if (fDeltaAlpha3Dq < 10) {cout << "qVector * PnVector = " << qVector * PnVector << ",  qMag = " << qMag << ", fPn = " << fPn << ", fDeltaAlpha3Dq = " << fDeltaAlpha3Dq << endl;}

  fDeltaAlpha3DMu = TMath::ACos( -(MuonVector * PnVector) / ( MuonVector.Mag() * fPn ) ) * 180./TMath::Pi();
  if (fDeltaAlpha3DMu > 180.) { fDeltaAlpha3DMu -= 180.; }
  if (fDeltaAlpha3DMu < 0.) { fDeltaAlpha3DMu += 180.; }

  fDeltaPhi3D = TMath::ACos( (qVector * ProtonVector) / ( qMag * ProtonVector.Mag() ) ) * 180./TMath::Pi();
  if (fDeltaPhi3D > 180.) { fDeltaPhi3D -= 180.; }
  if (fDeltaPhi3D < 0.) { fDeltaPhi3D += 180.; }

  // Magnitudes
  fPnPerp = fPn * sin(fDeltaAlpha3Dq * TMath::Pi() / 180.);
  fPnPar = fPn * cos(fDeltaAlpha3Dq * TMath::Pi() / 180.);

  // fPnPerp = - (zUnit.Cross(qVectorUnit) ).Dot(PnVector);
  // fPnPar = qVectorUnit.Dot(PnVector);
  // fPnPar = fECalMB - MuonVector.Z() - ProtonVector.Z();

  fPnPerpx = ( qTVectorUnit.Cross(zUnit) ).Dot(PnVector);
  fPnPerpy = ( qVectorUnit.Cross( (qTVectorUnit.Cross(zUnit) ) ) ).Dot(PnVector);
}
