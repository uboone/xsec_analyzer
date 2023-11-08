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
		STV_Tools(TVector3 MuonVector,TVector3 ProtonVector, double MuonEnergy, double ProtonEnergy);

		// Default destructor
		~STV_Tools(){}

		double ReturnkMiss();
		double ReturnEMiss();
		double ReturnPMissMinus();
		double ReturnPMiss();
		double ReturnPt();
		double ReturnPL();
		double ReturnPn();
		double ReturnDeltaAlphaT();
		double ReturnDeltaAlpha3Dq();
		double ReturnDeltaAlpha3DMu();			
		double ReturnDeltaPhiT();
		double ReturnDeltaPhi3D();		
		double ReturnECal();
		double ReturnECalMB();		
		double ReturnEQE();
		double ReturnQ2();		
		double ReturnA();
		double ReturnPtx();
		double ReturnPty();	
		double ReturnPnPerp();
		double ReturnPnPerpx();
		double ReturnPnPerpy();				
		double ReturnPnPar();			

};

#endif
