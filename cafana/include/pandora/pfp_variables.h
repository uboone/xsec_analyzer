/**
 * @file pfp_variables.h
 * @brief Header file for definitions of variables which act on single
 * pfps.
 * @details This file contains definitions of variables which act on single
 * pfps. Each variable is implemented as a function which takes a pfp
 * object as an argument and returns a double. These variables are intended to
 * be used to define more complex variables which act on slices.
 * @author mueller@fnal.gov
*/
#ifndef PFP_VARIABLES_H
#define PFP_VARIABLES_H
#define ELECTRON_MASS 0.5109989461
#define MUON_MASS 105.6583745
#define PION_MASS 139.57039
#define PROTON_MASS 938.2720813

#include "sbnana/SBNAna/Cuts/VolumeDefinitions.h"

#include <TMath.h>

/**
 * @namespace pvars
 * @brief Namespace for organizing generic variables which act on single
 * pfps.
 * @details This namespace is intended to be used for organizing variables which
 * act on single pfps. Each variable is implemented as a function which
 * takes a pfp object as an argument and returns a double. The function
 * should be templated on the type of pfp object if the variable is
 * intended to be used on both true and reconstructed pfps.
 * @note The namespace is intended to be used in conjunction with the
 * vars namespace, which is used for organizing variables which act on
 * interactions.
 */
namespace pfpvars
{
    //Shower variables

    // -direction
    template<class T>
        double kShowerDirx(const T & p) {return p.shw.dir.x;}
    template<class T>
        double kShowerDiry(const T & p) {return p.shw.dir.y;}
    template<class T>
        double kShowerDirz(const T & p) {return p.shw.dir.z;}

    // -other properties (still important)
    template<class T>
        double kShowerEnergy(const T & p) {return p.shw.bestplane_energy;}
    template<class T>
        double kShowerdEdx(const T & p) {return p.shw.bestplane_dEdx;}
    template<class T>
        double kShowerConversionGap(const T & p) {return p.shw.conversion_gap;}
    template<class T>
        double kShowerLength(const T & p) {return p.shw.len;}

    //Track variables
    
    // -direction
    template<class T>
        double kTrackDirx(const T & p) {return p.trk.dir.x;}
    template<class T>
        double kTrackDiry(const T & p) {return p.trk.dir.y;}
    template<class T>
        double kTrackDirz(const T & p) {return p.trk.dir.z;}
    
    // -other properties (still important)
    template<class T>
        bool kTrackCont(const T & p) {return (PtInVolAbsX(p.trk.start,fvndExit) && PtInVolAbsX(p.trk.end,fvndExit));}

    template<class T>
        int kTrackBestPlane(const T & p) {return p.trk.bestplane;}

    template<class T>
        double kTrackLength(const T & p) {return p.trk.len;}

    // -momentum
    template<class T>
        double kTrackMomPion(const T & p){
            if (kTrackCont(p)) return p.trk.rangeP.p_pion;
            else return p.trk.mcsP.fwdP_pion;
        }
    template<class T>
        double kTrackMomMuon(const T & p){
            if (kTrackCont(p)) return p.trk.rangeP.p_muon;
            else return p.trk.mcsP.fwdP_muon;
        }
    template<class T>
        double kTrackMomProton(const T & p){
            if (kTrackCont(p)) return p.trk.rangeP.p_proton;
            else return p.trk.mcsP.fwdP_proton;
        }

    // -chi2 scores
    template<class T>
        double kTrackChi2Muon(const T & p) {
            return p.trk.chi2pid[kTrackBestPlane(p)].chi2_muon;
        }
    template<class T>
        double kTrackChi2Pion(const T & p) {
            return p.trk.chi2pid[kTrackBestPlane(p)].chi2_pion;
        }
    template<class T>
        double kTrackChi2Proton(const T & p) {
            return p.trk.chi2pid[kTrackBestPlane(p)].chi2_proton;
        }

    /**
    * @brief Find the PDG value using chi2.
    * @details This function finds the PDG value using the chi2 method.
    * @param p the pfp to apply the variable on.
    * @return the PDG value (13 for muon, 211 for pion, 2212 for proton).
     */
    template<class T>
        double kTrackPDG(const T & p) {
            if ((kTrackLength(p) > 32) && (kTrackChi2Muon(p) < 18) && (kTrackChi2Proton(p) > 87)){
                return 13;
            }
            else if (kTrackChi2Proton(p) < 97){
                return 2212;
            }
            else if ((!std::isnan(kTrackChi2Muon(p))) && (!std::isnan(kTrackChi2Pion(p))) && (!std::isnan(kTrackChi2Proton(p)))){
                return 211;
            }
            else return 0;
        }

    //PFP variables
    template<class T>
        double kPFPTrackScore(const T & p) {return p.trackScore;}
        
    template<class T>
        double kPFPEnergy(const T & p) {
            if (kPFPTrackScore(p) < 0.5) return kShowerEnergy(p);
            else if(!std::isnan(kPFPTrackScore(p))){
                if (kTrackPDG(p) == 13) return TMath::Sqrt(TMath::Power(MUON_MASS,2) + TMath::Power(kTrackMomMuon(p),2));
                else if (kTrackPDG(p) == 211) return TMath::Sqrt(TMath::Power(PION_MASS,2) + TMath::Power(kTrackMomPion(p),2));
                else if (kTrackPDG(p) == 2212) return TMath::Sqrt(TMath::Power(PROTON_MASS,2) + TMath::Power(kTrackMomProton(p),2));
            }
            else return 0.0;
        }
    
    template<class T>
        double kPFPKineticEnergy(const T & p){
            if (kPFPTrackScore(p) < 0.5) return kShowerEnergy(p);
            else if(!std::isnan(kPFPTrackScore(p))){
                if (kTrackPDG(p) == 13) return kPFPEnergy(p) - MUON_MASS;
                else if (kTrackPDG(p) == 211) return kPFPEnergy(p) - PION_MASS;
                else if (kTrackPDG(p) == 2212) return kPFPEnergy(p) - PROTON_MASS;
            }
            else return 0.0;
        }
    

}
#endif // PFP_VARIABLES_H

