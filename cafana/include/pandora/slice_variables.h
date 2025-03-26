/**
 * @file variables.h
 * @brief Header file for definitions of analysis variables.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from interactions. Each variable is implemented
 * as a function which takes an object as an argument and returns a
 * double. Building blocks for constructing the xsec_analyzer tree.
 * @author brindenc@fnal.gov
*/

#ifndef SLICE_VARIABLES_H
#define SLICE_VARIABLES_H
#define ELECTRON_MASS 0.5109989461
#define MUON_MASS 105.6583745
#define PION_MASS 139.57039
#define PROTON_MASS 938.2720813

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRSlice.h"
#include "sbnanaobj/StandardRecord/SRSliceRecoBranch.h"
#include "sbnanaobj/StandardRecord/Proxy/EpilogFwd.h"

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/SBNAna/Vars/TruthVars.h"

#include "pfp_variables.h"
//#include "slice_cuts.h"

#include <vector>
#include <iostream>

/**
 * @namespace slicevars
 * @brief Namespace for organizing generic variables which act on interactions and
 * slices.
 * @details This namespace is intended to be used for organizing variables which
 * act on interactions and slices. Each variable is implemented as a function which takes
 * an object as an argument and returns a double. The function
 * should be templated on the type of object if the variable is
 * intended to be used on both true and reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the
 * pfpvars namespace, which is used for organizing variables which act on single
 * particles and pfps.
 */
namespace slicevars
{
    int kSliceIndex(const caf::SRSliceProxy & slc){
        return slc.self;
    }
    
    /**
    * @brief Finds the index corresponding to the leading particle of the specified
    * particle type.
    * @details The leading particle is defined as the particle with the highest
    * kinetic energy. If the slice is a true slice, the initial kinetic
    * energy is used instead of the CSDA kinetic energy.
    * @tparam T the type of slice (true or reco).
    * @param obj the slice to operate on.
    * @param pgd of the particle type.
    * @return the index of the leading particle (highest KE).
    */
    size_t kLeadingTrackIndex(const caf::SRSliceProxy & slc, uint16_t pdg)
    {
        double leading_ke(0);
        size_t index(-1);
        for(size_t i(0); i < slc.reco.npfp; ++i)
        {
            const auto & p = slc.reco.pfp[i];
            const int _pdg = pfpvars::kTrackPDG(p);
            const double trk_score = pfpvars::kPFPTrackScore(p);
            if(_pdg == pdg && trk_score > 0.5)
            {
                const double energy = pfpvars::kPFPKineticEnergy(p);
                if(energy > leading_ke)
                {
                    leading_ke = energy;
                    index = i;
                }
            }
        }
        return index;
    }

    /**
     * @brief Score of how neutrino-like the slice is according to pandora
     * @details This variable is based on the neutrino score assigned by pandora.
     * @param obj the slice to apply the variable on.
     * @return the neutrino score.
     */
    double kNuScore(const caf::SRSliceProxy & slc){ return slc.nu_score; };

    double kOpt0Score(const caf::SRSliceProxy & slc){ return slc.opt0.score; };

    double kNuVx(const caf::SRSliceProxy & slc){return slc.vertex.x;};
    double kNuVy(const caf::SRSliceProxy & slc){return slc.vertex.y;};
    double kNuVz(const caf::SRSliceProxy & slc){return slc.vertex.z;};

    double kNPFP(const caf::SRSliceProxy & slc){return slc.reco.npfp;};

    double kNTracks(const caf::SRSliceProxy & slc)
    {
        int n_tracks = 0;
        for(size_t i(0); i < slc.reco.npfp; ++i)
        {
            const auto & p = slc.reco.pfp[i];
            if(pfpvars::kPFPTrackScore(p) > 0.6) n_tracks++;
        }
        return n_tracks;
    }

    // Leading muon variables
    double kLeadingMuonMom(const caf::SRSliceProxy & slc)
    {
        const size_t index = kLeadingTrackIndex(slc, 13);
        if(index == -1) return 0;

        const auto & p = slc.reco.pfp[index];
        return pfpvars::kTrackMomMuon(p);
    }

    // Leading muon costheta
    double kLeadingMuonCostheta(const caf::SRSliceProxy & slc)
    {
        const size_t index = kLeadingTrackIndex(slc, 13);
        if(index == -1) return 9999;

        const auto & p = slc.reco.pfp[index];
        return pfpvars::kTrackDirz(p);
    }

    // True leading muon variables
    double kTrueLeadingMuonMom(const caf::SRSliceProxy & slc)
    {
        double mom(0);
        for(size_t i(0); i < slc.truth.prim.size(); ++i)
        {
            const auto & p = slc.truth.prim[i];
            if(std::abs(p.pdg) == 13)
            {
                double _mom = TMath::Sqrt(TMath::Power(p.genp.x, 2) + TMath::Power(p.genp.y, 2) + TMath::Power(p.genp.z, 2));
                if(_mom > mom) mom = _mom;
            }
        }
        return mom;
    }

    // True leading muon costheta
    double kTrueLeadingMuonCostheta(const caf::SRSliceProxy & slc)
    {
        double costheta(9999);
        for(size_t i(0); i < slc.truth.prim.size(); ++i)
        {
            const auto & p = slc.truth.prim[i];
            if(std::abs(p.pdg) == 13)
            {
                costheta = p.genp.z / TMath::Sqrt(TMath::Power(p.genp.x, 2) + TMath::Power(p.genp.y, 2) + TMath::Power(p.genp.z, 2));
            }
        }
        return costheta;
    }

} // namespace vars
#endif // SLICE_VARIABLES_H