/**
 * @file cuts.h
 * @brief Header file for definitions of analysis cuts.
 * @details This file contains definitions of analysis cuts which can be used
 * to select interactions/slices. Each cut is implemented as a function which takes an
 * interaction/slice object as an argument and returns a boolean. These are the
 * building blocks for defining more complex selections.
 * @author brindenc@fnal.gov
*/
#ifndef SLICE_CUTS_H
#define SLICE_CUTS_H
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "pfp_variables.h"
#include "slice_variables.h"

const FidVol fvnumucc = {0, 190.,
                         -190., 190.,
                         10., 450.};

namespace slicecuts
{   
     /**
     * @brief Apply no cut; all interactions/slices passed.
     * @details This is a placeholder function for a cut which does not apply
     * any selection criteria. It is intended to be used in cases where a cut
     * function is required, but no selection is desired.
     * @tparam T the type of object (true or reco).
     * @param obj the interaction/slice to select on.
     * @return true (always).
     */
    template<class T>
        bool kNoCut(const T & obj) { return true; }

    bool kIsFV(const caf::SRSliceProxy & slc)
    {
        return PtInVolAbsX(slc.vertex, fvnumucc);
    }

    bool kPrecut(const caf::SRSliceProxy & slc)
    {
        return ((slicevars::kNuScore(slc) > 0.5) && kIsFV(slc));
    }

    bool kCosmicCut(const caf::SRSliceProxy & slc)
    {
        return slicevars::kOpt0Score(slc) > 320;
    }

    bool kTrackCut(const caf::SRSliceProxy & slc)
    {
        return slicevars::kNTracks(slc) > 0;
    }

    bool kMuonCut(const caf::SRSliceProxy & slc)
    {
        return slicevars::kLeadingMuonMom(slc) > 0;
    }

    bool kIsSignal(const caf::SRSliceProxy & slc)
    {
        return kPrecut(slc) & kCosmicCut(slc) & kTrackCut(slc) & kMuonCut(slc);
    }


    
}
#endif // SLICE_CUTS_H