/**
 * @file  xsec_analyzer/pandora/definitions.h
  * @brief Definitions of analysis variables.
  * @details This file contains definitions of analysis variables which can be
  * used to extract information from slices. Each variable is implemented
  * as a function which takes a slice object as an argument and returns a
  * double. These variables are intended to define slice definitions.
  * @author brindenc@fnal.gov
 */

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/SBNAna/Cuts/TruthCuts.h"

#include "slice_cuts.h"

using namespace ana;

namespace defs
{
    int kEventType(const caf::SRSliceProxy & slc)
    {
        if (kTrueActiveVolumeND(&slc) && kIsNumu(&slc) && kIsCC(&slc)) return 0; //numucc
        else if (kTrueActiveVolumeND(&slc) & kIsNC(&slc)) return 1; //NC
        else if (kTrueActiveVolumeND(&slc) & kIsNue(&slc) & kIsCC(&slc)) return 2; //nuecc
        else if (!kIsNu(&slc)) return 3; //Cosmic
        else if (!kTrueActiveVolumeND(&slc) & kIsNu(&slc)) return 4; //Out of volume
        else return -1; //Other 
    }
}
#endif // DEFINITIONS_H