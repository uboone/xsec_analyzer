// ROOT includes
#include "TH1.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/EventCategoriesNeutron.hh"

std::map< int, std::pair< std::string, int > > CC1muNnXp_MAP = {
  { kUnknown, { "Unknown", kGray } },
  { kNuMuCC1n0pi0p_CCQE, { "CCmu1n0pi0p (CCQE)", kBlue - 2 } },		//  numuCC1n0pi0p
  { kNuMuCC1n0pi0p_CCMEC, { "CCmu1n0pi0p (CCMEC)", kBlue - 6 } },	// '		 '
  { kNuMuCC1n0pi0p_CCRES, { "CCmu1n0pi0p (CCRES)", kBlue - 9 } },	// '             '
  { kNuMuCC1n0pi0p_Other, { "CCmu1n0pi0p (Other)", kBlue - 10 } },	// '             '
  { kNuMuCCNn0pi0p_CCQE, { "CCmuNn0pi0p (CCQE)", kCyan - 3 } },		//  numuCCNn0pi0p
  { kNuMuCCNn0pi0p_CCMEC, { "CCmuNn0pi0p (CCMEC)", kCyan - 6 } },	// '             '
  { kNuMuCCNn0pi0p_CCRES, { "CCmuNn0pi0p (CCRES)", kCyan - 4 } },	// '             '
  { kNuMuCCNn0pi0p_Other, { "CCmuNn0pi0p (Other)", kCyan - 9 } },	// '             '
  { kNuMuCC1n0piNp_CCQE, { "CCmu1n0piNp (CCQE)",kGreen - 2 } },		//  numuCC1n0piNp
  { kNuMuCC1n0piNp_CCMEC, { "CCmu1n0piNp (CCMEC)",kGreen - 6 } },	// '             '
  { kNuMuCC1n0piNp_CCRES, { "CCmu1n0piNp (CCRES)",kGreen - 9 } },	// '             '
  { kNuMuCC1n0piNp_Other, { "CCmu1n0piNp (Other)",kGreen - 10 } },	// '             '
  { kNuMuCCNn0piNp_CCQE, { "CCmuNn0piNp (CCQE)",kYellow - 3 } },	//  numuCCNn0piNp
  { kNuMuCCNn0piNp_CCMEC, { "CCmuNn0piNp (CCMEC)",kYellow - 6 } },	// '             '
  { kNuMuCCNn0piNp_CCRES, { "CCmuNn0piNp (CCRES)",kYellow - 4 } },	// '             '
  { kNuMuCCNn0piNp_Other, { "CCmuNn0piNp (Other)",kYellow - 9 } },	// '             '
  { kNuMuCC0n0piNp, {"CCmu0n0piNp", kOrange + 4 } },			//  numuCC0n0piNp
  { kNuMuCCNpi, { "#nu_{#mu} CCN#pi", kOrange + 2 } },			//  numuCCNpi
  { kNuMuCCOther, { "Other #nu_{#mu} CC", kOrange - 1 } },		//  numuCCother
  { kNuECC, { "#nu_{e} CC", kViolet } },				//  nuECC
  { kNC, { "NC", kOrange } },						//  Neutral Current
  { kOOFV, {"Out FV", kRed + 3 } },					//  Out of fiducial volume
  { kOther, { "Other", kRed + 1 } }					//  Other (probably cosmic)
};
