#pragma once

// ROOT includes
#include "TF1.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class CC1mu1p0pi : public SelectionBase {

public:

  CC1mu1p0pi();

  virtual int categorize_event( AnalysisEvent* event ) override final;
  virtual void compute_reco_observables( AnalysisEvent* event ) override final;
  virtual void compute_true_observables( AnalysisEvent* event ) override final;
  virtual void define_category_map() override final;
  virtual void define_constants() override final;
  virtual void define_output_branches() override final;
  virtual bool define_signal( AnalysisEvent* event ) override final;
  virtual void reset() override final;
  virtual bool selection( AnalysisEvent* event ) override final;

private:

  bool sel_nslice_eq_1_;
  bool sel_nshower_eq_0_;
  bool sel_ntrack_eq_2_;
  bool sel_muoncandidate_tracklike_;
  bool sel_protoncandidate_tracklike_;
  bool sel_nuvertex_contained_;
  bool sel_muoncandidate_above_p_thresh;
  bool sel_protoncandidate_above_p_thresh;
  bool sel_muoncandidate_contained;
  bool sel_protoncandidate_contained;
  bool sel_muon_momentum_quality;
  bool sel_no_flipped_tracks_;
  bool sel_proton_cand_passed_LLRCut;
  bool sel_muon_momentum_in_range;
  bool sel_muon_costheta_in_range;
  bool sel_muon_phi_in_range;
  bool sel_proton_momentum_in_range;
  bool sel_proton_costheta_in_range;
  bool sel_proton_phi_in_range;

  bool sig_truevertex_in_fv_;
  bool sig_ccnc_;
  bool sig_is_numu_;
  bool sig_one_muon_above_thresh_;
  bool sig_one_proton_above_thresh_;
  bool sig_no_pions_;
  bool sig_no_heavy_mesons_;
  int sig_mc_n_threshold_muon;
  int sig_mc_n_threshold_proton;
  int sig_mc_n_threshold_pion0;
  int sig_mc_n_threshold_pionpm;
  int sig_mc_n_heaviermeson;

  int CandidateMuonIndex;
  int CandidateProtonIndex;

  double Reco_Pt;
  double Reco_Ptx;
  double Reco_Pty;
  double Reco_PL;
  double Reco_Pn;
  double Reco_PnPerp;
  double Reco_PnPerpx;
  double Reco_PnPerpy;
  double Reco_PnPar;
  double Reco_DeltaAlphaT;
  double Reco_DeltaAlpha3Dq;
  double Reco_DeltaAlpha3DMu;
  double Reco_DeltaPhiT;
  double Reco_DeltaPhi3D;
  double Reco_ECal;
  double Reco_EQE;
  double Reco_Q2;
  double Reco_A;
  double Reco_EMiss;
  double Reco_kMiss;
  double Reco_PMiss;
  double Reco_PMissMinus;
  double BackTrack_Pt;
  double BackTrack_Ptx;
  double BackTrack_Pty;
  double BackTrack_PL;
  double BackTrack_Pn;
  double BackTrack_PnPerp;
  double BackTrack_PnPerpx;
  double BackTrack_PnPerpy;
  double BackTrack_PnPar;
  double BackTrack_DeltaAlphaT;
  double BackTrack_DeltaAlpha3Dq;
  double BackTrack_DeltaAlpha3DMu;
  double BackTrack_DeltaPhiT;
  double BackTrack_DeltaPhi3D;
  double BackTrack_ECal;
  double BackTrack_EQE;
  double BackTrack_Q2;
  double BackTrack_A;
  double BackTrack_EMiss;
  double BackTrack_kMiss;
  double BackTrack_PMiss;
  double BackTrack_PMissMinus;
  double True_Pt;
  double True_Ptx;
  double True_Pty;
  double True_PL;
  double True_Pn;
  double True_PnPerp;
  double True_PnPerpx;
  double True_PnPerpy;
  double True_PnPar;
  double True_DeltaAlphaT;
  double True_DeltaAlpha3Dq;
  double True_DeltaAlpha3DMu;
  double True_DeltaPhiT;
  double True_DeltaPhi3D;
  double True_ECal;
  double True_EQE;
  double True_Q2;
  double True_A;
  double True_EMiss;
  double True_kMiss;
  double True_PMiss;
  double True_PMissMinus;

  STVCalcType CalcType;
  TF1* fPP;

  int truemuonindex;
  int trueprotonindex;
};
