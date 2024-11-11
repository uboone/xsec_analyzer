#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class CC1mu2p0pi : public SelectionBase {

public:

  CC1mu2p0pi();

  virtual bool selection(AnalysisEvent* Event) override final;
  virtual int categorize_event(AnalysisEvent* Event) override final;
  virtual void compute_reco_observables(AnalysisEvent* Event) override final;
  virtual void compute_true_observables(AnalysisEvent* Event) override final;
  virtual void define_output_branches() override final;
  virtual bool define_signal(AnalysisEvent* Event) override final;
  virtual void define_constants() override final;
  virtual void define_category_map() override final;
  virtual void reset() override final;

private:

  bool sel_reco_vertex_in_FV_;
  bool sel_has_muon_candidate_;
  bool sel_nu_mu_cc_;
  bool sel_npfps_eq_3;
  bool sel_ntracks_eq_3;
  bool sel_containedparticles;
  bool sel_correctparticles;
  bool sel_momentum_threshold_passed_;

  bool sig_ccnc_;
  bool sig_is_numu_;
  bool sig_two_protons_above_thresh_;
  bool sig_one_muon_above_thres_;
  bool sig_no_pions_;
  bool sig_truevertex_in_fv_;
  int sig_mc_n_threshold_muon;
  int sig_mc_n_threshold_proton;
  int sig_mc_n_threshold_pion0;
  int sig_mc_n_threshold_pionpm;

  int muon_candidate_idx_;
  int LeadingProtonIndex;
  int RecoilProtonIndex;

  double Reco_CosPlPr;
  double Reco_CosMuPsum;

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

  int TrueLeadingProtonIndex;
  int TrueRecoilProtonIndex;
  int TrueMuonIndex;

  STVCalcType CalcType;
};
