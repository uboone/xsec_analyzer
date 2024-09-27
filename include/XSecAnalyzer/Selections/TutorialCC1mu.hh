#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class TutorialCC1mu : public SelectionBase {

  public:

    TutorialCC1mu();

    int CategorizeEvent(AnalysisEvent* Event) override;
    bool Selection(AnalysisEvent* Event) override;
    bool DefineSignal(AnalysisEvent* Event) override;
    void ComputeRecoObservables(AnalysisEvent* Event) override;
    void ComputeTrueObservables(AnalysisEvent* Event) override;
    void DefineOutputBranches() override;
    void DefineConstants() override;
    void DefineCategoryMap() override;

  private:

    bool sig_isNuMu_;
    bool sig_inFV_;
    bool sig_has_fs_muon_;
    bool sig_is_signal_;

    bool sel_reco_vertex_in_FV_;
    bool sel_pfp_starts_in_PCV_;
    bool sel_has_muon_candidate_;
    bool sel_topo_cut_passed_;
    bool sel_nu_mu_cc_;
    bool sel_muon_contained_;

    int muon_candidate_idx_;

    MyPointer< TVector3 > reco_p3mu_;

    MyPointer< TVector3 > mc_p3mu_;

};
