#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class NuMICC1e : public SelectionBase {

    public:

    NuMICC1e();

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

    bool sig_isNuE_;
    bool sig_isCC_;
    bool sig_inFV_;
    bool sig_has_fs_electron_;
    bool sig_is_signal_;

    bool sel_pass_preselection_;
    bool sel_pass_cosmic_rejection_;
    bool sel_pass_shower_identification_;
    bool sel_pass_electron_identification_;
    bool sel_nu_e_cc_;

    // observables
    float reco_electron_energy_;
    float mc_electron_energy_;
};