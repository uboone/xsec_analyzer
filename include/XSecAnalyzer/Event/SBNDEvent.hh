/*
Description: Implementation of AnalysisEvent class for SBND
Date: 2024-10-16
Authors: Brinden Carlson (bcarlson1@ufl.edu)
*/

#include "AnalysisEvent.hh"
#include "XSecAnalyzer/SBND/Constants.hh"

#include "TLeaf.h"

class SBNDEvent : public AnalysisEvent {
public:
  SBNDEvent() {}
  ~SBNDEvent() {}

  // Declare member variables for SBND
  //FIXME: Load MC from CAFAna even though it's shorter than number of slices
  bool is_mc_ = true; //default to MC

  // Flux systematic weights
  MyPointer<std::vector<double>> expskin_Flux_;
  MyPointer<std::vector<double>> horncurrent_Flux_;
  MyPointer<std::vector<double>> nucleoninexsec_Flux_;
  MyPointer<std::vector<double>> nucleonqexsec_Flux_;
  MyPointer<std::vector<double>> nucleontotxsec_Flux_;
  MyPointer<std::vector<double>> pioninexsec_Flux_;
  MyPointer<std::vector<double>> pionqexsec_Flux_;
  MyPointer<std::vector<double>> piontotxsec_Flux_;
  MyPointer<std::vector<double>> piplus_Flux_;
  MyPointer<std::vector<double>> piminus_Flux_;
  MyPointer<std::vector<double>> kplus_Flux_;
  MyPointer<std::vector<double>> kminus_Flux_;
  MyPointer<std::vector<double>> kzero_Flux_;

  // xsec systematic weights
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_RPA_CCQE_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_CoulombCCQE_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NormCCMEC_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NormNCMEC_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NCELVariationResponse_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_COHVariationResponse_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse_;
  MyPointer<std::vector<double>> GENIEReWeight_SBN_v1_multisim_FSI_N_VariationRespons_;

  // Helper function to set branch addresses for reading information
  // from the Event TTree
  void set_event_branch_addresses( TTree& etree) override
  {
    //Header info
    //SetBranchAddress(etree, "is_mc", &this->is_mc_);

    //SetBranchAddress(etree, "nu_score", &this->topological_score_ );
    //Muon selection variables
    SetBranchAddress(etree, "leading_muon_costheta", &this->leading_muon_costheta_);
    SetBranchAddress(etree, "leading_muon_momentum", &this->leading_muon_momentum_);
    SetBranchAddress(etree, "is_signal", &this->is_signal_);
    
    // Set the branch addresses for the truth information
    SetBranchAddress(etree, "true_leading_muon_costheta", &this->leading_muon_costheta_truth_);
    SetBranchAddress(etree, "true_leading_muon_momentum", &this->leading_muon_momentum_truth_);
    SetBranchAddress(etree, "event_type", &this->event_type_);

    bool has_flux_weights = etree.GetBranch("expskin_Flux") != nullptr;
    //std::cout << "has_flux_weights: " << has_flux_weights << std::endl;
    if ( has_flux_weights ) {
      //std::cout << "Setting branch address for multisim weights" << std::endl;
      set_object_input_branch_address( etree, "expskin_Flux", this->expskin_Flux_ );
      set_object_input_branch_address( etree, "horncurrent_Flux", this->horncurrent_Flux_ );
      set_object_input_branch_address( etree, "nucleoninexsec_Flux", this->nucleoninexsec_Flux_ );
      set_object_input_branch_address( etree, "nucleonqexsec_Flux", this->nucleonqexsec_Flux_ );
      set_object_input_branch_address( etree, "nucleontotxsec_Flux", this->nucleontotxsec_Flux_ );
      set_object_input_branch_address( etree, "pioninexsec_Flux", this->pioninexsec_Flux_ );
      set_object_input_branch_address( etree, "pionqexsec_Flux", this->pionqexsec_Flux_ );
      set_object_input_branch_address( etree, "piontotxsec_Flux", this->piontotxsec_Flux_ );
      set_object_input_branch_address( etree, "piplus_Flux", this->piplus_Flux_ );
      set_object_input_branch_address( etree, "piminus_Flux", this->piminus_Flux_ );
      set_object_input_branch_address( etree, "kplus_Flux", this->kplus_Flux_ );
      set_object_input_branch_address( etree, "kminus_Flux", this->kminus_Flux_ );
      set_object_input_branch_address( etree, "kzero_Flux", this->kzero_Flux_ );
    }

    bool has_xsec_weights = etree.GetBranch("GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse") != nullptr;
    //std::cout << "has_xsec_weights: " << has_xsec_weights << std::endl;
    if (has_xsec_weights) {
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse", this->GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_RPA_CCQE", this->GENIEReWeight_SBN_v1_multisim_RPA_CCQE_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_CoulombCCQE", this->GENIEReWeight_SBN_v1_multisim_CoulombCCQE_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NormCCMEC", this->GENIEReWeight_SBN_v1_multisim_NormCCMEC_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NormNCMEC", this->GENIEReWeight_SBN_v1_multisim_NormNCMEC_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NCELVariationResponse", this->GENIEReWeight_SBN_v1_multisim_NCELVariationResponse_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse", this->GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse", this->GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_COHVariationResponse", this->GENIEReWeight_SBN_v1_multisim_COHVariationResponse_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse", this->GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse", this->GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse_);
      set_object_input_branch_address(etree, "GENIEReWeight_SBN_v1_multisim_FSI_N_VariationRespons", this->GENIEReWeight_SBN_v1_multisim_FSI_N_VariationRespons_);
    }

    // Set branch addresses for spline and tuned CV weights
    bool has_genie_mc_weights = etree.GetBranch("weightSpline") != nullptr;
    if ( has_genie_mc_weights ) {
      SetBranchAddress(etree, "weightSpline", &this->spline_weight_ );
      SetBranchAddress(etree, "weightTune", &this->tuned_cv_weight_ );
    }

  }

  // Helper function to set branch addresses for the output TTree
  void set_event_output_branch_addresses(TTree &out_tree, bool create = false) override
  {
    // Header info
    set_output_branch_address( out_tree, "is_mc", &this->is_mc_, create, "is_mc/O" );

    //Muon selection variables
    set_output_branch_address( out_tree, "leading_muon_costheta", &this->leading_muon_costheta_, create, "leading_muon_costheta/D" );
    set_output_branch_address( out_tree, "leading_muon_momentum", &this->leading_muon_momentum_, create, "leading_muon_momentum/D" );
    set_output_branch_address( out_tree, "is_signal", &this->is_signal_, create, "is_signal/D" );

    // Set the branch addresses for the truth information
    set_output_branch_address( out_tree, "true_leading_muon_costheta", &this->leading_muon_costheta_truth_, create, "true_leading_muon_costheta/D" );
    set_output_branch_address( out_tree, "true_leading_muon_momentum", &this->leading_muon_momentum_truth_, create, "true_leading_muon_momentum/D" );
    set_output_branch_address( out_tree, "event_type", &this->event_type_, create, "event_type/D" );

    // Set the branch addresses for the flux weights
    if ( !this->expskin_Flux_.empty() ) {
      set_object_output_branch_address< std::vector<double> >( out_tree, "weight_expskin_Flux", this->expskin_Flux_, create );
      set_object_output_branch_address< std::vector<double> >( out_tree, "weight_horncurrent_Flux", this->horncurrent_Flux_, create );
      set_object_output_branch_address< std::vector<double> >( out_tree, "weight_nucleoninexsec_Flux", this->nucleoninexsec_Flux_, create );
      set_object_output_branch_address< std::vector<double> >( out_tree, "weight_nucleonqexsec_Flux", this->nucleonqexsec_Flux_, create );
      set_object_output_branch_address< std::vector<double> >( out_tree, "weight_nucleontotxsec_Flux", this->nucleontotxsec_Flux_, create );
      set_object_output_branch_address< std::vector<double> >( out_tree, "weight_pioninexsec_Flux", this->pioninexsec_Flux_, create );
      set_object_output_branch_address< std::vector<double> >( out_tree, "weight_pionqexsec_Flux", this->pionqexsec_Flux_, create );
      set_object_output_branch_address< std::vector<double> >( out_tree, "weight_piontotxsec_Flux", this->piontotxsec_Flux_, create );
      set_object_output_branch_address< std::vector<double> >( out_tree, "weight_piplus_Flux", this->piplus_Flux_, create );
      set_object_output_branch_address< std::vector<double> >( out_tree, "weight_piminus_Flux", this->piminus_Flux_, create );
      set_object_output_branch_address< std::vector<double> >( out_tree, "weight_kplus_Flux", this->kplus_Flux_, create );
      set_object_output_branch_address< std::vector<double> >( out_tree, "weight_kminus_Flux", this->kminus_Flux_, create );
      set_object_output_branch_address< std::vector<double> >( out_tree, "weight_kzero_Flux", this->kzero_Flux_, create );
    }

    // Set the branch addresses for the xsec weights
    if ( !this->GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse_.empty()){
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse", this->GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE", this->GENIEReWeight_SBN_v1_multisim_RPA_CCQE_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE", this->GENIEReWeight_SBN_v1_multisim_CoulombCCQE_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC", this->GENIEReWeight_SBN_v1_multisim_NormCCMEC_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC", this->GENIEReWeight_SBN_v1_multisim_NormNCMEC_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse", this->GENIEReWeight_SBN_v1_multisim_NCELVariationResponse_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse", this->GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse", this->GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi", this->GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse", this->GENIEReWeight_SBN_v1_multisim_COHVariationResponse_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse", this->GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse", this->GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse_, create);
      set_object_output_branch_address<std::vector<double>>(out_tree, "weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationRespons", this->GENIEReWeight_SBN_v1_multisim_FSI_N_VariationRespons_, create);
    }

    // Set the branch addresses for the spline and tuned CV weights
    set_output_branch_address( out_tree, "spline_weight", &this->spline_weight_, create, "spline_weight/F" );
    set_output_branch_address( out_tree, "tuned_cv_weight", &this->tuned_cv_weight_, create, "tuned_cv_weight/F" );

  }
};