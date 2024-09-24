// Standard library includes
#include <iostream>

// XSecAnalyzer includes
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/Selections/SelectionBase.hh"

SelectionBase::SelectionBase(std::string fSelectionName_) {
  fSelectionName = fSelectionName_;
  nPassedEvents = 0;

  eventNumber = 0;

  TrueFV = {BOGUS,BOGUS,BOGUS,BOGUS,BOGUS,BOGUS};
  RecoFV = {BOGUS,BOGUS,BOGUS,BOGUS,BOGUS,BOGUS};

  STVTools = STV_Tools();

}

void SelectionBase::Setup(TTree* Tree_, bool Create_) {
  SetupTree(Tree_, Create_);
  DefineCategoryMap();
  DefineConstants();
}

void SelectionBase::ApplySelection(AnalysisEvent* Event) {
  Reset();

  MC_Signal = DefineSignal(Event);
  Selected = Selection(Event);
  EvtCategory = CategorizeEvent(Event);

  ComputeRecoObservables(Event);
  if (Event->is_mc_) {   //Event->is_mc_ is set in CategorizeEvent
    ComputeTrueObservables(Event);
  }

  if (Selected) {
    nPassedEvents++;
  }
  eventNumber++;
}

void SelectionBase::Summary() {
  std::cout << fSelectionName << " has " << nPassedEvents << " events which passed" << std::endl;
}

void SelectionBase::SetupTree(TTree* Tree_, bool Create_) {
  Tree = Tree_;
  Create = Create_;

  std::string BranchName;

  BranchName = "Selected";
  SetBranch(&Selected,BranchName,kBool);

  BranchName = "MC_Signal";
  SetBranch(&MC_Signal,BranchName,kBool);

  BranchName = fSelectionName+"Category";
  SetBranch(&EvtCategory,"EventCategory",kInteger);

  DefineAdditionalInputBranches();
  DefineOutputBranches();
}

void SelectionBase::SetBranch(void* Variable, std::string VariableName, VarType VariableType) {
  SaveVariablePointer(Variable,VariableType);

  VariableName = fSelectionName+"_"+VariableName;
  std::string Leaflist = VariableName;

  switch (VariableType) {
  case kBool:
    Leaflist += "/O";
    break;
  case kDouble:
    Leaflist += "/D";
    break;
  case kFloat:
    Leaflist += "/F";
    break;
  case kInteger:
    Leaflist += "/I";
    break;
  case kTVector:
    //set_object_output_branch_address< TVector >(*Tree,VariableName,Variable,Create);
    break;
  case kSTDVector:
    //set_object_output_branch_address< std::vector<double> >(*Tree,VariableName,Variable,Create);
    break;
  default:
    std::cerr << "Unexpected variable type:" << VariableType << std::endl;
    throw;
  }

  if (Leaflist!="") {
    set_output_branch_address(*Tree,VariableName,Variable,Create,Leaflist);
  } else {
    set_output_branch_address(*Tree,VariableName,Variable,Create);
  }
}

void SelectionBase::SaveVariablePointer(void* Variable, VarType VariableType) {
  switch (VariableType) {
  case kBool:
    Pointer_Bool.push_back((bool*)Variable);
    break;
  case kDouble:
    Pointer_Double.push_back((double*)Variable);
    break;
  case kFloat:
    Pointer_Float.push_back((float*)Variable);
    break;
  case kInteger:
    Pointer_Integer.push_back((int*)Variable);
    break;
  case kTVector:
    Pointer_TVector.push_back((TVector3*)Variable);
    break;
  case kSTDVector:
    Pointer_STDVector.push_back((std::vector<double>*)Variable);
    break;
  default:
    std::cerr << "Unexpected variable type:" << VariableType << std::endl;
    throw;
  }
}

void SelectionBase::Reset() {
  for (size_t i=0;i<Pointer_Bool.size();i++) {
    *(Pointer_Bool[i]) = false;
  }
  for (size_t i=0;i<Pointer_Double.size();i++) {
    *(Pointer_Double[i]) = BOGUS;
  }
  for (size_t i=0;i<Pointer_Float.size();i++) {
    *(Pointer_Float[i]) = BOGUS;
  }
  for (size_t i=0;i<Pointer_Integer.size();i++) {
    *(Pointer_Integer[i]) = BOGUS_INDEX;
  }
  for (size_t i=0;i<Pointer_TVector.size();i++) {
    /*
    for (size_t j=0;j<(*(Pointer_TVector[i])).GetNrows();j++) {
      (*(Pointer_TVector[i]))[j] = 0.;
    }
    */
    *(Pointer_TVector[i]) = TVector3(0,0,0);
  }
  for (size_t i=0;i<Pointer_STDVector.size();i++) {
    (*(Pointer_STDVector[i])).clear();
  }
}
