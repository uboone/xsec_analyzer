#ifndef __SELECTION_BASE_H__
#define __SELECTION_BASE_H__

#include "AnalysisEvent.h"
#include "TTree.h"
#include <string>
#include "Constants.h"
#include "TVector3.h"

#include "STV_Tools.h"

#include "FiducialVolume.hh"

class SelectionBase {
public:
  SelectionBase(std::string fSelectionName_);

  void Setup(TTree* Tree_, bool Create_=true);
  void ApplySelection(AnalysisEvent* Event);
  void Summary();

  bool IsEventMCSignal() {return MC_Signal;}
  bool IsEventSelected() {return Selected;}
  
protected:
  void SetBranch(void* Variable, std::string VariableName, VarType VariableType);

  template <typename T> void SetBranch(MyPointer<T>& u_ptr, std::string VariableName, VarType VariableType) {
    T*& address = u_ptr.get_bare_ptr();

    VariableName = fSelectionName+"_"+VariableName;
    SaveVariablePointer(address,VariableType);
    set_object_output_branch_address<T>(*Tree,VariableName,address,Create);
  }
    
  void SaveVariablePointer(void* Variable, VarType VariableType);
  void SetupTree(TTree* Tree_, bool Create_=true);
  void Reset();
  int GetEventNumber() {return eventNumber;}

  inline void DefineTrueFV(double XMin, double XMax, double YMin, double YMax, double ZMin, double ZMax) {
    TrueFV = {XMin, XMax, YMin, YMax, ZMin, ZMax};
    TrueFVSet = true;
  }
  inline FiducialVolume ReturnTrueFV() {
    if (!TrueFVSet) {std::cerr << "True Fiducial volume has not been defined for selection:" << fSelectionName << std::endl; throw;}
    return TrueFV;
  }

  inline void DefineRecoFV(double XMin, double XMax, double YMin, double YMax, double ZMin, double ZMax) {
    RecoFV = {XMin, XMax, YMin, YMax, ZMin, ZMax};
    RecoFVSet = true;
  }
  inline FiducialVolume ReturnRecoFV() {
    if (!RecoFVSet) {std::cerr << "Reco Fiducial volume has not been defined for selection:" << fSelectionName << std::endl; throw;}
    return RecoFV;
  }
  
  virtual bool Selection(AnalysisEvent* Event) = 0;
  virtual EventCategory CategorizeEvent(AnalysisEvent* Event) = 0;
  virtual void ComputeRecoObservables(AnalysisEvent* Event) = 0;
  virtual void ComputeTrueObservables(AnalysisEvent* Event) = 0;
  virtual void DefineOutputBranches() = 0;
  virtual bool DefineSignal(AnalysisEvent* Event) = 0;
  virtual void DefineConstants() = 0;
  void DefineAdditionalInputBranches() {};

  TTree* Tree;
  bool Create;

  STV_Tools STVTools;
  
private:
  std::string fSelectionName;
  int nPassedEvents;

  std::vector<bool*> Pointer_Bool;
  std::vector<double*> Pointer_Double;
  std::vector<float*> Pointer_Float;
  std::vector<int*> Pointer_Integer;
  std::vector<TVector3*> Pointer_TVector;
  std::vector<std::vector<double>*> Pointer_STDVector;

  EventCategory EvtCategory;
  
  bool Selected;
  bool MC_Signal;

  FiducialVolume TrueFV;
  FiducialVolume RecoFV;
  bool TrueFVSet = false;
  bool RecoFVSet = false;

  int eventNumber;
};

#endif
