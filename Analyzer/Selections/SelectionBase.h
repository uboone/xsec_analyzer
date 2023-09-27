#ifndef __SELECTION_BASE_H__
#define __SELECTION_BASE_H__

#include "AnalysisEvent.h"
#include "TTree.h"
#include <string>
#include "Constants.h"
#include "TVector.h"

class SelectionBase {
public:
  SelectionBase(std::string fSelectionName_);

  void Setup(TTree* Tree_, bool Create_=true);
  void ApplySelection(AnalysisEvent* Event);
  void Summary();

  bool IsMCSignal() {return MC_Signal;}
  
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
  
  virtual bool Selection(AnalysisEvent* Event) = 0;
  virtual void ComputeRecoObservables(AnalysisEvent* Event) = 0;
  virtual void ComputeTrueObservables(AnalysisEvent* Event) = 0;
  virtual void DefineBranches() = 0;
  virtual bool DefineSignal(AnalysisEvent* Event) = 0;
  virtual void DefineConstants() = 0;  

  FiducialVolume TrueFV;
  FiducialVolume RecoFV;

  TTree* Tree;
  bool Create;
  
private:
  std::string fSelectionName;
  int nPassedEvents;

  std::vector<bool*> Pointer_Bool;
  std::vector<double*> Pointer_Double;
  std::vector<float*> Pointer_Float;
  std::vector<int*> Pointer_Integer;
  std::vector<TVector*> Pointer_TVector;
  std::vector<std::vector<double>*> Pointer_STDVector;

  EventCategory EventCategory;
  
  bool Selected;
  bool MC_Signal;

  int eventNumber;
};

#endif
