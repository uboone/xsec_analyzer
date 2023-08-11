#ifndef __SELECTION_BASE_H__
#define __SELECTION_BASE_H__

#include "AnalysisEvent.h"
#include "TTree.h"
#include <string>
#include "Constants.h"

class SelectionBase {
public:
  SelectionBase(std::string fSelectionName_);
  
  void Print();
  void ApplySelection(AnalysisEvent* Event);
  void SetupTree(TTree* Tree_, bool Create_=true);
  
protected:
  void SetBranch(void* Variable, std::string VariableName, VarType VariableType);
  void SaveVariablePointer(void* Variable, VarType VariableType);
  void Reset();
  
  virtual bool Selection(AnalysisEvent* Event) = 0;
  virtual void DefineBranches() = 0;
  virtual void ComputeObservables() = 0;

private:
  std::string fSelectionName;
  int nPassedEvents;
  
  TTree* Tree;
  bool Create;

  std::vector<bool*> Pointer_Bool;
  std::vector<double*> Pointer_Double;
  std::vector<float*> Pointer_Float;
  std::vector<int*> Pointer_Integer;

  bool Selected;
};

#endif
