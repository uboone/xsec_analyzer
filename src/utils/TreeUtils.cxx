/*
Description: Implementation of TTree utility functions
Date: 2024-10-16
Authors: Steven Gardiner and Brinden Carlson
*/

#include "XSecAnalyzer/TreeUtils.hh"

// Helper function to set the branch addresses for the AnalysisEvent
void SetBranchAddress(TTree& etree, std::string BranchName, void* Variable) {
  if (etree.GetBranch(BranchName.c_str()) == nullptr) {
    std::cerr << "Branch " << BranchName << " not found in TTree" << std::endl;
    return;
  }
  if (Variable == nullptr) {
    std::cerr << "Variable is nullptr for branch " << BranchName << std::endl;
    return;
  }
  etree.SetBranchAddress(BranchName.c_str(),Variable);
  //std::cout << "Set branch address for " << BranchName << std::endl;
}