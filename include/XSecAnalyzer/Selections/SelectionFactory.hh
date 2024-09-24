#pragma once

// Standard library includes
#include <string>

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"

class SelectionFactory{
 public:
  SelectionFactory();
  SelectionBase* CreateSelection(std::string SelectionName);
};
