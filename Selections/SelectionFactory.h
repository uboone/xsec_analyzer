#include "SelectionBase.h"
#include <string>

class SelectionFactory{
 public:
  SelectionFactory();
  SelectionBase* CreateSelection(std::string SelectionName);
};
