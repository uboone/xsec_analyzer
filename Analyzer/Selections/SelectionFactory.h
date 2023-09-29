#include "CC1mu1p0pi.h"
#include "CC1mu2p0pi.h"
#include "CC1muNp0pi.h"

#include "SelectionBase.h"
#include <string>

class SelectionFactory{
 public:
  SelectionFactory() {};
  SelectionBase* CreateSelection(std::string SelectionName) {

    SelectionBase* Selection;
    if (SelectionName == "CC1mu1p0pi") {
      CC1mu1p0pi* CC1mu1p0piSel = new CC1mu1p0pi();
      Selection = (SelectionBase*)CC1mu1p0piSel;
    } else if (SelectionName == "CC1mu2p0pi") {
      CC1mu2p0pi* CC1mu2p0piSel = new CC1mu2p0pi();
      Selection = (SelectionBase*)CC1mu2p0piSel;
    } else if (SelectionName == "CC1muNp0pi") {
      CC1muNp0pi* CC1muNp0piSel = new CC1muNp0pi();
      Selection = (SelectionBase*)CC1muNp0piSel;
    } else {
      std::cerr << "Selection name requested: " << SelectionName << " is not implemented in " << __FILE__ << std::endl;
      throw;
    }

    return Selection;
  }  
};
