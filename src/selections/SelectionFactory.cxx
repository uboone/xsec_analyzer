
// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/CC1mu1p0pi.hh"
#include "XSecAnalyzer/Selections/CC1mu2p0pi.hh"
#include "XSecAnalyzer/Selections/CC1muNp0pi.hh"
#include "XSecAnalyzer/Selections/CC1muXp0pi.h"
#include "XSecAnalyzer/Selections/CC1mu1p0piXn.hh"
#include "XSecAnalyzer/Selections/CC1muNnXp0pi.h"
#include "XSecAnalyzer/Selections/CC1muNp0piNn.hh"
#include "XSecAnalyzer/Selections/DummySelection.hh"
#include "XSecAnalyzer/Selections/SelectionFactory.hh"

SelectionFactory::SelectionFactory() {
}

SelectionBase* SelectionFactory::CreateSelection(std::string SelectionName) {
  SelectionBase* Selection;
  if ( SelectionName == "CC1mu1p0pi" ) {
    CC1mu1p0pi* CC1mu1p0piSel = new CC1mu1p0pi();
    Selection = (SelectionBase*)CC1mu1p0piSel;
  }
  else if ( SelectionName == "CC1mu2p0pi" ) {
    CC1mu2p0pi* CC1mu2p0piSel = new CC1mu2p0pi();
    Selection = (SelectionBase*)CC1mu2p0piSel;
  }
  else if ( SelectionName == "CC1muNp0pi" ) {
    CC1muNp0pi* CC1muNp0piSel = new CC1muNp0pi();
    Selection = (SelectionBase*)CC1muNp0piSel;
  }
  else if ( SelectionName == "CC1muXp0pi" ) {
    CC1muXp0pi* CC1muXp0piSel = new CC1muXp0pi();
    Selection = (SelectionBase*)CC1muXp0piSel;
  }
  else if ( SelectionName == "CC1muNnXp0pi" ) {
    CC1muNnXp0pi* CC1muNnXp0piSel = new CC1muNnXp0pi();
    Selection = (SelectionBase*)CC1muNnXp0piSel;
  }
  else if ( SelectionName == "CC1mu1p0piXn" ) {
    CC1mu1p0piXn* CC1mu1p0piXnSel = new CC1mu1p0piXn();
    Selection = (SelectionBase*)CC1mu1p0piXnSel;
  }
  else if ( SelectionName == "CC1muNp0piNn" ) {
    CC1muNp0piNn* CC1muNp0piNnSel = new CC1muNp0piNn();
    Selection = (SelectionBase*)CC1muNp0piNnSel;
  }
  else if ( SelectionName == "Dummy" ) {
    DummySelection* DummySel = new DummySelection();
    Selection = (SelectionBase*)DummySel;
  }
  else {
    std::cerr << "Selection name requested: " << SelectionName
      << " is not implemented in " << __FILE__ << '\n';
    throw;
  }

  return Selection;
}
