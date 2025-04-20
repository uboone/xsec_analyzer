#pragma once

// Standard library includes
#include <map>
#include <string>

// ROOT includes
#include "TTree.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/TreeUtils.hh"

using TreeAndTreeMap = std::pair< TTree*, TreeMap >;

// Class that automatically manages storage for TTree branches
class TreeHandler {

  public:

    TreeHandler() {};

    void add_input_tree( TTree* in_tree, const std::string& name = "" );

    // Access the owned map
    inline const auto& tree_maps() const { return tree_maps_; }
    inline auto& tree_maps() { return tree_maps_; }

    // Call TTree::GetEntry() for all input trees
    void get_entry( long long entry );

    // Access a TTree from the map with a given name
    TTree* tree( const std::string& name );

    // Access a TreeMap from the map with a given name
    TreeMap& map( const std::string& name );

  protected:

    // Helper function for looking up map elements
    TreeAndTreeMap* find_element( const std::string& name );

    // Keys are names of the TTrees provided to add_input_tree(),
    // values are pairs in which the first element is a pointer
    // to the original TTree and the second is a TreeMap used
    // to access its branches during event processing.
    std::map< std::string, TreeAndTreeMap > tree_maps_;
};
