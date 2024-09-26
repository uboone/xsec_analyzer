#pragma once

// Standard library includes
#include <algorithm>
#include <array>
#include <map>
#include <string>

// ROOT includes
#include "TTree.h"

// STV analysis includes
#include "TreeUtils.hh"

// Class that provides temporary storage for event weights being processed by a
// UniverseMaker object
class WeightHandler {
  public:

    WeightHandler() {};

    // Configure storage and set input branch addresses for processing event
    // weights from the input TTree. If a non-null pointer to a vector of
    // branch names is supplied, use it to decide which branches to include.
    // Otherwise, auto-detect appropriate branches using the TTree itself.
    void set_branch_addresses( TTree& in_tree,
      const std::vector<std::string>* branch_names = nullptr );

    // Add a single branch from the input TTree with the specified name.
    // This function also sets the branch address appropriately.
    void add_branch( TTree& in_tree, const std::string& branch_name,
      bool throw_when_missing = true );

    // Overloaded version that allows for easy configuration of a single input
    // branch
    void set_branch_addresses( TTree& in_tree, const std::string& branch_name );

    // Access the owned map
    inline const auto& weight_map() const { return weight_map_; }
    inline auto& weight_map() { return weight_map_; }

  protected:

    // Keys are branch names in the input TTree, values point to vectors of
    // event weights
    std::map< std::string, MyPointer< std::vector<double> > > weight_map_;
};
