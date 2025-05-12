#pragma once

// Standard library includes
#include <map>
#include <set>
#include <string>

// ROOT includes
#include "TTree.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/TreeMap.hh"

// Class that provides access to a TreeHandler without allowing the user
// to change the current entry or otherwise directly manipulate the
// owned TTrees. The input TTrees are read-only, while a single output TTree
// is available with write access.
class AnalysisEvent {

  friend class TreeHandler;

  public:

    // Access an input TreeMap by name
    const TreeMap::Wrapper& in( const std::string& tree_name ) const;

    // Access an input TreeMap by index
    const TreeMap::Wrapper& in( size_t index = 0u ) const;

    // Acces the output TreeMap
    inline TreeMap::Wrapper& out() { return out_tree_; }

  protected:

    AnalysisEvent( const TreeMap::Wrapper& out_tmw ) : out_tree_( out_tmw ) {}

    inline void add_input( const std::string& name,
      const TreeMap::Wrapper& tmw )
    {
      in_trees_.emplace( name,  tmw );
    }

    std::map< std::string, const TreeMap::Wrapper > in_trees_;
    TreeMap::Wrapper out_tree_;
};

// Class that automatically manages storage for TTree branches
class TreeHandler {

  public:

    TreeHandler() {};

    void add_input_tree( TTree* in_tree, const std::string& name = "",
      bool load_all_branches = false );

    TreeMap::Variant* add_input_branch( const std::string& in_tree_key,
      const std::string& branch_name );

    void add_output_tree( TDirectory* out_dir, const std::string& name,
      const std::string& title = "" );

    // Access the owned maps
    inline const auto& input_tree_maps() const { return in_tree_maps_; }
    inline auto& input_tree_maps() { return in_tree_maps_; }

    inline const auto& output_tree_maps() const { return out_tree_maps_; }
    inline auto& output_tree_maps() { return out_tree_maps_; }

    // Call TTree::GetEntry() for all input trees
    void get_entry( long long entry );

    // Call TTree::LoadTree() for all input trees
    long long load_tree( long long entry );

    // Call TTree::Fill() for all output trees
    void fill();

    // Call TTree::Write() for all output trees
    void write();

    // Access an TTree from the maps with a given name
    TTree* in_tree( const std::string& name );
    TTree* out_tree( const std::string& name );

    // Access a TreeMap from the input map with a given name
    // using a wrapper class to allow us to overload
    // some unary operators
    TreeMap::Wrapper in_map( const std::string& name );
    TreeMap::Wrapper out_map( const std::string& name );

    // Get direct access to a TreeMap by name
    TreeMap& access_in_map( const std::string& name );
    TreeMap& access_out_map( const std::string& name );

    AnalysisEvent get_event( const std::set< std::string >& input_names,
      const std::string& output_name );

  protected:

    // Helper function for looking up map elements
    TreeMap* find_element( const std::string& name,
      bool input = true );

    // Keys are names of the TTrees provided to add_input_tree(), values are
    // TreeMap objects used to access the TTree branches during event
    // processing.
    std::map< std::string, TreeMap > in_tree_maps_;

    // Keys are names of the TTrees provided to add_output_tree(), values are
    // TreeMap objects used to access the TTree branches during event
    // processing.
    std::map< std::string, TreeMap > out_tree_maps_;

    // Outer keys are names of the TTrees provided to add_input_tree(),
    // inner keys are the leaf names storing sizes of variable-length
    // arrays, values are sets of branch names for which that size is
    // used.
    std::map< std::string, std::map< std::string,
      std::set< std::string > > > in_var_size_map_;

    // Output TTrees become "locked" after the first call to fill(),
    // ensuring that the branch structure remains the same for all entries.
    // This flag indicates whether the output TTrees are currently locked or
    // not.
    bool output_locked_ = false;
};
