#pragma once

// Standard library includes
#include <map>
#include <set>
#include <string>

// ROOT includes
#include "TTree.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/TreeMap.hh"

class MyVariantWrapper {

  public:

    MyVariantWrapper( MyVariant* v ) : variant_( v ) {}

    // Overload the = operator to allow setting new values of the variant
    // easily even when working with types that are internally wrapped
    // in a MyPointer.
    template < typename T > MyVariantWrapper& operator=( const T& in ) {

      // Check that the variant either already holds a value of the
      // requested type or has not previously been assigned a value
      // and thus has type std::monostate. If neither of these things are
      // true, then we are attempting to switch types. To prevent
      // issues with inconsistent types in a TTree branch, throw an exception.
      bool is_monostate = std::holds_alternative< std::monostate >( *variant_ );

      bool type_changed = false;
      // For fundamental types, we can use the type directly
      if constexpr ( std::is_fundamental_v< T > ) {
        bool is_T = std::holds_alternative< T >( *variant_ );
        type_changed = !is_T && !is_monostate;
      }
      // For MyPointer types, we need to access the type of the pointed-to
      // object
      else {
        bool is_MyPointer_to_T = std::holds_alternative< MyPointer< T > >(
          *variant_ );
        type_changed = !is_MyPointer_to_T && !is_monostate;
      }

      if ( type_changed ) {
        throw std::runtime_error( "Type switch detected in assignment to"
          " previously initialized MyVariant object" );
      }

      // For simple types, just directly assign to the variant
      if constexpr ( std::is_fundamental_v< T > ) {
        *variant_ = in;
      }
      // Otherwise, assign using a MyPointer to wrap the input type
      else {
        // Check if the active variant is already a MyPointer< T >
        auto* my_ptr = std::get_if< MyPointer< T > >( variant_ );

        // If so, then access the bare pointer inside and do the assignment
        if ( my_ptr ) {
          T* bare_ptr = my_ptr->get();
          if ( !bare_ptr ) {
            throw std::runtime_error( "Invalid dereference in MyVariantWrapper"
              " assignment" );
          }
          *bare_ptr = in;
        }
        // If not, create a new MyPointer externally and move it into the
        // variant, overwriting any prior content
        else {
          MyPointer< T > temp_ptr;
          *temp_ptr = in;
          *variant_ = std::move( temp_ptr );
        }
      }

      return *this;
    }

    // Overload the >> operator to allow easy loading of the active variant
    // into a target variable. The type is inferred from the target without the
    // need for explicit use of a template parameter.
    template < typename T > const MyVariantWrapper& operator>>( T& out ) const
    {
      if constexpr ( std::is_pointer_v< T > ) {
        get_my_variant( *variant_, out );
      }
      else {
        copy_my_variant( *variant_, out );
      }

      return *this;
    }

  protected:

    MyVariant* variant_ = nullptr;
};

class TreeMapWrapper {

  public:

    TreeMapWrapper( TreeMap* tm, const std::string& tree_name )
      : tm_( tm ), tree_name_( tree_name ) {}

    MyVariantWrapper at( const std::string& name ) {
      return this->helper_for_at( name );
    }

    const MyVariantWrapper at( const std::string& name ) const {
      return this->helper_for_at( name );
    }

    MyVariantWrapper operator[]( const std::string& name ) {
      MyVariant& var = tm_->operator[]( name );
      return MyVariantWrapper( &var );
    }

  protected:

    MyVariantWrapper helper_for_at( const std::string& name ) const {
      // Check if we have an existing entry in the map
      MyVariant* var = nullptr;
      auto it = tm_->find( name );
      // If we do, then just retrieve the stored variant
      if ( it != tm_->end() ) var = &it->second;
      else {
        // Set up storage for the input branch and add
        // a corresponding entry to the map
        TBranch* added_branch = nullptr;
        var = tm_->add_input_branch( name, added_branch );
        // Find the owned TTree's current entry number
        long long cur_entry = tm_->get_read_entry();
        // Store the value of the branch in the current
        // entry in the variant
        added_branch->GetEntry( cur_entry );
      }
      return MyVariantWrapper( var );
    }

    TreeMap* tm_ = nullptr;
    const std::string tree_name_;
};

// Class that provides access to a TreeHandler without allowing the user
// to change the current entry or otherwise directly manipulate the
// owned TTrees. The input TTrees are read-only, while a single output TTree
// is available with write access.
class AnalysisEvent {

  friend class TreeHandler;

  public:

    // Access an input TreeMap by name
    const TreeMapWrapper& in( const std::string& tree_name ) const;

    // Access an input TreeMap by index
    const TreeMapWrapper& in( size_t index = 0u ) const;

    // Acces the output TreeMap
    inline TreeMapWrapper& out() { return out_tree_; }

  protected:

    AnalysisEvent( const TreeMapWrapper& out_tmw ) : out_tree_( out_tmw ) {}

    inline void add_input( const std::string& name, const TreeMapWrapper& tmw )
      { in_trees_.emplace( name,  tmw ); }

    std::map< std::string, const TreeMapWrapper > in_trees_;
    TreeMapWrapper out_tree_;
};

// Class that automatically manages storage for TTree branches
class TreeHandler {

  public:

    TreeHandler() {};

    void add_input_tree( TTree* in_tree, const std::string& name = "",
      bool load_all_branches = false );

    MyVariant* add_input_branch( const std::string& in_tree_key,
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
    TreeMapWrapper in_map( const std::string& name );
    TreeMapWrapper out_map( const std::string& name );

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
