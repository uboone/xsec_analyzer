#pragma once

// Standard library includes
#include <map>
#include <set>
#include <string>

// ROOT includes
#include "TTree.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/TreeUtils.hh"

using TreeAndTreeMap = std::pair< TTree*, TreeMap >;

class MyVariantWrapper {

  public:

    MyVariantWrapper( MyVariant* v ) : variant_( v ) {}

    // Overload the = operator to allow setting new values of the variant
    // easily even when working with types that are internally wrapped
    // in a MyPointer.
    template < typename T > MyVariantWrapper& operator=( const T& in ) {
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

    TreeMapWrapper( TreeMap* tm ) : tm_( tm ) {}

    MyVariantWrapper at( const std::string& name ) const {
      MyVariant& var = tm_->at( name );
      return MyVariantWrapper( &var );
    }

    MyVariantWrapper operator[]( const std::string& name ) {
      MyVariant& var = tm_->operator[]( name );
      return MyVariantWrapper( &var );
    }

  protected:

    TreeMap* tm_ = nullptr;
};

// Class that automatically manages storage for TTree branches
class TreeHandler {

  public:

    TreeHandler() {};

    void add_input_tree( TTree* in_tree, const std::string& name = "" );

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

  protected:

    // Helper function for looking up map elements
    TreeAndTreeMap* find_element( const std::string& name,
      bool input = true );

    // Keys are names of the TTrees provided to add_input_tree(),
    // values are pairs in which the first element is a pointer
    // to the original TTree and the second is a TreeMap used
    // to access its branches during event processing.
    std::map< std::string, TreeAndTreeMap > in_tree_maps_;

    // Keys are names of the TTrees provided to add_output_tree(),
    // values are pairs in which the first element is a pointer
    // to the output TTree and the second is a TreeMap used
    // to manage its branches during event processing.
    std::map< std::string, TreeAndTreeMap > out_tree_maps_;

    // Outer keys are names of the TTrees provided to add_input_tree(),
    // inner keys are the leaf names storing sizes of variable-length
    // arrays, values are sets of branch names for which that size is
    // used.
    std::map< std::string, std::map< std::string,
      std::set< std::string > > > in_var_size_map_;
};
