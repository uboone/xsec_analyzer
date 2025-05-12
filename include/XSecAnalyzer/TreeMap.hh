#pragma once

// Standard library includes
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include <vector>

// ROOT includes
#include "TTree.h"
#include "TVector3.h"

// A std::unique_ptr with redundant storage of a bare pointer to the managed
// object. This is a hacky workaround for setting TTree branch addresses for
// objects managed by a std::unique_ptr.
template < typename T > class MyPointer : public std::unique_ptr< T > {
  public:

    MyPointer() : std::unique_ptr< T >( new T ) {}

    T*& get_bare_ptr() {
      bare_ptr_ = this->get();
      return bare_ptr_;
    }

  protected:

    T* bare_ptr_ = nullptr;
};

// Create a type trait to allow "constexpr if" to distinguish between
// instantiations of the MyPointer class template and other types
template < typename T >
  struct is_MyPointer : std::false_type {};

template < typename T >
  struct is_MyPointer< MyPointer< T > > : std::true_type {};

template < typename T >
  constexpr bool is_MyPointer_v = is_MyPointer< T >::value;

// Create a type trait to allow "constexpr if" to distinguish between
// MyPointer< std::vector< T > > and other types
template < typename T >
  struct is_MyPointerToVector : std::false_type {};

template < typename T >
  struct is_MyPointerToVector< MyPointer< std::vector< T > > >
    : std::true_type {};

template < typename T >
  constexpr bool is_MyPointerToVector_v = is_MyPointerToVector< T >::value;

// Manages a map used to automatically allocate storage for managing TTree
// branch variables using TreeMap::Variant objects. Keys in this map are strings
// giving the branch names (no duplicates allowed). Values in the map are
// TreeMap::Variant objects that will be used to store the branch data.
class TreeMap {
  public:

    // A std::variant that can store any of the branch types of interest for
    // the TTrees managed by a TreeMap. New branches with a different type
    // will require extending the template arguments given here.
    using Variant = std::variant<
      // ** IMPORTANT: Keep std::monostate first in the template arguments **
      // This is needed to ensure that default construction of Variant leads
      // to it holding the std::monostate type, which is used as a marker
      // that another type has not been assigned to the variant.
      std::monostate,
      // Now list the fundamental types that are allowed here before moving
      // on to objects
      unsigned char, char, unsigned short, short, unsigned int, int,
      unsigned long, long, unsigned long long, long long, bool, float, double,
      // Classes, structs, etc. that are allowed should be wrapped in
      // a MyPointer here
      MyPointer< std::string >, MyPointer< TVector3 >,
      MyPointer< std::vector< bool > >,
      MyPointer< std::vector< unsigned char > >,
      MyPointer< std::vector< char > >,
      MyPointer< std::vector< unsigned short > >,
      MyPointer< std::vector< short > >,
      MyPointer< std::vector< unsigned int > >,
      MyPointer< std::vector< int > >,
      MyPointer< std::vector< unsigned long > >,
      MyPointer< std::vector< long > >,
      MyPointer< std::vector< unsigned long long > >,
      MyPointer< std::vector< long long > >,
      MyPointer< std::vector< float > >,
      MyPointer< std::vector< double > >,
      MyPointer< std::vector< std::vector< double > > >,
      MyPointer< std::vector< TVector3> >,
      MyPointer< std::map< std::string, std::vector< double > > >
    >;

    using iterator = std::map< std::string, Variant >::iterator;
    using const_iterator = std::map< std::string, Variant >::const_iterator;

    class VariantWrapper {

      public:

        VariantWrapper( Variant* v ) : variant_( v ) {}

        // Overload the = operator to allow setting new values of the variant
        // easily even when working with types that are internally wrapped
        // in a MyPointer.
        template < typename T > VariantWrapper& operator=( const T& in ) {

          // Check that the variant either already holds a value of the
          // requested type or has not previously been assigned a value
          // and thus has type std::monostate. If neither of these things are
          // true, then we are attempting to switch types. To prevent
          // issues with inconsistent types in a TTree branch, throw an
          // exception.
          bool is_monostate = std::holds_alternative<
            std::monostate >( *variant_ );

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
              " previously initialized TreeMap::Variant object" );
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
                throw std::runtime_error( "Invalid dereference in"
                  " TreeMap::VariantWrapper assignment" );
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
        // into a target variable. The type is inferred from the target without
        // the need for explicit use of a template parameter.
        template < typename T > const
          VariantWrapper& operator>>( T& out ) const
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

        Variant* variant_ = nullptr;
    };

    TreeMap() = default;
    TreeMap( TTree* tree ) : tree_( tree ) {}

    Variant& at( const std::string& name ) { return map_.at( name ); }

    const Variant& at( const std::string& name ) const
      { return map_.at( name ); }

    Variant& operator[]( const std::string& name )
      { return map_[ name ]; }

    TTree* tree() { return tree_; }

    const TTree* tree() const { return tree_; }

    template< class... Args > auto emplace( Args&&... args )
      { return map_.emplace( std::forward< Args >( args )... ); }

    auto begin() { return map_.begin(); }
    auto cbegin() const { return map_.cbegin(); }

    auto end() { return map_.end(); }
    auto cend() const { return map_.end(); }

    auto find( const std::string& key ) { return map_.find( key ); }

    Variant* add_input_branch( const std::string& name,
      TBranch*& added_branch );

    inline Variant* add_input_branch( const std::string& name ) {
      TBranch* dummy = nullptr;
      return add_input_branch( name, dummy );
    }

    void set_tree( TTree* tree ) { tree_ = tree; }

    const auto& var_size_map() const { return var_size_map_; }

    // Calls TTree::Write() for the owned TTree
    inline void write() { if ( tree_ ) tree_->Write(); }

    // Calls TTree::LoadTree() for the owned TTree
    inline long long load_tree( long long entry ) {
      if ( !tree_ ) return -1;
      return tree_->LoadTree( entry );
    }

    // Calls TTree::GetReadEntry() for the owned TTree
    inline long long get_read_entry() const {
      if ( !tree_ ) return -1;
      else return tree_->GetReadEntry();
    }

    // Call TTree::GetEntry() for the owned TTree
    void get_entry( long long entry );

    // Call TTree::Fill() for the owned TTree. Before doing so, create
    // any missing branches and set branch addresses
    void fill();

  protected:

    // Flag that prevents adding new branches after fill() has been called
    // at least once
    bool output_locked_ = false;

    // Pointer to the owned TTree
    TTree* tree_ = nullptr;

    // Storage for branches of various types
    std::map< std::string, Variant > map_;

    // Storage for the names of branches storing sizes for variable-size
    // C-style arrays (needed for dynamic resizing of vectors used to
    // handle the input)
    std::map< std::string, std::set< std::string > > var_size_map_;
};

// Returns a T* to the content of a TreeMap::Variant if T is the actively-stored
// type, or nullptr otherwise
template < typename T > T* get_active_ptr_as( TreeMap::Variant& v ) {
  // For fundamental types, the TreeMap::Variant uses them directly, so just get
  // a pointer via std::get_if()
  if constexpr ( std::is_fundamental_v< T > ) {
    return std::get_if< T >( &v );
  }
  else {
    // For classes, structs, etc., we need to go inside the owned MyPointer
    auto* my_ptr = std::get_if< MyPointer< T > >( &v );
    if ( my_ptr ) {
      return my_ptr->get();
    }
  }
  return nullptr;
}

// Returns a T& to the content of a TreeMap::Variant if T is the actively-stored
// type. Throws an expression if it is not.
template < typename T > T& get_active_ref_as( TreeMap::Variant& v ) {
  auto* temp_ptr = get_active_ptr_as< T >( v );
  if ( !temp_ptr ) {
    throw std::runtime_error( "Invalid TreeMap::Variant dereference" );
  }
  return *temp_ptr;
}

// Determines T from the type of the second argument and
// copies the contents of the branch into the variable "out"
template< typename T > void copy_my_variant( TreeMap::Variant& v, T& out ) {
  out = get_active_ref_as< T >( v );
}

// Determines T from the type of the second argument and
// gets a bare T* to the branch variable owned by the variant
template< typename T > void get_my_variant( TreeMap::Variant& v, T*& out_ptr ) {
  out_ptr = get_active_ptr_as< T >( v );
}

// Overload the >> operator to allow easy loading of the active variant into a
// target variable. The type is inferred from the target without the
// need for explicit use of a template parameter.
template < typename T > void operator>>( TreeMap::Variant& v, T& out ) {
  if constexpr ( std::is_pointer_v< T > ) {
    get_my_variant( v, out );
  }
  else {
    copy_my_variant( v, out );
  }
}

// Helper function used to set branch addresses when populating a TreeMap
template< typename T > void set_variant_input_branch_address(
  const std::pair< TreeMap::iterator, bool >& emplace_result,
  TTree& in_tree, const std::string& branch_name )
{
  const bool& was_emplaced = emplace_result.second;
  if ( !was_emplaced ) {
    std::cerr << "WARNING: Duplicate branch name \""
      << branch_name << "\"" << " in the \""
      << in_tree.GetName() << "\" tree will be ignored.\n";
    return;
  }
  const auto& iter = emplace_result.first;
  auto& var = iter->second;
  auto* var_ptr = std::get_if< T >( &var );

  if ( !var_ptr ) {
    throw std::runtime_error( "Unsuccessful type retrieval for variant"
      " representing the branch \"" + branch_name + "\"" );
    return;
  }

  // If the branch type of interest is a MyPointer object, then call the helper
  // function for setting the pointed-to object's branch address
  if constexpr ( is_MyPointer_v< T > ) {
    typename T::element_type*& address = var_ptr->get_bare_ptr();
    in_tree.SetBranchAddress( branch_name.c_str(), &address );
  }
  // Otherwise, the branch is a simple type, and we can use
  // TTree::SetBranchAddress() with the bare pointer directly
  else {
    in_tree.SetBranchAddress( branch_name.c_str(), var_ptr );
  }
}

// Creates a new TreeMap::Variant in a TreeMap and sets the corresponding
// branch address in an associated TTree
template < typename T > std::pair< TreeMap::iterator, bool >
  emplace_variant_and_set_input_address( TreeMap& branch_map,
    TTree& in_tree, const std::string& branch_name )
{
  // If we're working with a fundamental type, then just use a new instance
  // for branch data storage. Otherwise, construct a MyPointer object for easy
  // memory management.
  using StorageType = std::conditional_t<
    std::is_fundamental< T >::value,
    T,  // type T is fundamental (int, char, double, bool, etc.)
    MyPointer< T > // type T is a class, struct, etc.
  >;

  TBranch* br = in_tree.GetBranch( branch_name.c_str() );
  if ( !br ) {
    throw std::runtime_error( "Branch \"" + branch_name + "\" is not"
      " present in the TTree " + in_tree.GetName() );
  }

  auto er = branch_map.emplace( branch_name, StorageType() );
  set_variant_input_branch_address< StorageType >( er, in_tree, branch_name );
  return er;
}
