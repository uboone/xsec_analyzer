#pragma once

// Standard library includes
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include <vector>

// ROOT includes
#include "TClass.h"
#include "TDataType.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TTree.h"

//#include "AnalysisEvent.hh"

// **** Helper code to facilitate setting TTree branch addresses, etc. ****

// A std::unique_ptr with redundant storage of a bare pointer to the managed
// object. This is a hacky workaround for setting TTree branch addresses for
// objects managed by a std::unique_ptr.
template <typename T> class MyPointer : public std::unique_ptr<T> {
  public:

    MyPointer() : std::unique_ptr<T>( new T ) {}

    T*& get_bare_ptr() {
      bare_ptr_ = this->get();
      return bare_ptr_;
    }

  protected:

    T* bare_ptr_ = nullptr;
};

// Create a type trait to allow "constexpr if" to distinguish between
// instantiations of the MyPointer class template and other types
template <typename T>
  struct is_MyPointer : std::false_type {};

template < typename T >
  struct is_MyPointer< MyPointer< T > > : std::true_type {};

template < typename T >
  constexpr bool is_MyPointer_v = is_MyPointer< T >::value;

// A std::variant that can store any of the branch types of interest for
// the TTrees used by xsec_analyzer. New branches with a different type
// will require extending the template arguments given here.
using MyVariant = std::variant<
  // List the simple types that are allowed here first
  unsigned char,
  unsigned int,
  bool,
  float,
  int,
  double,
  // Classes, structs, etc. that are allowed should be wrapped in
  // a MyPointer here
  MyPointer< std::string >,
  MyPointer< std::vector< bool > >,
  MyPointer< std::vector< unsigned short > >,
  MyPointer< std::vector< unsigned int > >,
  MyPointer< std::vector< unsigned long > >,
  MyPointer< std::vector< int > >,
  MyPointer< std::vector< float > >,
  MyPointer< std::vector< double > >,
  MyPointer< std::vector< std::vector< double > > >,
  MyPointer< std::map< std::string, std::vector< double > > >
>;

// A map used to automatically allocate storage for managing TTree branch
// variables using MyVariant objects. Keys are strings giving the branch
// names (no duplicates allowed). Values are MyVariant objects that will
// be used to store the branch data.
using TreeMap = std::map< std::string, MyVariant >;

// Returns a T* to the content of a MyVariant if T is the actively-stored type,
// or nullptr otherwise
template < typename T > T* get_active_ptr_as( MyVariant& v ) {
  // For fundamental types, the MyVariant uses them directly, so just get a
  // pointer via std::get_if()
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

// Returns a T& to the content of a MyVariant if T is the actively-stored type.
// Throws an expression if it is not.
template < typename T > T& get_active_ref_as( MyVariant& v ) {
  auto* temp_ptr = get_active_ptr_as< T >( v );
  if ( !temp_ptr ) {
    throw std::runtime_error( "Invalid MyVariant dereference" );
  }
  return *temp_ptr;
}

// Determines T from the type of the second argument and
// copies the contents of the branch into the variable "out"
template< typename T > void copy_my_variant( MyVariant& v, T& out ) {
  out = get_active_ref_as< T >( v );
}

// Determines T from the type of the second argument and
// gets a bare T* to the branch variable owned by the variant
template< typename T > void get_my_variant( MyVariant& v, T*& out_ptr ) {
  out_ptr = get_active_ptr_as< T >( v );
}

// Helper function template that sets a new address for a pointer to an object
// in an input TTree
template <typename T> void set_object_input_branch_address( TTree& in_tree,
  const std::string& branch_name, T*& address )
{
  in_tree.SetBranchAddress( branch_name.c_str(), &address );
}

// Overloaded version that uses a MyPointer argument instead
template <typename T> void set_object_input_branch_address( TTree& in_tree,
  const std::string& branch_name, MyPointer<T>& u_ptr )
{
  T*& address = u_ptr.get_bare_ptr();
  set_object_input_branch_address( in_tree, branch_name, address );
}

// Helper function that creates a branch (or just sets a new address) for a
// simple variable in an output TTree
inline void set_output_branch_address( TTree& out_tree,
  const std::string& branch_name, void* address, bool create = false,
  const std::string& leaf_spec = "" )
{
  if ( create ) {
    if ( leaf_spec != "" ) {
      out_tree.Branch( branch_name.c_str(), address, leaf_spec.c_str() );
    }
    else {
      out_tree.Branch( branch_name.c_str(), address );
    }
  }
  else {
    out_tree.SetBranchAddress( branch_name.c_str(), address );
  }
}

// Helper function template that creates a branch (or just sets a new address)
// for a pointer to an object in an output TTree
template <typename T> void set_object_output_branch_address( TTree& out_tree,
  const std::string& branch_name, T*& address, bool create = false )
{
  if ( create ) out_tree.Branch( branch_name.c_str(), &address );
  else out_tree.SetBranchAddress( branch_name.c_str(), &address );
}

// Overloaded version that uses a MyPointer argument instead
template <typename T> void set_object_output_branch_address( TTree& out_tree,
  const std::string& branch_name, MyPointer<T>& u_ptr, bool create = false )
{
  T*& address = u_ptr.get_bare_ptr();
  set_object_output_branch_address( out_tree, branch_name, address, create );
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
    set_object_input_branch_address( in_tree, branch_name, *var_ptr );
  }
  // Otherwise, the branch is a simple type, and we can use
  // TTree::SetBranchAddress() with the bare pointer directly
  else {
    in_tree.SetBranchAddress( branch_name.c_str(), var_ptr );
  }
}

// Creates a new MyVariant in a TreeMap and sets the corresponding
// branch address in an associated TTree
template < typename T > void emplace_variant_and_set_input_address(
  TreeMap& branch_map, TTree& in_tree, const std::string& branch_name )
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
}
