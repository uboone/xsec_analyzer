#pragma once

// Standard library includes
#include <iostream>

// XSecAnalyzer includes
#include "XSecAnalyzer/TreeMap.hh"

// **** Helper code to facilitate setting TTree branch addresses, etc. ****

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
