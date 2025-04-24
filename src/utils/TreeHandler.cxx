// ROOT includes
#include "TLeaf.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/TreeHandler.hh"
#include "XSecAnalyzer/TreeUtils.hh"

namespace {

  // Emplaces either a scalar type or a std::vector for an array of that
  // type based on the boolean input from the user
  template < typename T > std::pair< TreeMap::iterator, bool >
    emplace_variant_helper( TreeMap& tree_map, TTree& in_tree,
      const std::string& branch_name, bool is_scalar_branch )
  {
    if ( is_scalar_branch ) {
      return emplace_variant_and_set_input_address< T >( tree_map,
        in_tree, branch_name );
    }

    // The branch represents an array, so we handle that special case here
    TBranch* br = in_tree.GetBranch( branch_name.c_str() );
    if ( !br ) {
      throw std::runtime_error( "Array branch \"" + branch_name + "\" is not"
        " present in the TTree " + in_tree.GetName() );
    }

    // Use a MyPointer< std::vector< T > > as the storage type in the variant
    auto er = tree_map.emplace( branch_name,
      MyPointer< std::vector< T > >() );

    // Check for a duplicate branch indicated by a failed emplace
    const bool& was_emplaced = er.second;
    if ( !was_emplaced ) {
      std::cerr << "WARNING: Duplicate array branch name \""
        << branch_name << "\"" << " in the \""
        << in_tree.GetName() << "\" tree will be ignored.\n";
      return er;
    }
    const auto& iter = er.first;
    auto& var = iter->second;
    auto* var_ptr = std::get_if< MyPointer< std::vector< T > > >( &var );

    if ( !var_ptr ) {
      throw std::runtime_error( "Unsuccessful type retrieval for array variant"
        " representing the branch \"" + branch_name + "\"" );
      return er;
    }

    // Here's the trick we use for the C-style arrays: set the branch address
    // using a pointer to the first element in the array owned by the
    // std::vector in the variant. Since a std::vector has elements that are
    // guaranteed by the C++ standard to be contiguous in memory, this trick
    // will work as long as we check the array size before calling
    // TTree::GetEntry to make sure that the vector is large enough to hold the
    // array.
    in_tree.SetBranchAddress( branch_name.c_str(),
      ( *var_ptr )->data() );

    return er;
  }

  // Sets the branch address for output TTree variables in a unified way
  template < typename T, typename S > TBranch* set_branch(
    TTree& out_tree, S var_name, T*& var )
  {
    std::string var_name_str( var_name );

    // Check if the branch already exists. If it doesn't, then create it.
    TBranch* br = out_tree.GetBranch( var_name_str.c_str() );
    if ( !br ) {

      std::string leaf_list( var_name );

      // The use of if constexpr here is a C++17 feature
      if constexpr ( std::is_same_v< T, bool > ) leaf_list += "/O";
      else if constexpr ( std::is_same_v< T, double > ) leaf_list += "/D";
      else if constexpr ( std::is_same_v< T, float > ) leaf_list += "/F";
      else if constexpr ( std::is_same_v< T, int > ) leaf_list += "/I";
      else if constexpr ( std::is_same_v< T, unsigned int > ) leaf_list += "/i";
      else leaf_list = "";

      // Branches for objects do not use a leaf list and use a
      // pointer-to-a-pointer to set the address
      if ( leaf_list.empty() ) {
        br = out_tree.Branch( var_name_str.c_str(), &var );
      }
      // Branches for simple types need to specify a leaf list upon creation
      // and use a regular pointer to set the address
      else {
        br = out_tree.Branch( var_name_str.c_str(), var, leaf_list.c_str() );
      }
    }

    // Check the branch again. If it's nullptr here, we failed to create it.
    if ( !br ) {
      throw std::runtime_error( "Failed to create output branch \""
        + var_name_str + "\"" );
      return nullptr;
    }

    // Now set the branch address. Note that fundamental types use a pointer
    // while objects use a pointer-to-a-pointer.
    if constexpr ( std::is_fundamental_v< T > ) {
      // The branch represents a fundamental type
      br->SetAddress( var );
    }
    else {
      // The branch represents a class, struct, etc.
      br->SetAddress( &var );
    }

    return br;
  }

  // Overloaded version that takes an rvalue reference to a pointer argument
  // instead (see, e.g., https://stackoverflow.com/a/5465371/4081973)
  template < typename T, typename S >
    void set_branch( TTree& out_tree, S var_name, T*&& rval_ptr )
  {
    T* address = rval_ptr;
    set_branch( out_tree, var_name, address );
  }

  // Overloaded version that takes a reference to the type T instead. This
  // version auto-detects a MyPointer argument and adjusts the approach
  // accordingly
  template < typename T, typename S >
    void set_branch( TTree& out_tree, S var_name, T& ref )
  {
    if constexpr ( is_MyPointer_v< T > ) {
      auto*& address = ref.get_bare_ptr();
      set_branch( out_tree, var_name, address );
    }
    else {
      T* address = &ref;
      set_branch( out_tree, var_name, address );
    }
  }

}

// Uses an input TTree to build a map of branch names each associated with a
// MyVariant object. These objects are linked to the TTree using
// TTree::SetBranchAddress() to provide temporary storage for event-by-event
// data processing.
void TreeHandler::add_input_tree( TTree* in_tree, const std::string& name ) {

  // Complain if a null pointer is passed to this function
  if ( !in_tree ) {
    throw std::runtime_error( "nullptr passed to TreeHandler"
      "::add_input_tree()" );
    return;
  }

  // If the input key is empty, default to using the name reported
  // by the TTree itself
  std::string key = name;
  if ( key.empty() ) {
    key = in_tree->GetName();
  }

  // If other trees have already been loaded, double-check that the
  // new one has the same number of entries. If there is a mismatch,
  // the trees will not be compatible. Complain if a problem is found.
  if ( !in_tree_maps_.empty() ) {
    long long num_events = in_tree->GetEntries();
    for ( auto pair_iter = in_tree_maps_.cbegin();
      pair_iter != in_tree_maps_.cend(); ++pair_iter )
    {
      TTree* old_tree = pair_iter->second.first;
      long long old_num_events = old_tree->GetEntries();
      if ( num_events != old_num_events ) {
        throw std::runtime_error( "Mismatch of event entry counts"
          " between the \"" + key + "\" and \"" + old_tree->GetName()
          + "\" TTree objects" );
      }
    }
  }

  // Double-check that the new TTree will have a unique key in the map
  if ( in_tree_maps_.count( key ) ) {
    throw std::runtime_error( "Duplicate TTree name \"" + key
      + "\" encountered in TreeHandler::add_input_tree()" );
  }

  // Create a new entry in the map with a default-constructed value
  auto& map_value = in_tree_maps_[ key ];

  // Set the TTree pointer in the new map entry
  map_value.first = in_tree;

  // Get access to the TreeMap object owned by the new map entry
  auto& tree_map = map_value.second;

  auto branch_list = in_tree->GetListOfBranches();
  int num_branches = branch_list->GetEntries();

  // Temporary variables used to query each branch about its data type
  // via the TBranch::GetExpectedType() function
  TClass* temp_class = nullptr;
  EDataType temp_type;

  for ( int b = 0; b < num_branches; ++b ) {

    // Retrieve data type information for the current branch
    auto* br = dynamic_cast< TBranch* >( branch_list->At(b) );
    std::string branch_name = br->GetName();
    br->GetExpectedType( temp_class, temp_type );

    // If the TClass pointer is non-null, the branch corresponds to an object
    if ( temp_class ) {
      // Get the name of the class stored in the current branch
      std::string class_name = temp_class->GetName();

      // Initialize an appropriate variant in the map matching the class of
      // interest
      // TODO: add more cases as needed below
      if ( class_name == "string" ) {
        emplace_variant_and_set_input_address< std::string >(
          tree_map, *in_tree, branch_name );
      }
      else if ( class_name == "TVector3" ) {
        emplace_variant_and_set_input_address< TVector3 >(
          tree_map, *in_tree, branch_name );
      }
      else if ( class_name == "vector<bool>" ) {
        emplace_variant_and_set_input_address< std::vector< bool > >(
          tree_map, *in_tree, branch_name );
      }
      else if ( class_name == "vector<double>" ) {
        emplace_variant_and_set_input_address< std::vector< double > >(
          tree_map, *in_tree, branch_name );
      }
      else if ( class_name == "vector<float>" ) {
        emplace_variant_and_set_input_address< std::vector< float > >(
          tree_map, *in_tree, branch_name );
      }
      else if ( class_name == "vector<int>" ) {
        emplace_variant_and_set_input_address< std::vector< int > >(
          tree_map, *in_tree, branch_name );
      }
      else if ( class_name == "vector<unsigned int>" ) {
        emplace_variant_and_set_input_address< std::vector< unsigned int > >(
          tree_map, *in_tree, branch_name );
      }
      else if ( class_name == "vector<unsigned long>" ) {
        emplace_variant_and_set_input_address< std::vector< unsigned long > >(
          tree_map, *in_tree, branch_name );
      }
      else if ( class_name == "vector<unsigned short>" ) {
        emplace_variant_and_set_input_address< std::vector< unsigned short > >(
          tree_map, *in_tree, branch_name );
      }
      else if ( class_name == "vector<vector<double> >" ) {
        emplace_variant_and_set_input_address< std::vector<
          std::vector< double > > >( tree_map, *in_tree, branch_name );
      }
      else if ( class_name == "vector<TVector3>" ) {
        emplace_variant_and_set_input_address< std::vector< TVector3 > >(
          tree_map, *in_tree, branch_name );
      }
      else if ( class_name == "map<string,vector<double> >" ) {
        emplace_variant_and_set_input_address< std::map<
          std::string, std::vector< double > > >( tree_map, *in_tree,
          branch_name );
      }
      else throw std::runtime_error( "Unrecognized class type \""
        + class_name + "\"" + " for branch \"" + branch_name + "\"" );
    }
    // Otherwise the branch has a simple data type indicated by an enum value.
    // However, it could still be a fixed- or variable-size array of the simple
    // data type. Check this first before instantiating an appropriate variant
    // in the map.
    else {

      // Get the first leaf associated with the current branch
      TLeaf* lf = dynamic_cast< TLeaf* >( br->GetListOfLeaves()->First() );
      if ( !lf ) {
        throw std::runtime_error( "Missing leaf for branch \""
          + branch_name + "\"" );
      }

      // TLeaf::GetLeafCounter() will store a pointer to the TLeaf that is used
      // to store the size if this is indeed a variable-length array. If it is
      // nullptr, you can use count_result to get the fixed size. A nonpositive
      // value of count_result represents an error condition.
      int count_result;
      TLeaf* lfc = lf->GetLeafCounter( count_result );

      if ( count_result <= 0 ) {
        throw std::runtime_error( "Call to TLeaf::GetLeafCounter() failed" );
        return;
      }

      // A null leaf counter pointer and a size of one indicates that
      // we are working with a branch that represents a single instance
      // of a fundamental type
      bool is_scalar_branch = ( !lfc && count_result == 1 );

      // A null leaf counter pointer and a size greater than one indicates
      // that we are working with a fixed size array
      bool is_fixed_size_array = ( !lfc && count_result > 1 );

      // Based on the reported type, allocate storage for the branch with
      // a variant and set the branch address
      std::pair< TreeMap::iterator, bool > er;
      if ( temp_type == kBool_t ) {
        // Handle only the scalar case for now since std::vector< bool >
        // is not guaranteed to store its elements in a contiguous array.

        if ( !is_scalar_branch ) {
          std::cerr << "WARNING: Handling of C-style arrays of bool is not"
            << " yet implemented in TreeHandler. Branch \"" << branch_name
            << "\" will be ignored.\n";
          continue;
        }

        emplace_variant_and_set_input_address< bool >( tree_map, *in_tree,
          branch_name );

        //er = emplace_variant_helper< bool >( tree_map, *in_tree, branch_name,
        //  is_scalar_branch );
      }
      else if ( temp_type == kUChar_t ) {
        er = emplace_variant_helper< unsigned char >( tree_map, *in_tree,
          branch_name, is_scalar_branch );
      }
      else if ( temp_type == kChar_t ) {
        er = emplace_variant_helper< char >( tree_map, *in_tree,
          branch_name, is_scalar_branch );
      }
      else if ( temp_type == kUShort_t ) {
        er = emplace_variant_helper< unsigned short >( tree_map, *in_tree,
          branch_name, is_scalar_branch );
      }
      else if ( temp_type == kShort_t ) {
        er = emplace_variant_helper< short >( tree_map, *in_tree,
          branch_name, is_scalar_branch );
      }
      else if ( temp_type == kUInt_t ) {
        er = emplace_variant_helper< unsigned int >( tree_map, *in_tree,
          branch_name, is_scalar_branch );
      }
      else if ( temp_type == kInt_t ) {
        er = emplace_variant_helper< int >( tree_map, *in_tree,
          branch_name, is_scalar_branch );
      }
      else if ( temp_type == kULong_t ) {
        er = emplace_variant_helper< unsigned long >( tree_map, *in_tree,
          branch_name, is_scalar_branch );
      }
      else if ( temp_type == kLong_t ) {
        er = emplace_variant_helper< long >( tree_map, *in_tree,
          branch_name, is_scalar_branch );
      }
      else if ( temp_type == kULong64_t ) {
        er = emplace_variant_helper< unsigned long long >( tree_map, *in_tree,
          branch_name, is_scalar_branch );
      }
      else if ( temp_type == kLong64_t ) {
        er = emplace_variant_helper< long long >( tree_map, *in_tree,
          branch_name, is_scalar_branch );
      }
      else if ( temp_type == kFloat_t
        || temp_type == kFloat16_t )
      {
        er = emplace_variant_helper< float >( tree_map, *in_tree,
          branch_name, is_scalar_branch );
      }
      else if ( temp_type == kDouble_t
        || temp_type == kDouble32_t )
      {
        er = emplace_variant_helper< double >( tree_map, *in_tree,
          branch_name, is_scalar_branch );
      }
      else throw std::runtime_error( "Unrecognized EDataType value "
        + std::to_string(temp_type) + " for branch \""
        + branch_name + "\"" );

      // For non-scalar branches (i.e., C-style arrays), manage the
      // storage depending on whether it is a fixed-size array or not.
      if ( !is_scalar_branch ) {
        if ( is_fixed_size_array ) {
          // For a fixed-size array, use std::visit to adjust the
          // wrapper vector owned by the variant to the proper size
          auto& temp_variant = er.first->second;
          std::visit( [ &in_tree, &branch_name, count_result ](
            auto& var_val ) -> void
            {
              using T = std::decay_t< decltype( var_val ) >;
              // We can't currently handle fixed-length arrays of bool
              // using the vector strategy because std::vector< bool >
              // is not guaranteed by the C++ standard to be stored
              // as a contiguous array in memory.
              if constexpr ( is_MyPointerToVector_v< T >
                && !std::is_same_v< T, MyPointer< std::vector< bool > > > )
              {
                var_val->resize( count_result );
                // std::vector::resize() can invalidate pointers to vector
                // elements, so update the branch address accordingly
                in_tree->SetBranchAddress( branch_name.c_str(),
                  var_val->data() );
              }
            }, temp_variant );
        }
        else {
	  // For variable-size arrays, we resort to adjusting the vector size
	  // during every call to get_entry(). Record the name of the current
	  // TTree, size leaf, and the branch to be automatically sized in the
	  // map used for this purpose.
          const std::string& size_leaf_name = lfc->GetName();
          in_var_size_map_[ key ][ size_leaf_name ].insert( branch_name );
        }
      }
    }
  }
}

// Call TTree::GetEntry() for each input TTree
void TreeHandler::get_entry( long long entry ) {
  for ( auto pair_iter = in_tree_maps_.begin();
    pair_iter != in_tree_maps_.end(); ++pair_iter )
  {
    // First get access to the current TTree and
    // a map storing information about its variable-size arrays
    TTree* temp_tree = pair_iter->second.first;
    if ( !temp_tree ) continue;

    TreeMap& tm = pair_iter->second.second;

    const std::string& tree_name = pair_iter->first;
    auto& tree_size_map = in_var_size_map_[ tree_name ];

    // Pre-load all branches that contain the variable sizes
    for ( const auto& size_pair : tree_size_map ) {
      const std::string& size_leaf_name = size_pair.first;
      TBranch* size_br = temp_tree->GetBranch( size_leaf_name.c_str() );
      if ( !size_br ) throw std::runtime_error( "Missing array size branch \""
        + size_leaf_name + "\" encountered" );
      else size_br->GetEntry( entry );
    }

    // For each variable size, resize the vectors corresponding
    // to the branches that are controlled by it
    for ( const auto& size_pair : tree_size_map ) {
      const std::string& size_leaf_name = size_pair.first;

      // Retrieve the stored size from the TreeMap and save it
      // to an integer for later use
      int vec_size;
      tm.at( size_leaf_name ) >> vec_size;
      if ( vec_size < 0 ) {
        throw std::runtime_error( "Invalid array size encountered in"
          " TreeHandler::get_entry()" );
      }

      // Update the sizes of all vector wrappers for the array. Since the call
      // to resize() can invalidate pointers to vector elements, update the
      // branch address immediately after resizing.
      const auto& branch_name_set = size_pair.second;
      for ( const auto& br_name : branch_name_set ) {
        auto& var = tm.at( br_name );
        std::visit( [ &temp_tree, &br_name, vec_size ]( auto& var_val ) -> void
          {
            using T = std::decay_t< decltype( var_val ) >;
            if constexpr ( is_MyPointerToVector_v< T >
              && !std::is_same_v< T, MyPointer< std::vector< bool > > > )
            {
              var_val->resize( vec_size );
              temp_tree->SetBranchAddress( br_name.c_str(), var_val->data() );
            }
          }, var );
      }
    }

    // We're ready. Now retrieve the current entry reading from all enabled
    // TTree branches
    temp_tree->GetEntry( entry );
  }
}

// Call TTree::Fill() for each output TTree. Before doing so, create
// any missing branches and set branch addresses
void TreeHandler::fill() {
  for ( auto pair_iter = out_tree_maps_.begin();
    pair_iter != out_tree_maps_.end(); ++pair_iter )
  {
    // Get access to the current output TTree. If we get
    // a nullptr, just move on to the next one.
    TTree* temp_tree = pair_iter->second.first;
    if ( !temp_tree ) {
      continue;
    }

    // Now get access to the TreeMap used to manage the branches
    TreeMap& tm = pair_iter->second.second;

    // Iterate over the TreeMap to set the branch addresses before filling.
    for ( auto& tm_pair : tm ) {
      const std::string& br_name = tm_pair.first;
      MyVariant& var = tm_pair.second;

      // Set the branch address (and create a new branch if needed).  Use
      // std::visit to automatically handle the type of the active variant.
      std::visit( [ &temp_tree, &br_name ]( auto& var_val )
        -> void { set_branch( *temp_tree, br_name, var_val ); }, var );
    }

    // We're done preparing everything. Fill the tree.
    temp_tree->Fill();
  }
}

// Call TTree::Write() for each output TTree
// TODO: reduce code duplication here
void TreeHandler::write() {
  for ( auto pair_iter = out_tree_maps_.cbegin();
    pair_iter != out_tree_maps_.cend(); ++pair_iter )
  {
    TTree* temp_tree = pair_iter->second.first;
    if ( temp_tree ) {
      temp_tree->Write();
    }
  }
}

TreeAndTreeMap* TreeHandler::find_element( const std::string& name,
  bool input )
{
  auto end = in_tree_maps_.end();
  if ( !input ) end = out_tree_maps_.end();

  decltype( end ) iter;
  if ( input ) iter = in_tree_maps_.find( name );
  else iter = out_tree_maps_.find( name );

  if ( iter == end ) {
    std::string err_message = "No loaded ";
    if ( input ) err_message += "input";
    else err_message += "output";

    throw std::runtime_error( err_message + " has the name \""
      + name + "\"" );
    return nullptr;
  }
  return &( iter->second );
}

TTree* TreeHandler::in_tree( const std::string& name ) {
  auto* pair = this->find_element( name, true );
  return pair->first;
}

TreeMapWrapper TreeHandler::in_map( const std::string& name ) {
  auto* pair = this->find_element( name, true );
  return TreeMapWrapper( &(pair->second) );
}

TreeMap& TreeHandler::access_in_map( const std::string& name ) {
  auto* pair = this->find_element( name, true );
  return pair->second;
}

TTree* TreeHandler::out_tree( const std::string& name ) {
  auto* pair = this->find_element( name, false );
  return pair->first;
}

TreeMapWrapper TreeHandler::out_map( const std::string& name ) {
  auto* pair = this->find_element( name, false );
  return TreeMapWrapper( &(pair->second) );
}

TreeMap& TreeHandler::access_out_map( const std::string& name ) {
  auto* pair = this->find_element( name, false );
  return pair->second;
}

// Defines an output TTree to be written to a given TDirectory.
void TreeHandler::add_output_tree( TDirectory* out_dir,
  const std::string& name, const std::string& title )
{
  // Complain if a null pointer is passed to this function
  if ( !out_dir ) {
    throw std::runtime_error( "nullptr passed to TreeHandler"
      "::add_output_tree()" );
    return;
  }

  // Complain if the tree name is empty
  if ( name.empty() ) {
    throw std::runtime_error( "Empty tree name passed to TreeHandler"
      "::add_output_tree()" );
  }

  // Double-check that the new TTree will have a unique key in the map
  if ( out_tree_maps_.count( name ) ) {
    throw std::runtime_error( "Duplicate TTree name \"" + name
      + "\" encountered in TreeHandler::add_output_tree()" );
  }

  // Create a new entry in the map with a default-constructed value
  auto& map_value = out_tree_maps_[ name ];

  // Initialize the new TTree and associate it with the output TDirectory
  TTree*& out_tree = map_value.first;
  out_tree = new TTree( name.c_str(), title.c_str() );
  out_tree->SetDirectory( out_dir );
}
