// XSecAnalyzer includes
#include "XSecAnalyzer/TreeHandler.hh"
#include "XSecAnalyzer/TreeUtils.hh"

namespace {

  // Sets the branch address for output TTree variables managed by this
  // SelectionBase object
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
      else if ( class_name == "map<string,vector<double> >" ) {
        emplace_variant_and_set_input_address< std::map<
          std::string, std::vector< double > > >( tree_map, *in_tree,
          branch_name );
      }
      else throw std::runtime_error( "Unrecognized class type \""
        + class_name + "\"" + " for branch \"" + branch_name + "\"" );
    }
    // Otherwise the branch has a simple data type indicated by an enum value.
    // Instantiate an appropriate variant for this type in the map.
    else {
      if ( temp_type == kUChar_t ) {
        emplace_variant_and_set_input_address< unsigned char >(
          tree_map, *in_tree, branch_name );  
      }
      else if ( temp_type == kUInt_t ) {
        emplace_variant_and_set_input_address< unsigned int >(
          tree_map, *in_tree, branch_name );  
      }
      else if ( temp_type == kBool_t ) {
        emplace_variant_and_set_input_address< bool >(
          tree_map, *in_tree, branch_name );  
      }
      else if ( temp_type == kFloat_t
        || temp_type == kFloat16_t )
      {
       emplace_variant_and_set_input_address< float >(
          tree_map, *in_tree, branch_name );  
      }
      else if ( temp_type == kInt_t ) {
        emplace_variant_and_set_input_address< int >(
          tree_map, *in_tree, branch_name );  
      }
      else if ( temp_type == kDouble_t ) {
        emplace_variant_and_set_input_address< double >(
          tree_map, *in_tree, branch_name );  
      }
      else throw std::runtime_error( "Unrecognized EDataType value "
        + std::to_string(temp_type) + " for branch \""
        + branch_name + "\"" );
    }
  }
}

// Call TTree::GetEntry() for each input TTree
void TreeHandler::get_entry( long long entry ) {
  for ( auto pair_iter = in_tree_maps_.cbegin();
    pair_iter != in_tree_maps_.cend(); ++pair_iter )
  {
    TTree* temp_tree = pair_iter->second.first;
    if ( temp_tree ) {
      temp_tree->GetEntry( entry );
    }
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
