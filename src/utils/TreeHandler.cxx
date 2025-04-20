// XSecAnalyzer includes
#include "XSecAnalyzer/TreeHandler.hh"
#include "XSecAnalyzer/TreeUtils.hh"

// Uses an input TTree to build a map of branch names each associated with a
// MyVariant object. These objects are linked to the TTree using
// TTree::SetBranchAddress() to provide temporary storage for event-by-event
// data processing.
void TreeHandler::add_input_tree( TTree* in_tree, const std::string& name ) {

  // Complain if a null pointer is passed to this function
  if ( !in_tree ) {
    throw std::runtime_error( "nullptr passed to TreeHandler"
      "::add_input_tree()" );
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
  if ( !tree_maps_.empty() ) {
    long long num_events = in_tree->GetEntries();
    for ( auto pair_iter = tree_maps_.cbegin();
      pair_iter != tree_maps_.cend(); ++pair_iter )
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
  if ( tree_maps_.count( key ) ) {
    throw std::runtime_error( "Duplicate TTree name \"" + key
      + "\" encountered in TreeHandler::add_input_tree()" );
  }

  // Create a new entry in the map with a default-constructed value
  auto& map_value = tree_maps_[ key ];

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
  for ( auto pair_iter = tree_maps_.cbegin();
    pair_iter != tree_maps_.cend(); ++pair_iter )
  {
    TTree* temp_tree = pair_iter->second.first;
    if ( temp_tree ) {
      temp_tree->GetEntry( entry );
    }
  }
}

TreeAndTreeMap* TreeHandler::find_element( const std::string& name ) {
  auto iter = tree_maps_.find( name );
  if ( iter == tree_maps_.end() ) {
    throw std::runtime_error( "No loaded tree has the name \""
      + name + "\"" );
    return nullptr;
  }
  return &( iter->second );
}

TTree* TreeHandler::tree( const std::string& name ) {
  auto* pair = this->find_element( name );
  return pair->first;
}

TreeMapWrapper TreeHandler::map( const std::string& name ) {
  auto* pair = this->find_element( name );
  return TreeMapWrapper( &(pair->second) );
}

TreeMap& TreeHandler::access_map( const std::string& name ) {
  auto* pair = this->find_element( name );
  return pair->second;
}
