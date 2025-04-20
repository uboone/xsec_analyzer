// XSecAnalyzer includes
#include "XSecAnalyzer/TreeUtils.hh"

// Uses an input TTree to build a map of branch names each associated with a
// MyVariant object. These objects are linked to the TTree using
// TTree::SetBranchAddress() to provide temporary storage for event-by-event
// data processing.
TreeMap make_tree_map( TTree& p_tree ) {
  TreeMap branch_map;

  auto branch_list = p_tree.GetListOfBranches();
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
  
      // Initialize an appropriate variant in the map matching the class of interest
      // TODO: add more cases as needed below
      if ( class_name == "string" ) {
        emplace_variant_and_set_input_address< std::string >(
          branch_map, p_tree, branch_name );
      }
      else if ( class_name == "vector<bool>" ) {
        emplace_variant_and_set_input_address< std::vector< bool > >(
          branch_map, p_tree, branch_name );
      }
      else if ( class_name == "vector<double>" ) {
        emplace_variant_and_set_input_address< std::vector< double > >(
          branch_map, p_tree, branch_name );
      }
      else if ( class_name == "vector<float>" ) {
        emplace_variant_and_set_input_address< std::vector< float > >(
          branch_map, p_tree, branch_name );
      }
      else if ( class_name == "vector<int>" ) {
        emplace_variant_and_set_input_address< std::vector< int > >(
          branch_map, p_tree, branch_name );
      }
      else if ( class_name == "vector<unsigned int>" ) {
        emplace_variant_and_set_input_address< std::vector< unsigned int > >(
          branch_map, p_tree, branch_name );
      }
      else if ( class_name == "vector<unsigned long>" ) {
        emplace_variant_and_set_input_address< std::vector< unsigned long > >(
          branch_map, p_tree, branch_name );
      }
      else if ( class_name == "vector<unsigned short>" ) {
        emplace_variant_and_set_input_address< std::vector< unsigned short > >(
          branch_map, p_tree, branch_name );
      }
      else if ( class_name == "vector<vector<double> >" ) {
        emplace_variant_and_set_input_address< std::vector< std::vector< double > > >(
          branch_map, p_tree, branch_name );
      }
      else if ( class_name == "map<string,vector<double> >" ) {
        emplace_variant_and_set_input_address< std::map< std::string, std::vector< double > > >(
          branch_map, p_tree, branch_name );
      }
      else throw std::runtime_error( "Unrecognized class type \""
        + class_name + "\"" + " for branch \"" + branch_name + "\"" );
    }
    // Otherwise the branch has a simple data type indicated by an enum value.
    // Instantiate an appropriate variant for this type in the map.
    else {
      if ( temp_type == kUChar_t ) {
        emplace_variant_and_set_input_address< unsigned char >(
          branch_map, p_tree, branch_name );  
      }
      else if ( temp_type == kUInt_t ) {
        emplace_variant_and_set_input_address< unsigned int >(
          branch_map, p_tree, branch_name );  
      }
      else if ( temp_type == kBool_t ) {
        emplace_variant_and_set_input_address< bool >(
          branch_map, p_tree, branch_name );  
      }
      else if ( temp_type == kFloat_t
        || temp_type == kFloat16_t )
      {
       emplace_variant_and_set_input_address< float >(
          branch_map, p_tree, branch_name );  
      }
      else if ( temp_type == kInt_t ) {
        emplace_variant_and_set_input_address< int >(
          branch_map, p_tree, branch_name );  
      }
      else if ( temp_type == kDouble_t ) {
        emplace_variant_and_set_input_address< double >(
          branch_map, p_tree, branch_name );  
      }
      else throw std::runtime_error( "Unrecognized EDataType value "
        + std::to_string(temp_type) + " for branch \""
        + branch_name + "\"" );
    }
  }
    
  return branch_map;
}
