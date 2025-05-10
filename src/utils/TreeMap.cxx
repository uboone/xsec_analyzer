// ROOT includes
#include "TLeaf.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/TreeMap.hh"
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

}

// Associates a branch name from an input TTree with a MyVariant object.
// This object is linked to the TTree using TTree::SetBranchAddress() to
// provide temporary storage for event-by-event data processing. Returns
// a pointer to the MyVariant object associated with the branch.
MyVariant* TreeMap::add_input_branch( const std::string& branch_name,
  TBranch*& added_branch )
{

  // Check whether the input branch has already been added to the TreeMap.
  // If it has, just return rather than trying to add it again.
  auto iter = this->find( branch_name );
  if ( iter != this->end() ) {
    return &iter->second;
  }

  // Complain if the tree pointer is null
  if ( !tree_ ) {
    throw std::runtime_error( "Null tree pointer encountered in"
      " TreeMap::add_input_branch()" );
    return nullptr;
  }

  // Temporary variables used to query the branch about its data type
  // via the TBranch::GetExpectedType() function
  TClass* temp_class = nullptr;
  EDataType temp_type;

  // Retrieve data type information for the current branch
  added_branch = tree_->GetBranch( branch_name.c_str() );
  if ( !added_branch ) {
    std::string tree_name = tree_->GetName();
    throw std::runtime_error( "Missing branch \"" + branch_name + "\" in"
      " the \"" + tree_name + "\" tree" );
  }
  added_branch->GetExpectedType( temp_class, temp_type );

  // Get the result of the emplace operation, which will include an
  // iterator to the relevant TreeMap element and a bool indicating whether
  // a new element was added
  std::pair< TreeMap::iterator, bool > er;

  // If the TClass pointer is non-null, the branch corresponds to an object
  if ( temp_class ) {
    // Get the name of the class stored in the current branch
    std::string class_name = temp_class->GetName();

    // Initialize an appropriate variant in the map matching the class of
    // interest
    // TODO: add more cases as needed below
    if ( class_name == "string" ) {
      er = emplace_variant_and_set_input_address< std::string >(
        *this, *tree_, branch_name );
    }
    else if ( class_name == "TVector3" ) {
      er = emplace_variant_and_set_input_address< TVector3 >(
        *this, *tree_, branch_name );
    }
    else if ( class_name == "vector<bool>" ) {
      er = emplace_variant_and_set_input_address< std::vector< bool > >(
        *this, *tree_, branch_name );
    }
    else if ( class_name == "vector<double>" ) {
      er = emplace_variant_and_set_input_address< std::vector< double > >(
        *this, *tree_, branch_name );
    }
    else if ( class_name == "vector<float>" ) {
      er = emplace_variant_and_set_input_address< std::vector< float > >(
        *this, *tree_, branch_name );
    }
    else if ( class_name == "vector<int>" ) {
      er = emplace_variant_and_set_input_address< std::vector< int > >(
        *this, *tree_, branch_name );
    }
    else if ( class_name == "vector<unsigned int>" ) {
      er = emplace_variant_and_set_input_address<
        std::vector< unsigned int > >( *this, *tree_, branch_name );
    }
    else if ( class_name == "vector<unsigned long>" ) {
      er = emplace_variant_and_set_input_address<
        std::vector< unsigned long > >( *this, *tree_, branch_name );
    }
    else if ( class_name == "vector<unsigned short>" ) {
      er = emplace_variant_and_set_input_address<
        std::vector< unsigned short > >( *this, *tree_, branch_name );
    }
    else if ( class_name == "vector<vector<double> >" ) {
      er = emplace_variant_and_set_input_address< std::vector<
        std::vector< double > > >( *this, *tree_, branch_name );
    }
    else if ( class_name == "vector<TVector3>" ) {
      er = emplace_variant_and_set_input_address< std::vector< TVector3 > >(
        *this, *tree_, branch_name );
    }
    else if ( class_name == "map<string,vector<double> >" ) {
      er = emplace_variant_and_set_input_address< std::map<
        std::string, std::vector< double > > >( *this, *tree_,
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
    TLeaf* lf = dynamic_cast< TLeaf* >(
      added_branch->GetListOfLeaves()->First() );

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
      return nullptr;
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
    if ( temp_type == kBool_t ) {
      // Handle only the scalar case for now since std::vector< bool >
      // is not guaranteed to store its elements in a contiguous array.

      if ( !is_scalar_branch ) {
        std::cerr << "WARNING: Handling of C-style arrays of bool is not"
          << " yet implemented in TreeHandler. Branch \"" << branch_name
          << "\" will be ignored.\n";
        return nullptr;
      }

      er = emplace_variant_and_set_input_address< bool >( *this, *tree_,
        branch_name );

      //er = emplace_variant_helper< bool >( *this, *tree_, branch_name,
      //  is_scalar_branch );
    }
    else if ( temp_type == kUChar_t ) {
      er = emplace_variant_helper< unsigned char >( *this, *tree_,
        branch_name, is_scalar_branch );
    }
    else if ( temp_type == kChar_t ) {
      er = emplace_variant_helper< char >( *this, *tree_,
        branch_name, is_scalar_branch );
    }
    else if ( temp_type == kUShort_t ) {
      er = emplace_variant_helper< unsigned short >( *this, *tree_,
        branch_name, is_scalar_branch );
    }
    else if ( temp_type == kShort_t ) {
      er = emplace_variant_helper< short >( *this, *tree_,
        branch_name, is_scalar_branch );
    }
    else if ( temp_type == kUInt_t ) {
      er = emplace_variant_helper< unsigned int >( *this, *tree_,
        branch_name, is_scalar_branch );
    }
    else if ( temp_type == kInt_t ) {
      er = emplace_variant_helper< int >( *this, *tree_,
        branch_name, is_scalar_branch );
    }
    else if ( temp_type == kULong_t ) {
      er = emplace_variant_helper< unsigned long >( *this, *tree_,
        branch_name, is_scalar_branch );
    }
    else if ( temp_type == kLong_t ) {
      er = emplace_variant_helper< long >( *this, *tree_,
        branch_name, is_scalar_branch );
    }
    else if ( temp_type == kULong64_t ) {
      er = emplace_variant_helper< unsigned long long >( *this, *tree_,
        branch_name, is_scalar_branch );
    }
    else if ( temp_type == kLong64_t ) {
      er = emplace_variant_helper< long long >( *this, *tree_,
        branch_name, is_scalar_branch );
    }
    else if ( temp_type == kFloat_t
      || temp_type == kFloat16_t )
    {
      er = emplace_variant_helper< float >( *this, *tree_,
        branch_name, is_scalar_branch );
    }
    else if ( temp_type == kDouble_t
      || temp_type == kDouble32_t )
    {
      er = emplace_variant_helper< double >( *this, *tree_,
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
        std::visit( [ this, &branch_name, count_result ](
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
              tree_->SetBranchAddress( branch_name.c_str(),
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

        // If this is our first time encountering this array size branch,
        // then set up a MyVariant object to store its value using recursion.
        auto size_iter = var_size_map_.find( size_leaf_name );
        if ( size_iter == var_size_map_.end() ) {
          this->add_input_branch( size_leaf_name );
        }

        // Store the name of the current branch that needs to be dynamically
        // resized for later use in get_entry().
        var_size_map_[ size_leaf_name ].insert( branch_name );
      }
    }
  }

  // Return a pointer to the MyVariant object storing the current branch
  return &er.first->second;
}
