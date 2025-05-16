// Standard library includes
#include <regex>
#include <type_traits>

// XSecAnalyzer includes
#include "XSecAnalyzer/TreeMap.hh"

namespace {

  // Helper function for TreeMap::Formula::Error() that takes variadic arguments
  // and formats a std::string using them as input
  std::string format_variadic_string( const char* fmt, va_list args ) {
    // Make a copy of the arguments since vsnprintf consumes them
    std::va_list args_copy;
    va_copy( args_copy, args );

    // Determine the length of the string we need to format all of the
    // variadic arguments
    int length = std::vsnprintf( nullptr, 0, fmt, args_copy );
    va_end( args_copy );

    if ( length < 0 ) {
      throw std::runtime_error( "vsnprintf failed inside"
        " format_variadic_string()" );
    }

    std::string result( length, '\0' );
    // Length + 1 used in this function call includes the null terminator (+1)
    std::vsnprintf( result.data(), length + 1, fmt, args );
    return result;
  }

  // Emplaces either a scalar type or a std::vector for an array of that
  // type based on the boolean input from the user
  template < typename T > std::pair< TreeMap::iterator, bool >
    emplace_variant_helper( TreeMap& tree_map, TTree& in_tree,
      const std::string& branch_name, const Dimensions& dims )
  {
    bool is_scalar_branch = dims.empty();
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

    // Build a vector of size_t values to specify the Array dimensions in the
    // constructor below. Use a placeholder of 1 for variable-size dimensions.
    // Upon reading data from the TTree branch, these will be dynamically
    // resized as needed.
    std::vector< size_t > temp_dims;
    constexpr size_t PLACEHOLDER_DIM_SIZE = 1;
    for ( size_t d = 0u; d < dims.size(); ++d ) {
      if ( dims.is_fixed( d ) ) {
        temp_dims.push_back( dims.fixed_dim( d ) );
      }
      else {
        temp_dims.push_back( PLACEHOLDER_DIM_SIZE );
      }
    }

    // Use an Array< T > as the storage type in the variant
    auto er = tree_map.emplace( branch_name, Array< T >( temp_dims ) );

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
    auto* var_ptr = std::get_if< Array< T > >( &var );

    if ( !var_ptr ) {
      throw std::runtime_error( "Unsuccessful type retrieval for array variant"
        " representing the branch \"" + branch_name + "\"" );
      return er;
    }

    // Here's the trick we use for the C-style arrays: set the branch address
    // using a pointer to the first element in the std:vector owned by the
    // Array object in the variant. Since a std::vector has elements that are
    // guaranteed by the C++ standard to be contiguous in memory, this trick
    // will work as long as we check the array size before calling
    // TTree::GetEntry to make sure that the vector is large enough to hold the
    // array.
    in_tree.SetBranchAddress( branch_name.c_str(), var_ptr->data() );

    return er;
  }

  // Sets the branch address for output TTree variables in a unified way
  template < typename T, typename S > TBranch* set_branch(
    TTree& out_tree, S var_name, T*& var, bool output_locked )
  {
    std::string var_name_str( var_name );

    // Check if the branch already exists. If it doesn't, then try to create it.
    TBranch* br = out_tree.GetBranch( var_name_str.c_str() );
    if ( !br ) {

      // If the output TTrees are "locked" due to previously calling fill(),
      // then we should not add any more branches. Otherwise, we will end up
      // with some TTree entries having more branches than others, rendering
      // the output invalid. We therefore check here for the locked state.
      // If we find it, throw an exception.
      if ( output_locked ) {
        throw std::runtime_error( "Attempted to add branch \"" + var_name_str
          + "\" to the output TTree \"" + out_tree.GetName()
          + "\" after calling TreeHandler::fill()" );
      }

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
    void set_branch( TTree& out_tree, S var_name, T*&& rval_ptr,
      bool output_locked )
  {
    T* address = rval_ptr;
    set_branch( out_tree, var_name, address, output_locked );
  }

  // Overloaded version that takes a reference to the type T instead. This
  // version auto-detects a MyPointer argument and adjusts the approach
  // accordingly
  template < typename T, typename S >
    void set_branch( TTree& out_tree, S var_name, T& ref, bool output_locked )
  {
    if constexpr ( is_MyPointer_v< T > ) {
      auto*& address = ref.get_bare_ptr();
      set_branch( out_tree, var_name, address, output_locked );
    }
    // For a wrapped C-style array, set the branch to use the owned std::vector
    else if constexpr ( is_Array_v< T > ) {
      auto*& address = ref.vector();
      set_branch( out_tree, var_name, address, output_locked );
    }
    else {
      T* address = &ref;
      set_branch( out_tree, var_name, address, output_locked );
    }
  }

}

// Associates a branch name from an input TTree with a TreeMap::Variant object.
// This object is linked to the TTree using TTree::SetBranchAddress() to
// provide temporary storage for event-by-event data processing. Returns
// a pointer to the TreeMap::Variant object associated with the branch.
TreeMap::Variant* TreeMap::add_input_branch( const std::string& branch_name,
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

    // Parse the branch title to check for C-style array dimensions. These
    // are indicated in ROOT using square brackets [] in the title.
    auto dims = get_branch_dimensions( *added_branch );

    // Based on the reported type, allocate storage for the branch with
    // a variant and set the branch address
    if ( temp_type == kBool_t ) {

      bool is_scalar_branch = dims.empty();

      // Handle only the scalar case for now since std::vector< bool >
      // is not guaranteed to store its elements in a contiguous array.
      if ( !is_scalar_branch ) {
        std::cerr << "WARNING: Handling of C-style arrays of bool is not"
          << " yet implemented in TreeMap. Branch \"" << branch_name
          << "\" will be ignored.\n";
        return nullptr;
      }

      er = emplace_variant_and_set_input_address< bool >( *this, *tree_,
        branch_name );

      //er = emplace_variant_helper< bool >( *this, *tree_, branch_name,
      //  dims );
    }
    else if ( temp_type == kUChar_t ) {
      er = emplace_variant_helper< unsigned char >( *this, *tree_,
        branch_name, dims );
    }
    else if ( temp_type == kChar_t ) {
      er = emplace_variant_helper< char >( *this, *tree_,
        branch_name, dims );
    }
    else if ( temp_type == kUShort_t ) {
      er = emplace_variant_helper< unsigned short >( *this, *tree_,
        branch_name, dims );
    }
    else if ( temp_type == kShort_t ) {
      er = emplace_variant_helper< short >( *this, *tree_,
        branch_name, dims );
    }
    else if ( temp_type == kUInt_t ) {
      er = emplace_variant_helper< unsigned int >( *this, *tree_,
        branch_name, dims );
    }
    else if ( temp_type == kInt_t ) {
      er = emplace_variant_helper< int >( *this, *tree_,
        branch_name, dims );
    }
    else if ( temp_type == kULong_t ) {
      er = emplace_variant_helper< unsigned long >( *this, *tree_,
        branch_name, dims );
    }
    else if ( temp_type == kLong_t ) {
      er = emplace_variant_helper< long >( *this, *tree_,
        branch_name, dims );
    }
    else if ( temp_type == kULong64_t ) {
      er = emplace_variant_helper< unsigned long long >( *this, *tree_,
        branch_name, dims );
    }
    else if ( temp_type == kLong64_t ) {
      er = emplace_variant_helper< long long >( *this, *tree_,
        branch_name, dims );
    }
    else if ( temp_type == kFloat_t
      || temp_type == kFloat16_t )
    {
      er = emplace_variant_helper< float >( *this, *tree_,
        branch_name, dims );
    }
    else if ( temp_type == kDouble_t
      || temp_type == kDouble32_t )
    {
      er = emplace_variant_helper< double >( *this, *tree_,
        branch_name, dims );
    }
    else throw std::runtime_error( "Unrecognized EDataType value "
      + std::to_string(temp_type) + " for branch \""
      + branch_name + "\"" );

    // Check whether the current branch has at least one variable-size dimension
    if ( dims.dynamic() ) {
      // If it does, then save its dimensions to a map for later reuse in
      // dynamic resizing
      array_dim_map_[ branch_name ] = dims;

      // Ensure that we have access to all dynamic sizes by checking
      // for the presence of the corresponding branch names in the TreeMap.
      size_t num_dims = dims.size();
      for ( size_t d = 0u; d < num_dims; ++d ) {
        // For each variable array size, check if we already have storage for
        // its TTree branch in this TreeMap
        if ( !dims.is_fixed( d ) ) {
          const auto& d_name = dims.dim_name( d );

          // Cache the branch name in the set of branch names representing
          // dynamic array sizes
          dynamic_size_branch_names_.insert( d_name );

          auto iter = map_.find( d_name );
          // For the sizes without a loaded branch, use recursion to load them
          if ( iter == map_.end() ) this->add_input_branch( d_name );
        }
      }
    }

  }

  // Return a pointer to the TreeMap::Variant object storing the current branch
  return &er.first->second;
}

void TreeMap::get_entry( long long entry ) {
  // If we're working with a nullptr, just return without doing anything
  if ( !tree_ ) return;

  // Pre-load all branches that contain variable C-style array sizes
  for ( const auto& size_br_name : dynamic_size_branch_names_ ) {
    this->get_size_entry( entry, size_br_name );
  }

  for ( const auto& arr_dim_pair : array_dim_map_ ) {
    const std::string& arr_br_name = arr_dim_pair.first;
    const Dimensions& dims = arr_dim_pair.second;

    // Evaluate the current array dimensions, retrieving any dynamic sizes
    // from the TTree branches that store them
    std::vector< size_t > size_vec = this->get_dimensions( arr_br_name, dims );

    this->resize_variant_array( arr_br_name, size_vec );
  }

  // We're ready. Now retrieve the current entry reading from all enabled
  // TTree branches
  tree_->GetEntry( entry );
}

// Call TTree::Fill() for the owned TTree. Before doing so, create
// any missing branches and set branch addresses
void TreeMap::fill() {

  // If the pointer to the owned TTree is null, then just return without
  // doing anything
  if ( !tree_ ) return;

  // Iterate over the TreeMap to set the branch addresses before filling.
  for ( auto& tm_pair : map_ ) {
    const std::string& br_name = tm_pair.first;
    TreeMap::Variant& var = tm_pair.second;

    // Set the branch address (and create a new branch if needed).  Use
    // std::visit to automatically handle the type of the active variant.
    std::visit( [ this, &br_name ]( auto& var_val ) -> void
      {
        // Valid branches can be handled for all TreeMap::Variant types except
        // std::monostate, which is used to represent an uninitialized
        // variant
        using T = std::decay_t< decltype( var_val ) >;
        if constexpr ( !std::is_same_v< std::monostate, T > ) {
          set_branch( *tree_, br_name, var_val, output_locked_ );
        }
      }, var );
  }

  // We're done preparing everything. Fill the tree.
  tree_->Fill();

  // Lock the TTree to avoid adding new branches now that fill()
  // has been called at least once. This ensures a consistent branch
  // structure across all TTree entries;
  output_locked_ = true;
}

void TreeMap::Formula::Error( const char* location,
  const char* fmt, ... ) const
{
  std::va_list args;
  va_start( args, fmt );
  std::string formatted_args = format_variadic_string( fmt, args );
  va_end( args );

  // Format the final message with location
  std::string message( "Invalid TreeMapFormula led to error in"
    " TTreeFormula::" );
  message += location;
  message += ": " + formatted_args;

  throw std::runtime_error( message );
}

TreeMap::FormulaWrapper TreeMap::formula( const std::string& expr ) {
  // Search for the formula in the map
  auto it = formulas_.find( expr );
  // If there is no match, then cache a new formula object for later re-use
  if ( it == formulas_.end() ) {
    auto f = std::make_shared< Formula >( expr, *tree_ );
    it = formulas_.emplace( expr, f ).first;
  }
  // Provide access to the formula via a wrapper
  return FormulaWrapper( it->second.get() );
}

Dimensions TreeMap::get_branch_dimensions( const TBranch& br ) const {

  // General array dimension regex (i.e., square brackets [] with something
  // inside)
  const std::regex bracket_re( R"(\[([^\]]+)\])" );

  // Matches only digits, used to test whether the array dimension is a
  // fixed or variable value
  const std::regex digits_only_re( R"(^\d+$)" );

  std::string title = br.GetTitle();

  auto begin = title.cbegin();
  auto end   = title.cend();

  Dimensions dims;

  std::smatch match;
  while ( std::regex_search( begin, end, match, bracket_re ) ) {

    // Capture the next regex group inside square brackets []
    std::string token = match[ 1 ];

    if ( std::regex_match( token, digits_only_re ) ) {
      // The current dimension is a fixed size
      int cur_size = std::stoi( token );
      dims.add_dimension( cur_size );
    }
    else {
      // The current dimension is a variable size
      // TODO: add handling of variable sizes
      dims.add_dimension( token );
    }

    begin = match.suffix().first;  // move past this match
  }

  return dims;
}

// Build a vector of the current sizes along each dimension of an Array< T >.
// All of these are cast to size_t for compatibility with Array< T >::resize().
// If the last argument is set to a value other than DUMMY_ENTRY_NUMBER,
// then TBranch::GetEntry() is called to update any variable array sizes
// found before they are used. The default behavior skips this step
// and assumes that the stored values are already up-to-date.
std::vector< size_t > TreeMap::get_dimensions(
  const std::string& array_branch_name, const Dimensions& dims,
  long long entry )
{
  std::vector< size_t > size_vec;

  size_t num_dimensions = dims.size();

  for ( size_t d = 0u; d < num_dimensions; ++d ) {

    // For fixed dimensions, just add the size directly
    if ( dims.is_fixed( d ) ) {
      size_vec.push_back( dims.fixed_dim( d ) );
      continue;
    }

    // For variable dimensions, look up the stored Variant to retrieve
    // its current value. Cast the result to size_t to ensure compatibility
    // with a later call to Array< T >::resize()
    const std::string& size_br_name = dims.dim_name( d );

    // If requested, call TBranch::GetEntry() to update the dynamic array size
    // before retrieving it from the owned variant
    if ( entry != DUMMY_ENTRY_NUMBER ) {
      this->get_size_entry( entry, size_br_name );
    }
    auto& var = this->at( size_br_name );

    size_t dyn_size = std::visit( [ this,
      &size_br_name ]( auto& var_val ) -> size_t
      {
        using T = std::decay_t< decltype( var_val ) >;
        if constexpr ( std::is_convertible_v< T, size_t > ) {
          return static_cast< size_t >( var_val );
        }
        throw std::runtime_error( "Variant storing a dynamic array size"
          " from the tree branch \"" + size_br_name + "\" is not convertible"
          " to size_t" );
        return 0u;
      }, var );

    size_vec.push_back( dyn_size );
  }

  return size_vec;
}

void TreeMap::get_size_entry( long long entry,
  const std::string& size_br_name )
{
  TBranch* size_br = tree_->GetBranch( size_br_name.c_str() );
  if ( !size_br ) throw std::runtime_error( "Missing array size branch \""
    + size_br_name + "\" encountered" );
  else size_br->GetEntry( entry );
}

void TreeMap::resize_variant_array( const std::string& arr_br_name,
  const std::vector< size_t >& size_vec )
{
  // Using the input vector of array sizes along each dimension, resize
  // an Array< T > used for storage of a C-style array. Use std::visit
  // to expose the contents of the Variant with the correct type.
  auto& array_variant = this->at( arr_br_name );
  std::visit( [ this, &size_vec, &arr_br_name ]( auto& var_val ) -> void
    {
      using T = std::decay_t< decltype( var_val ) >;
      if constexpr ( is_Array_v< T > ) {
        var_val.resize( size_vec );

        // Resizing can invalidate pointers to std::vector elements,
        // so update the address for the array branch
        tree_->SetBranchAddress( arr_br_name.c_str(), var_val.data() );

        return;
      }
      throw std::runtime_error( "Variant storing a dynamic array"
        " from the tree branch \"" + arr_br_name
        + "\" could not be resized" );
    }, array_variant );
}
