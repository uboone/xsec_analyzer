// Standard library includes
#include <iostream>
#include <iterator>

// XSecAnalyzer includes
#include "XSecAnalyzer/TreeHandler.hh"

namespace {

  // Sets the branch address for output TTree variables in a unified way
  template < typename T, typename S > TBranch* set_branch(
    TTree& out_tree, S var_name, T*& var, bool output_locked )
  {
    std::string var_name_str( var_name );

    // Check if the branch already exists. If it doesn't, then try to create
    // it.
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
    else {
      T* address = &ref;
      set_branch( out_tree, var_name, address, output_locked );
    }
  }

}

// Uses an input TTree to build a map of branch names each associated with a
// MyVariant object. These objects are linked to the TTree using
// TTree::SetBranchAddress() to provide temporary storage for event-by-event
// data processing.
void TreeHandler::add_input_tree( TTree* in_tree, const std::string& name,
  bool load_all_branches )
{
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

  // Check for an exact match for this TTree (based on pointer value)
  // among the already-loaded TTree objects. If one is found, this function
  // will just return while noting the duplicate.
  for ( const auto& tree_pair : in_tree_maps_ ) {
    const TTree* old_tree = tree_pair.second.tree();
    if ( in_tree == old_tree ) {
      std::cout << "Skipping already-loaded TTree \"" + key + "\"\n";
      return;
    }
  }

  // If other trees have already been loaded, double-check that the
  // new one has the same number of entries. If there is a mismatch,
  // the trees will not be compatible. Complain if a problem is found.
  if ( !in_tree_maps_.empty() ) {
    long long num_events = in_tree->GetEntries();
    for ( auto pair_iter = in_tree_maps_.cbegin();
      pair_iter != in_tree_maps_.cend(); ++pair_iter )
    {
      const TTree* old_tree = pair_iter->second.tree();
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
  auto& tree_map = in_tree_maps_[ key ];

  // Set the TTree pointer in the new map entry
  tree_map.set_tree( in_tree );

  // If the user hasn't requested to automatically load all branches,
  // then we're done
  if ( !load_all_branches ) return;

  auto branch_list = in_tree->GetListOfBranches();
  int num_branches = branch_list->GetEntries();

  for ( int b = 0; b < num_branches; ++b ) {
    auto* br = dynamic_cast< TBranch* >( branch_list->At(b) );
    std::string branch_name = br->GetName();
    tree_map.add_input_branch( branch_name );
  }
}

// Call TTree::GetEntry() for each input TTree
void TreeHandler::get_entry( long long entry ) {
  for ( auto pair_iter = in_tree_maps_.begin();
    pair_iter != in_tree_maps_.end(); ++pair_iter )
  {
    // First get access to the current TTree and
    // a map storing information about its variable-size arrays
    TreeMap& tm = pair_iter->second;
    TTree* temp_tree = tm.tree();
    if ( !temp_tree ) continue;

    const auto& tree_size_map = tm.var_size_map();

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

    // Get access to the TreeMap used to manage the current output TTree
    TreeMap& tm = pair_iter->second;

    // If the pointer to the tree it owns is null, then just move on to
    // the next one.
    TTree* temp_tree = tm.tree();
    if ( !temp_tree ) {
      continue;
    }

    // Iterate over the TreeMap to set the branch addresses before filling.
    for ( auto& tm_pair : tm ) {
      const std::string& br_name = tm_pair.first;
      MyVariant& var = tm_pair.second;

      // Set the branch address (and create a new branch if needed).  Use
      // std::visit to automatically handle the type of the active variant.
      std::visit( [ &temp_tree, &br_name, this ]( auto& var_val ) -> void
        {
          // Valid branches can be handled for all MyVariant types except
          // std::monostate, which is used to represent an uninitialized
          // variant
          using T = std::decay_t< decltype( var_val ) >;
          if constexpr ( !std::is_same_v< std::monostate, T > ) {
            set_branch( *temp_tree, br_name, var_val, output_locked_ );
          }
        }, var );
    }

    // We're done preparing everything. Fill the tree.
    temp_tree->Fill();
  }

  // Lock the output TTrees to avoid adding new branches now that fill()
  // has been called at least once. This ensures a consistent branch
  // structure across all output TTree entries.
  output_locked_ = true;
}

// Call TTree::Write() for each output TTree
void TreeHandler::write() {
  for ( auto& pair : out_tree_maps_ ) pair.second.write();
}

// Call TTree::LoadTree() for each input TTree
long long TreeHandler::load_tree( long long entry ) {
  long long result = -1;
  for ( auto& pair : in_tree_maps_ ) result = pair.second.load_tree( entry );
  return result;
}

TreeMap* TreeHandler::find_element( const std::string& name,
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

    throw std::runtime_error( err_message + " tree has the name \""
      + name + "\"" );
    return nullptr;
  }
  return &( iter->second );
}

TTree* TreeHandler::in_tree( const std::string& name ) {
  auto* tm = this->find_element( name, true );
  return tm->tree();
}

TreeMapWrapper TreeHandler::in_map( const std::string& name ) {
  auto* tm = this->find_element( name, true );
  return TreeMapWrapper( tm, name );
}

TreeMap& TreeHandler::access_in_map( const std::string& name ) {
  auto* tm = this->find_element( name, true );
  return *tm;
}

TTree* TreeHandler::out_tree( const std::string& name ) {
  auto* tm = this->find_element( name, false );
  return tm->tree();
}

TreeMapWrapper TreeHandler::out_map( const std::string& name ) {
  auto* tm = this->find_element( name, false );
  return TreeMapWrapper( tm, name );
}

TreeMap& TreeHandler::access_out_map( const std::string& name ) {
  auto* tm = this->find_element( name, false );
  return *tm;
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
  auto& tree_map = out_tree_maps_[ name ];

  // Initialize the new TTree and associate it with the output TDirectory
  TTree* out_tree = new TTree( name.c_str(), title.c_str() );
  out_tree->SetDirectory( out_dir );

  // Associate the TTree with the TreeMap
  tree_map.set_tree( out_tree );
}

AnalysisEvent TreeHandler::get_event(
  const std::set< std::string >& input_names, const std::string& output_name )
{
  auto out_tmw = this->out_map( output_name );
  AnalysisEvent event( out_tmw );

  for ( const auto& in_name : input_names ) {
    event.add_input( in_name, this->in_map( in_name ) );
  }

  return event;
}

const TreeMapWrapper& AnalysisEvent::in( const std::string& tree_name ) const
{
  const auto it = in_trees_.find( tree_name );
  if ( it != in_trees_.cend() ) return it->second;
  throw std::runtime_error( "Could not find input TTree named \"" + tree_name
    + "\" in call to AnalysisEvent::in()" );
}

// Access an input TreeMap by index
const TreeMapWrapper& AnalysisEvent::in( size_t index ) const
{
  if ( index >= in_trees_.size() ) {
    throw std::runtime_error( "AnalysisEvent::in() called with an out-of-range"
      " numerical index" );
  }
  auto it = in_trees_.begin();
  std::advance( it, index );
  return it->second;
}
