// Standard library includes
#include <iostream>
#include <iterator>

// XSecAnalyzer includes
#include "XSecAnalyzer/TreeHandler.hh"

// Creates a TreeMap object to manage the branches of the input TTree,
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
  for ( auto& tm_pair : in_tree_maps_ ) {
    auto& tm = tm_pair.second;
    tm.get_entry( entry );
  }
}

// Call TTree::Fill() for each output TTree
void TreeHandler::fill() {
  for ( auto& tm_pair : out_tree_maps_ ) {
    TreeMap& tm = tm_pair.second;
    tm.fill();
  }
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

TreeMap::Wrapper TreeHandler::in_map( const std::string& name ) {
  auto* tm = this->find_element( name, true );
  return TreeMap::Wrapper( tm, name );
}

TreeMap& TreeHandler::access_in_map( const std::string& name ) {
  auto* tm = this->find_element( name, true );
  return *tm;
}

TTree* TreeHandler::out_tree( const std::string& name ) {
  auto* tm = this->find_element( name, false );
  return tm->tree();
}

TreeMap::Wrapper TreeHandler::out_map( const std::string& name ) {
  auto* tm = this->find_element( name, false );
  return TreeMap::Wrapper( tm, name );
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

const TreeMap::Wrapper& AnalysisEvent::in( const std::string& tree_name ) const
{
  const auto it = in_trees_.find( tree_name );
  if ( it != in_trees_.cend() ) return it->second;
  throw std::runtime_error( "Could not find input TTree named \"" + tree_name
    + "\" in call to AnalysisEvent::in()" );
}

// Access an input TreeMap by index
const TreeMap::Wrapper& AnalysisEvent::in( size_t index ) const
{
  if ( index >= in_trees_.size() ) {
    throw std::runtime_error( "AnalysisEvent::in() called with an out-of-range"
      " numerical index" );
  }
  auto it = in_trees_.begin();
  std::advance( it, index );
  return it->second;
}
