// ROOT includes
#include "TCollection.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/DirectoryMap.hh"

void DirectoryMap::save() const {
  for ( const auto& pair : params_ ) {
    const std::string& key_name = pair.first;
    const auto& variant = pair.second;

    std::visit( [ this, &key_name ]( auto& var ) -> void {
      // Don't do anything if the owned TDirectory pointer is null
      if ( !dir_ ) return;

      using T = std::decay_t< decltype( var ) >;
      // Ignore the parameter if the active variant is std::monostate.
      // This corresponds to an uninitialized entry in the parameter map.
      if constexpr ( std::is_same_v< T, std::monostate > ) return;
      // If the active variant is a string, handle it with TNamed
      else if constexpr ( std::is_same_v< T, std::string > ) {
        TNamed temp_named( key_name.c_str(), var.c_str() );
        dir_->WriteObject( &temp_named, key_name.c_str(), "WriteDelete" );
      }
      // Otherwise, use a TParameter of the appropriate type
      else {
        TParameter< T > temp_param( key_name.c_str(), var );
        dir_->WriteObject( &temp_param, key_name.c_str(), "WriteDelete" );
      }
    }, variant );
  }
}


void DirectoryMap::load() {
  if ( !dir_ ) return;
  params_.clear();
  TIter next_key( dir_->GetListOfKeys() );
  while ( auto key = static_cast< TKey* >(next_key()) ) {
    std::string key_name = key->GetName();
    std::string class_name = key->GetClassName();
    auto iter = loaders_.find( class_name );
    if ( iter != loaders_.end() ) iter->second( key_name );
  }
}
