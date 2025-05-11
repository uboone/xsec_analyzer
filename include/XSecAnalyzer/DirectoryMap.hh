#pragma once

// Standard library includes
#include <functional>
#include <map>
#include <string>
#include <unordered_map>
#include <variant>

// ROOT includes
#include "TDirectory.h"
#include "TKey.h"
#include "TNamed.h"
#include "TParameter.h"

class DirectoryMap {
  public:

    using Variant = std::variant< std::monostate, bool, int, long, long long,
      float, double, std::string >;

    using Parameters = std::map< std::string, Variant >;

    // Expose a lot of the functionality from the owned map
    using iterator = Parameters::iterator;
    using const_iterator = Parameters::const_iterator;

    template< class... Args > auto emplace( Args&&... args )
      { return params_.emplace( std::forward< Args >( args )... ); }

    auto begin() { return params_.begin(); }
    auto cbegin() const { return params_.cbegin(); }

    auto end() { return params_.end(); }
    auto cend() const { return params_.end(); }

    auto find( const std::string& key ) { return params_.find( key ); }

    const auto find( const std::string& key ) const
      { return params_.find( key ); }

    // Wrapper class exposing owned Variant objects in a way that allows
    // for easy assignment and retrieval of values
    class VariantWrapper {
    
      public:
    
        VariantWrapper( Variant* v ) : variant_( v ) {}
    
        // Overload the = operator to allow setting new values of the variant
        template < typename T > VariantWrapper& operator=( const T& in ) {
          *variant_ = in;
          return *this;
        }

        // Overload the >> operator to allow easy loading of the active variant
        // into a target variable. The type is inferred from the target without
        // the need for explicit use of a template parameter.
        template < typename T > const VariantWrapper &operator>>( T& out )
          const
        {
          // Set a pointer argument to point to the active variant
          if constexpr ( std::is_pointer_v< T > ) {
            out = std::get_if< T >( variant_ );
          }
          // Otherwise copy the value of the active variant
          else {
            auto* temp_ptr = std::get_if< T >( variant_ );
            if ( !temp_ptr ) {
              throw std::runtime_error( "Invalid DirectoryMap"
                "::Variant dereference" );
            }
            else out = *temp_ptr;
          }
          return *this;
        }

      protected:
    
        Variant* variant_ = nullptr;
    };

    VariantWrapper at( const std::string& name )
      { return VariantWrapper( &params_.at( name ) ); }

    const VariantWrapper at( const std::string& name ) const {
      return VariantWrapper( const_cast< Variant* >( &params_.at( name ) ) );
    }

    VariantWrapper operator[]( const std::string& name )
      { return VariantWrapper( &params_[ name ] ); }

    // Save all defined parameter values to the TDirectory
    void save() const;

    // General template for loading parameters stored in the TDirectory
    template < typename T > bool load_param( const std::string& key_name ) {
      if ( !dir_ ) return false;
      auto* param = dir_->Get< TParameter< T > >( key_name.c_str() );
      if ( !param ) return false;
    
      params_[ key_name ] = param->GetVal();
      return true;
    }

    // Specialization for std::string, which uses TNamed rather than TParameter
    template <> bool load_param< std::string >( const std::string& key_name ) {
      if ( !dir_ ) return false;
      auto* param = dir_->Get< TNamed >( key_name.c_str() );
      if ( !param ) return false;

      params_[ key_name ] = param->GetTitle();
      return true;
    }

    // Type used to associate lambda-wrapped specializations of load_param()
    // with class names for easy loading logic.
    using Loaders = std::unordered_map< std::string,
      std::function< bool( const std::string& key_name ) > >;

    DirectoryMap( TDirectory* dir ) : dir_( dir ) {
      this->load();
    }

    // Clear any existing parameter definitions and load parameter
    // values from the owned TDirectory
    void load();

  protected:

    TDirectory* dir_ = nullptr;
    Parameters params_;

    // Convenience macro for writing the lambda functions used to load
    // parameter values from the owned TDirectory. Templated lambdas
    // could also work, but that feature is not available until C++20
    #define LOAD(type) [ this ]( const std::string& key_name ) -> bool \
      { return load_param< type >( key_name ); }

    Loaders loaders_ = {
      { "TParameter<bool>", LOAD(bool) }, { "TParameter<int>", LOAD(int) },
      { "TParameter<long>", LOAD(long) },
      { "TParameter<Long64_t>", LOAD(long long) },
      { "TParameter<float>", LOAD(float) },
      { "TParameter<double>", LOAD(double) },
      { "TNamed", LOAD(std::string) },
    };

    #undef LOAD
};
