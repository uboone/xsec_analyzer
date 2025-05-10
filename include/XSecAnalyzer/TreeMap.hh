#pragma once

// Standard library includes
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include <vector>

// ROOT includes
#include "TTree.h"
#include "TVector3.h"

// A std::unique_ptr with redundant storage of a bare pointer to the managed
// object. This is a hacky workaround for setting TTree branch addresses for
// objects managed by a std::unique_ptr.
template <typename T> class MyPointer : public std::unique_ptr<T> {
  public:

    MyPointer() : std::unique_ptr<T>( new T ) {}

    T*& get_bare_ptr() {
      bare_ptr_ = this->get();
      return bare_ptr_;
    }

  protected:

    T* bare_ptr_ = nullptr;
};

// Create a type trait to allow "constexpr if" to distinguish between
// instantiations of the MyPointer class template and other types
template <typename T>
  struct is_MyPointer : std::false_type {};

template < typename T >
  struct is_MyPointer< MyPointer< T > > : std::true_type {};

template < typename T >
  constexpr bool is_MyPointer_v = is_MyPointer< T >::value;

// Create a type trait to allow "constexpr if" to distinguish between
// MyPointer< std::vector< T > > and other types
template <typename T>
  struct is_MyPointerToVector : std::false_type {};

template < typename T >
  struct is_MyPointerToVector< MyPointer< std::vector< T > > >
    : std::true_type {};

template < typename T >
  constexpr bool is_MyPointerToVector_v = is_MyPointerToVector< T >::value;

// A std::variant that can store any of the branch types of interest for
// the TTrees used by xsec_analyzer. New branches with a different type
// will require extending the template arguments given here.
using MyVariant = std::variant<
  // **** IMPORTANT: Keep std::monostate first in the template arguments ****
  // This is needed to ensure that default construction of MyVariant
  // leads to it holding the std::monostate type, which is used as a marker
  // in the framework that another type has not been assigned to the variant.
  std::monostate,
  // Now list the fundamental types that are allowed here before moving on to
  // objects
  unsigned char,
  char,
  unsigned short,
  short,
  unsigned int,
  int,
  unsigned long,
  long,
  unsigned long long,
  long long,
  bool,
  float,
  double,
  // Classes, structs, etc. that are allowed should be wrapped in
  // a MyPointer here
  MyPointer< std::string >,
  MyPointer< TVector3 >,
  MyPointer< std::vector< bool > >,
  MyPointer< std::vector< unsigned char > >,
  MyPointer< std::vector< char > >,
  MyPointer< std::vector< unsigned short > >,
  MyPointer< std::vector< short > >,
  MyPointer< std::vector< unsigned int > >,
  MyPointer< std::vector< int > >,
  MyPointer< std::vector< unsigned long > >,
  MyPointer< std::vector< long > >,
  MyPointer< std::vector< unsigned long long > >,
  MyPointer< std::vector< long long > >,
  MyPointer< std::vector< float > >,
  MyPointer< std::vector< double > >,
  MyPointer< std::vector< std::vector< double > > >,
  MyPointer< std::vector< TVector3 > >,
  MyPointer< std::map< std::string, std::vector< double > > >
>;



// Returns a T* to the content of a MyVariant if T is the actively-stored type,
// or nullptr otherwise
template < typename T > T* get_active_ptr_as( MyVariant& v ) {
  // For fundamental types, the MyVariant uses them directly, so just get a
  // pointer via std::get_if()
  if constexpr ( std::is_fundamental_v< T > ) {
    return std::get_if< T >( &v );
  }
  else {
    // For classes, structs, etc., we need to go inside the owned MyPointer
    auto* my_ptr = std::get_if< MyPointer< T > >( &v );
    if ( my_ptr ) {
      return my_ptr->get();
    }
  }
  return nullptr;
}

// Returns a T& to the content of a MyVariant if T is the actively-stored type.
// Throws an expression if it is not.
template < typename T > T& get_active_ref_as( MyVariant& v ) {
  auto* temp_ptr = get_active_ptr_as< T >( v );
  if ( !temp_ptr ) {
    throw std::runtime_error( "Invalid MyVariant dereference" );
  }
  return *temp_ptr;
}

// Determines T from the type of the second argument and
// copies the contents of the branch into the variable "out"
template< typename T > void copy_my_variant( MyVariant& v, T& out ) {
  out = get_active_ref_as< T >( v );
}

// Determines T from the type of the second argument and
// gets a bare T* to the branch variable owned by the variant
template< typename T > void get_my_variant( MyVariant& v, T*& out_ptr ) {
  out_ptr = get_active_ptr_as< T >( v );
}

// Overload the >> operator to allow easy loading of the active variant into a
// target variable. The type is inferred from the target without the
// need for explicit use of a template parameter.
template < typename T > void operator>>( MyVariant& v, T& out ) {
  if constexpr ( std::is_pointer_v< T > ) {
    get_my_variant( v, out );
  }
  else {
    copy_my_variant( v, out );
  }
}

// Manages a map used to automatically allocate storage for managing TTree
// branch variables using MyVariant objects. Keys in this map are strings giving
// the branch names (no duplicates allowed). Values in the map are MyVariant
// objects that will be used to store the branch data.
class TreeMap {
  public:

    using iterator = std::map< std::string, MyVariant >::iterator;
    using const_iterator = std::map< std::string, MyVariant >::const_iterator;

    TreeMap() = default;
    TreeMap( TTree* tree ) : tree_( tree ) {}

    MyVariant& at( const std::string& name ) { return map_.at( name ); }

    const MyVariant& at( const std::string& name ) const
      { return map_.at( name ); }

    MyVariant& operator[]( const std::string& name ) { return map_[ name ]; }

    TTree* tree() { return tree_; }

    const TTree* tree() const { return tree_; }

    template< class... Args > auto emplace( Args&&... args )
      { return map_.emplace( std::forward< Args >( args )... ); }

    auto begin() { return map_.begin(); }
    auto cbegin() const { return map_.cbegin(); }

    auto end() { return map_.end(); }
    auto cend() const { return map_.end(); }

    auto find( const std::string& key ) { return map_.find( key ); }

    MyVariant* add_input_branch( const std::string& name,
      TBranch*& added_branch );

    inline MyVariant* add_input_branch( const std::string& name ) {
      TBranch* dummy = nullptr;
      return add_input_branch( name, dummy );
    }

    void set_tree( TTree* tree ) { tree_ = tree; }

    const auto& var_size_map() const { return var_size_map_; }

    // Calls TTree::Write() for the owned TTree
    inline void write() { if ( tree_ ) tree_->Write(); }

    // Calls TTree::LoadTree() for the owned TTree
    inline long long load_tree( long long entry ) {
      if ( !tree_ ) return -1;
      return tree_->LoadTree( entry );
    }

    // Calls TTree::GetReadEntry() for the owned TTree
    inline long long get_read_entry() const {
      if ( !tree_ ) return -1;
      else return tree_->GetReadEntry();
    }

  protected:

    TTree* tree_ = nullptr;
    std::map< std::string, MyVariant > map_;
    std::map< std::string, std::set< std::string > > var_size_map_;
};
