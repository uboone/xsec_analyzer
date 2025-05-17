#pragma once

#include <vector>
#include <stdexcept>

// A dynamically-sized array with one or more dimensions. A one-dimensional
// std::vector is used to store the array elements
template< typename T > class Array {

  public:

    using value_type = T;

    Array( const std::vector< size_t >& dims ) {
      data_ptr_ = &data_;
      this->resize( dims );
    }

    void resize( const std::vector< size_t >& dims )
    {
      if ( dims.empty() ) {
        throw std::invalid_argument( "Empty vector of dimensions passed"
          " to Array< T >::resize()" );
      }
      dims_ = dims;

      // Compute strides for each dimension as well as the total number of
      // array elements
      strides_.resize( dims_.size() );
      size_t total_size = 1u;
      for ( size_t i = dims_.size(); i-- > 0; ) {
        strides_[ i ] = ( i + 1 < dims_.size() )
          ? strides_[ i + 1 ] * dims_[ i + 1 ] : 1;
        total_size *= dims_[ i ];
      }
      data_.resize( total_size );
    }

    // Returns a vector of the array sizes for each dimension
    const std::vector< size_t >& dims() const { return dims_; }

    // Returns the total number of elements in the array
    const size_t size() const { return data_.size(); }

    // Get a raw pointer to the one-dimensional array storing all elements
    T* data() { return data_.data(); }
    const T* data() const { return data_.data(); }

    // Get a raw pointer to the owned vector object
    std::vector< T >*& vector() { return data_ptr_; }
    const std::vector< T >*& vector() const { return data_ptr_; }

    // Directly access an element of the array by providing all needed indices
    template< typename... Indices > T& operator()( Indices... indices ) {

      // Check that all indices passed in can be interpreted as size_t
      static_assert( ( std::is_convertible<Indices, size_t>::value && ... ),
        "All indices passed to ArrayWrapper::operator() must be convertible"
        " to size_t" );

      // Build an array of the indices, converting via static_cast to have
      // size_t type in all cases
      std::array< size_t, sizeof...( Indices ) >
        idx = { static_cast< size_t >( indices )... };

      // Check that the number of indices provided matches the expectation
      if ( idx.size() != dims_.size() ) {
        size_t i_size = idx.size();
        size_t d_size = dims_.size();
        throw std::runtime_error( "Incorrect number of indices passed"
          " to ArrayWrapper::operator(). " + std::to_string( d_size )
          + " indices were expected but " + std::to_string( i_size )
          + " were provided." );
      }

      // Check bounds across all dimensions
      for ( size_t d = 0; d < idx.size(); ++d ) {
        size_t temp_idx = idx[ d ];
        size_t temp_dim = dims_[ d ];
        if ( temp_idx >= temp_dim ) {
          throw std::out_of_range( "Index " + std::to_string( temp_idx )
            + " passed to ArrayWrapper::operator()() is out of bounds in"
            + " dimension " + std::to_string( d ) + ", which has size "
            + std::to_string( temp_dim ) );
        }
      }

      // Compute the index for the desired element in the one-dimensional
      // vector used to store the multidimensional data. Row-major ordering
      // is assumed, which is the same as multidimensional C-style arrays.
      size_t global_idx = 0u, stride = 1u;
      // Loop backwards over the dimensions. Note that the postfix decrement
      // operator is needed here (a prefix would yield different behavior).
      for ( size_t i = dims_.size(); i-- > 0u; ) {
        global_idx += idx[ i ] * stride;
        stride *= dims_[ i ];
      }
      return data_[ global_idx ];
    }

    // Build a top‐level view of the data to allow C-style array access
    class View;

    View view() { return { this->data(), dims_.data(),
      strides_.data(), 0, dims_.size() }; }

    View view() const { return { this->data(), dims_.data(),
      strides_.data(), 0, dims_.size() }; }

    // If you use square brackets to access a sub-array, build a top-level view
    // and then call its own operator[]()
    View operator[]( size_t i ) {
      return this->view()[ i ];
    }

    View at( size_t i ) { return this->operator[]( i ); }

  protected:

    std::vector< size_t > dims_, strides_;
    std::vector< T > data_;
    std::vector< T >* data_ptr_; // Bare pointer for ROOT TTree branch handling
};

// Alias for the nested View class template
template< typename T > using ArrayView = typename Array< T >::View;

// Provides non-owning access to sub-arrays and individual elements via
// C-style use of operator[]()
template< typename T > class Array< T >::View {

public:

  using value_type = T;

  View() {}

  View( T* data, const size_t* dims, const size_t* strides,
    size_t offset, size_t dims_left ) : data_( data ), dims_( dims ),
    strides_( strides ), offset_( offset ), dims_left_( dims_left ) {}

  View at( size_t i ) const {
    return this->operator[]( i );
  }

  View operator[]( size_t i ) const {
    if ( dims_left_ == 0 ) {
      throw std::out_of_range( "Array dimensions exceeded in call to"
        " Array< T >::View::operator[]()" );
    }
    if ( i >= dims_[ 0 ] ) {
      throw std::out_of_range( "Index " + std::to_string( i )
        + " passed to Array< T >::View::operator[]() exceeds "
        + " the size " + std::to_string( dims_[ 0 ] ) + " of the"
        + " current dimension" );
    }

    return {
      data_,
      dims_ + 1,
      strides_ + 1,
      offset_ + i * strides_[ 0 ],
      dims_left_ - 1
    };
  }

  operator T&() const {
    if ( dims_left_ != 0 ) {
      throw std::logic_error( "Cannot convert a multi-element"
        " Array< T >::View to T&" );
    }
    return *( data_ + offset_ );
  }

  View& operator=( const T& value ) {
    if ( dims_left_ != 0 ) {
      throw std::logic_error("Cannot assign a single T value to a multi-element"
        " Array< T >::View" );
    }
    *( data_ + offset_ ) = value;
    return *this;
  }

  // Returns a vector of the array sizes for each dimension accessible
  // to the view
  const std::vector< size_t > dims() const {
    std::vector< size_t > temp_dims;
    for ( size_t d = 0u; d < dims_left_; ++d ) {
      temp_dims.push_back( dims_[ d ] );
    }
    return temp_dims;
  }

  // Returns the total number of elements in the array
  const size_t size() const {
    size_t result = 1u;
    for ( size_t d = 0u; d < dims_left_; ++d ) result *= dims_[ d ];
    return result;
  }

protected:

  T* data_ = nullptr;
  const size_t* dims_ = nullptr;
  const size_t* strides_ = nullptr;
  size_t offset_ = 0u;
  size_t dims_left_ = 0u;
};

// Container used to store dimension information for managing TTree branches
// that store C-style arrays wrapped by an Array< T >
class Dimensions {
  public:

    static constexpr int DUMMY_DIM_SIZE = -1;

    Dimensions() {}

    bool empty() const { return fixed_dims_.empty(); }
    size_t size() const { return fixed_dims_.size(); }

    void add_dimension( int size, const std::string& name = "" ) {
      if ( size != DUMMY_DIM_SIZE && size <= 0 ) {
        throw std::runtime_error( "Invalid size " + std::to_string( size )
          + " passed to Array< T >::Dimensions::add_dimension()" );
        return;
      }
      fixed_dims_.push_back( size );
      dim_names_.push_back( name );
    }

    void add_dimension( const std::string& name ) {
      this->add_dimension( DUMMY_DIM_SIZE, name );
    }

    // Checks whether the i-th dimension has a fixed or variable size. Throws
    // an exception if the value of i is equal to or greater than the current
    // number of dimensions
    bool is_fixed( size_t i ) const {
      this->check_bounds( i );
      return fixed_dims_.at( i ) != DUMMY_DIM_SIZE;
    }

    // Returns true if at least one dimension of the array has a variable size
    // that requires dynamic resizing
    bool dynamic() const {
      for ( const auto& d : fixed_dims_ ) {
        if ( d == DUMMY_DIM_SIZE ) return true;
      }
      return false;
    }

    int fixed_dim( size_t i ) const {
      this->check_bounds( i );
      return fixed_dims_.at( i );
    }

    const std::string& dim_name( size_t i ) const {
      this->check_bounds( i );
      return dim_names_.at( i );
    }

    std::pair< int, std::string > dim( size_t i ) const {
      this->check_bounds( i );
      return std::make_pair( fixed_dims_.at( i ), dim_names_.at( i ) );
    }

  protected:

    // Helper function for vector bounds checking
    void check_bounds( size_t i ) const {
      if ( i >= fixed_dims_.size() ) throw std::runtime_error( "Requested"
        " index " + std::to_string( i ) + " exceeds the current number of"
        " array dimensions in a Array< T >::Dimensions object" );
    }

    // Vector giving numerical sizes of fixed dimensions. Variable-sized
    // dimensions are marked by the special value DUMMY_DIM_SIZE.
    std::vector< int > fixed_dims_;

    // Vector containing names for each dimension. For variable-sized
    // dimensions, this name is interpreted as the name of a TTree branch
    // storing the current value of the dynamic size.
    std::vector< std::string > dim_names_;
};
