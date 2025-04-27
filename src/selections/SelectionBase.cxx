// Standard library includes
#include <fstream>
#include <iostream>
#include <regex>
#include <set>
#include <sstream>
#include <vector>

// ROOT includes
#include "TColor.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/Selections/SelectionBase.hh"
#include "XSecAnalyzer/TreeHandler.hh"

// Define this static member of the SelectionBase class
std::unique_ptr< SelectionBase::SelectionConfig >
  SelectionBase::selection_config_;

SelectionBase::SelectionBase( const std::string& sel_name ) {

  selection_name_ = sel_name;
  num_passed_events_ = 0;

  if ( !selection_config_ ) this->load_selection_config();
}

const std::pair< std::shared_ptr< SelectionBase::TreeNameSet >,
  std::shared_ptr< SelectionBase::CategoryMap > >&
  SelectionBase::load_info() const
{
  if ( !selection_config_ ) {
    throw std::runtime_error( "Missing selection configurations encountered"
      " in SelectionBase::load_info()" );
  }

  if ( !selection_config_->count( selection_name_ ) ) {
    throw std::runtime_error( "Missing configuration for selection \""
      + selection_name_ + "\" encountered in SelectionBase::load_info()" );
  }

  const auto& result = selection_config_->at( selection_name_ );
  return result;
}

const SelectionBase::CategoryMap& SelectionBase::category_map() const {
  const auto& info = this->load_info();
  const auto& categ_map = info.second;
  if ( !categ_map ) throw std::runtime_error( "Null category map"
    " for selection \"" + selection_name_ + '\"' );
  return *categ_map;
}

const SelectionBase::TreeNameSet& SelectionBase::input_tree_names() const {
  const auto& info = this->load_info();
  const auto& tree_names = info.first;
  if ( !tree_names ) {
    throw std::runtime_error( "Null input TTree names for"
      " selection \"" + selection_name_ + '\"' );
  }
  return *tree_names;
}

void SelectionBase::apply_selection( bool is_mc, AnalysisEvent& event ) {

  // The operations in this block rely on the presence of MC truth information.
  // Set some variables to default values in case we are working with real
  // data.
  bool mc_signal = false;
  int event_category = SelectionBase::UNKNOWN_CATEGORY_CODE;
  if ( is_mc ) {
    const std::string& categ_name = this->categorize_event( event );
    event_category = this->get_category_code( categ_name, mc_signal );
  }

  bool selected = this->is_selected( event );

  if ( selected ) {
    ++num_passed_events_;
  }

  // Always store the results below in the output TTree
  auto out_tree = event.out();
  out_tree[ selection_name_ + "_MC_Signal" ] = mc_signal;
  out_tree[ selection_name_ + "_Selected" ] = selected;
  out_tree[ selection_name_ + "_Category" ] = event_category;

}

void SelectionBase::load_selection_config() {

  selection_config_ = std::make_unique< SelectionConfig >();

  // Get the path to the xsec_analyzer source code folder
  const char* path = std::getenv( "XSEC_ANALYZER_DIR" );
  if ( !path ) {
    throw std::runtime_error( "The environment variable XSEC_ANALYZER_DIR"
      " is not set. Please set it and try again." );
  }

  // Access the configs/categories.conf file
  std::string in_file_name = std::string( path ) + "/configs/selections.conf";
  std::ifstream config_file( in_file_name );

  // This regular expression matches all valid non-comment tokens in the
  // category configuration file, namely the "trees" keyword (without
  // double quotes), double-quote (") delimited strings, and hex color
  // codes (e.g., #33ccff)
  const std::regex token_rx(
    R"!((trees)|(signal)|"([^"]*)"|(#([0-9A-Fa-f]{6})))!" );

  // Stores the selection name(s) whose configuration options are currently
  // being processed. This set is initially empty and should remain so
  // until the default event category is set.
  std::set< std::string > sel_names;

  // Storage for the event category names and integer ROOT color codes
  // associated with the current selection name(s)
  std::shared_ptr< CategoryMap > categ_map;

  // Parse the configuration file line-by-line. Keep track of the line number
  // to allow for its inclusion in parsing error messages
  std::string line;
  size_t line_num = 0u;
  while ( std::getline( config_file, line ) ) {

    ++line_num;

    // Skip empty lines and comment lines (those for which the first
    // non-whitespace character is '#')
    size_t first_non_space = line.find_first_not_of( " \t" );
    if ( first_non_space == std::string::npos
      || line.at( first_non_space ) == '#' ) continue;

    // Flag that indicates whether TTree names should be parsed from the
    // current line. This is set to true below in response to the presence
    // of the "trees" keyword.
    bool tree_name_mode = false;

    // Flag that indicates whether the "signal" keyword occurred as the
    // first token on the current line. For categories, this indicates
    // that it should be considered part of the signal definition.
    bool signal_line = false;

    // Parse allowed tokens in the line. Keep track of the offset from the
    // start of the line to detect unmatched (and thus unallowed) text
    size_t offset = 0;
    auto begin = std::sregex_iterator( line.begin(), line.end(), token_rx );
    auto end = std::sregex_iterator();

    // Temporary storage for the allowed token types that need it
    std::vector< std::string > string_tokens;
    std::vector< std::string > hex_color_tokens;

    // Loop over the matched tokens, also keeping track of the current index
    size_t token_index = 0u;
    for ( auto it = begin; it != end; ++it, ++token_index )
    {
      const std::smatch& m = *it;

      // If there was any leading non-token and non-whitespace text before
      // the current regex match, complain
      if ( m.position() > offset ) {
        std::string junk = line.substr( offset, m.position() - offset );
        if ( junk.find_first_not_of(" \t") != std::string::npos ) {
          throw std::runtime_error( "Unexpected text \"" + junk
            + "\" encountered on line " + std::to_string( line_num )
            + " of the selection configuration file" );
        }
      }

      // Classify the current matched token based on the regex options.
      // The first option is the "trees" keyword, which is only allowed to
      // occur as the first token on a line.
      if ( m[ 1 ].matched ) {
        if ( token_index != 0 ) {
          throw std::runtime_error( "The \"trees\" keyword must be the first"
            " token on a line of the selection configuration file. Please"
            " check line " + std::to_string( line_num ) + " and try again." );
        }
        tree_name_mode = true;
      }
      // The second option is the "signal" keyword, which likewise may only
      // appear as the first token on a line
      else if ( m[ 2 ].matched ) {
        if ( token_index != 0 ) {
          throw std::runtime_error( "The \"signal\" keyword must be the first"
            " token on a line of the selection configuration file. Please"
            " check line " + std::to_string( line_num ) + " and try again." );
        }
        signal_line = true;
      }
      // The third option is a double-quoted string. Parse it and store
      // it for later.
      else if ( m[ 3 ].matched ) {
        string_tokens.push_back( m[ 3 ].str() );
      }
      // The third and final option is a hex color code (e.g., #cc33ff).
      // Store it as a string in a separate vector for later processing.
      else if ( m[ 4 ].matched ) {
        hex_color_tokens.push_back( m[ 4 ].str() );
      }

      offset = m.position() + m.length();
    }

    // Check for any unallowed trailing text after the last matched token
    if ( offset < line.size() ) {
      std::string junk = line.substr( offset );
      if ( junk.find_first_not_of(" \t") != std::string::npos ) {
        throw std::runtime_error( "Unexpected text \"" + junk
          + "\" encountered at the end of line " + std::to_string( line_num )
          + " of the selection configuration file" );
      }
    }

    // Handle the two valid line formats (neglecting the "trees" keyword,
    // which has already been processed). The two formats are
    //
    // 1. Zero or more double-quoted strings without any hex color code:
    //
    //   "CC1mu1p0pi" "CC1mu2p0pi" "CC1muNp0pi"
    //
    // 2. A single double-quoted category name followed by a single hex
    //    color code:
    //
    //    "CCQE" #ccee99

    size_t num_strings = string_tokens.size();
    size_t num_hex_colors = hex_color_tokens.size();

    if ( num_hex_colors == 0u ) {

      // The "signal" keyword is only valid for category definition lines.
      // If it occurred on a line with zero hex colors, the input is
      // incorrect, and we need to complain about it.
      if ( signal_line ) throw std::runtime_error( "The \"signal\" keyword"
        " was used on line " + std::to_string( line_num ) + " of the"
        " selections configuration file, but this line does not contain"
        " a valid category definition" );

      if ( tree_name_mode ) {
        // If the current line represents TTree names, populate a new set of
        // names by inserting all the string tokens from the line
        auto tree_names = std::make_shared< TreeNameSet >(
          string_tokens.begin(), string_tokens.end() );

        if ( sel_names.empty() ) {
          throw std::runtime_error( "\"trees\" keyword used for an undefined"
            " selection in the selections configuration file" );
        }

        // Associate these TTree names with the names of all active selections
        for ( const auto& sn : sel_names ) {
          auto& config_pair = selection_config_->at( sn );
          config_pair.first = tree_names;
        }
      }
      else {
        // The current line gives a list of selection names, so set these
        // to be the active ones. Also check first whether the list is
        // empty and complain if it is.
        if ( string_tokens.empty() ) {
          throw std::runtime_error( "Missing selection names on line "
            + std::to_string( line_num ) + " of the selections configuration"
            + " file" );
        }

        // Note that this replaces the existing set with new contents
        sel_names = std::set< std::string >( string_tokens.begin(),
          string_tokens.end() );

	// Create new entries for these selections in the configuration map.
	// First check if the selection already exists and complain about
	// duplicates.
        for ( const auto& sn : sel_names ) {
          if ( selection_config_->count( sn ) ) {
            throw std::runtime_error( "Duplicate selection name \"" + sn
              + "\" encountered on line " + std::to_string( line_num )
              + " of the selections configuration file" );
          }

          // Accessing the map using the [] operator will create
          // a new value by calling the default constructor.
          auto& sel_pair = selection_config_->operator[]( sn );

          // Insert the default category with integer code 0 at the start
          // of the category map for the current selection
          auto& categ_map = sel_pair.second;
          categ_map = std::make_shared< CategoryMap >();
          categ_map->emplace( UNKNOWN_CATEGORY_CODE,
            CategoryDefinition( UNKNOWN_CATEGORY_NAME,
              UNKNOWN_CATEGORY_COLOR, false )
          );
        }

      }

    }
    else if ( num_hex_colors == 1u && num_strings == 1u ) {

      // The "trees" keyword is only valid for lines giving input TTree
      // names. If we have encountered it on a line with a hex color
      // (i.e., a category definition line), then complain about it.
      if ( tree_name_mode ) throw std::runtime_error( "The \"trees\" keyword"
        " was used on line " + std::to_string( line_num ) + " of the"
        " selections configuration file, but this line does not contain"
        " a valid set of input TTree names" );

      // We are defining a new category in terms of a name and hex color code.
      // If we do not have any active selections defined, then this is an error
      // condition.
      if ( sel_names.empty() ) {
        throw std::runtime_error( "Category definition provided for an"
          " undefined selection on line " + std::to_string( line_num )
          + " of the selections configuration file" );
      }

      // Add this category definition to the maps for all active selections
      for ( const auto& sn : sel_names ) {
	// If this selection name hasn't been previously added to the
	// configuration map, complain
        if ( !selection_config_->count( sn ) ) {
          throw std::runtime_error( "Undefined selection name \"" + sn
            + "\" encountered when processing line "
            + std::to_string( line_num ) + " of the selections configuration"
            + " file" );
        }

        // Make sure that the category map for the current selection has
        // already been initialized
        auto& categ_map = selection_config_->at( sn ).second;
        if ( !categ_map ) {
          throw std::runtime_error( "Uninitialized category map encountered"
            " for selection \"" + sn + "\" while parsing the selections"
            " configuration file" );
        }

	// Retrieve the ROOT integer color code for the requested hex color,
	// defining a new TColor object as needed in the global system
        const std::string& categ_name = string_tokens.front();
        const auto& hex_color_str = hex_color_tokens.front();
        int color_code = TColor::GetColor( hex_color_str.c_str() );

        // Check for duplication of a category already present in the map.
        // If we find one, complain.
        auto it = std::find_if( categ_map->cbegin(), categ_map->cend(),
          [&]( const auto& c_pair ) -> bool {
            return c_pair.second.name_ == categ_name;
          }
        );

        if ( it != categ_map->cend() ) {
          throw std::runtime_error( "Category name \"" + categ_name
            + "\" found on line " + std::to_string( line_num ) + " of"
            + " the selections configuration file is already defined for"
            + " the \"" + sn + "\" selection" );
        }

        // Insert a new category entry into the map using the current
        // number of categories pre-insertion as the integer code for the
        // new category
        size_t num_existing_categories = categ_map->size();
        categ_map->emplace( num_existing_categories,
          CategoryDefinition( categ_name, color_code, signal_line )
        );
      }

    }
    else {
      throw std::runtime_error( "Invalid format on line "
        + std::to_string( line_num ) + " of the selections"
        + " configuration file" );
    }

  } // loop over selection configuration file lines

}

const FiducialVolume& SelectionBase::get_fv() const {
  if ( !fv_ ) {
    throw std::runtime_error( "Fiducial volume has not been defined"
      " for selection " + selection_name_ + '\n' );
  }
  return *fv_;
}

void SelectionBase::define_FV( double x_min, double x_max, double y_min,
  double y_max, double z_min, double z_max )
{
  fv_ = std::make_unique< FiducialVolume >( x_min, x_max, y_min, y_max,
    z_min, z_max );
}

int SelectionBase::get_category_code( const std::string& categ_name,
  bool& is_signal_category ) const
{
  // Grab the category map and search for the integer key corresponding
  // to the category name
  const auto& categ_map = this->category_map();
  auto it = std::find_if( categ_map.cbegin(), categ_map.cend(),
    [&]( const auto& c_pair ) -> bool {
      return c_pair.second.name_ == categ_name;
    }
  );

  // Obtain the integer key from the lookup, or complain if it failed
  int result = DUMMY_CATEGORY_CODE;
  if ( it == categ_map.cend() ) throw std::runtime_error( "Undefined category"
    " \"" + categ_name + "\" encountered for the selection \""
    + selection_name_ + '\"' );
  else result = it->first;

  // Also retrieve the flag indicating whether this category corresponds
  // to signal (true) or background (false)
  is_signal_category = it->second.signal_;

  std::cout << "Looked up category \"" << categ_name << "\": " << result
    << " is signal = " << is_signal_category << '\n';

  return result;
}
