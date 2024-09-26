// Standard library includes
#include <fstream>
#include <sstream>
#include <vector>

// ROOT includes
#include "TFile.h"
#include "TMatrixD.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/StandaloneUnfolding.hh"
#include "XSecAnalyzer/DAgostiniUnfolder.hh"
#include "XSecAnalyzer/WienerSVDUnfolder.hh"

using DCC = DAgostiniUnfolder::ConvergenceCriterion;
using RMT = WienerSVDUnfolder::RegularizationMatrixType;

StandaloneUnfolding::StandaloneUnfolding( const std::string& config_file_name )
{
  std::ifstream config_file( config_file_name );
  if ( !config_file.good() ) {
    throw std::runtime_error( "Could not read StandaloneUnfolding"
      " configuration from the file \"" + config_file_name + '\"' );
  }

  // Strings that will be initialized using the contents of the configuration
  // file
  std::string input_root_file_name;
  std::string output_root_file_name;
  std::string unfolding_tech;
  std::string unfolding_opt;

  std::string line;
  while ( std::getline(config_file, line) ) {
    // Skip lines starting with a '#' character
    if ( line.front() == '#' ) continue;

    // Prepare a stringstream to parse the line contents
    std::istringstream iss( line );
    std::string first_word;
    iss >> first_word;

    if ( first_word == "InputFile" ) {
      iss >> input_root_file_name;
      input_root_file_.reset(
        new TFile( input_root_file_name.c_str(), "read" )
      );
    }
    else if ( first_word == "OutputFile" ) {
      // Read in the non-default setting for the name of the configuration
      // file for the FilePropertiesManager
      iss >> output_root_file_name;
      output_root_file_.reset(
        new TFile( output_root_file_name.c_str(), "recreate" )
      );
    }
    else if ( first_word == "BlocksFile" ) {
      iss >> blocks_file_;
    }
    else if ( first_word == "Unfold" ) {
      // Get the string indicating which unfolding method should be used
      std::string unf_type;
      iss >> unf_type;

      // Construct the appropriate Unfolder object
      Unfolder* temp_unfolder = nullptr;
      if ( unf_type == "DAgostini" ) {
	unfolding_tech = "DAgostini";

        // Determine how to configure the D'Agostini unfolding algorithm
        std::string dagost_mode;
        iss >> dagost_mode;

        if ( dagost_mode == "fm" ) {
          double fig_of_merit;
          iss >> fig_of_merit;
          temp_unfolder = new DAgostiniUnfolder( DAgostiniUnfolder
            ::ConvergenceCriterion::FigureOfMerit, fig_of_merit );

	  unfolding_opt = Form("Figure of Merit : %4.2f", fig_of_merit);
        }
        else if ( dagost_mode == "iter" ) {
          int iterations;
          iss >> iterations;
          temp_unfolder = new DAgostiniUnfolder( iterations );

	  unfolding_opt = Form("Iterations : %i",iterations);
        }
        else {
          throw std::runtime_error( "Unrecognized D'Agostini unfolding"
            " mode \"" + dagost_mode + '\"' );
        }
      }
      else if ( unf_type == "WienerSVD" ) {
        unfolding_tech = "WienerSVD";

        bool use_filter;
        iss >> use_filter;

	unfolding_opt = Form("use_filter : %i",use_filter);

        // Default to regularizing using the second derivative. The choice
        // actually doesn't matter unless the Wiener filter is used.
        RMT reg_type = RMT::kSecondDeriv;

        if ( use_filter ) {
          // Determine the regularization recipe to use with the WSVD
          // method
          std::string reg_mode;
          iss >> reg_mode;

          if ( reg_mode == "identity" ) {
            reg_type = RMT::kIdentity;
	    unfolding_opt = "Identity";
          }
          else if ( reg_mode == "first-deriv" ) {
            reg_type = RMT::kFirstDeriv;
	    unfolding_opt = "First Derivative";
          }
          else if ( reg_mode == "second-deriv" ) {
            reg_type = RMT::kSecondDeriv;
	    unfolding_opt = "Second Derivative";
          }
          else {
            throw std::runtime_error( "Unrecognized Wiener-SVD"
              " regularization mode \"" + reg_mode + '\"' );
          }
        }

        temp_unfolder = new WienerSVDUnfolder( use_filter, reg_type );
      }
      else {
        throw std::runtime_error( "Unrecognized unfolder type \""
          + unf_type + '\"' );
      }

      // Store the complete Unfolder object for later use
      unfolder_.reset( temp_unfolder );
    }
    else {
      throw std::runtime_error( "Unrecognized StandaloneUnfolding "
        " configuration file command \"" + first_word + '\"' );
    }
  }

  std::cout << "StandaloneUnfolding initialized with options:" << '\n';
  std::cout << "\tinput ROOT file: " << input_root_file_name << '\n';
  std::cout << "\toutput ROOT file: " << output_root_file_name << '\n';
  std::cout << "\tunfolding_tech: " << unfolding_tech << '\n';
  std::cout << "\tblocks file: " << blocks_file_ << '\n';
  std::cout << "\t\tOption: " << unfolding_opt << "\n\n";

  // We've finished parsing the configuration file. Check that we have the
  // required information.
  if ( !unfolder_ ) {
    throw std::runtime_error( "Missing \"Unfold\" command in the"
      " StandaloneUnfolding configuration file" );
  }
  if ( !input_root_file_ ) {
    throw std::runtime_error( "Missing \"InputFile\" command in the"
      " StandaloneUnfolding configuration file" );
  }
  if ( !output_root_file_ ) {
    throw std::runtime_error( "Missing \"OutputFile\" command in the"
      " StandaloneUnfolding configuration file" );
  }
  if ( blocks_file_.empty() ) {
    throw std::runtime_error( "Missing \"BlocksFile\" command in the"
      " StandaloneUnfolding configuration file" );
  }

}

void StandaloneUnfolding::run_unfolding() const {

  // Retrieve the pre-calculated unfolding inputs
  TMatrixD* data_signal = nullptr;
  TMatrixD* data_covmat = nullptr;
  TMatrixD* smearcept = nullptr;
  TMatrixD* prior_true_signal = nullptr;

  input_root_file_->GetObject( "data_signal", data_signal );
  input_root_file_->GetObject( "data_covmat", data_covmat );
  input_root_file_->GetObject( "smearcept", smearcept );
  input_root_file_->GetObject( "prior_true_signal", prior_true_signal );

  // Perform the unfolding
  UnfoldedMeasurement result = unfolder_->unfold( *data_signal, *data_covmat,
    *smearcept, *prior_true_signal, blocks_file_ );

  TMatrixD* unf_sig = result.unfolded_signal_.get();
  output_root_file_->WriteObject( unf_sig, "unfolded_signal" );

  TMatrixD* cov = result.cov_matrix_.get();
  output_root_file_->WriteObject( cov, "cov_matrix" );

  TMatrixD* unf_mat = result.unfolding_matrix_.get();
  output_root_file_->WriteObject( unf_mat, "unfolding_matrix" );

  TMatrixD* Ac = result.unfolding_matrix_.get();
  output_root_file_->WriteObject( Ac, "add_smear_matrix" );

}
