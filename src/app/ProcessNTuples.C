// Post-processing program for the MicroBooNE xsec_analyzer framework
//
// Updated 26 April 2025
// Steven Gardiner <gardiner@fnal.gov>
// Daniel Barrow <daniel.barrow@physics.ox.ac.uk>

// Standard library includes
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// ROOT includes
#include "TFile.h"
#include "TBranch.h"
#include "TParameter.h"
#include "TSystem.h"
#include "TTree.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Constants.hh"
#include "XSecAnalyzer/Functions.hh"

#include "XSecAnalyzer/Selections/SelectionBase.hh"
#include "XSecAnalyzer/Selections/SelectionFactory.hh"
#include "XSecAnalyzer/TreeHandler.hh"

int main( int argc, char* argv[] ) {

  if ( argc != 4 ) {
    std::cout << "Usage: " << argv[0]
      << " INPUT_NTUPLE_FILE SELECTION_NAMES OUTPUT_FILE\n";
    return 1;
  }

  std::string input_file_name( argv[ 1 ] );
  std::string output_file_name( argv[ 3 ] );

  std::vector< std::string > selection_names;

  std::stringstream sel_ss( argv[ 2 ] );
  std::string sel_name;
  while ( std::getline( sel_ss, sel_name, ',' ) ) {
    selection_names.push_back( sel_name );
  }

  std::cout << "\nRunning ProcessNTuples with options:\n";
  std::cout << "\tinput_file_name: " << input_file_name << '\n';
  std::cout << "\toutput_file_name: " << output_file_name << '\n';
  std::cout << "\n\nselection names:\n";
  for ( const auto& sel_name : selection_names ) {
    std::cout << "\t\t- " << sel_name << '\n';
  }

  // Copy the input ntuple file to the output ntuple file
  int copy_result = gSystem->CopyFile( input_file_name.c_str(),
    output_file_name.c_str(), true );

  // Complain if we encountered a problem during the copying process
  if ( copy_result != 0 ) {
    throw std::runtime_error( "Unsuccessful attempt to duplicate"
      " the input file." );
  }

  // Open the output file and create a new subfolder to hold the
  // post-processing results
  TFile out_file( output_file_name.c_str(), "update" );
  TDirectory* sub_dir = out_file.mkdir( "XSecAnalyzer" );
  if ( !sub_dir ) {
    throw std::runtime_error( "Failed to create XSecAnalyzer subdirectory" );
  }
  sub_dir->cd();

  // Get the PeLEE TTree containing subrun POT information. Note that real data
  // does not have this TTree, so just skip summing the POT if it is missing.
  TTree* subruns_tree = nullptr;
  out_file.GetObject( "nuselection/SubRun", subruns_tree );

  // Get the total POT from the subruns TTree
  float summed_pot = 0.;
  if ( subruns_tree ) {
    TBranch* pot_br = subruns_tree->GetBranch( "pot" );
    if ( pot_br ) {
      float pot;
      pot_br->SetAddress( &pot );
      for ( int se = 0; se < pot_br->GetEntries(); ++se ) {
        pot_br->GetEntry( se );
        summed_pot += pot;
      }
    }
  }

  // Save the total POT in the output subdirectory as a TParameter< float >
  auto summed_pot_param = std::make_unique< TParameter< float > >(
    "summed_pot", summed_pot );

  sub_dir->cd();
  summed_pot_param->Write();

  // Determine whether we are working with data or MC using the presence
  // or absence of the SubRun TTree
  bool is_mc = ( subruns_tree != nullptr );

  // Also save the result as a TParameter< bool >
  auto is_mc_param = std::make_unique< TParameter< bool > >( "is_mc", is_mc );
  is_mc_param->Write();

  // Instantiate the requested selections to use for post-processing
  std::vector< std::unique_ptr< SelectionBase > > selections;

  SelectionFactory sf;
  for ( const auto& sel_name : selection_names ) {
    selections.emplace_back().reset( sf.CreateSelection( sel_name ) );
  }

  // Create a TreeHandler object to manage TTree input/output. Apply
  // setup procedures to it from each of the configured selections
  TreeHandler th;

  for ( auto& sel : selections ) {

    // Configure an output TTree to store the selection results, and give
    // it the same name as the selection object itself
    const std::string& sel_name = sel->name();
    th.add_output_tree( sub_dir, sel_name, sel_name + " selection results" );

    // Load all input TTrees required by the current selection object
    const auto& in_tree_names = sel->input_tree_names();
    for ( const auto& tn : in_tree_names ) {
      TTree* temp_tree = nullptr;
      out_file.GetObject( tn.c_str(), temp_tree );

      if ( !temp_tree ) throw std::runtime_error( "Could not load TTree \""
        + tn + "\" from the file" + output_file_name + " (this TTree"
        + " is required to apply the \"" + sel_name + "\" selection)" );

      th.add_input_tree( temp_tree, tn );
    }
  }

  // EVENT LOOP
  long long events_entry = 0;
  while ( true ) {

    if ( events_entry % 1000 == 0 ) {
      std::cout << "Processing event #" << events_entry << '\n';
    }

    // If the return value is negative, there was an I/O error, or we've
    // attempted to read past the end of a TTree
    int local_entry = th.load_tree( events_entry );

    // If we've reached the end of the TTrees (or encountered an I/O error),
    // then terminate the event loop
    if ( local_entry < 0 ) break;

    // Load the current event using branches from the input TTree(s)
    th.get_entry( events_entry );

    for ( auto& sel : selections ) {
      sel->apply_selection( is_mc, th );
    }

    // We're done. Save the results and move on to the next event.
    th.fill();
    ++events_entry;
  }

  for ( auto& sel : selections ) {
    std::cout << sel->name() << " has " << sel->passed_events()
      << " events which passed\n";
  }
  std::cout << "Wrote output to:" << output_file_name << '\n';

  th.write();
}
