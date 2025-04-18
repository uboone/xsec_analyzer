// Post-processing program for the MicroBooNE xsec_analyzer framework. This is
// currently designed for use with the PeLEE group's "searchingfornues" ntuples
//
// Updated 24 September 2024
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
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/AnalysisEvent.hh"
#include "XSecAnalyzer/Branches.hh"
#include "XSecAnalyzer/Constants.hh"
#include "XSecAnalyzer/Functions.hh"

#include "XSecAnalyzer/Selections/SelectionBase.hh"
#include "XSecAnalyzer/Selections/SelectionFactory.hh"

struct Arguments {
  std::vector<std::string> input_files;
  std::string output_file;
  std::vector<std::string> selection_names;
  int nevents = -1;
};

void analyze( const Arguments & arguments )
{
  std::cout << "\nRunning ProcessNTuples with options:\n";
  std::cout << "\toutput_filename: " << arguments.output_file << '\n';
  std::cout << "\tinput_file_names:\n";
  for (const auto & name : arguments.input_files) {
    std::cout << "\t\t- " << name << '\n';
  }
  std::cout << "\n\nselection names:\n";
  for ( const auto& sel_name : arguments.selection_names ) {
    std::cout << "\t\t- " << sel_name << '\n';
  }

  // Get the TTrees containing the event ntuples and subrun POT information
  // Use TChain objects for simplicity in manipulating multiple files
  TChain events_ch( "nuselection/NeutrinoSelectionFilter" );
  TChain subruns_ch( "nuselection/SubRun" );

  for ( const auto& f_name : arguments.input_files ) {
    events_ch.Add( f_name.c_str() );
    subruns_ch.Add( f_name.c_str() );
  }

  // OUTPUT TTREE
  // Make an output TTree for plotting (one entry per event)
  TFile* out_file = new TFile( arguments.output_file.c_str(), "recreate" );
  out_file->cd();
  TTree* out_tree = new TTree( "stv_tree", "STV analysis tree" );

  // Get the total POT from the subruns TTree. Save it in the output
  // TFile as a TParameter<float>. Real data doesn't have this TTree,
  // so check that it exists first.
  float pot;
  float summed_pot = 0.;
  bool has_pot_branch = ( subruns_ch.GetBranch("pot") != nullptr );
  if ( has_pot_branch ) {
    subruns_ch.SetBranchAddress( "pot", &pot );
    for ( int se = 0; se < subruns_ch.GetEntries(); ++se ) {
      subruns_ch.GetEntry( se );
      summed_pot += pot;
    }
  }

  TParameter<float>* summed_pot_param = new TParameter<float>( "summed_pot",
    summed_pot );

  summed_pot_param->Write();

  std::vector< std::unique_ptr<SelectionBase> > selections;

  SelectionFactory sf;
  for ( const auto& sel_name : arguments.selection_names ) {
    selections.emplace_back().reset( sf.CreateSelection(sel_name) );
  }

  out_file->cd();
  for ( auto& sel : selections ) {
    sel->setup( out_tree );
  }

  // EVENT LOOP
  // TChains can potentially be really big (and spread out over multiple
  // files). When that's the case, calling TChain::GetEntries() can be very
  // slow. I get around this by using a while loop instead of a for loop.
  bool created_output_branches = false;
  long events_entry = 0;

  while ( true ) {

    //If not doing all events (-1), break after nevents
    if ( (arguments.nevents != -1) && (events_entry > arguments.nevents) )
      break;

    if ( events_entry % 1000 == 0 ) {
      std::cout << "Processing event #" << events_entry << '\n';
    }

    // Create a new AnalysisEvent object. This will reset all analysis
    // variables for the current event.
    AnalysisEvent cur_event;

    // Set branch addresses for the member variables that will be read
    // directly from the Event TTree.
    set_event_branch_addresses( events_ch, cur_event );

    // TChain::LoadTree() returns the entry number that should be used with
    // the current TTree object, which (together with the TBranch objects
    // that it owns) doesn't know about the other TTrees in the TChain.
    // If the return value is negative, there was an I/O error, or we've
    // attempted to read past the end of the TChain.
    int local_entry = events_ch.LoadTree( events_entry );

    // If we've reached the end of the TChain (or encountered an I/O error),
    // then terminate the event loop
    if ( local_entry < 0 ) break;

    // Load all of the branches for which we've called
    // TChain::SetBranchAddress() above
    events_ch.GetEntry( events_entry );

    // Set the output TTree branch addresses, creating the branches if needed
    // (during the first event loop iteration)
    bool create_them = false;
    if ( !created_output_branches ) {
      create_them = true;
      created_output_branches = true;
    }
    set_event_output_branch_addresses(*out_tree, cur_event, create_them );

    for ( auto& sel : selections ) {
      sel->apply_selection( &cur_event );
    }

    // We're done. Save the results and move on to the next event.
    out_tree->Fill();
    ++events_entry;
  }

  for ( auto& sel : selections ) {
    sel->summary();
  }
  std::cout << "Wrote output to:" << arguments.output_file << std::endl;

  for ( auto& sel : selections ) {
    sel->final_tasks();
  }

  out_tree->Write();
  out_file->Close();
  delete out_file;
}

void analyzer( const Arguments & arguments )
{
  analyze( arguments );
}

bool parse_args( int argc, char* argv[], Arguments & arg_results ) {
  for (int iArg = 1; iArg < argc; iArg++) {
    if (!strcasecmp(argv[iArg],"-o")) {
      arg_results.output_file = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-i")) {
      arg_results.input_files.push_back(argv[++iArg]);
    }
    if (!strcasecmp(argv[iArg],"-s")) {
      std::stringstream sel_ss(argv[++iArg]);
      std::string sel_name;
      while ( std::getline(sel_ss, sel_name, ',') ) {
        arg_results.selection_names.push_back( sel_name );
      }
    }
    if (!strcasecmp(argv[iArg],"-n")) {
      arg_results.nevents = std::atoi(argv[++iArg]);
      if (arg_results.nevents < -1) {
        std::cerr << "Error: Must provide -1 for all events or positive number" <<
                     std::endl;
        return false;
      }
    }
    if (!strcasecmp(argv[iArg],"-h")) {
      std::cout << argv[0] <<
          "-i <input_pelee_file> -o <output_file> " <<
          "-s <comma-separated selection names list> " <<
          "-n <nevents: default -1 for all> " <<
          std::endl;
      return false;
    }
  }
  return true;
}

int main( int argc, char* argv[] ) {

  Arguments arguments;
  if (!parse_args(argc, argv, arguments)) {
    return 1;
  }

  analyzer( arguments );

  return 0;
}
