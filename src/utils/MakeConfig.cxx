// Standard library includes
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

// XSecAnalyzer includes
#include "XSecAnalyzer/ConfigMakerUtils.hh"
#include "XSecAnalyzer/HistUtils.hh"
#include "XSecAnalyzer/SliceBinning.hh"
#include "XSecAnalyzer/UniverseMaker.hh"

#include "XSecAnalyzer/Binning/MakeConfig.hh"

// Lazy way to get compiled code for all the functions defined in headers
#include "XSecAnalyzer/AnalysisEvent.hh"
#include "XSecAnalyzer/Branches.hh"
#include "XSecAnalyzer/Constants.hh"
#include "XSecAnalyzer/ConstrainedCalculator.hh"
#include "XSecAnalyzer/CovMatUtils.hh"
#include "XSecAnalyzer/CrossSectionExtractor.hh"
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/FilePropertiesManager.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/IntegratedFluxUniverseManager.hh"
#include "XSecAnalyzer/MatrixUtils.hh"
#include "XSecAnalyzer/MCC9SystematicsCalculator.hh"
#include "XSecAnalyzer/NormShapeCovMatrix.hh"
#include "XSecAnalyzer/PGFPlotsDumpUtils.hh"
#include "XSecAnalyzer/PlotUtils.hh"
#include "XSecAnalyzer/SliceHistogram.hh"
#include "XSecAnalyzer/SystematicsCalculator.hh"
#include "XSecAnalyzer/TreeUtils.hh"
#include "XSecAnalyzer/TruthSystematicsCalculator.hh"
#include "XSecAnalyzer/Unfolder.hh"
#include "XSecAnalyzer/WeightHandler.hh"

// Initialize the static dummy counter owned by the MakeConfig class. This
// counter is used to ensure that each automatically-generated histogram has a
// unique name to use with TTree::Draw()
int MakeConfig::hist_count = 0;

MakeConfig::MakeConfig( const std::string& bin_scheme_name ) {
  RUNS = {1};
  bin_scheme_name_ = bin_scheme_name;
}

MakeConfig::~MakeConfig(){
  std::cout << "free make config\n";
  for(int i = 0; i < vect_block->size(); i++){
    delete vect_block->at(i).block_true_;
    delete vect_block->at(i).block_reco_;
  }
}

void MakeConfig::ResPlots() {

  for(int i = 0; i < vect_block->size(); i++){
    ++hist_count;
    if(vect_block->at(i).block_true_->Is1D()){
      make_res_plots( vect_block->at(i).block_reco_->GetXName(),
          vect_block->at(i).block_reco_->GetXTitle(),
          vect_block->at(i).block_reco_->GetSelection(),
          RUNS,
          vect_block->at(i).block_reco_->GetVector(),
          false, false,
          DEFAULT_TRUE_BINS,
          vect_block->at(i).block_true_->GetXName(),
          vect_block->at(i).block_true_->GetSelection(),
          DEFAULT_MC_EVENT_WEIGHT );
    }
    else{

      std::vector< TrueBin > true_bins;
      std::vector< RecoBin > reco_bins;
      for(int j = 0; j < vect_block->at(i).block_true_->GetNBinsX(); j++){
        double xlow = vect_block->at(i).block_true_->GetBinXLow(j);
        double xhigh = vect_block->at(i).block_true_->GetBinXHigh(j);
        for(int k = 0; k < vect_block->at(i).block_true_->GetNBinsY(j); k++){
          true_bins.emplace_back(
            vect_block->at(i).block_true_->GetBinDef(j, k),
            TrueBinType(vect_block->at(i).block_true_->GetBinType()), i );
          reco_bins.emplace_back(
            vect_block->at(i).block_reco_->GetBinDef(j, k),
            RecoBinType(vect_block->at(i).block_reco_->GetBinType()), i );
        }
      }

      std::stringstream temp_ss;
      temp_ss << "temp\n";
      temp_ss << "stv_tree\n";
      temp_ss << true_bins.size() << '\n';
      for ( const auto& tb : true_bins ) temp_ss << tb << '\n';

      temp_ss << reco_bins.size() << '\n';
      for ( const auto& rb : reco_bins ) temp_ss << rb << '\n';
      true_bins.clear();
      reco_bins.clear();

      // Move the input position to the start of the stream
      temp_ss.seekg( 0 );

      make_res_plots( temp_ss, RUNS );
    }
  }
}

void MakeConfig::Print(){
  constexpr int DUMMY_BLOCK_INDEX = -1;

  SliceBinning sb;
  std::vector< TrueBin > true_bins;
  std::vector< RecoBin > reco_bins;


  for(int i = 0; i < vect_block->size(); i++){

    if(vect_block->at(i).block_true_->Is1D()){

      sb.slice_vars_.emplace_back( vect_block->at(i).block_true_->GetXTitle(),
          vect_block->at(i).block_true_->GetXTexTitle(),
          vect_block->at(i).block_true_->GetXTitleUnit(),
          vect_block->at(i).block_true_->GetXTexTitleUnit() );
      int var_idx = find_slice_var_index(
        vect_block->at(i).block_true_->GetXTitle(), sb.slice_vars_ );

      size_t num_bins = vect_block->at(i).block_true_->GetNBinsX();

      auto& slice = add_slice(sb, vect_block->at(i).block_true_->GetVector(), var_idx);

      for(int j = 0; j < vect_block->at(i).block_true_->GetNBinsX(); j++){
        slice.bin_map_[ j + 1 ].insert( true_bins.size() );
        true_bins.emplace_back( vect_block->at(i).block_true_->GetBinDef(j),
          TrueBinType(vect_block->at(i).block_true_->GetBinType()), i );

      }
      for(int j = 0; j < vect_block->at(i).block_reco_->GetNBinsX(); j++){
        reco_bins.emplace_back( vect_block->at(i).block_reco_->GetBinDef(j),
          RecoBinType(vect_block->at(i).block_reco_->GetBinType()), i );
      }
    }
    else{
      sb.slice_vars_.emplace_back( vect_block->at(i).block_true_->GetXTitle(),
          vect_block->at(i).block_true_->GetXTexTitle(),
          vect_block->at(i).block_true_->GetXTitleUnit(),
          vect_block->at(i).block_true_->GetXTexTitleUnit() );
      sb.slice_vars_.emplace_back(
        vect_block->at(i).block_true_->GetYTitle(),
        vect_block->at(i).block_true_->GetYTexTitle(),
        vect_block->at(i).block_true_->GetYTitleUnit(),
        vect_block->at(i).block_true_->GetYTexTitleUnit() );


      int xvar_idx = find_slice_var_index(
        vect_block->at(i).block_true_->GetXTitle(), sb.slice_vars_ );
      int yvar_idx = find_slice_var_index(
        vect_block->at(i).block_true_->GetYTitle(), sb.slice_vars_ );

      for( int j = 0; j < vect_block->at(i).block_true_->GetNBinsX(); j++ ){
        double xlow = vect_block->at(i).block_true_->GetBinXLow(j);
        double xhigh = vect_block->at(i).block_true_->GetBinXHigh(j);
        auto& slice = add_slice( sb, vect_block->at(i).block_true_->GetVector(j),
            xvar_idx, yvar_idx, xlow, xhigh );
        for( int k = 0; k < vect_block->at(i).block_true_->GetNBinsY(j); k++ ){
          slice.bin_map_[ k + 1 ].insert( true_bins.size() );
          true_bins.emplace_back(vect_block->at(i).block_true_->GetBinDef(j, k),
            TrueBinType(vect_block->at(i).block_true_->GetBinType()), i );
          reco_bins.emplace_back(vect_block->at(i).block_reco_->GetBinDef(j, k),
            RecoBinType(vect_block->at(i).block_reco_->GetBinType()), i );
        }
      }
    }
  }

  sb.slice_vars_.emplace_back( "bin number",
      "",
      "bin number",
      "" );
  int bin_number_var_idx = find_slice_var_index( "bin number", sb.slice_vars_ );

  // Create a slice showing all results together as a function of bin number
  auto& bin_num_slice = add_slice( sb, reco_bins.size(), 0, reco_bins.size(),
      bin_number_var_idx );
  for ( int ab = 0; ab < reco_bins.size(); ++ab ) {
    // The ROOT histogram bins are one-based, so we correct for this here
    bin_num_slice.bin_map_[ ab + 1 ].insert( ab );
  }

  // Add a single true bin to collect background events by inverting the
  // signal definition for the input selection
  std::string bkgd_bdef = "!" + SELECTION + "_MC_Signal";
  true_bins.emplace_back( bkgd_bdef, kBackgroundTrueBin, DUMMY_BLOCK_INDEX );

  std::cout << DIRECTORY << '\n';
  std::cout << TREE << '\n';
  std::cout << SELECTION << '\n';
  std::cout << true_bins.size() << '\n';
  for ( const auto& tb : true_bins ) std::cout  << tb << '\n';
  std::cout << reco_bins.size() << '\n';
  for ( const auto& tb : reco_bins ) std::cout  << tb << '\n';

  std::string bin_config_output = std::getenv( "XSEC_ANALYZER_DIR" )
    + std::string( "/configs/" ) + BIN_CONFIG + "bin_config.txt";

  std::ofstream out_file( bin_config_output );
  out_file <<  DIRECTORY << '\n';
  out_file << TREE << '\n';
  out_file << SELECTION << '\n';
  out_file << true_bins.size() << '\n';
  for ( const auto& tb : true_bins ) out_file << tb << '\n';

  out_file << reco_bins.size() << '\n';
  for ( const auto& rb : reco_bins ) out_file << rb << '\n';
  out_file.close();

  std::string slice_config_output = std::getenv( "XSEC_ANALYZER_DIR" )
    + std::string( "/configs/" ) + BIN_CONFIG + "slice_config.txt";

  std::cout << sb << '\n';
  std::ofstream slice_out_file( slice_config_output );
  slice_out_file << sb;
  slice_out_file.close();

  std::cout << "Save universes bin configuration into => "
    << bin_config_output << '\n';
  std::cout << "Save slice configuration into         => "
   << slice_config_output << '\n';
}

void MakeConfig::make_res_plots( const std::string& branchexpr,
    const std::string& variable_title, const std::string& selection,
    const std::set<int>& runs, std::vector<double> bin_low_edges,
    bool show_bin_plots,
    bool show_smear_numbers,
    int num_true_bins,
    const std::string& mc_branchexpr,
    const std::string& signal_cuts,
    const std::string& mc_event_weight )
{
  // Get the outer edges of the reco-space bins. This will be used to set the
  // plot range for the true-space histograms.
  double xmin = bin_low_edges.front();
  double xmax = bin_low_edges.back();

  // If the user hasn't explicitly specified a branch expression for the
  // true quantity, assume that it's the same as the reco quantity but
  // with the prefix "mc_" added.
  std::string true_branchexpr = mc_branchexpr;
  //if ( true_branchexpr.empty() ) {
  //  true_branchexpr = "mc_" + branchexpr;
  //}

  // Get access to the singleton utility class that manages the processed
  // ntuple files
  const FilePropertiesManager& fpm = FilePropertiesManager::Instance();

  // Make a TChain to process the CV numu ntuples from the requested run(s).
  // For the resolution studies, this is all we really need. Add the
  // appropriate ntuples to the TChain. Also tally the total simulated
  // POT for later scaling purposes.
  TChain chain( "stv_tree" );
  double total_simulated_POT = 0.;

  const auto& ntuple_map = fpm.ntuple_file_map();
  for ( const auto& run : runs ) {
    const auto& ntuple_files = ntuple_map.at( run )
      .at( NtupleFileType::kNumuMC );
    for ( const auto& file_name : ntuple_files ) {
      chain.Add( file_name.c_str() );

      TFile temp_file( file_name.c_str(), "read" );
      TParameter<float>* temp_pot = nullptr;
      temp_file.GetObject( "summed_pot", temp_pot );
      double pot = temp_pot->GetVal();
      total_simulated_POT += pot;
    }
  }

  // Dummy counter used to ensure that each histogram generated by this
  // function has a unique name to use with TTree::Draw()
  // static int hist_count = 0; put it in header file

  if ( show_bin_plots ) {
    for ( size_t b = 1u; b < bin_low_edges.size(); ++b ) {
      ++hist_count;

      TCanvas* c = new TCanvas;
      std::string true_hist_name = "true_hist" + std::to_string( hist_count );

      TH1D* true_hist = new TH1D( true_hist_name.c_str(),
          ("true events in " + variable_title + " reco bin "
           + std::to_string(b) + "; " + variable_title + "; events").c_str(),
          num_true_bins, xmin, xmax );

      double reco_bin_min = bin_low_edges.at( b - 1 );
      double reco_bin_max = bin_low_edges.at( b );
      std::string cuts = mc_event_weight + " * (is_mc && " + signal_cuts
        + " && " + selection + " && " + branchexpr + " >= "
        + std::to_string( reco_bin_min ) + " && " + branchexpr
        + " < " + std::to_string( reco_bin_max ) + ')';

      chain.Draw( (true_branchexpr + " >> " + true_hist_name).c_str(),
          cuts.c_str(), "goff" );

      true_hist->SetStats( false );
      true_hist->SetLineWidth( 2 );
      true_hist->SetLineColor( kBlack );

      true_hist->Draw( "hist pe" );

      // Prepare vertical lines to draw on the plot. These will show the
      // boundaries of the reco bin in true space
      double max_for_lines = std::numeric_limits<double>::max();

      TLine* line_bin_min = new TLine( reco_bin_min, 0.,
          reco_bin_min, max_for_lines );

      TLine* line_bin_max = new TLine( reco_bin_max, 0.,
          reco_bin_max, max_for_lines );

      line_bin_min->SetLineColor( kRed );
      line_bin_min->SetLineWidth( 2 );
      line_bin_min->SetLineStyle( 1 );
      line_bin_min->Draw( "same" );

      line_bin_max->SetLineColor( kRed );
      line_bin_max->SetLineWidth( 2 );
      line_bin_max->SetLineStyle( 1 );
      line_bin_max->Draw( "same" );

    } // loop over reco bins

  } // show bin plots

  // Also get the total number of reco bins for the 2D smearing plot
  int num_reco_bins = bin_low_edges.size() - 1u;

  // Compute the smearing matrix for a choice of true bins that exactly
  // match the ones in reco space.
  std::string smear_hist_name = "smear_hist" + std::to_string( hist_count );
  TH2D* smear_hist = new TH2D( smear_hist_name.c_str(),
      ("smearing matrix for " + variable_title + "; true " + variable_title
       + "; reco " + variable_title).c_str(), num_reco_bins, bin_low_edges.data(),
      num_reco_bins, bin_low_edges.data() );

  std::string smear_expr = branchexpr + " : " + true_branchexpr
    + " >> " + smear_hist_name;

  std::string smear_cuts = mc_event_weight + " * (is_mc && " + signal_cuts
    + " && " + selection + ')';

  chain.Draw( smear_expr.c_str(), smear_cuts.c_str(), "goff" );

  // Before renormalizing the smearing matrix histogram, take a projection
  // along the reco (y) axis. This will show the expected number of signal
  // events in each reco bin according to our central value MC model. Reco bins
  // should be chosen to have sufficient expected statistics in addition to
  // small smearing.
  TH1D* expected_reco_hist = smear_hist->ProjectionY();

  // Scale the expected reco bin counts to the POT analyzed for the full
  // dataset. Also set the bin stat uncertainties to the square root of their
  // contents. This is not correct for getting the MC statistical uncertainties
  // (which should use the sum of the squares of the weights to get the
  // variance), but we're less interested in those. Primarily we'd like to know
  // what the anticipated statistical uncertainties on the *measurement* will
  // be. We can estimate that by choosing the bin errors in this way. This will
  // help in the effort to choose suitable bins for reporting the final result.
  expected_reco_hist->Scale( EXPECTED_POT / total_simulated_POT );
  for ( int eb = 0; eb <= num_reco_bins + 1; ++eb ) {
    double bin_events = expected_reco_hist->GetBinContent( eb );
    double bin_stat_err = std::sqrt( std::max(0., bin_events) );
    expected_reco_hist->SetBinError( eb, bin_stat_err );
  }

  expected_reco_hist->SetStats( false );
  expected_reco_hist->SetLineColor( kBlack );
  expected_reco_hist->SetLineWidth( 2 );

  std::stringstream temp_ss;
  temp_ss << "expected reco bin counts (" << EXPECTED_POT << " POT);"
    << " reco " << variable_title << "; events";

  expected_reco_hist->SetTitle( temp_ss.str().c_str() );

  TCanvas* c_expected = new TCanvas;
  expected_reco_hist->Draw( "hist e" );

  // Normalize the smearing matrix elements so that a sum over all reco bins
  // (including the under/overflow bins) yields a value of one. This means that
  // every selected signal event must end up somewhere in reco space.
  int num_bins_x = smear_hist->GetXaxis()->GetNbins();
  int num_bins_y = smear_hist->GetYaxis()->GetNbins();

  // Loop over the true (x) bins. Include the underflow (index zero) and
  // overflow (index num_bins_x + 1) bins.
  for ( int bx = 0; bx <= num_bins_x + 1; ++bx ) {

    // For the current true (x) bin, compute the sum of all reco (y) bins.
    double y_sum = 0.;
    for ( int by = 0; by <= num_bins_y + 1; ++by ) {
      y_sum += smear_hist->GetBinContent( bx, by );
    }

    // Normalize each of the reco (y) bins so that the sum over y is unity.
    for ( int by = 0; by <= num_bins_y + 1; ++by ) {

      // To avoid dividing by zero, set the bin content to zero if the sum of
      // the reco (y) bins is not positive.
      if ( y_sum <= 0. ) {
        //smear_hist->SetBinContent( bx, by, REALLY_SMALL );
        smear_hist->SetBinContent( bx, by, 0. );
      }
      else {
        // Otherwise, normalize in the usual way
        double bc = smear_hist->GetBinContent( bx, by );

        double content = std::max( bc / y_sum, REALLY_SMALL );

        smear_hist->SetBinContent( bx, by, content );
      }
    } // loop over reco (y) bins

  } // loop over true (x) bins

  // Smearing matrix histogram style options
  smear_hist->GetXaxis()->SetTitleFont( FONT_STYLE);
  smear_hist->GetYaxis()->SetTitleFont( FONT_STYLE );
  smear_hist->GetXaxis()->SetTitleSize( 0.05 );
  smear_hist->GetYaxis()->SetTitleSize( 0.05 );
  smear_hist->GetXaxis()->SetLabelFont( FONT_STYLE );
  smear_hist->GetYaxis()->SetLabelFont( FONT_STYLE );
  smear_hist->GetZaxis()->SetLabelFont( FONT_STYLE );
  smear_hist->GetZaxis()->SetLabelSize( 0.03 );
  smear_hist->GetXaxis()->CenterTitle();
  smear_hist->GetYaxis()->CenterTitle();
  smear_hist->GetXaxis()->SetTitleOffset( 1.2 );
  smear_hist->GetYaxis()->SetTitleOffset( 1.1 );
  smear_hist->SetStats( false );
  smear_hist->SetMarkerSize( 1.8 ); // text size
  smear_hist->SetMarkerColor( kWhite ); // text color

  // Draw the smearing matrix plot
  TCanvas* c_smear = new TCanvas;
  c_smear->SetBottomMargin( 0.15 );
  c_smear->SetLeftMargin( 0.13 );

  if ( show_smear_numbers ) {
    // Round all numbers to this precision when rendering them
    gStyle->SetPaintTextFormat( "4.2f" );

    smear_hist->Draw("text colz");
  }
  else {
    smear_hist->Draw( "colz" );
  }

  // For each true bin, print the fraction of events that are reconstructed
  // in the correct corresponding reco bin.
  printf("Diagonal of 1D %dx%d smear matrix of ", num_reco_bins, num_reco_bins);
  std::cout << variable_title;
  printf(" is as follows.\n");
  printf("index  | Low Edge | Diag \n");
  for ( int bb = 1; bb <= num_reco_bins; ++bb ) {
    printf("bin #%2d: %6.3f, %8.3f \n", bb, expected_reco_hist->GetBinLowEdge( bb ), smear_hist->GetBinContent(bb, bb));
  }

}

// Overloaded version that uses a fixed number of equal-width bins
void MakeConfig::make_res_plots( const std::string& branchexpr,
    const std::string& variable_title, const std::string& selection,
    const std::set<int>& runs,
    double xmin, double xmax, int Nbins,
    bool show_bin_plots,
    bool show_smear_numbers,
    int num_true_bins,
    const std::string& mc_branchexpr,
    const std::string& signal_cuts,
    const std::string& mc_event_weight)
{
  auto low_edges = get_bin_low_edges( xmin, xmax, Nbins );
  return make_res_plots( branchexpr, variable_title, selection, runs,
      low_edges, show_bin_plots, show_smear_numbers, num_true_bins,
      mc_branchexpr, signal_cuts, mc_event_weight );
}

void MakeConfig::make_res_plots( std::istream& in_stream,
  const std::set<int>& runs, const std::string& universe_branch_name,
  size_t universe_index, bool show_smear_numbers )
{
  const std::string variable_title = "bin";

  // Create a UniverseMaker object that will handle the actual
  // calculation of the smearing matrix
  UniverseMaker um( in_stream );

  // Get access to the singleton utility class that manages the processed
  // ntuple files
  const FilePropertiesManager& fpm = FilePropertiesManager::Instance();

  // TODO: Reduce code duplication for the POT tallying

  // Add the appropriate CV numu ntuples from the requested run(s) to the
  // TChain owned by the UniverseMaker object. For the resolution
  // studies, this is all we really need. Also tally the total simulated POT
  // for later scaling purposes.
  double total_simulated_POT = 0.;

  const auto& ntuple_map = fpm.ntuple_file_map();
  for ( const auto& run : runs ) {
    const auto& ntuple_files = ntuple_map.at( run )
      .at( NtupleFileType::kNumuMC );
    for ( const auto& file_name : ntuple_files ) {
      um.add_input_file( file_name );
      std::cout << file_name << '\n';

      TFile temp_file( file_name.c_str(), "read" );
      TParameter<float>* temp_pot = nullptr;
      temp_file.GetObject( "summed_pot", temp_pot );
      double pot = temp_pot->GetVal();
      total_simulated_POT += pot;
    }
  }


  // Look up the MC event weights from the input files and construct the
  // response matrices in the usual way. For speed, restrict the calculation to
  // just the universe branch requested by the user (typically the CV branch).
  um.build_universes( { universe_branch_name } );

  // For all but the "unweighted" universe, the key used to look up
  // the map entry is "weight_" prepended to the original ntuple branch name.
  std::string universe_key = "unweighted";
  if ( universe_key != universe_branch_name ) {
    universe_key = "weight_" + universe_branch_name;
  }

  std::cout << universe_key << " : " << universe_index << '\n';
  // Get access to the Universe object that stores the histograms of summed MC
  // event weights that we need
  const auto& universe = um.universe_map().at( universe_key )
    .at( universe_index );

  TH2D* smear_hist = dynamic_cast< TH2D* >(
      universe.hist_2d_->Clone("smear_hist") );

  // TODO: also reduce code duplication here

  // Before renormalizing the smearing matrix histogram, take a projection
  // along the reco (y) axis. This will show the expected number of signal
  // events in each reco bin according to our central value MC model. Reco bins
  // should be chosen to have sufficient expected statistics in addition to
  // small smearing.
  TH1D* expected_reco_hist = smear_hist->ProjectionY();
  int num_reco_bins = expected_reco_hist->GetNbinsX();

  // Scale the expected reco bin counts to the POT analyzed for the full
  // dataset. Also set the bin stat uncertainties to the square root of their
  // contents. This is not correct for getting the MC statistical uncertainties
  // (which should use the sum of the squares of the weights to get the
  // variance), but we're less interested in those. Primarily we'd like to know
  // what the anticipated statistical uncertainties on the *measurement* will
  // be. We can estimate that by choosing the bin errors in this way. This will
  // help in the effort to choose suitable bins for reporting the final result.
  expected_reco_hist->Scale( EXPECTED_POT / total_simulated_POT );
  for ( int eb = 0; eb <= num_reco_bins + 1; ++eb ) {
    double bin_events = expected_reco_hist->GetBinContent( eb );
    double bin_stat_err = std::sqrt( std::max(0., bin_events) );
    expected_reco_hist->SetBinError( eb, bin_stat_err );
  }

  expected_reco_hist->SetStats( false );
  expected_reco_hist->SetLineColor( kBlack );
  expected_reco_hist->SetLineWidth( 2 );

  std::stringstream temp_ss;
  temp_ss << "expected reco bin counts (" << EXPECTED_POT << " POT);"
    << " reco " << variable_title << "; events";

  expected_reco_hist->SetTitle( temp_ss.str().c_str() );

  TCanvas* c_expected = new TCanvas;
  expected_reco_hist->Draw( "hist e" );

  // Normalize the smearing matrix elements so that a sum over all reco bins
  // (including the under/overflow bins) yields a value of one. This means that
  // every selected signal event must end up somewhere in reco space.
  int num_bins_x = smear_hist->GetXaxis()->GetNbins();
  int num_bins_y = smear_hist->GetYaxis()->GetNbins();

  // Loop over the true (x) bins. Include the underflow (index zero) and
  // overflow (index num_bins_x + 1) bins.
  for ( int bx = 0; bx <= num_bins_x + 1; ++bx ) {

    // For the current true (x) bin, compute the sum of all reco (y) bins.
    double y_sum = 0.;
    for ( int by = 0; by <= num_bins_y + 1; ++by ) {
      y_sum += smear_hist->GetBinContent( bx, by );
    }

    // Normalize each of the reco (y) bins so that the sum over y is unity.
    for ( int by = 0; by <= num_bins_y + 1; ++by ) {

      // To avoid dividing by zero, set the bin content to zero if the sum of
      // the reco (y) bins is not positive.
      if ( y_sum <= 0. ) {
        //smear_hist->SetBinContent( bx, by, REALLY_SMALL );
        smear_hist->SetBinContent( bx, by, 0. );
      }
      else {
        // Otherwise, normalize in the usual way
        double bc = smear_hist->GetBinContent( bx, by );

        double content = std::max( bc / y_sum, REALLY_SMALL );

        smear_hist->SetBinContent( bx, by, content );
      }
    } // loop over reco (y) bins

  } // loop over true (x) bins

  // Smearing matrix histogram style options
  smear_hist->GetXaxis()->SetTitleFont( FONT_STYLE);
  smear_hist->GetYaxis()->SetTitleFont( FONT_STYLE );
  smear_hist->GetXaxis()->SetTitleSize( 0.05 );
  smear_hist->GetYaxis()->SetTitleSize( 0.05 );
  smear_hist->GetXaxis()->SetLabelFont( FONT_STYLE );
  smear_hist->GetYaxis()->SetLabelFont( FONT_STYLE );
  smear_hist->GetZaxis()->SetLabelFont( FONT_STYLE );
  smear_hist->GetZaxis()->SetLabelSize( 0.03 );
  smear_hist->GetXaxis()->CenterTitle();
  smear_hist->GetYaxis()->CenterTitle();
  smear_hist->GetXaxis()->SetTitleOffset( 1.2 );
  smear_hist->GetYaxis()->SetTitleOffset( 1.1 );
  smear_hist->SetStats( false );
  smear_hist->SetMarkerSize( 1.8 ); // text size
  smear_hist->SetMarkerColor( kWhite ); // text color

  // Draw the smearing matrix plot
  TCanvas* c_smear = new TCanvas;
  c_smear->SetBottomMargin( 0.15 );
  c_smear->SetLeftMargin( 0.13 );

  if ( show_smear_numbers ) {
    // Round all numbers to this precision when rendering them
    gStyle->SetPaintTextFormat( "4.2f" );

    smear_hist->Draw("text colz");
  }
  else {
    smear_hist->Draw( "colz" );
  }

  // Draw a vertical red line to separate the signal true bins from the
  // background true bins in the smearing matrix plot. Start by finding
  // where the first background bin is. Here we assume that the full set
  // of signal bins comes before any background bins.
  const auto& true_bins = um.true_bins();
  size_t num_true_bins = true_bins.size();
  size_t first_bkgd_bin_idx = num_true_bins;
  for ( size_t t = 0u; t < num_true_bins; ++t ) {
    const auto& tbin = true_bins.at( t );
    if ( tbin.type_ == TrueBinType::kBackgroundTrueBin ) {
      first_bkgd_bin_idx = t;
      break;
    }
  }

  // We've found the bin index. Now draw the line to indicate the
  // signal/background boundary in true space
  TLine* bkgd_line = new TLine( first_bkgd_bin_idx, 0.,
      first_bkgd_bin_idx, num_reco_bins );

  bkgd_line->SetLineColor( kRed );
  bkgd_line->SetLineWidth( 3 );
  bkgd_line->SetLineStyle( 2 );
  bkgd_line->Draw( "same" );

  // For each true bin, print the fraction of events that are reconstructed
  // in the correct corresponding reco bin.
  //for ( int bb = 1; bb <= num_reco_bins; ++bb ) {
  //  std::cout << "bin #" << bb << ": "
  //    << expected_reco_hist->GetBinLowEdge( bb ) << ", "
  //    << smear_hist->GetBinContent(bb, bb) << '\n';
  //}

  // For each true bin, print the fraction of events that are reconstructed
  // in the correct corresponding reco bin.
  printf("Diagonal of 2D %dx%d smear matrix of ", num_reco_bins, num_reco_bins);
  std::cout << variable_title;
  printf(" is as follows.\n");
  printf("index  | Low Edge | Diag \n");
  for ( int bb = 1; bb <= num_reco_bins; ++bb ) {
    printf("bin #%2d: %6f, %8.3f \n", bb, expected_reco_hist->GetBinLowEdge( bb ), smear_hist->GetBinContent(bb, bb));
  }


}




void Block1D::Init(){
  this->SetTitle( fTitle );
  this->SetName( fName );
  this->SetTexTitle( fTexTitle );

  for ( size_t b = 0u; b < fblock.size() - 1; ++b ) {
    double low = fblock.at( b );
    double high = fblock.at( b + 1u );
    std::string bin_def = "";
    if ( fselection.size() != 0 ) {
      bin_def = fselection + " && ";
    }
    bin_def += xName + Form(" >= %f && ", low) + xName + Form(" < %f", high);
    binDef.push_back( bin_def );
  }
}

void Block1D::SetTitle( const std::string& title ) {
  TString temp_title = title;
  temp_title.ReplaceAll("#;",2,"#semicolon",10);
  fTitle = temp_title;

  std::vector< TString > str_container;
  TString str1 = fTitle;
  Int_t isc = str1.Index(";");
  Int_t lns = str1.Length();
  while(isc >=0){
    str_container.push_back(str1(0,isc));
    str1 = str1(isc+1, lns);
    isc = str1.Index(";");
    lns = str1.Length();
  }
  str_container.push_back(str1);
  if(str_container.size() == 2u){
    xTitle = str_container[0];
    xTitleUnit = str_container[1];
    yTitle = "Events";
    yTitleUnit = "";
    str_container.clear();
  }
  else if(str_container.size() == 1u){
    xTitle = str_container[0];
    xTitleUnit = "";
    yTitle = "Events";
    yTitleUnit = "";
    str_container.clear();
  }
  else throw std::runtime_error( "Wrong title -> " + fTitle
      + ". The format of the title of 1D block must be <branch title>;"
      + "<unit>." );
}



void Block1D::SetTexTitle( const std::string& textitle ) {
  TString temp_textitle = textitle;
  temp_textitle.ReplaceAll("#;",2,"#semicolon",10);
  fTexTitle = temp_textitle;
  std::vector<TString> str_container;
  TString str1 = fTexTitle;
  Int_t isc = str1.Index(";");
  Int_t lns = str1.Length();
  while(isc >=0){
    str_container.push_back(str1(0,isc));
    str1 = str1(isc+1, lns);
    isc = str1.Index(";");
    lns = str1.Length();
  }
  str_container.push_back(str1);
  if(str_container.size() == 2u){
    xTexTitle = str_container[0];
    xTexTitleUnit = str_container[1];
    yTexTitle = "Events";
    yTexTitleUnit = "";
    str_container.clear();
  }
  else if(str_container.size() == 1u){
    xTexTitle = str_container[0];
    xTexTitleUnit = "";
    yTexTitle = "Events";
    yTexTitleUnit = "";
    str_container.clear();
  }
  else
    throw std::runtime_error("Wrong title -> " + fTexTitle + ". The format of the title of 1D block must be <branch title>; <unit>.");

}



void Block1D::SetName( const std::string& name ) {
  TString temp_name = name;
  temp_name.ReplaceAll("#;",2,"#semicolon",10);
  fName = temp_name;

  std::vector<TString> str_container;
  TString str1 = fName;
  Int_t isc = str1.Index(";");
  Int_t lns = str1.Length();
  while(isc >=0){
    str_container.push_back(str1(0,isc));
    str1 = str1(isc+1, lns);
    isc = str1.Index(";");
    lns = str1.Length();
  }
  str_container.push_back(str1);
  if(str_container.size() == 2u){
    xName = str_container[0];
    xNameUnit = str_container[1];
    yName = "Events";
    yNameUnit = "";
    str_container.clear();
  }
  else if(str_container.size() == 1u){
    xName = str_container[0];
    xNameUnit = "";
    yName = "Events";
    yNameUnit = "";
    str_container.clear();
  }
  else
    throw std::runtime_error("Wrong name -> " + fName + ". The format of the name of 1D block must be <branch name>; <unit>.");
}

void Block2D::Init(){
  this->SetTitle( fTitle );
  this->SetName( fName );
  this->SetTexTitle( fTexTitle );

  for (auto iter = fblock.cbegin(); iter != fblock.cend(); ++iter ) {
    // Get an iterator to the map element after the current one. Due to the
    // automatic sorting, this is guaranteed to contain the upper edge of the
    // current delta_pT bin
    xbin.push_back(iter->first);
    auto next = iter;
    ++next;
    if(next == fblock.cend()) continue;
    double slice_low = iter->first;
    double slice_high = next->first;

    fblock_vv_.push_back(iter->second);

    for(size_t b = 0u; b < iter->second.size() - 1; b++){
      double bin_low = iter->second.at(b);
      double bin_high = iter->second.at(b + 1u);
      std::string bin_def = "";
      if ( fselection.size() != 0 ) {
        bin_def = fselection + " && ";
      }
      bin_def += xName + Form(" >= %f && ", slice_low)
        + xName + Form(" < %f && ", slice_high) + yName
        + Form(" >= %f && ", bin_low) + yName + Form(" < %f ", bin_high);
      binDef.push_back( bin_def );
    }
  }
}


void Block2D::SetName( const std::string& name ) {
  TString temp_name = name;
  temp_name.ReplaceAll("#;",2,"#semicolon",10);
  fName = temp_name;

  std::vector<TString> str_container;

  TString str1 = fName;
  Int_t isc = str1.Index(";");
  Int_t lns = str1.Length();
  if(isc < 0)
    throw std::runtime_error("Wrong name -> " + fName + ". The format of the name of 2D block must be <branch name>; <unit>; <y branch name>; <y unit>.");
  while(isc >=0){
    str_container.push_back(str1(0,isc));
    str1 = str1(isc+1, lns);
    isc = str1.Index(";");
    lns = str1.Length();
  }
  str_container.push_back(str1);
  if(str_container.size() == 4u){
    xName = str_container[0];
    xNameUnit = str_container[1];
    yName =  str_container[2];
    yNameUnit =  str_container[3];
    str_container.clear();
  }
  else if(str_container.size() == 2u){
    xName = str_container[0];
    yName =  str_container[1];
    xNameUnit = "";
    yNameUnit = "";
    str_container.clear();
  }
  else
    throw std::runtime_error("Wrong name -> " + fName + ". The format of the name of 2D block must be <branch name>; <unit>; <y branch name>; <y unit>.");
}



void Block2D::SetTitle( const std::string& title ) {
  TString temp_title = title;
  temp_title.ReplaceAll("#;",2,"#semicolon",10);
  fTitle = temp_title;

  std::vector<TString> str_container;

  TString str1 = fTitle;
  Int_t isc = str1.Index(";");
  Int_t lns = str1.Length();
  if ( isc < 0 ) throw std::runtime_error("Wrong title -> " + fTitle
      + ". The format of the title of 2D block must be <branch title>;"
      + " <unit>; <y branch title>; <y unit>." );
  while ( isc >=0 ) {
    str_container.push_back(str1(0,isc));
    str1 = str1(isc+1, lns);
    isc = str1.Index(";");
    lns = str1.Length();
  }
  str_container.push_back(str1);
  if(str_container.size() == 4u){
    xTitle = str_container[0];
    xTitleUnit = str_container[1];
    yTitle =  str_container[2];
    yTitleUnit =  str_container[3];
    str_container.clear();
  }
  else if(str_container.size() == 2u){
    xTitle = str_container[0];
    yTitle =  str_container[1];
    xTitleUnit = "";
    yTitleUnit = "";
    str_container.clear();
  }
  else throw std::runtime_error( "Wrong title -> " + fTitle
    + ". The format of the title of 2D block must be <branch title>;"
    + " <unit>; <y branch title>; <y unit>." );
}


void Block2D::SetTexTitle( const std::string& textitle ) {
  TString temp_textitle = textitle;
  temp_textitle.ReplaceAll("#;",2,"#semicolon",10);
  fTexTitle = temp_textitle;

  std::vector<TString> str_container;

  TString str1 = fTexTitle;
  Int_t isc = str1.Index(";");
  Int_t lns = str1.Length();
  if( isc < 0 ) return;

  while ( isc >=0 ) {
    str_container.push_back(str1(0,isc));
    str1 = str1(isc+1, lns);
    isc = str1.Index(";");
    lns = str1.Length();
  }
  str_container.push_back(str1);
  if(str_container.size() == 4u){
    xTexTitle = str_container[0];
    xTexTitleUnit = str_container[1];
    yTexTitle =  str_container[2];
    yTexTitleUnit =  str_container[3];
    str_container.clear();
  }
  else if(str_container.size() == 2u){
    xTexTitle = str_container[0];
    yTexTitle =  str_container[1];
    xTexTitleUnit = "";
    yTexTitleUnit = "";
    str_container.clear();
  }
  else throw std::runtime_error( "Wrong title -> " + fTexTitle
    + ". The format of the title of 2D block must be <branch title>;"
    + " <unit>; <y branch title>; <y unit>." );
}

void MakeConfig::BinScheme() {

  // Instantiate the request binning scheme using the factory
  BinSchemeFactory bsf;
  BinSchemeBase* bs = bsf.CreateBinScheme( bin_scheme_name_ );
  bin_scheme_.reset( bs );

  // The name of a TDirectoryFile which will store all of the histograms within
  // the output ROOT file
  DIRECTORY = bin_scheme_->out_tdir_name_;

  // Prefix of output bin configure file and slice configure file
  BIN_CONFIG = bin_scheme_->out_config_prefix_;

  SELECTION = bin_scheme_->selection_name_;

  TREE = bin_scheme_->ntuple_ttree_name_;

  // Runs used to plot smearing matrix
  RUNS = bin_scheme_->runs_to_use_;

  vect_block = &bin_scheme_->vect_block;
}
